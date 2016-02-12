package org.cbio.mutex;

import org.cbio.causality.data.GeneCards;
import org.cbio.causality.util.FDR;
import org.cbio.causality.util.Kronometre;

import java.io.*;
import java.util.*;

/**
 * This is the main executor class of the mutex search.
 *
 * @author Ozgun Babur
 */
public class Main
{
	/**
	 * False discovery rate threshold.
	 */
	public static double fdrThr = -1;

	/**
	 * Maximum group size to use during the searches.
	 */
	public static int maxGroupSize = 5;

	/**
	 * Maximum number of iterations for estimating the null distribution of initial p-values.
	 */
	public static int randIter1 = 10000;

	/**
	 * Number of iterations to estimate the null distribution of final group scores.
	 */
	public static int randIter2 = 100;

	/**
	 * Directory that contains parameters file.
	 */
	public static String dir;

	/**
	 * Whether to reduce the search space using signaling networks.
	 */
	private static boolean useGraph = true;

	/**
	 * The name of the tab-delimited file containing gene alterations.
	 */
	public static String dataFileName;

	/**
	 * Name of the signaling network file. No need to specify this to use the default network.
	 */
	private static String networkFilename;

	/**
	 * The signaling network.
	 */
	private static Network network;

	/**
	 * Users can limit the search to certain genes using this file.
	 */
	private static String symbolsFile;

	/**
	 * While producing results, if there is an available sub-typing information, nodes are
	 * highlighted according to the subtype.
	 */
	private static String subtypeDatasetName;

	/**
	 * Keywords associated with the dataset, separated by comma.
	 */
	private static String literatureKeywords;

	/**
	 * Minimum nnumber of altered samples for a gene to be included in the study.
	 */
	private static Integer minAltCntThr;

	/**
	 * Number of genes to limit the study. This is useful to restrict the analysis to top X most significantly
	 * altered genes.
	 */
	private static Integer geneLimit;

	/**
	 * A file to provide ordering of genes, higher priority first.
	 */
	private static String geneRankingFile;

	/**
	 * Tells even if there is need of random runs to estimate FDR, don't do it and just write result
	 * without FDR estimation.
	 */
	private static boolean noRandomRun = false;

	/**
	 * Parameters for auto-downloading data matrix from cBioPortal.
	 */
	private static String portalStudyID;
	private static String portalCaseListID;
	private static String portalExpProfileID;
	private static String portalCNAProfileID;
	private static String portalMutProfileID;
	private static Double minAltRatio;

	/**
	 * Parameter to run the analysis on a randomized set of alterations.
	 */
	private static boolean randomizeDataMatrix = false;

	public static void main(String[] args) throws IOException, ClassNotFoundException
	{
		Kronometre kron = new Kronometre();

		if (args.length < 1)
		{
			displayHelp();
			return;
		}
		dir = args[0];
		if (!dir.endsWith(File.separator)) dir += File.separator;

		if (!loadParameters()) System.exit(1);

		if (useGraph) network = networkFilename == null ?
			new Network() : new Network(networkFilename);

		if (args.length > 1 && args[1].equals("random"))
		{
			int howMany = 1;

			if (args.length > 2)
			{
				howMany = Integer.parseInt(args[2]);
				System.out.println("howMany = " + howMany);
			}

			generateRamdomRun(howMany);
		}
		else
		{
			noRandomRun = args.length > 1 && args[1].equals("no-random-run");
			search();
//			searchOnRandomized();
		}

		kron.stop();
		kron.print();
	}

	public static void reset()
	{
		fdrThr = -1;
		maxGroupSize = 5;
		randIter1 = 10000;
		randIter2 = 100;
		dir = null;
		useGraph = true;
		dataFileName = null;
		networkFilename = null;
		network = null;
		symbolsFile = null;
		subtypeDatasetName = null;
		literatureKeywords = null;
		portalStudyID = null;
		portalCaseListID = null;
		portalExpProfileID = null;
		portalCNAProfileID = null;
		portalMutProfileID = null;
		randomizeDataMatrix = false;
		minAltCntThr = null;
		geneLimit = null;
		geneRankingFile = null;
	}

	/**
	 * Makes a run for generating the null distribution of final group scores.
	 * @param howMany number of iterations for this run
	 * @throws IOException
	 */
	public static void generateRamdomRun(int howMany) throws IOException
	{
		// load the alteration data
		Map<String, GeneAlt> genesMap = loadAlterations();

		MutexGreedySearcher searcher = new MutexGreedySearcher(genesMap, network);
		Set<String> symbols = genesMap.keySet();
		Set<String> noShuffle = loadHighlySignificantGenes();
		generateRandomPvals(searcher, symbols, noShuffle, null, howMany);
	}

	/**
	 *
	 * @throws IOException
	 * @throws ClassNotFoundException
	 */
	public static void search() throws IOException, ClassNotFoundException
	{
		System.out.println("----------------------------------------\n");
		System.out.println("Directory = " + dir);

		// load the alteration data
		Map<String, GeneAlt> genesMap = readCache("data-cache");
		if (genesMap == null) genesMap = loadAlterations();

		// if the data needs to be downloaded from cBioPortal, do it
		if (genesMap == null && portalStudyID != null)
		{
			downloadPortalData();
			loadParameters();
			genesMap = loadAlterations();
		}

		if (randomizeDataMatrix)
		{
			System.out.print("Randomizing data matrix ... ");
			for (GeneAlt gene : genesMap.values())
			{
				gene.shufflePermanent();
			}
			System.out.println("done");
		}

		System.out.println("Number of genes = " + genesMap.size());
		System.out.println("Number of samples = " + genesMap.values().iterator().next().size());

		MutexGreedySearcher searcher = new MutexGreedySearcher(genesMap, network);

		Set<String> symbols = genesMap.keySet();
		if (network != null) symbols.retainAll(network.getSymbols());

		Map<String, Group> groupsOfSeeds = searcher.getGroupsOfSeeds(symbols, maxGroupSize,
			randIter1);
		cacheData(genesMap, "data-cache");

		writeRankedGroups(groupsOfSeeds, null, "ranked-groups.txt");

		// we are done if we won't cutoff from an fdr
		if (randIter2 <= 0) return;

		// Load and/or generate final scores null distribution
		List<Double> nullDist = new ArrayList<Double>();
		int cnt = readRandomPvals(nullDist);
		if (cnt < randIter2)
		{
			if (noRandomRun) return;

			generateRandomPvals(searcher, symbols, loadHighlySignificantGenes(), nullDist,
				randIter2 - cnt);
		}

		writeRankedGroups(groupsOfSeeds, nullDist, "ranked-groups.txt");

		// Apply FDR cutoff
		Map<String, Double> resultScores = new HashMap<String, Double>();
		for (String id : groupsOfSeeds.keySet())
		{
			resultScores.put(id, groupsOfSeeds.get(id).calcFinalScore());
		}

		if (dir.contains("simulation1")) Simulation.plotEstimatedVsActualFDR(groupsOfSeeds, resultScores, nullDist, randIter2, 830, 150);
		else if (dir.contains("simulation2")) Simulation.plotEstimatedVsActualFDR(groupsOfSeeds, resultScores, nullDist, randIter2, 156, 60);

		if (fdrThr < 0)
		{
			fdrThr = decideBestFDR(resultScores, nullDist, randIter2);
		}

		System.out.println("Selected FDR = " + fdrThr);
		if (fdrThr < 0) return;

		List<String> selectedSeeds = FDR.select(resultScores, fdrThr, nullDist, randIter2);

		// Find threshold score
		double thrScore = 0;
		for (String seed : selectedSeeds)
		{
			if (resultScores.get(seed) > thrScore) thrScore = resultScores.get(seed);
		}

		List<Group> groups = new ArrayList<Group>(selectedSeeds.size());
		for (String seed : selectedSeeds)
		{
			Group group = groupsOfSeeds.get(seed);
			if (useGraph) group.fetchTragets(network, genesMap);
			groups.add(group);
		}

		System.out.println("Number of mutex groups in results = " + groups.size());

		// clean subsets in the result
		Group.removeSubsets(groups);

		System.out.println("Groups after removing subsets = " + groups.size());

		// sort the list favoring high-coverage
		Group.sortToCoverage(groups);

		SubtypeAligner sa = null;
//		SubtypeAligner sa = subtypeDatasetName == null ? null :
//			new SubtypeAligner(PortalDatasetEnum.find(subtypeDatasetName), Group.collectGenes(groups));

		// Write the output graph to visualize in ChiBE
		if (useGraph) GraphWriter.write(groups, genesMap, network, dir, dir, sa);


		BufferedWriter writer = new BufferedWriter(new FileWriter(dir + "oncoprint.txt"));

		// Print textual results
		for (Group group : groups)
		{
			writer.write(group.getPrint(sa, null, true, useGraph) + "\n\n");
		}
		writer.close();

//		if (literatureKeywords != null) printAnnotations(getGenes(groups, genesMap));
	}

	public static void searchOnRandomized() throws IOException, ClassNotFoundException
	{
		System.out.println("----------------------------------------\n");
		System.out.println("Directory = " + dir);

		// load the alteration data
		Map<String, GeneAlt> genesMap = loadAlterations();

		System.out.print("Randomizing data matrix ... ");
		for (GeneAlt gene : genesMap.values())
		{
			gene.shufflePermanent();
		}
		System.out.println("done");

		System.out.println("Number of genes = " + genesMap.size());
		System.out.println("Number of samples = " + genesMap.values().iterator().next().size());

		MutexGreedySearcher searcher = new MutexGreedySearcher(genesMap, network);

		Set<String> symbols = genesMap.keySet();

		Map<String, Group> groupsOfSeeds = searcher.getGroupsOfSeeds(symbols, maxGroupSize,
			randIter1);

		cacheData(genesMap, "random-cache");
		writeRankedGroups(groupsOfSeeds, null, "ranked-groups-random.txt");
	}

		private static List<String> getGenes(List<Group> groups,
		final Map<String, GeneAlt> genesMap)
	{
		Set<String> genes = new HashSet<String>();
		for (Group group : groups)
		{
			genes.addAll(group.getGeneNames());
		}
		List<String> list = new ArrayList<String>(genes);

		Collections.sort(list, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return new Integer(genesMap.get(o2).getAltCnt()).compareTo(
					genesMap.get(o1).getAltCnt());
			}
		});

		return list;
	}

	public static void printAnnotations(List<String> genes)
	{
		System.out.println("Number genes in mutex groups = " + genes.size());
		System.out.println(genes + "\n");

		String[] kywd = literatureKeywords.split(",");
		for (int i = 0; i < kywd.length; i++)
		{
			kywd[i] = kywd[i].trim();
		}

		List<String> withKw = new ArrayList<String>();
		List<String> withOtherCancer = new ArrayList<String>(genes);
		List<String> notAssoc = new ArrayList<String>();

		System.out.println("Reported cancer associations in GeneCards database\n" +
			"--------------------------------------------------\n");

		for (String gene : genes)
		{
			System.out.print(gene + ": ");
			Set<String> relatedCancers = GeneCards.getRelatedCancers(gene);
			for (String cancer : relatedCancers)
			{
				System.out.print(cancer + ", ");
				if (contains(cancer, kywd) && !withKw.contains(gene)) withKw.add(gene);
			}
			System.out.println();
			if (relatedCancers.isEmpty()) notAssoc.add(gene);
		}
		withOtherCancer.removeAll(withKw);
		withOtherCancer.removeAll(notAssoc);

		System.out.println("\nClassification summary\n----------------------");

		System.out.println("\nAlready associated with \"" + literatureKeywords + "\": " + withKw.size());
		System.out.println(withKw);
		System.out.println("\nAssociated with other types of cancer: " + withOtherCancer.size());
		System.out.println(withOtherCancer);
		System.out.println("\nNot associated with any cancer: " + notAssoc.size());
		System.out.println(notAssoc);
		System.out.println("\n\n");
	}

	private static boolean contains(String s, String[] q)
	{
		for (String x : q)
		{
			if (!s.isEmpty() && s.contains(x)) return true;
		}
		return false;
	}

	private static void displayHelp()
	{
		System.out.println("Please provide the directory that contains the parameters.txt file " +
			"as the first program argument. In this file, the parameter keys and the values " +
			"should be separated with the = sign. parameters.txt file can contain the below " +
			"parameters:\n\n");

		System.out.println(
			"data-file: Name of the data file. Mandatory.\n\n" +
			"max-group-size: The maximum size of a result mutex group. Integer value. Default is 5.\n\n" +
			"first-level-random-iteration: Number of randomization to estimate null distribution of member p-values in mutex groups. Integer. Default is 10000.\n\n" +
			"second-level-random-iteration: Number of runs to estimate the null distribution of final scores. Integer. Default is 100. If FDR control on results is not required and only the ranking of the result groups is sufficient, set this parameter to 0.\n\n" +
			"fdr-cutoff: Users can select a specific FDR cutoff. When not provided, or when set to a negative value, the FDR cutoff that maximizes the expected value of true positives - false positives is used.\n\n" +
			"search-on-signaling-network: Whether to reduce the search space using the signaling network. true or false. Default is true.\n\n" +
			"genes-file: This parameter can be used to limit the search to a subset of genes. The file should contain a gene symbol per line.\n\n" +
			"network-file: To customize the signaling network, users can use this parameter. The tab-delimited network file should contain 3 columns (Gene Symbol 1, interaction-type, Gene Symbol 2).");
	}

	private static void cacheData(Map<String, GeneAlt> geneMap, String filename) throws IOException
	{
		ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(dir + filename));
		out.writeObject(geneMap);
		out.close();
	}

	public static Map<String, GeneAlt> readCache(String filename) throws IOException, ClassNotFoundException
	{
		File file = new File((dir == null ? "" : dir + File.separator) + filename);
		if (file.exists())
		{
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(file));
			Map<String, GeneAlt> map = (Map<String, GeneAlt>) in.readObject();
			in.close();
			return map;
		}
		return null;
	}

	public static Map<String, GeneAlt> loadAlterations() throws IOException
	{
		if (dataFileName == null) return null;
		if (!new File(dataFileName).exists()) return null;

		Map<String, GeneAlt> map = new HashMap<String, GeneAlt>();
		BufferedReader reader = new BufferedReader(new FileReader(dataFileName));

		// skip header
		reader.readLine();

		for (String line = reader.readLine(); line != null; line = reader.readLine())
		{
			String[] token = line.split("\t");
			GeneAlt gene = new GeneAlt(token);

			if (minAltCntThr != null)
			{
				int altCnt = gene.countAltered();
				if (altCnt < minAltCntThr) continue;
			}

			map.put(gene.id, gene);
		}

		reader.close();
		crop(map);

		if (geneLimit != null && map.size() > geneLimit)
		{
			List<String> ranking = readGeneRanking();
			if (ranking != null)
			{
				Collections.reverse(ranking);
				Iterator<String> iter = ranking.iterator();
				while (map.size() > geneLimit && iter.hasNext())
				{
					map.remove(iter.next());
				}
			}
			else
			{
				List<String> names = new ArrayList<String>(map.keySet());
				final Map<String, Integer> cnt = new HashMap<String, Integer>();
				for (String name : map.keySet())
				{
					cnt.put(name, map.get(name).countAltered());
				}
				Collections.sort(names, new Comparator<String>()
				{
					@Override
					public int compare(String o1, String o2)
					{
						return cnt.get(o2).compareTo(cnt.get(o1));
					}
				});
				int thr = cnt.get(names.get(geneLimit + 1));
				System.out.println("Alteration cnt threshold > " + thr);
				for (String name : cnt.keySet())
				{
					if (cnt.get(name) <= thr) map.remove(name);
				}
			}
		}

		return map;
	}

	/**
	 * Crops data to the user-provided symbols.
	 */
	private static void crop(Map<String, GeneAlt> map) throws FileNotFoundException
	{
		if (symbolsFile == null) return;

		Set<String> symbols = readSymbolsFile();

		Set<String> remove = new HashSet<String>(map.keySet());
		remove.removeAll(symbols);
		for (String symbol : remove)
		{
			map.remove(symbol);
		}
	}

	private static List<String> readGeneRanking() throws FileNotFoundException
	{
		if (geneRankingFile == null) return null;
		if (!new File(geneRankingFile).exists()) return null;

		List<String> list = new ArrayList<String>();
		Scanner sc = new Scanner(new File(geneRankingFile));
		sc.nextLine();
		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			if (token.length > 0) list.add(token[0]);
		}
		return list;
	}

	private static Set<String> readSymbolsFile() throws FileNotFoundException
	{
		Set<String> symbols = new HashSet<String>();

		Scanner sc = new Scanner(new File(symbolsFile));
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			if (line.startsWith("#")) continue;
			for (String s : line.split("\\s+"))
			{
				if (!s.isEmpty()) symbols.add(s);
			}
		}
		return symbols;
	}

	/**
	 * Writes ordered mutex groups in a file.
	 */
	private static void writeRankedGroups(Map<String, Group> groupMap, List<Double> nullDist, String filename) throws IOException
	{
		final Map<String, Double> scoreMap = new HashMap<String, Double>();
		for (String id : groupMap.keySet())
		{
			scoreMap.put(id, groupMap.get(id).calcFinalScore());
		}

		List<String> sorted = new ArrayList<String>(groupMap.keySet());
		Collections.sort(sorted, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return scoreMap.get(o1).compareTo(scoreMap.get(o2));
			}
		});

		BufferedWriter writer = new BufferedWriter(new FileWriter(dir + filename));

		writer.write(nullDist == null? "Score\tMembers" : "Score\tq-val\tMembers");
		int i = 0;
		for (String seed : sorted)
		{
			Group group = groupMap.get(seed);
			i++;

			double score = group.calcFinalScore();
			writer.write("\n" + score);

			if (nullDist != null)
			{
				double qval = (countFalsePositive(nullDist, score) / (double) randIter2) / i;
				writer.write("\t" + qval);
			}

			for (String name : group.getGeneNames())
			{
				writer.write("\t" + name);
			}
		}

		writer.close();
	}

	private static int countFalsePositive(List<Double> nullDist, double thr)
	{
		int cnt = 0;
		for (Double v : nullDist)
		{
			if (v <= thr) cnt++;
		}
		return cnt;
	}

	private static int readRandomPvals(List<Double> vals) throws FileNotFoundException
	{
		String directory = dir + "randscores/";
		File d = new File(directory);
		if (!d.exists()) d.mkdirs();

		int cnt = 0;
		for (File file : new File(directory).listFiles())
		{
			if (file.getName().endsWith(".txt"))
			{
				Scanner sc = new Scanner(file);
				while (sc.hasNextLine())
				{
					String line = sc.nextLine();
					if (!line.isEmpty()) vals.add(new Double(line));
				}
			}
			cnt++;

			if (cnt == randIter2) break;
		}
		return cnt;
	}

	private static void generateRandomPvals(MutexGreedySearcher searcher, Set<String> genes,
		Set<String> noShuffle, List<Double> vals, int howMany) throws IOException
	{
		String directory = dir + "randscores/";
		File d = new File(directory);
		if (!d.exists()) d.mkdirs();

		for (int i = 0; i < howMany; i++)
		{
			System.out.println("iteration = " + (i + 1));
			List<Double> list = searcher.generateRandPvals(genes, noShuffle, maxGroupSize, randIter1);
			if (vals != null) vals.addAll(list);

			Collections.sort(list);
			BufferedWriter writer = new BufferedWriter(new FileWriter(directory + "/randfile-" +
				System.currentTimeMillis() + "-" + new Random().nextInt(1000) + ".txt"));

			for (Double v : list)
			{
				writer.write(v + "\n");
			}

			writer.close();
		}
	}

	public static double decideBestFDR(final  Map<String,  Double>  results,
		List<Double>  randomized, int randMultiplier) throws IOException
	{
		double bestFDR = -1;
		double maxScore = -Double.MAX_VALUE;

		BufferedWriter writer = new BufferedWriter(new FileWriter(dir + "fdr-guide.txt"));

		writer.write("cutoff-val\tResult size\tFDR\tExpected true positives (tp)\tExpected false positives (fp)\ttp-fp");
		for (int i = 1; i <= 50; i++)
		{
			double fdr = i / 100D;
			List<String> select = FDR.select(results, fdr, randomized, randMultiplier);
			double tp = select.size() * (1 - fdr);
			double fp = select.size() * fdr;
			double score = tp - fp;

			double cutoffVal = -1;
			for (String s : select)
			{
				if (results.get(s) > cutoffVal) cutoffVal = results.get(s);
			}

			if (score >= maxScore)
			{
				maxScore = score;
				bestFDR = fdr;
			}

			writer.write("\n" + cutoffVal + "\t" + select.size() + "\t" + fdr + "\t" + tp + "\t" +
				fp + "\t" + (tp - fp));
		}
		writer.close();

		if (bestFDR == 0.5) return -1;

		return bestFDR;
	}

	private static Set<String> loadHighlySignificantGenes() throws FileNotFoundException
	{
		Set<String> set = new HashSet<String>();
		File f = new File(dir + "ranked-groups.txt");
		if (!f.exists()) return null;

		Scanner sc = new Scanner(f);
		boolean hasQval = sc.nextLine().split("\t").length > 2;

		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			double pval = Double.parseDouble(token[0]);
			if (pval >= 0.01) break;
			set.addAll(Arrays.asList(token).subList(hasQval ? 2 : 1, token.length));
		}
		sc.close();
		return set;
	}

	private static void downloadPortalData() throws IOException
	{
		System.out.println("Downloading alterations from cBioPortal for study " + portalStudyID);
		PortalDataset dataset = PortalDatasetEnum.findByStudyID(portalStudyID);
		if (dataset == null)
		{
			dataset = new PortalDataset("noname", minAltRatio, portalStudyID, portalCaseListID,
				new String[]{portalMutProfileID, portalCNAProfileID, portalExpProfileID},
				null, null);
		}
		else
		{
			if (portalCaseListID != null) dataset.caseList = portalCaseListID;
			if (portalMutProfileID != null) dataset.profile[0] = portalMutProfileID;
			if (portalCNAProfileID != null) dataset.profile[1] = portalCNAProfileID;
			if (portalExpProfileID != null) dataset.profile[2] = portalExpProfileID;
			if (minAltRatio != null) dataset.minAltThr = minAltRatio;
		}
		PortalReader reader = new PortalReader(dataset);
		if (symbolsFile != null) reader.setGeneSymbols(readSymbolsFile());
		reader.setUseNetwork(useGraph);
		reader.setOutputDir(dir);
		reader.setOutputFileName("DataMatrix.txt");
		reader.loadDataset();
		reader.updateParametersFile(new File(dir + "parameters.txt"));
	}

	private static boolean loadParameters() throws FileNotFoundException
	{
		File file = new File(dir + "parameters.txt");
		if (!file.exists() || file.isDirectory())
		{
			System.err.println("The file \"parameters.txt\" does not exist.");
			return false;
		}

		try{
		Scanner sc = new Scanner(file);
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			if (line.startsWith("#")) continue;
			String[] token = line.split("=");

			if (token.length < 2) continue;

			for (int i = 0; i < token.length; i++) token[i] = token[i].trim();

			if (token[0].equals("fdr-cutoff"))
			{
				fdrThr = Double.parseDouble(token[1]);
			}
			else if (token[0].equals("max-group-size"))
			{
				maxGroupSize = Integer.parseInt(token[1]);
			}
			else if (token[0].equals("first-level-random-iteration"))
			{
				randIter1 = Integer.parseInt(token[1]);
			}
			else if (token[0].equals("second-level-random-iteration"))
			{
				randIter2 = Integer.parseInt(token[1]);
			}
			else if (token[0].equals("data-file"))
			{
				dataFileName = dir + token[1];
			}
			else if (token[0].equals("search-on-signaling-network"))
			{
				useGraph = Boolean.parseBoolean(token[1]);
			}
			else if (token[0].equals("genes-file"))
			{
				symbolsFile = dir + token[1];
			}
			else if (token[0].equals("network-file"))
			{
				networkFilename = dir + token[1];
			}
			else if (token[0].equals("dataset-for-subtype"))
			{
				subtypeDatasetName = token[1];
			}
			else if (token[0].equals("keywords"))
			{
				literatureKeywords = token[1];
			}
			else if (token[0].equals("portal-study-id"))
			{
				portalStudyID = token[1];
			}
			else if (token[0].equals("portal-caselist-id"))
			{
				portalCaseListID = token[1];
			}
			else if (token[0].equals("portal-expression-profile-id"))
			{
				portalExpProfileID = token[1];
			}
			else if (token[0].equals("portal-cna-profile-id"))
			{
				portalCNAProfileID = token[1];
			}
			else if (token[0].equals("portal-mutation-profile-id"))
			{
				portalMutProfileID = token[1];
			}
			else if (token[0].equals("minimum-alteration-ratio"))
			{
				minAltRatio = new Double(token[1]);
			}
			else if (token[0].equals("randomize-data-matrix"))
			{
				randomizeDataMatrix = Boolean.parseBoolean(token[1]);
			}
			else if (token[0].equals("minimum-alteration-count-threshold"))
			{
				minAltCntThr = Integer.parseInt(token[1]);
			}
			else if (token[0].equals("gene-limit"))
			{
				geneLimit = Integer.parseInt(token[1]);
			}
			else if (token[0].equals("gene-ranking-file"))
			{
				geneRankingFile = dir + token[1];
			}
		}
		return true;
		} catch (Exception e)
		{
			System.err.println("Error while reading the parameters file.");
			e.printStackTrace();
			return false;
		}
	}
}
