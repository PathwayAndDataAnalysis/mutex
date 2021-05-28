package org.cbio.mutex;

import org.panda.utility.ArrayUtil;
import org.panda.utility.Kronometre;
import org.panda.utility.statistics.FDR;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

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
	public static double fdrThr;

	/**
	 * Score threshold is used as an alternative to FDR threshold.
	 */
	public static double scoreThr;

	/**
	 * Maximum group size to use during the searches.
	 */
	public static int maxGroupSize;

	/**
	 * Maximum number of iterations for estimating the null distribution of initial p-values.
	 */
	public static int randIter1;

	/**
	 * Number of iterations to estimate the null distribution of final group scores.
	 */
	public static int randIter2;

	/**
	 * Directory that contains parameters file.
	 */
	public static String dir;

	/**
	 * Whether to reduce the search space using signaling networks.
	 */
	private static boolean useGraph;

	/**
	 * The name of the tab-delimited file containing gene alterations.
	 */
	public static String dataFileName;

	/**
	 * Name of the signaling network file. No need to specify this to use the default network.
	 */
	private static Set<String> networkFilename;

	/**
	 * The signaling network.
	 */
	private static Network network;

	/**
	 * Users can limit the search to certain genes using this file.
	 */
	private static String symbolsFile;

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
	 * A file to provide mapping between tissue types that exist in the dataset to the sample names.
	 */
	private static String sampleToTissueMappingFile;

	/**
	 * Tells even if there is need of random runs to estimate FDR, don't do it and just write result
	 * without FDR estimation.
	 */
	private static boolean noRandomRun = false;

	/**
	 * Parameter to run the analysis on a randomized set of alterations.
	 */
	private static boolean randomizeDataMatrix;

	public static void main(String[] args) throws IOException, ClassNotFoundException
	{
		Kronometre kron = new Kronometre();

		if (args.length < 1)
		{
			displayHelp();
			return;
		}

		reset();

		dir = args[0];
		if (!dir.endsWith(File.separator)) dir += File.separator;

		if (!loadParameters()) System.exit(1);

		network = null;
		if (useGraph)
		{
			if (networkFilename == null || networkFilename.isEmpty())
			{
				network = new Network();
			}
			else
			{
				for (String name : networkFilename)
				{
					if (network == null) network = new Network(name);
					else network.addResource(name);
				}
			}
		}

		if (args.length > 1 && args[1].equals("random"))
		{
			int howMany = 1;

			if (args.length > 2)
			{
				howMany = Integer.parseInt(args[2]);
				System.out.println("howMany = " + howMany);
			}

			generateRandomRun(howMany);
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
		scoreThr = -1;
		maxGroupSize = 5;
		randIter1 = 10000;
		randIter2 = 0;
		dir = null;
		useGraph = false;
		dataFileName = null;
		networkFilename = null;
		network = null;
		symbolsFile = null;
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
	public static void generateRandomRun(int howMany) throws IOException
	{
		// load the alteration data
		Map<String, GeneAlt> genesMap = loadAlterations();

		MutexGreedySearcher searcher = new MutexGreedySearcher(genesMap, network);
		Set<String> symbols = genesMap.keySet();
		if (network != null) symbols.retainAll(network.getSymbols());
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
		Map<String, GeneAlt> genesMap = loadAlterations();

		// if the data needs to be downloaded from cBioPortal, do it
		if (genesMap == null)
		{
			System.err.println("Cannot load alterations.");
			return;
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
		searcher.setTypeToInds(readTissueToSampleMapping());

		Set<String> symbols = genesMap.keySet();
		if (network != null) symbols.retainAll(network.getSymbols());

		Map<String, Group> groupsOfSeeds = searcher.getGroupsOfSeeds(symbols, maxGroupSize, randIter1);

		writeRankedGroups(groupsOfSeeds, null, "ranked-groups.txt");

		// we are done if we won't cutoff from an fdr or a score
		if (randIter2 <= 0 && scoreThr < 0) return;

		// Load and/or generate final scores null distribution
		List<Double> nullDist = new ArrayList<Double>();

		if (!noRandomRun && randIter2 > 0)
		{
			int cnt = readRandomPvals(nullDist);
			if (cnt < randIter2)
			{
				generateRandomPvals(searcher, symbols, loadHighlySignificantGenes(), nullDist,
					randIter2 - cnt);
			}
			writeRankedGroups(groupsOfSeeds, nullDist, "ranked-groups.txt");
		}

		// Apply FDR cutoff
		Map<String, Double> resultScores = new HashMap<String, Double>();
		for (String id : groupsOfSeeds.keySet())
		{
			resultScores.put(id, groupsOfSeeds.get(id).calcFinalScore());
		}

//		if (dir.contains("simulation1")) Simulation.plotEstimatedVsActualFDR(groupsOfSeeds, resultScores, nullDist, randIter2, 830, 150);
//		else if (dir.contains("simulation2")) Simulation.plotEstimatedVsActualFDR(groupsOfSeeds, resultScores, nullDist, randIter2, 156, 60);

		if (fdrThr < 0 && scoreThr < 0)
		{
			fdrThr = decideBestFDR(resultScores, nullDist, randIter2);
		}

		System.out.println("Selected FDR = " + fdrThr);
		if (fdrThr < 0 && scoreThr < 0) return;

		List<String> selectedSeeds = fdrThr >= 0 ? FDR.select(resultScores, fdrThr, nullDist, randIter2) :
			selectWithScore(resultScores, scoreThr);

		List<Group> groups = new ArrayList<>(selectedSeeds.size());
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

		// Write the output graph to visualize in ChiBE
		if (useGraph) GraphWriter.write(groups, genesMap, network, dir, dir);


		BufferedWriter writer = new BufferedWriter(new FileWriter(dir + "oncoprint.txt"));

		// Print textual results
		for (Group group : groups)
		{
			writer.write(group.getPrint(null, true, useGraph) + "\n\n");
		}
		writer.close();
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
		Set<String> genes = new HashSet<>();
		for (Group group : groups)
		{
			genes.addAll(group.getGeneNames());
		}
		List<String> list = new ArrayList<>(genes);

		list.sort((o1, o2) -> Integer.compare(genesMap.get(o2).getAltCnt(), genesMap.get(o1).getAltCnt()));

		return list;
	}

	private static List<String> selectWithScore(final Map<String, Double> resultScores, double thr)
	{
		List<String> list = new ArrayList<>();
		for (String s : resultScores.keySet())
		{
			if (resultScores.get(s) <= thr) list.add(s);
		}
		list.sort(Comparator.comparing(resultScores::get));
		return list;
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

	/**
	 * @deprecated
	 */
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

				for (String gene : new HashSet<>(map.keySet()))
				{
					if (!ranking.contains(gene)) map.remove(gene);
				}

				Iterator<String> iter = ranking.iterator();
				while (map.size() > geneLimit && iter.hasNext())
				{
					map.remove(iter.next());
				}
			}
			else
			{
				List<String> names = new ArrayList<>(map.keySet());
				final Map<String, Integer> cnt = new HashMap<>();
				for (String name : map.keySet())
				{
					cnt.put(name, map.get(name).countAltered());
				}
				Collections.sort(names, (o1, o2) -> cnt.get(o2).compareTo(cnt.get(o1)));
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

		if (!new File(geneRankingFile).exists())
		{
			throw new IllegalArgumentException("File does not exists: " + geneRankingFile);
		}

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

	private static Map<String, int[]> readTissueToSampleMapping() throws IOException
	{
		if (sampleToTissueMappingFile != null)
		{
			String[] header = Files.lines(Paths.get(dataFileName)).findFirst().get().split("\t");

			Map<String, String> sampleToType = Files.lines(Paths.get(sampleToTissueMappingFile)).skip(1)
				.map(l -> l.split("\t")).collect(Collectors.toMap(t -> t[0], t -> t[1], (s, s2) -> s));

			Map<String, List<Integer>> map = new HashMap<>();
			for (int i = 1; i < header.length; i++)
			{
				String type = sampleToType.get(header[i]);
				if (!map.containsKey(type)) map.put(type, new ArrayList<>());
				map.get(type).add(i - 1);
			}

			return map.keySet().stream()
				.collect(Collectors.toMap(t -> t, t -> ArrayUtil.convertToBasicIntArray(map.get(t))));
		}
		return null;
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
		Collections.sort(sorted, (o1, o2) -> scoreMap.get(o1).compareTo(scoreMap.get(o2)));

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
		if (!d.exists()) return 0;

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
			else if (token[0].equals("score-cutoff"))
			{
				scoreThr = Double.parseDouble(token[1]);
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
				if (networkFilename == null) networkFilename = new HashSet<>();
				networkFilename.add(dir + token[1]);
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
			else if (token[0].equals("sample-to-tissue-mapping-file"))
			{
				sampleToTissueMappingFile = dir + token[1];
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
