package org.cbio.mutex;

import org.cbio.causality.analysis.Graph;
import org.cbio.causality.data.GeneCards;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;
import org.cbio.causality.util.*;

import java.io.*;
import java.util.*;

/**
 * Searcher class for mutex groups on the network.
 * @author Ozgun Babur
 */
public class MutexGreedySearcher implements Serializable
{
	/**
	 * Gene alterations.
	 */
	private Map<String, GeneAlt> genes;

	/**
	 * Gene alterations with low alteration rate.
	 */
	private Map<String, GeneAlt> omitted;

	/**
	 * Network provider.
	 */
	private Graph graph;

	private boolean[] hyper;

	private static final long serialVersionUID = 5998727214775118493L;

	/**
	 * A null distribution is sampled until the number of values smaller than the current compared
	 * value is equal to this number, or up to iteration limit.
	 */
	private static final int LOW_ACCURACY = 10;
	private static final int HIGH_ACCURACY = 100;
	private static final double ACCURACY_SWITCH = 0.2;

	/**
	 * Constructor with network and alterations.
	 * @param graph the network helper
	 * @param packs alterations
	 * @param alterationThreshold minimum ratio of altered samples, should be between 0 and 1.
	 */
	public MutexGreedySearcher(Graph graph, Map<String, AlterationPack> packs,
		double alterationThreshold, boolean removeHypermutated)
	{
		this.graph = graph;
		this.genes = new HashMap<String, GeneAlt>();
		this.omitted = new HashMap<String, GeneAlt>();

//		AlterationPack sampleGene = packs.get("TP53");
//		if (sampleGene == null) sampleGene = packs.get("PIK3CA");
//		if (sampleGene == null) sampleGene = packs.get("BRAF");
//		boolean[] missing = getSamplesWithMissingData(sampleGene);
//		int mis = ArrayUtil.countValue(missing, true);
//		System.out.println("Samples missing cn or exp = " + mis);

		hyper = removeHypermutated ? AltDistr.getOutlierAltered(packs.values()) : null;

//		if (hyper != null) ArrayUtil.ORWith(hyper, missing);
//		else if (mis > 0) hyper = missing;

		System.out.println("unfiltered genes size = " + packs.size());
		System.out.println("alterationThreshold = " + alterationThreshold);
		System.out.println("original sample size = " + packs.values().iterator().next().getSize());

		for (String s : new HashSet<String>(packs.keySet()))
		{
			AlterationPack pack = packs.get(s);

			GeneAlt gene = new GeneAlt(pack, Alteration.GENOMIC, hyper);

			// Filter out genes altered less than the alteration
			if (gene.getAlteredRatio() >= alterationThreshold)
				this.genes.put(gene.getId(), gene);
			else omitted.put(gene.getId(), gene);
		}

		System.out.println("filtered to min alt = " + this.genes.size());

		// Remove disconnected genes
		for (String s : new HashSet<String>(this.genes.keySet()))
		{
			Set<String> related = graph.getGenesWithCommonDownstream(s);
			related.retainAll(this.genes.keySet());
			if (related.isEmpty())
			{
				omitted.put(s, this.genes.remove(s));
			}
		}

		System.out.println("filtered disconnected  = " + this.genes.size());

//		printCancerAssociationRatio(genes.keySet());

//		checkAGroup();
//		System.exit(0);
	}

	public boolean[] getHyper()
	{
		return hyper;
	}

	public Graph getGraph()
	{
		return graph;
	}

	private static void printCancerAssociationRatio(Set<String> genes)
	{
		int related = 0;
		for (String gene : genes)
		{
			Set<String> set = GeneCards.getRelatedCancers(gene);
			if (!set.isEmpty()) related++;
//			System.out.println(gene + "\t" + set);
		}
		System.out.println("related = " + related);
		System.out.println("genes.size() = " + genes.size());
		System.out.println("cancer related ratio = " + related / (double) genes.size());
		System.exit(0);
	}

	/**
	 * Adds the omitted genes back to the genes map.
	 */
	public void addOmitted()
	{
		genes.putAll(omitted);
	}

	public Map<String, GeneAlt> getGenes()
	{
		return genes;
	}

	private boolean[] getSamplesWithMissingData(AlterationPack pack)
	{
		boolean[] b = new boolean[pack.getSize()];

		for (int i = 0; i < b.length; i++)
		{
			if (pack.getChange(Alteration.COPY_NUMBER, i).equals(Change.NO_DATA))
			{
				b[i] = true;
			}
		}
		return b;
	}

	/**
	 * Writes ordered mutex groups in a file.
	 */
	public void writeRankedGroups(int maxGroupSize, int randMult1, Map<String, String> labelMap, String filename) throws IOException
	{
		System.out.println("\nCalculating seed scores");

		Set<String> seedCandidates = genes.keySet();
		final Map<String, Double> seedScores = getGenePvals(seedCandidates, maxGroupSize, randMult1);

		List<String> seeds = new ArrayList<String>(seedScores.keySet());
		Collections.sort(seeds, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return seedScores.get(o1).compareTo(seedScores.get(o2));
			}
		});

		Map<String, Group> groupMap = getGroupsOfSeeds(seeds, maxGroupSize, -1, randMult1, genes);

		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

		for (String seed : seeds)
		{
			Group group = groupMap.get(seed);

			if (group != null)
			{
				String s = "";

				for (String name : group.getGeneNames())
				{
					s += (labelMap == null ? name : labelMap.get(name)) + "\t";
				}
				writer.write(s.trim() + "\n");
			}
		}

		writer.close();
	}

	/**
	 * Performs a greedy search for each altered gene, assuming it is the common downstream. At each
	 * step the qualifying upstream of members are considered for expansion.
	 * @return significant mutex groups
	 */
	public List<Group> search(double fdrThr, int maxGroupSize, int randMult1, int randMult2,
		String studyName, String fileCache) throws IOException
	{
		System.out.println("\nCalculating seed scores");

		Set<String> seedCandidates = genes.keySet();
//		Map<String, Double> seedScores = null;
		Map<String, Double> seedScores = getGenePvals(seedCandidates, maxGroupSize, randMult1);
		List<Double> randScores = getRandPvals(seedCandidates, maxGroupSize, randMult1, randMult2,
			studyName);

		if (fileCache != null) serialize(fileCache);

		List<String> seeds = FDR.select(seedScores, fdrThr, randScores, randMult2);
		double fsThr = findThrScore(seeds, seedScores);
		System.out.println("Thr final score = " + fsThr);

		Map<String, Group> groupMap = getGroupsOfSeeds(seeds, maxGroupSize, fsThr, randMult1, genes);

		System.out.println("selected seed size = " + seeds.size());
		System.out.println("selected seeds = " + seeds);

		List<Group> groups = new ArrayList<Group>();

		Set<String> totalGenes = new HashSet<String>();

		for (String seed : seeds)
		{
			Group group = groupMap.get(seed);
			if (group == null)
			{
				System.out.println("Group for seed " + seed + " is null!");
				continue;
			}
			group.fetchTragets(graph, genes);
			groups.add(group);
			totalGenes.addAll(group.getGeneNames());
		}

		System.out.println("genes in groups = " + totalGenes.size());

		// clean subsets in the result
		Group.removeSubsets(groups);

		return groups;
	}

	private double findThrScore(List<String> seeds, Map<String, Double> seedScores)
	{
		double max = 0;
		for (String seed : seeds)
		{
			if (seedScores.get(seed) > max) max = seedScores.get(seed);
		}
		return max;
	}

	private Map<String, Group> getGroupsOfSeeds(Collection<String> seeds, int maxGroupSize,
		double thr, int randIter, Map<String, GeneAlt> genes)
	{
		Map<String, Group> s2g = new HashMap<String, Group>();
		for (String seed : seeds)
		{
			Group group = getGroupOfSeed(seed, maxGroupSize, thr, randIter, genes);
			if (group == null)
			{
				System.out.println("group for seed " + seed + " is null!");
				group = getGroupOfSeed(seed, maxGroupSize, thr, randIter, genes);
			}
			if (group != null) s2g.put(seed, group);
		}

		return s2g;
	}

	private Group getGroupOfSeed(String seed, int maxGroupSize, double thr,
		int randIter, Map<String, GeneAlt> genes)
	{
		GeneAlt gene = genes.get(seed);
		Group group = new Group(gene);

		Set<GeneAlt> candidates;

		do
		{
			candidates = determineCandidates(group, genes);
			if (!expandGroup(group, candidates, false, thr, maxGroupSize, randIter, genes))
				break;
		}
		while(group.size() < maxGroupSize);


		if (group.size() > 1)
		{
			return group;
		}
		return null;
	}

	private Map<String, Double> getGenePvals(Set<String> names, int maxGroupSize,
		int randIteration)
	{
		Map<String, Double> map = new HashMap<String, Double>();

		Progress prg = new Progress(names.size());
		for (String seed : names)
		{
			Group group = getGroupOfSeed(seed, maxGroupSize, -1, randIteration, genes);
			if (group != null) map.put(seed, group.calcFinalScore());
			prg.tick();
		}

		return map;
	}

	private List<Double> getRandPvals(final Set<String> names, final int maxGroupSize,
		final int randIter1, int randIter2, String studyName)
	{
		final List<Double> list = new ArrayList<Double>(genes.size() * randIter2);

		String dir = "data/randscores/" + studyName + "/" + randIter1 + "/";
		File file = new File(dir);
		if (!file.exists()) file.mkdirs();

		int read = readRandomPvals(dir, list, randIter2);

		for (int i = read; i < randIter2; i++)
		{
			List<Double> ll = getRandPvals(names, maxGroupSize, randIter1);

			writeRandomScores(ll, dir + System.currentTimeMillis() + ".txt");

			list.addAll(ll);
		}

		return list;
	}

	private List<Double> getRandPvals(Set<String> names, int maxGroupSize, int randIter1)
	{
		List<Double> ll = new ArrayList<Double>(genes.size());

		for (GeneAlt gene : genes.values()) gene.shuffleSticky();

		Map<String, Double> map = new HashMap<String, Double>();

		final Progress prg = new Progress(names.size());
		for (String seed : names)
		{
			Group group = getGroupOfSeed(seed, maxGroupSize, -1, randIter1, genes);
			if (group != null) map.put(seed, group.calcFinalScore());
			prg.tick();
		}

		for (GeneAlt gene : genes.values()) gene.unshuffleSticky();

		for (String s : map.keySet())
		{
			ll.add(map.get(s));
		}
		return ll;
	}

	private void writeRandomScores(List<Double> list, String filename)
	{
		try
		{
			BufferedWriter writer = new BufferedWriter(new FileWriter(filename));
			for (Double val : list)
			{
				writer.write(val + "\n");
			}
			writer.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	private int readRandomPvals(String dir, List<Double> list, int count)
	{
		int cnt = 0;
		try
		{
			for (File file : new File(dir).listFiles())
			{
				if (file.getName().endsWith(".txt"))
				{
					Scanner sc = new Scanner(file);
					while (sc.hasNextLine())
					{
						String line = sc.nextLine();
						if (!line.isEmpty()) list.add(new Double(line));
					}
				}
				cnt++;

				if (cnt == count) break;
			}
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		return cnt;
	}

//	private double calcGeneVal(GeneAlt gene, int maxGroupSize, int randIter, Map<String, GeneAlt> genes)
//	{
//		Group group = new Group(gene);
//
//		double bestVal = 1;
//		boolean expanded;
//
//		do
//		{
//			Set<GeneAlt> candidates = determineCandidates(group, genes);
//			expanded = expandGroup(group, candidates, true, -1, maxGroupSize, randIter, genes);
//			if (expanded)
//			{
//				double val = group.calcScore();
//
//				if (val < bestVal) bestVal = val;
//			}
//		}
//		while(expanded && group.size() < maxGroupSize);
//		return bestVal;
//	}

	private double calcGeneVal(GeneAlt gene, int maxGroupSize, int randIter, Map<String, GeneAlt> genes)
	{
		Group group = new Group(gene);

		boolean expanded;

		do
		{
			Set<GeneAlt> candidates = determineCandidates(group, genes);
			expanded = expandGroup(group, candidates, true, -1, maxGroupSize, randIter, genes);
		}
		while(expanded && group.size() < maxGroupSize);

		return group.calcPVals1().get(gene.getId());
	}

	private void assignNullScoreDistr(GeneAlt gene, int maxGroupSize, int randomIteration,
		double score, Map<String, GeneAlt> genes)
	{
		if (gene.randScores != null &&
			(gene.randScores.size() == randomIteration || // already met highest
				(gene.randScores.size() >= HIGH_ACCURACY && gene.randScores.get(HIGH_ACCURACY - 1) <= score) || // already high accuracy
				((gene.randScores.get((int) (gene.randScores.size() * ACCURACY_SWITCH)) < score) && // only need low accuracy ..
					(gene.randScores.size() >= LOW_ACCURACY && gene.randScores.get(LOW_ACCURACY - 1) <= score)))) // but already low accuracy
			return;

		List<Double> dist = getNullDist(gene, maxGroupSize, randomIteration, genes, gene.randScores, score);
		Collections.sort(dist);
		gene.setRandScores(dist);
		gene.unshuffle();
	}

	private List<Double> getNullDist(GeneAlt gene, int maxGroupSize, int randomIteration,
		Map<String, GeneAlt> genes, List<Double> startWith, double forValue)
	{
		int cnt = 0;
		if (startWith == null) startWith = new ArrayList<Double>();
		else cnt = countLessThanOrEqual(startWith, forValue);

		double p = cnt / (double) startWith.size();

		while ((cnt < LOW_ACCURACY || (cnt < HIGH_ACCURACY && p < ACCURACY_SWITCH)) &&
			 startWith.size() < randomIteration)
		{
			gene.shuffle();
			double val = calcGeneVal(gene, maxGroupSize, randomIteration, genes);
			startWith.add(val);
			if (val <= forValue)
			{
				cnt++;
				p = cnt / (double) startWith.size();
			}
		}
		return startWith;
	}

	private int countLessThanOrEqual(List<Double> randScore, double val)
	{
		int cnt  = 0;
		for (Double v : randScore)
		{
			if (v <= val) cnt++;
			else break;
		}
		return cnt;
	}

	/**
	 * Expands the group with the best candidate in its candidates list. Also adds upstream of the
	 * added gene among candidates of the group to use in next expansion.
	 * @param group groups to expand
	 * @return true if expanded
	 */
	private boolean expandGroup(Group group, Set<GeneAlt> candidates, boolean useInitialScore,
		double thr, int maxGroupSize, int randIter, Map<String, GeneAlt> genes)
	{
		assert thr < 0 || !useInitialScore;

		if (candidates.isEmpty()) return false;

		// Choose the best candidate

		GeneAlt best = null;
		double bestVal = 1;

		double currentScore = group.calcScore();

		if (!useInitialScore)
		{
			for (GeneAlt member : group.members)
			{
				if (member.randScores == null)
					assignNullScoreDistr(member, maxGroupSize, randIter, currentScore, genes);
			}
		}

		double currentVal = useInitialScore ? currentScore : group.calcFinalScore();

		for (GeneAlt cand : candidates)
		{
			if (!useInitialScore)
			{
				double futScore = group.calcFutureScore(cand);
				assignNullScoreDistr(cand, maxGroupSize, randIter, futScore, genes);

				for (GeneAlt member : group.members)
				{
					assignNullScoreDistr(member, maxGroupSize, randIter, futScore, genes);
				}
			}

			double val = useInitialScore ?
				group.calcFutureScore(cand) : group.calcFutureFinalScore(cand);

			if (val < bestVal && (val < currentVal || val <= thr))
			{
				bestVal = val;
				best = cand;
			}
		}

		if (best != null)
		{
			// Add the best gene
			group.addGene(best);

			return true;
		}
		else if (!useInitialScore && currentVal == 1)
		{
			return expandGroup(group, candidates, true, -1, maxGroupSize, randIter, genes);
		}
		else return false;
	}

	/**
	 * Gets the upstream of the given gene and adds the fitting one to the expansion candidates, and
	 * non-fitting ones to the black set.
	 * @param group group to update candidates
	 */
	private Set<GeneAlt> determineCandidates(Group group, Map<String, GeneAlt> genes)
	{
		Set<String> candNames = getCandidateNames(group);

		Set<GeneAlt> candidates = new HashSet<GeneAlt>();

		for (String cand : candNames)
		{
			// the upstream gene is either a candidate or we don't want to re-consider it

			GeneAlt candGene = genes.get(cand);
			if (group.isOKToConsider(candGene))
			{
				candidates.add(candGene);
			}
			else
			{
				group.black.add(candGene);
			}
		}

		return candidates;
	}

	private Set<String> getCandidateNames(Group group)
	{
		List<String> members = group.getGeneNames();
		HashSet<String> candNames = new HashSet<String>(members);
		Set<String> comm = graph.getLinkedCommonDownstream(candNames);
		candNames.addAll(comm);
		candNames.addAll(graph.getUpstream(candNames));
		candNames.removeAll(members);
		candNames.retainAll(genes.keySet());
		return candNames;
	}

	private Set<String> getCandidateNamesX(Group group)
	{
		List<String> members = group.getGeneNames();
		Set<String> comm = new HashSet<String>(genes.keySet());
		comm.removeAll(members);
		return comm;
	}

	public void serialize(String filename) throws IOException
	{
		ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(filename));
		out.writeObject(this);
		out.close();
		System.out.println("Wrote to file " + filename);
	}

	public static MutexGreedySearcher deserialize(String filename) throws IOException, ClassNotFoundException
	{
		ObjectInputStream in = new ObjectInputStream(new FileInputStream(filename));
		MutexGreedySearcher s = (MutexGreedySearcher) in.readObject();
		in.close();
		return s;
	}

	public static void main(String[] args) throws IOException, ClassNotFoundException
	{
		String study = args[0];
		MutexGreedySearcher searcher = deserialize("cache-" + study);

		String dir = "scores-" + study;
		File f = new File(dir);
		if (!f.exists()) f.mkdirs();

		List<Double> randPvals = searcher.getRandPvals(searcher.getGenes().keySet(), 10, 10000);
		searcher.writeRandomScores(randPvals, dir + "/randfile-" + System.currentTimeMillis() + "-" + new Random().nextInt(1000) + ".txt");
	}
}
