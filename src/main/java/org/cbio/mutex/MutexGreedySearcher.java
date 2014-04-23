package org.cbio.mutex;

import org.cbio.causality.analysis.Graph;
import org.cbio.causality.data.GeneCards;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;
import org.cbio.causality.util.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * Searcher class for mutex groups on the network.
 * @author Ozgun Babur
 */
public class MutexGreedySearcher
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

		AlterationPack sampleGene = packs.get("TP53");
		if (sampleGene == null) sampleGene = packs.get("PIK3CA");
		if (sampleGene == null) sampleGene = packs.get("BRAF");
//		boolean[] missing = getSamplesWithMissingData(sampleGene);
//		int mis = ArrayUtil.countValue(missing, true);
//		System.out.println("Samples missing cn or exp = " + mis);

		hyper = removeHypermutated ? AltDistr.getOutlierAltered(packs.values()) : null;

//		if (hyper != null) ArrayUtil.ORWith(hyper, missing);
//		else if (mis > 0) hyper = missing;

		System.out.println("unfiltered genes size = " + packs.size());
		System.out.println("alterationThreshold = " + alterationThreshold);

		for (String s : new HashSet<String>(packs.keySet()))
		{
			AlterationPack pack = packs.get(s);

			GeneAlt gene = new GeneAlt(pack, Alteration.GENOMIC, hyper);
//			gene.removeMinorCopyNumberAlts();

			// Filter out genes altered less than the alteration
			if (gene.getAlteredRatio() >= alterationThreshold)
				this.genes.put(gene.getId(), gene);
			else omitted.put(gene.getId(), gene);
		}

		System.out.println("filtered to min alt = " + this.genes.size());

		// Remove disconnected genes
		for (String s : new HashSet<String>(this.genes.keySet()))
		{
			Set<String> neigh = graph.getNeighbors(s);
			neigh.retainAll(this.genes.keySet());
			if (neigh.isEmpty())
			{
				omitted.put(s, this.genes.remove(s));
			}
		}

		// DEBUG CODE---------
//		Set<String> sss = new HashSet<String>(Arrays.asList("PIK3CA", "PIK3R1"));
//		for (String s : new HashSet<String>(this.genes.keySet()))
//		{
//			if (!sss.contains(s))
//			{
//				omitted.put(s, this.genes.remove(s));
//			}
//		}
		// DEBUG CODE---------

		System.out.println("filtered disconnected  = " + this.genes.size());

//		printCancerAssociationRatio(genes.keySet());

//		checkAGroup();
//		System.exit(0);
	}

	public boolean[] getHyper()
	{
		return hyper;
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
	 * Performs a greedy search for each altered gene, assuming it is the common downstream. At each
	 * step the qualifying upstream of members are considered for expansion.
	 * @return significant mutex groups
	 */
	public List<Group> search(double fdrThr, int maxGroupSize, int randMult1, int randMult2,
		String studyName)
	{
		System.out.println("\nCalculating seed scores");

		Set<String> seedCandidates = genes.keySet();
		Map<String, Double> seedScores = null;
//		Map<String, Double> seedScores = getGenePvals(seedCandidates, maxGroupSize, randMult1);
		List<Double> randScores = getRandPvals(seedCandidates, maxGroupSize, randMult1, randMult2,
			studyName);

		//-----------
//		List<Double> list = new ArrayList<Double>();
//		for (String key : seedScores.keySet())
//		{
//			list.add(seedScores.get(key));
//		}
//		DiscretePvalHisto h = new DiscretePvalHisto(list, 0.02);
//		h.plot();
		//-----------

//		Map<String, Double> limits = new HashMap<String, Double>();
//		for (String name : seedCandidates)
//		{
//			GeneAlt gene = genes.get(name);
//			limits.put(gene.getId(), gene.getMinPval());
//		}
//
//		List<String> seeds = FDR.select(seedScores, limits, fdrThr);
		List<String> seeds = FDR.select(seedScores, fdrThr, randScores, randMult2);
		double pvThr = findThrPval(seeds, seedScores);
		System.out.println("pvThr = " + pvThr);

		Map<String, Group> groupMap = getGroupsOfSeeds(seeds, maxGroupSize, pvThr, randMult1,
			genes);

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

		sanity(groups, randMult1);

		return groups;
	}

	private void sanity(List<Group> groups, int randMult1)
	{
		for (Group group : groups)
		{
			for (GeneAlt member : group.members)
			{
				if (member.randScores.length < randMult1)
				{
					System.out.println("Gene null distr has fever values: " + member.getId());
					System.out.println("size = " + member.randScores.length);
				}
			}
		}
	}

	private double findThrPval(List<String> seeds, Map<String, Double> seedScores)
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
			if (group != null) map.put(seed, group.calcPVal());
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
			List<Double> ll = new ArrayList<Double>(genes.size());

			for (GeneAlt gene : genes.values()) gene.shuffleSticky();

			Map<String, Double> map = new HashMap<String, Double>();

			final Progress prg = new Progress(names.size());
			for (String seed : names)
			{
				Group group = getGroupOfSeed(seed, maxGroupSize, -1, randIter1, genes);
				if (group != null) map.put(seed, group.calcPVal());
				prg.tick();
			}

			for (GeneAlt gene : genes.values()) gene.unshuffleSticky();

			for (String s : map.keySet())
			{
				ll.add(map.get(s));
			}

			writeRandomScores(ll, dir + System.currentTimeMillis() + ".txt");

			list.addAll(ll);
		}

		return list;
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

	private double calcGeneVal(GeneAlt gene, int maxGroupSize, boolean useScore,
		int randIter, Map<String, GeneAlt> genes)
	{
		Group group = new Group(gene);

		double bestVal = 1;
		boolean expanded;

		do
		{
			Set<GeneAlt> candidates = determineCandidates(group, genes);
			expanded = expandGroup(group, candidates, useScore, -1, maxGroupSize, randIter, genes);
			if (expanded)
			{
				double val = useScore ? group.calcScore() : group.calcPVal();

				if (val < bestVal) bestVal = val;
			}
		}
		while(expanded && group.size() < maxGroupSize);
		return bestVal;
	}

	private void assignNullScoreDistr(GeneAlt gene, int maxGroupSize, int randomIteration,
		double score, Map<String, GeneAlt> genes)
	{
		if (gene.randScores == null)
		{
			assignNullScores(gene, maxGroupSize,
				randomIteration > 100 ? randomIteration / 100 : randomIteration / 50, genes);
		}

		if (gene.randScores.length < randomIteration / 10 && gene.getPvalOfScore(score) < 0.5)
		{
			assignNullScores(gene, maxGroupSize, randomIteration / 10, genes);
		}

		if (gene.randScores.length < randomIteration && gene.getPvalOfScore(score) < 0.2)
		{
			assignNullScores(gene, maxGroupSize, randomIteration, genes);
		}

		gene.unshuffle();
	}

	private void assignNullScores(GeneAlt gene, int maxGroupSize, int randomIteration, Map<String, GeneAlt> genes)
	{
		double[] dist = getNullDist(gene, maxGroupSize, randomIteration, true, genes);
		Arrays.sort(dist);
		gene.setRandScores(dist);
	}

	private double[] getNullDist(GeneAlt gene, int maxGroupSize, int randomIteration,
		boolean score, Map<String, GeneAlt> genes)
	{
		double[] dist = new double[randomIteration];

		for (int i = 0; i < randomIteration; i++)
		{
			gene.shuffle();
			dist[i] = calcGeneVal(gene, maxGroupSize, score, randomIteration, genes);
		}
		return dist;
	}

	/**
	 * Expands the group with the best candidate in its candidates list. Also adds upstream of the
	 * added gene among candidates of the group to use in next expansion.
	 * @param group groups to expand
	 * @return true if expanded
	 */
	private boolean expandGroup(Group group, Set<GeneAlt> candidates, boolean useScore, double thr,
		int maxGroupSize, int randIter, Map<String, GeneAlt> genes)
	{
		if (candidates.isEmpty()) return false;

		// Choose the best candidate

		GeneAlt best = null;
		double bestVal = 1;

		double currentScore = group.calcScore();

		if (!useScore)
		{
			for (GeneAlt member : group.members)
			{
				if (member.randScores == null)
					assignNullScoreDistr(member, maxGroupSize, randIter, currentScore, genes);
			}
		}

		double currentVal = useScore ? currentScore : group.calcPVal();

		for (GeneAlt cand : candidates)
		{
			double futureScore = group.calcFutureScore(cand);

			if (!useScore)
			{
				assignNullScoreDistr(cand, maxGroupSize, randIter, futureScore, genes);

				for (GeneAlt member : group.members)
				{
					assignNullScoreDistr(member, maxGroupSize, randIter, futureScore, genes);
				}
			}

			double val = useScore ? futureScore : group.calcFuturePVal(cand);

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
}
