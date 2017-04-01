package org.cbio.mutex;

import org.panda.utility.Progress;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;

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
	 * Network provider.
	 */
	private DirectedGraph graph;

	/**
	 * Random number generator for shuffling operations.
	 */
	Random rand;

	/**
	 * A mapping from sample types to sample indices. Needed only if the dataset is heterogeneous, i.e. random
	 * alterations have unequal likelihood of distribution between types, independent from sample alteration rates.
	 */
	Map<String, int[]> typeToInds;

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
	 */
	public MutexGreedySearcher(Map<String, GeneAlt> geneAlts, DirectedGraph graph)
	{
		this.genes = geneAlts;
		this.graph = graph;
		this.rand = new Random();
	}

	public void setTypeToInds(Map<String, int[]> typeToInds)
	{
		this.typeToInds = typeToInds;
	}

	public Map<String, Group> getGroupsOfSeeds(Collection<String> seeds, int maxGroupSize,
		int randIter)
	{
		if (typeToInds != null) genes.values().forEach(g -> g.setTypeMap(typeToInds));

		Progress prg = new Progress(seeds.size(),
			"Searching for groups of " + seeds.size() + " seeds");

		Map<String, Group> s2g = new HashMap<>();
		for (String seed : seeds)
		{
			Group group = getGroupOfSeed(seed, maxGroupSize, randIter);
			if (group != null) s2g.put(seed, group);
			prg.tick();
		}

		return s2g;
	}

	private Group getGroupOfSeed(String seed, int maxGroupSize, int randIter)
	{
		GeneAlt gene = genes.get(seed);
		Group group = new Group(gene);

		Set<GeneAlt> candidates;

		do
		{
			candidates = determineCandidates(group);
			if (!expandGroup(group, candidates, true, maxGroupSize, randIter))
				break;
		}
		while(group.size() < maxGroupSize);


		if (group.size() > 1)
		{
			return group;
		}
		return null;
	}

	public void expandGroupIfPossible(Group group, double limitScore, int randIter)
	{
		Set<GeneAlt> candidates;

		do
		{
			candidates = determineCandidates(group);
			if (!expandGroupMore(group, candidates, Integer.MAX_VALUE, randIter, limitScore))
				break;
		}
		while(true);
	}

	public List<Double> generateRandPvals(Set<String> names, Set<String> noShuffle,
		int maxGroupSize, int randIter1)
	{
		Progress prog = new Progress(names.size(), "Generating a random run for final scores null distribution");
		List<Double> ll = new ArrayList<Double>(names.size());

		for (GeneAlt gene : genes.values())
		{
			if (noShuffle == null || !noShuffle.contains(gene.getId()))
			{
				gene.shuffleSticky();
			}
		}

		Map<String, Double> map = new HashMap<>();
		Map<String, Group> groups = new HashMap<>();

		for (String seed : names)
		{
			if (noShuffle == null || !noShuffle.contains(seed))
			{
				Group group = getGroupOfSeed(seed, maxGroupSize, randIter1);

				if (group != null) groups.put(seed, group);
			}
			prog.tick();
		}

		for (String seed : groups.keySet())
		{
			map.put(seed, groups.get(seed).calcFinalScore());
		}

		for (GeneAlt gene : genes.values()) gene.unshuffleSticky();

		for (String s : map.keySet())
		{
			ll.add(map.get(s));
		}
		return ll;
	}

	private double calcGeneVal(GeneAlt gene, int maxGroupSize, int randIter)
	{
		Group group = new Group(gene);

		boolean expanded;

		do
		{
			Set<GeneAlt> candidates = determineCandidates(group);
			expanded = expandGroup(group, candidates, false, maxGroupSize, randIter);
		}
		while(expanded && group.size() < maxGroupSize);

		return group.calcPVals1().get(gene.getId());
	}

	private void assignNullScoreDistr(GeneAlt gene, int maxGroupSize, int randomIteration,
		double score)
	{
		if (gene.randScores != null &&
			(gene.randScores.size() == randomIteration || // already met highest
				(gene.randScores.size() >= HIGH_ACCURACY && gene.randScores.get(HIGH_ACCURACY - 1) <= score) || // already high accuracy
				((gene.randScores.get((int) (gene.randScores.size() * ACCURACY_SWITCH)) < score) && // only need low accuracy ..
					(gene.randScores.size() >= LOW_ACCURACY && gene.randScores.get(LOW_ACCURACY - 1) <= score)))) // but already low accuracy
			return;

		List<Double> dist = getNullDist(gene, maxGroupSize, randomIteration, gene.randScores, score);
		Collections.sort(dist);
		gene.setRandScores(dist);
		gene.unshuffle();
	}

	private List<Double> getNullDist(GeneAlt gene, int maxGroupSize, int randomIteration,
		List<Double> startWith, double forValue)
	{
		int cnt = 0;
		if (startWith == null) startWith = new ArrayList<>();
		else cnt = countLessThanOrEqual(startWith, forValue);

		double p = cnt / (double) startWith.size();

		while ((cnt < LOW_ACCURACY || (cnt < HIGH_ACCURACY && p < ACCURACY_SWITCH)) &&
			 startWith.size() < randomIteration)
		{
			gene.shuffle(rand);
			double val = calcGeneVal(gene, maxGroupSize, randomIteration);
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
	 * Expands the group with the best candidate in its candidates list.
	 * @param group groups to expand
	 * @return true if expanded
	 */
	private boolean expandGroup(Group group, Set<GeneAlt> candidates, boolean useFinalScore,
		int maxGroupSize, int randIter)
	{
		if (candidates.isEmpty()) return false;

		// Choose the best candidate

		GeneAlt best = null;
		double bestVal = 1;

		if (useFinalScore)
		{
			Map<String, Double> pv = group.calcPVals1();
			for (GeneAlt member : group.members)
			{
				if (member.randScores == null)
					assignNullScoreDistr(member, maxGroupSize, randIter, pv.get(member.id));
			}
		}

		double currentVal = useFinalScore ? group.calcFinalScore() : group.calcScore();

		for (GeneAlt cand : candidates)
		{
			if (useFinalScore)
			{
				Map<String, Double> pv = group.calcFuturePvals1(cand);

				assignNullScoreDistr(cand, maxGroupSize, randIter, pv.get(cand.id));

				for (GeneAlt member : group.members)
				{
					assignNullScoreDistr(member, maxGroupSize, randIter, pv.get(member.id));
				}
			}

			double future = useFinalScore ? group.calcFutureFinalScore(cand) :
				group.calcFutureScore(cand);

			if (future < bestVal && future < currentVal)
			{
				bestVal = future;
				best = cand;
			}
		}

		if (best != null)
		{
			// Add the best gene
			group.addGene(best);

			return true;
		}
		else if (useFinalScore && group.members.size() == 1)
		{
			return expandGroup(group, candidates, false, maxGroupSize, randIter);
		}
		else return false;
	}

	/**
	 * Expands the group with the best candidate in its candidates list.
	 * @param group groups to expand
	 * @return true if expanded
	 */
	private boolean expandGroupMore(Group group, Set<GeneAlt> candidates, int maxGroupSize,
		int randIter, double scoreThr)
	{
		if (candidates.isEmpty()) return false;

		// Choose the best candidate

		GeneAlt best = null;
		double bestVal = 1;

		Map<String, Double> pv = group.calcPVals1();
		for (GeneAlt member : group.members)
		{
			if (member.randScores == null)
				assignNullScoreDistr(member, maxGroupSize, randIter, pv.get(member.id));
		}

		for (GeneAlt cand : candidates)
		{
			pv = group.calcFuturePvals1(cand);

			assignNullScoreDistr(cand, maxGroupSize, randIter, pv.get(cand.id));

			for (GeneAlt member : group.members)
			{
				assignNullScoreDistr(member, maxGroupSize, randIter, pv.get(member.id));
			}

			double future = group.calcFutureFinalScore(cand);

			if (future < bestVal && future < scoreThr)
			{
				bestVal = future;
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
	private Set<GeneAlt> determineCandidates(Group group)
	{
		Set<String> candNames = getCandidateNames(group);

		Set<GeneAlt> candidates = new HashSet<>();

		for (String cand : candNames)
		{
			// the upstream gene is either a candidate or we don't want to re-consider it

			GeneAlt candGene = genes.get(cand);
			if (group.isOKToConsider(candGene))
			{
				candidates.add(candGene);
			}
		}

		return candidates;
	}

	private Set<String> getCandidateNames(Group group)
	{
		if (graph == null) return getCandidateNamesWithoutUsingGraph(group);

		List<String> members = group.getGeneNames();
		HashSet<String> candNames = new HashSet<String>(members);
		Set<String> comm = graph.getLinkedCommonDownstream(candNames);
		candNames.addAll(comm);
		candNames.addAll(graph.getUpstream(candNames));
		candNames.removeAll(members);
		candNames.retainAll(genes.keySet());
		return candNames;
	}

	private Set<String> getCandidateNamesWithoutUsingGraph(Group group)
	{
		List<String> members = group.getGeneNames();
		Set<String> comm = new HashSet<String>(genes.keySet());
		comm.removeAll(members);
		return comm;
	}
}
