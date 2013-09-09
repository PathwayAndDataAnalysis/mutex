package org.cbio.mutex;

import org.cbio.causality.analysis.Traverse;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.util.FDR;
import org.cbio.causality.util.Histogram;
import org.cbio.causality.util.Progress;

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
	 * Network provider.
	 */
	private Traverse traverse;

	/**
	 * Constructor with network and alterations.
	 * @param traverse the network helper
	 * @param packs alterations
	 * @param alterationThreshold minimum ratio of altered samples, should be between 0 and 1.
	 */
	public MutexGreedySearcher(Traverse traverse, Map<String, AlterationPack> packs,
		double alterationThreshold)
	{
		this.traverse = traverse;
		this.genes = new HashMap<String, GeneAlt>();

		System.out.println("unfiltered genes size = " + packs.size());

		for (String s : new HashSet<String>(packs.keySet()))
		{
			AlterationPack pack = packs.get(s);

			GeneAlt gene = new GeneAlt(pack, Alteration.GENOMIC);
			gene.removeMinorCopyNumberAlts();

			// Filter out genes with les than 3% alteration
			if (gene.getAlteredRatio() >= alterationThreshold)
				this.genes.put(gene.getId(), gene);
		}

		System.out.println("filtered to min alt = " + this.genes.size());

		// Remove disconnected genes
		for (String s : new HashSet<String>(this.genes.keySet()))
		{
			Set<String> neigh = traverse.getNeighbors(s);
			neigh.retainAll(this.genes.keySet());
			if (neigh.isEmpty()) this.genes.remove(s);
		}

		System.out.println("filtered disconnected  = " + this.genes.size());

//		for (GeneAlt gene : genes.values())
//		{
//			gene.shuffle();
//		}
	}

	/**
	 * Performs a greedy search for each altered gene, assuming it is the common downstream. At each
	 * step the qualifying upstream of members are considered for expansion.
	 * @return significant mutex groups
	 */
	public List<Group> search(double fdrThr, int maxGroupSize, int randMult)
	{
		List<Group> groups = new ArrayList<Group>();

		Map<String, Double> seedScores = getGeneScores(genes.keySet(), maxGroupSize);

		List<Double> randList = new ArrayList<Double>();
		for (int i = 0; i < randMult; i++)
		{
			System.out.println("\nrandomization " + i);
			shuffleAlterations();
			Map<String, Double> randScores = getGeneScores(genes.keySet(), maxGroupSize);
			randList.addAll(randScores.values());
		}
		unshuffleAlterations();

		// Print distribution -------------------
//		List<Double> scoresList = new ArrayList<Double>();
//		for (String s : seedScores.keySet()) scoresList.add(seedScores.get(s));
//		printPvalDistr(scoresList, randList, randMult, 0.05);
		// end of print -------------------------

		List<String> seeds = FDR.select(seedScores, fdrThr, randList, randMult);

		System.out.println("selected seed size = " + seeds.size());

		// estimate the pval threshold

		double pvalThr = 0;
		for (String seed : seeds)
		{
			if (seedScores.get(seed) > pvalThr) pvalThr = seedScores.get(seed);
		}

		System.out.println("pvalThr = " + pvalThr);

		for (String seed : seeds)
		{
			GeneAlt gene = genes.get(seed);

			Group group = new Group(gene);

			Set<GeneAlt> candidates;
			do
			{
				candidates = determineCandidates(group);
			}
			while(expandGroup(group, candidates) && group.size() < maxGroupSize);

			assert group.size() > 1;

			group.shrinkToSignificantMembers(pvalThr);

			assert group.size() > 1;

			groups.add(group);
		}

		return groups;
	}

	private void shuffleAlterations()
	{
		for (GeneAlt gene : genes.values())
		{
			gene.shuffle();
		}
	}

	private void unshuffleAlterations()
	{
		for (GeneAlt gene : genes.values())
		{
			gene.unshuffle();
		}
	}

	private Map<String, Double> getGeneScores(Set<String> names, double maxGroupSize)
	{
		Map<String, Double> map = new HashMap<String, Double>();

		Progress prg = new Progress(names.size());

		for (String seed : names)
		{
			GeneAlt gene = genes.get(seed);

			Group group = new Group(gene);

			double bestScore = 1;
			boolean expanded;

			do
			{
				Set<GeneAlt> candidates = determineCandidates(group);
				expanded = expandGroup(group, candidates);
				if (expanded)
				{
					double score = group.calcOverallPVal();
					if (score < bestScore) bestScore = score;
				}
			}
			while(expanded && group.size() < maxGroupSize);

			map.put(seed, bestScore);
			prg.tick();
		}

		return map;
	}

	private void printPvalDistr(List<Double> result, List<Double> noise, int mult, double interval)
	{
		Histogram h1 = new Histogram(interval);
		for (Double pval : result)
		{
			h1.count(pval);
		}
		Histogram h2 = new Histogram(interval);
		for (Double pval : noise)
		{
			h2.count(pval);
		}

		h1.printTogether(h2, 1D / mult);
	}

	private List<String> selectToFDR(final Map<String, Double> pvMap, double thr)
	{
		List<String> names = new ArrayList<String>(pvMap.keySet());
		Collections.sort(names, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return pvMap.get(o1).compareTo(pvMap.get(o2));
			}
		});

		List<String> filtered = new ArrayList<String>();

		for (String name : names)
		{
			double pv = pvMap.get(name);

			double noise = pv * pvMap.size();

			if (noise / (filtered.size() + 1) <= thr)
			{
				filtered.add(name);
			}
			else break;
		}

		return filtered;
	}

	/**
	 * Expands the group with the best candidate in its candidates list. Also adds upstream of the
	 * added gene among candidates of the group to use in next expansion.
	 * @param group groups to expand
	 * @return true if expanded
	 */
	private boolean expandGroup(Group group, Set<GeneAlt> candidates)
	{
		if (candidates.isEmpty()) return false;

		// Choose the best candidate

		GeneAlt best = null;
		double bestPval = 1;
		double currentPval = group.calcOverallPVal();

		for (GeneAlt cand : candidates)
		{
			double pval = group.calcFuturePVal(cand, new HashSet<GeneAlt>(candidates));

			if (pval < bestPval)
			{
				bestPval = pval;
				best = cand;
			}
		}

		if (best != null && bestPval < currentPval)
		{
			// Add the best gene
			group.addGene(best, candidates, true);

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
		List<String> members = group.getGeneNames();
		HashSet<String> candNames = new HashSet<String>(members);
		Set<String> comm = traverse.getLinkedCommonDownstream(candNames);
		candNames.addAll(comm);
		candNames.addAll(traverse.goBFS(candNames, null, false));
		candNames.removeAll(members);
		candNames.retainAll(genes.keySet());

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
}
