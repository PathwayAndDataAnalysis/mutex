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
	 * @param pvalThr p-value threshold for member significance
	 * @return significant mutex groups
	 */
	public List<Group> search(double pvalThr, double fdrThr, int maxGroupSize, int randMult)
	{
		List<Group> groups = new ArrayList<Group>();

		Map<String, Double> seedScores = getGeneScores(genes.keySet(), pvalThr, maxGroupSize);

		List<Double> randList = new ArrayList<Double>();
		for (int i = 0; i < randMult; i++)
		{
			shuffleAlterations();
			Map<String, Double> randScores = getGeneScores(genes.keySet(), pvalThr, maxGroupSize);
			randList.addAll(randScores.values());
		}
		unshuffleAlterations();

		List<String> seeds = FDR.select(seedScores, fdrThr, randList, randMult);

		System.out.println("selected seed size = " + seeds.size());

		for (String seed : seeds)
		{
			GeneAlt gene = genes.get(seed);

			Group group = new Group(gene);

			do
			{
				determineCandidates(group);
			}
			while(expandGroup(group) && group.size() < maxGroupSize);

			if (group.size() > 1)
			{
				group.shrinkToSignificantMembers(pvalThr);
			}

			if (group.size() > 1) groups.add(group);
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

	private Map<String, Double> getGeneScores(Set<String> names, double pvalThr,
		double maxGroupSize)
	{
		Map<String, Double> map = new HashMap<String, Double>();

		Progress prg = new Progress(names.size());

		for (String seed : names)
		{
			GeneAlt gene = genes.get(seed);

			Group group = new Group(gene);

			double bestScore = 1;
			boolean expanded;

//			do
//			{
//				determineCandidates(group);
//			}
//			while(expandGroup(group) && group.size() < maxGroupSize);
//
//			if (group.size() > 1)
//			{
//				group.shrinkToSignificantMembers(pvalThr);
//			}
//
//			map.put(seed, group.calcOverallPVal());

			do
			{
				determineCandidates(group);
				expanded = expandGroup(group);
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

//		printPvalDistr(map.values());

		return map;
	}

	private void printPvalDistr(Collection<Double> pvals)
	{
		Histogram h = new Histogram(0.05);
		for (Double pval : pvals)
		{
			h.count(pval);
		}
		h.print();
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
	private boolean expandGroup(Group group)
	{
		// Choose the best candidate

		GeneAlt best = null;
		double bestPval = 1;
		for (GeneAlt cand : group.candidates)
		{
			double pval = group.calcFuturePVal(cand, group.candidates.size());

			if (pval < bestPval)
			{
				bestPval = pval;
				best = cand;
			}
		}

		if (best != null)
		{
			// Add the best gene
			group.addGene(best, group.candidates.size(), true);

			return true;
		}
		else return false;
	}

	/**
	 * Gets the upstream of the given gene and adds the fitting one to the expansion candidates, and
	 * non-fitting ones to the black set.
	 * @param group group to update candidates
	 */
	private int determineCandidates(Group group)
	{
		List<String> members = group.getGeneNames();
		HashSet<String> candidates = new HashSet<String>(members);
		Set<String> comm = traverse.getLinkedCommonDownstream(candidates);
		candidates.addAll(comm);
		candidates.addAll(traverse.goBFS(candidates, null, false));
		candidates.removeAll(members);
		candidates.retainAll(genes.keySet());

		group.candidates.clear();

		for (String cand : candidates)
		{
			// the upstream gene is either a candidate or we don't want to re-consider it

			GeneAlt candGene = genes.get(cand);
			if (group.isOKToConsider(candGene))
			{
				group.candidates.add(candGene);
			}
			else
			{
				group.black.add(candGene);
			}
		}

		return group.candidates.size();
	}
}
