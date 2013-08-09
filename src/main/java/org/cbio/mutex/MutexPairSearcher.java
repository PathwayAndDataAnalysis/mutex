package org.cbio.mutex;

import org.cbio.causality.analysis.Traverse;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.util.FDR;
import org.cbio.causality.util.Histogram;
import org.cbio.causality.util.Overlap;
import org.cbio.causality.util.Progress;

import java.util.*;

/**
 * Searcher class for mutex groups on the network.
 * @author Ozgun Babur
 */
public class MutexPairSearcher
{
	/**
	 * Gene alterations.
	 */
	private Map<String, GeneAlt> genes;

	/**
	 * Constructor with network and alterations.
	 * @param packs alterations
	 * @param alterationThreshold minimum ratio of altered samples, should be between 0 and 1.
	 */
	public MutexPairSearcher(Map<String, AlterationPack> packs,
		double alterationThreshold)
	{
		this.genes = new HashMap<String, GeneAlt>();

		System.out.println("unfiltered genes size = " + packs.size());

		for (String s : new HashSet<String>(packs.keySet()))
		{
			AlterationPack pack = packs.get(s);

			GeneAlt gene = new GeneAlt(pack, Alteration.GENOMIC);
//			gene.removeMinorCopyNumberAlts();

			// Filter out genes with les than 3% alteration
			double alteredRatio = gene.getAlteredRatio();
			if (alteredRatio > 0 && alteredRatio >= alterationThreshold)
				this.genes.put(gene.getId(), gene);
		}

		System.out.println("filtered to min alt = " + this.genes.size());
		System.out.println("Sample size = " + genes.values().iterator().next().size());
	}

	/**
	 * Performs a greedy search for each altered gene, assuming it is the common downstream. At each
	 * step the qualifying upstream of members are considered for expansion.
	 * @return significant mutex groups
	 */
	public List<Group> search(double fdrThr, int randMult)
	{
		Map<String, Double> seedScores = getGeneScores(genes.keySet());

		List<Double> randList = new ArrayList<Double>();
		for (int i = 0; i < randMult; i++)
		{
			shuffleAlterations();
			Map<String, Double> randScores = getGeneScores(genes.keySet());
			randList.addAll(randScores.values());
		}
		unshuffleAlterations();

		// Print distribution -------------------
		List<Double> scoresList = new ArrayList<Double>();
		for (String s : seedScores.keySet()) scoresList.add(seedScores.get(s));
		printPvalDistr(scoresList, randList, randMult, 0.05);
		// end of print -------------------------

		List<String> pairs = FDR.select(seedScores, fdrThr, randList, randMult);

		System.out.println("selected pairs size = " + pairs.size());

		for (String pair : pairs)
		{
			System.out.println(pair + "\t" + seedScores.get(pair));
		}

		return null;
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

	private Map<String, Double> getGeneScores(Set<String> names)
	{
		Map<String, Double> map = new HashMap<String, Double>();

		Progress prg = new Progress(names.size());

		for (String seed1 : names)
		{
			for (String seed2 : names)
			{
				if (seed1.compareTo(seed2) >= 0) continue;

				GeneAlt gene1 = genes.get(seed1);
				GeneAlt gene2 = genes.get(seed2);

				Group group = new Group(gene1);
				expandGroup(group, gene2);
//				double score = group.calcOverallPVal();
				double score = Overlap.calcCoocPval(gene1.getBooleanChanges(), gene2.getBooleanChanges());
				map.put(seed1 + "\t" + seed2, score);
			}

			prg.tick();
		}

		return map;
	}

	private void printPvalDistr(List<Double> result, List<Double> noise, int mult, double interval)
	{
		Histogram h1 = new Histogram(interval);
		h1.setBordered(true);
		for (Double pval : result)
		{
			h1.count(pval);
		}
		Histogram h2 = new Histogram(interval);
		h2.setBordered(true);
		for (Double pval : noise)
		{
			h2.count(pval);
		}

		h1.printTogether(h2, 1D / mult);
	}


	/**
	 * Expands the group with the best candidate in its candidates list. Also adds upstream of the
	 * added gene among candidates of the group to use in next expansion.
	 * @param group groups to expand
	 * @return true if expanded
	 */
	private void expandGroup(Group group, GeneAlt gene)
	{
		// Add the best gene
		group.addGene(gene, new HashSet<GeneAlt>(), true);
	}
}
