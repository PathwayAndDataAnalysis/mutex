package org.cbio.mutex;

import org.cbio.causality.analysis.Traverse;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;

import java.util.*;

/**
 * Searcher class for mutex groups on the network.
 * @author Ozgun Babur
 */
public class MutexExhaustiveSearcher
{
	/**
	 * Gene alterations.
	 */
	private Map<String, AlterationPack> alterations;

	/**
	 * Network provider.
	 */
	private Traverse traverse;

	/**
	 * Constructor with network and alterations.
	 * @param traverse the network helper
	 * @param alterations alterations
	 * @param alterationThreshold minimum ratio of altered samples, should be between 0 and 1.
	 */
	public MutexExhaustiveSearcher(Traverse traverse, Map<String, AlterationPack> alterations,
		double alterationThreshold)
	{
		this.traverse = traverse;
		this.alterations = alterations;

		if (alterationThreshold > 0 && alterationThreshold < 1)
		{
			for (String s : new HashSet<String>(alterations.keySet()))
			{
				AlterationPack pack = alterations.get(s);

				// Filter out genes with les than 3% alteration
				if (pack.getAlteredRatio(Alteration.GENOMIC) < alterationThreshold)
					alterations.remove(s);
			}
		}
	}

	/**
	 * Searches subsets of the upstream genes up to a threshold size.
	 * @param pvalThr p-value threshold for member significance
	 * @return significant mutex groups
	 */
	public List<Group> search(double pvalThr, double minAltRatio, double maxGroupSize)
	{
		List<Group> groups = new ArrayList<Group>();

		Alteration[] alts = new Alteration[]{Alteration.ACTIVATING, Alteration.INHIBITING};

		for (String seed : alterations.keySet())
		{
			Set<Group> seedGroups = new HashSet<Group>();

//			if (!seed.equals("NDRG1")) continue;
			for (Alteration alt : alts)
			{
				AlterationPack gene = alterations.get(seed);
				if (gene.getAlteredRatio(alt) < minAltRatio) continue;

				Group group = new Group(new GeneAlt(gene, alt));

				addUpstreamToCandidates(group, gene, alts, minAltRatio);

				// expand the group to maximum 10 members
				while(expandGroup(group, alts, minAltRatio) && group.size() < maxGroupSize);

				if (group.size() > 1)
				{
					updateSeedCandidateValue(group, minAltRatio, alts);
					group.shrinkToSignificantMembers(pvalThr);
				}

				if (group.size() > 1) seedGroups.add(group);
			}

			// Get only one result per seed
			if (!seedGroups.isEmpty()) groups.add(selectMostCovered(seedGroups));
		}

		return groups;
	}

	/**
	 * Selects the group with the most coverage.
	 * @param groups groups to select from
	 * @return group with most coverage
	 */
	private Group selectMostCovered(Collection<Group> groups)
	{
		if (groups.size() == 1) return groups.iterator().next();

		Group best = null;

		for (Group group : groups)
		{
			if (best == null || group.calcCoverage() > best.calcCoverage())
			{
				best = group;
			}
		}
		return best;
	}

	/**
	 * Expands the group with the best candidate in its candidates list. Also adds upstream of the
	 * added gene among candidates of the group to use in next expansion.
	 * @param group groups to expand
	 * @return true if expanded
	 */
	private boolean expandGroup(Group group, Alteration[] alts, double minAltRatio)
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

			// Filter out candidates which are not fit anymore after addition

			for (GeneAlt cand : new HashSet<GeneAlt>(group.candidates))
			{
				if (!group.isOKToConsider(cand))
				{
					group.candidates.remove(cand);
					group.black.add(cand);
				}
			}

			// Add the upstream of the gene to the candidates for the next expansion cycle
			addUpstreamToCandidates(group, best.gene, alts, minAltRatio);

			return true;
		}
		else return false;
	}

	/**
	 * Gets the upstream of the given gene and adds the fitting one to the expansion candidates, and
	 * non-fitting ones to the black set.
	 * @param group group to update candidates
	 * @param gene gene to consider upstream
	 */
	private void addUpstreamToCandidates(Group group, AlterationPack gene, Alteration[] alts,
		double minAltRatio)
	{
		Set<String> upstream = getUpstream(Collections.singleton(gene.getId()), null);

		for (String cand : upstream)
		{
			for (Alteration alt : alts)
			{
				// the upstream gene is either a candidate or we don't want to re-consider it

				GeneAlt consider = new GeneAlt(alterations.get(cand), alt);

				if (consider.getAlteredRatio() > minAltRatio && group.isOKToConsider(consider))
				{
					group.candidates.add(consider);
				}
				else
				{
					group.black.add(consider);
				}
			}
		}
	}

	private void updateSeedCandidateValue(Group group, double minAlt, Alteration... alts)
	{
		List<String> genes = group.getGeneNames();
		genes.remove(0);
		Set<String> common = getCommonDownstream(genes);
		common.retainAll(alterations.keySet());

		Group g = group.copy();
		g.unique = null;
		g.removeGene(g.members.get(0));
		g.initUniqueCoverageMap();

		int cnt = 0;
		for (String c : common)
		{
			for (Alteration alt : alts)
			{
				GeneAlt gene = new GeneAlt(alterations.get(c), alt);
				if (gene.getAlteredRatio() > minAlt && g.isOKToConsider(gene)) cnt++;
			}
		}
		group.alternatives.put(group.members.get(0).getId(), cnt);
	}

	/**
	 * Gets the common downstream of the given genes in a transitive way, i.e. common downstream may
	 * be distant but is reachable traversing only the seed genes.
	 * @param seed seed to query
	 * @return transitive common downstream
	 */
	private Set<String> getCommonDownstream(Collection<String> seed)
	{
		Map<String, Set<String>> downMap = new HashMap<String, Set<String>>();

		for (String s : seed)
		{
			downMap.put(s, traverse.goBFS(Collections.singleton(s), null, true));
		}

		boolean again = true;

		while (again)
		{
			again = false;

			for (String s : seed)
			{
				for (String d : new HashSet<String>(downMap.get(s)))
				{
					if (downMap.containsKey(d) && !downMap.get(s).containsAll(downMap.get(d)))
					{
						downMap.get(s).addAll(downMap.get(d));
						again = true;
					}
				}
			}
		}

		Set<String> common = null;

		for (Set<String> downs : downMap.values())
		{
			if (common == null) common = downs;
			else common.retainAll(downs);
		}

		if (common == null)
		{
			System.out.println();
		}

		return common;
	}

	/**
	 * Gets upstream of the given genes.
	 * @param seed seed of the query
	 * @param black black list. we don't want these genes in the result
	 * @return upstream genes
	 */
	private Set<String> getUpstream(Set<String> seed, Set<String> black)
	{
		Set<String> up = traverse.goBFS(seed, black, false);
		up.retainAll(alterations.keySet());
		keepOnesWithAlterations(up);
		return up;
	}

	/**
	 * Filters out the genes with no alteration data associated.
	 * @param genes gene names to filter
	 */
	private void keepOnesWithAlterations(Set<String> genes)
	{
		for (String gene : new HashSet<String>(genes))
		{
			if (!alterations.containsKey(gene)) genes.remove(gene);
		}
	}
}
