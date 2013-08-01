package org.cbio.mutex;

import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.Change;
import org.cbio.causality.util.ChiSquare;
import org.cbio.causality.util.Overlap;
import org.cbio.causality.util.Summary;

import java.text.DecimalFormat;
import java.util.*;

/**
 * A mutex group of genes.
 * @author Ozgun Babur
 */
public class Group
{
	/**
	 * Other candidates for this position.
	 */
	Map<String, Integer> alternatives;

	/**
	 * Genes in the group.
	 */
	List<GeneAlt> members;

	/**
	 * These are the genes that we do not want to consider for expanding.
	 */
	Set<GeneAlt> black;

	/**
	 * These are the genes that we consider for expanding in the next cycle.
	 */
	Set<GeneAlt> candidates;

	/**
	 * Boolean arrays showing unique coverage of each  member.
	 */
	Map<String, boolean[]> unique;

	/**
	 * Used to refer the merge of genes.
	 */
	public static final String MERGE = "MERGE";

	/**
	 * Used for formatting print output.
	 */
	static final DecimalFormat fmt = new DecimalFormat("0.0000");

	/**
	 * Constructor with the seed gene.
	 * @param seed initial gene alteration
	 */
	public Group(GeneAlt seed)
	{
		this();
		addGene(seed, 1, false);
		initUniqueCoverageMap();
	}

	/**
	 * Empty constructor that initializes sets and maps.
	 */
	public Group()
	{
		members = new ArrayList<GeneAlt>();
		alternatives = new HashMap<String, Integer>();
		black = new HashSet<GeneAlt>();
		candidates = new HashSet<GeneAlt>();
		unique = new HashMap<String, boolean[]>();
	}

	/**
	 * Calculates p-values for each gene in the group.
	 * @return p-values
	 */
	public Map<String, Double> calcPVals()
	{
		Map<String, Double> pvals = new HashMap<String, Double>();

		if (members.size() == 2)
		{
			double p = Overlap.calcMutexPval(
				members.get(0).getBooleanChanges(),
				members.get(1).getBooleanChanges());

			pvals.put(members.get(0).getId(), p);
			pvals.put(members.get(1).getId(), p);
		}
		else
		{
			for (int i = 0; i < members.size(); i++)
			{
				boolean[] others = getMergedAlterations(i);

				pvals.put(members.get(i).getId(), Overlap.calcMutexPval(
					members.get(i).getBooleanChanges(), others));
			}
		}
		adjustToMultipleHypothesisTesting(pvals);
		return pvals;
	}

	public double calcOverallPVal()
	{
		if (size() == 1) return 1;

		double[] pvals = calcPvalArray();

		double pval = 1;

		for (double v : pvals)
		{
			pval *= (1 - v);
		}

		pval = 1 - pval;

		return pval;
	}

	public double calcOverallPVal_old()
	{
		if (size() == 1) return 1;

		double[] pvals = calcIndependentPVals();

		if (pvals.length == 1) return pvals[0];

		double chisq = 0;

		for (double pval : pvals)
		{
			chisq += Math.log(pval);
		}

		chisq *= -2;

		return ChiSquare.pValue(chisq, 2 * pvals.length);
	}

	/**
	 * Calculates n-1 independent p-values for mutexness.
	 * @return p-values
	 */
	private double[] calcIndependentPVals()
	{
		boolean[][] alts = new boolean[members.size()][];
		for (int i = 0; i < alts.length; i++)
		{
			alts[i] = members.get(i).getBooleanChanges();
		}
		return calcIndependentPVals(alts);
	}

	/**
	 * Calculates n-1 independent p-values for mutexness.
	 * @param alts gene alterations
	 * @return p-values
	 */
	public double[] calcIndependentPVals(boolean[][] alts)
	{
		double[] pvals = new double[alts.length - 1];

		// prepare usage arrays

		boolean[] useGene = new boolean[alts.length];
		boolean[] useSample = new boolean[alts[0].length];

		for (int i = 0; i < useGene.length; i++) useGene[i] = true;
		for (int i = 0; i < useSample.length; i++) useSample[i] = true;


		for (int i = 0; i < pvals.length; i++)
		{
			useGene[i] = false;

			boolean[] single = alts[i];
			boolean[] merge = getMergedAlterations(alts, useGene);

			pvals[i] = Overlap.calcMutexPval(single, merge, useSample);

			for (int j = 0; j < useSample.length; j++)
			{
				if (single[j]) useSample[j] = false;
			}
		}

		return pvals;
	}

	/**
	 * Calculates the geometric mean of the p-values.
	 * @return geometric mean of p-values
	 */
	public double calcPvalGeomMean()
	{
		return Summary.geometricMean(calcPvalArray());
	}

	/**
	 * Calculates the highest of the p-values.
	 * @return highest of the p-values
	 */
	public double calcWorstPval()
	{
		return Summary.max(calcPvalArray());
	}

	/**
	 * Gets the p-values in an array.
	 * @return p-values
	 */
	public double[] calcPvalArray()
	{
		return getPvalArray(calcPVals());
	}

	/**
	 * Gets the p-values in an array.
	 * @return p-values
	 */
	public double[] getPvalArray(Map<String, Double> pvalMap)
	{
		double[] pv = new double[pvalMap.size()];
		int i = 0;
		for (Double d : pvalMap.values())
		{
			pv[i++] = d;
		}
		return pv;
	}

	/**
	 * Assuming the given gene is added to the group, calculates the new geometric mean of the
	 * p-values. Does not modify the group.
	 * @param gene gene alteration to consider
	 * @return geometric mean of the new p-values
	 */
	public double calcFuturePVal(GeneAlt gene, int candSize)
	{
		addGene(gene, candSize, false);
		double pval = calcOverallPVal();
		removeGene(gene);
		return pval;
	}

	/**
	 * Adds the given gene alteration to the group. Updates the unique coverage map only if the
	 * addition is permanent.
	 * @param gene gene alteration to add
	 * @param permanent whether the addition is final
	 */
	public void addGene(GeneAlt gene, int candSize, boolean permanent)
	{
		members.add(gene);
		alternatives.put(gene.getId(), candSize);

		if (permanent)
		{
			removeCandidate(gene.getId());
			updateUniqueCoverageMap(gene);
		}
	}

	/**
	 * Removes candidates with the given name.
	 * @param gene the gene name to remove
	 */
	public void removeCandidate(String gene)
	{
		for (GeneAlt cand : new HashSet<GeneAlt>(candidates))
		{
			if (cand.getId().equals(gene)) candidates.remove(cand);
		}
	}

	/**
	 * Removes a temporary gene from the group.
	 * @param gene gene alteration to remove
	 */
	public void removeGene(GeneAlt gene)
	{
		assert unique == null || !unique.containsKey(gene.getId());
		members.remove(gene);
		alternatives.remove(gene.getId());
	}

	/**
	 * Record unique coverage of each member, and the merge data.
	 */
	public void initUniqueCoverageMap()
	{
		unique = new HashMap<String, boolean[]>();

		boolean[][] a = new boolean[members.size()][];

		int k = 0;
		for (GeneAlt gene : members)
		{
			boolean[] b = getCoverage(gene);
			a[k++] = b;
			unique.put(gene.getId(), b);
		}

		boolean[] merge = new boolean[a[0].length];

		for (int i = 0; i < a[0].length; i++)
		{
			// record the merged coverage
			int n = 0;
			for (boolean[] b : a)
			{
				if (b[i])
				{
					merge[i] = true;
					n++;
				}
			}

			// if coverage at that sample is not unique then mark it
			if (n > 1)
			{
				for (boolean[] b : a)
				{
					b[i] = false;
				}
			}
		}

		unique.put(MERGE, merge);
	}

	/**
	 * Converts the change array to a boolean array where the value is true if altered.
	 * @param gene gene alteration to convert
	 * @return coverage array
	 */
	private boolean[] getCoverage(GeneAlt gene)
	{
		Change[] ch = gene.getChanges();
		boolean[] b = new boolean[ch.length];
		for (int i = 0; i < ch.length; i++)
		{
			if (ch[i].isAltered()) b[i] = true;
		}
		return b;
	}

	/**
	 * Updates the unique coverage map with the given gene alteration.
	 * @param gene gene alteration to add
	 */
	public void updateUniqueCoverageMap(GeneAlt gene)
	{
		assert !unique.containsKey(gene.getId());

		boolean[] x = getCoverage(gene);
		boolean[] merge = unique.get(MERGE);
		assert merge != null;

		for (int i = 0; i < x.length; i++)
		{
			if (x[i])
			{
				for (boolean[] b : unique.values())
				{
					if (b == merge) continue;

					if (b[i])
					{
						b[i] = false;
						x[i] = false;
						assert merge[i] : "merge[i] should be true if b[i] is true.";
						break;
					}
				}

				if (x[i]) merge[i] = true;
			}
		}

		unique.put(gene.getId(), x);
	}

	/**
	 * Checks if the given candidate gene alteration can be part of the group. The candidate has to
	 * contribute with alteration of unique samples, and should not completely cover the unique
	 * contributions of any existing member.
	 * @param gene candidate gene alteration
	 * @return true if the gene alteration can be considered for expansion
	 */
	public boolean isOKToConsider(GeneAlt gene)
	{
		// not ok if already a member
		if (getGeneNames().contains(gene.getId())) return false;

		// not ok if black-listed
		if (black.contains(gene)) return false;

		boolean[] x = getCoverage(gene);
		boolean[] merge = unique.get(MERGE);

		// if all alterations already covered by the merge of genes it is not ok
		if (secondIsSubset(merge, x)) return false;

		for (boolean[] b : unique.values())
		{
			if (b == merge) continue;

			// if a gene in the group becomes non-contributing after addition of that gene, then it
			// is not ok
			if (secondIsSubset(x, b)) return false;
		}
		return true;
	}

	/**
	 * Checks if the second array is a subset of the first one, i.e. if second[i] is true then
	 * first[i] is also true.
	 * @param first first array
	 * @param second second array
	 * @return true if second is subset
	 */
	private boolean secondIsSubset(boolean[] first, boolean[] second)
	{
		for (int i = 0; i < second.length; i++)
		{
			if (second[i] && !first[i]) return false;
		}
		return true;
	}

	/**
	 * Gets a merged change array for the genes in the group. Skips the gene with the given index.
	 * @param skipIndex index of the gene to skip. Use negative value if no skipping is required
	 * @return merged changes
	 */
	public boolean[] getMergedAlterations(int skipIndex)
	{
		boolean[] others = new boolean[members.get(0).size()];

		for (int k = 0; k < others.length; k++)
		{
			others[k] = false;

			for (int j = 0; j < members.size(); j++)
			{
				if (j == skipIndex) continue;

				if (members.get(j).getChanges()[k].isAltered())
				{
					others[k] = true;
					break;
				}
			}
		}
		return others;
	}

	/**
	 * Gets a merged change array for the subset of the genes in the group.
	 * @param alts gene alterations
	 * @param useGenes indices of the gene to merge alterations
	 * @return merged changes
	 */
	public boolean[] getMergedAlterations(boolean[][] alts, boolean[] useGenes)
	{
		boolean[] merge = new boolean[alts[0].length];
		for (int i = 0; i < useGenes.length; i++)
		{
			if (!useGenes[i]) continue;

			ORWith(merge, alts[i]);
		}

		return merge;
	}

	/**
	 * Merges the second array into the first one by OR operation.
	 * @param thisOne array to modify
	 * @param with array to add with OR
	 */
	private void ORWith(boolean[] thisOne, boolean[] with)
	{
		for (int i = 0; i < with.length; i++)
		{
			if (!thisOne[i] && with[i]) thisOne[i] = true;
		}
	}

	/**
	 * Adjusts the p-values to the multiple hypothesis testing performed.
	 * @param pvals p-values
	 */
	private void adjustToMultipleHypothesisTesting(Map<String, Double> pvals)
	{
		for (String gene : alternatives.keySet())
		{
			assert alternatives.containsKey(gene);

			int candidateSize = alternatives.get(gene);
			pvals.put(gene, 1 - pow(1 - pvals.get(gene), candidateSize));
		}
	}

	/**
	 * Takes power.
	 * @param v number to take its power
	 * @param exp exponent
	 * @return result
	 */
	private double pow(double v, int exp)
	{
		double e = 1;
		for (int i = 0; i < exp; i++)
		{
			e *= v;
		}
		return e;
	}

	/**
	 * Removes elements from the group in the reverse order of addition, until all the members are
	 * significant. Does nothing if all members are significant. It is assumed that this group won't
	 * re-expand after shrinking.
	 * @param pvalThr p-value threshold for significance
	 */
	public void shrinkToSignificantMembers(double pvalThr)
	{
		// Nullify the unique coverage map. This guarantees to get an error if we try to re-expand.
		unique = null;

		while(calcWorstPval() > pvalThr)
		{
			removeGene(members.get(members.size() - 1));

			if (size() == 1) return;
		}
	}

	/**
	 * Gets the member size of the group.
	 * @return the size
	 */
	public int size()
	{
		return members.size();
	}

	/**
	 * Checks if this group is a subset of the given group.
	 * @param g group to check
	 * @return true if this group is a subset of the given group
	 */
	public boolean isSubsetOf(Group g)
	{
		for (GeneAlt gene : members)
		{
			if (!g.members.contains(gene)) return false;
		}
		return true;
	}

	/**
	 * Calculates the coverage of the group.
	 * @return the coverage value between 0 and 1
	 */
	public double calcCoverage()
	{
		if (unique == null)
		{
			unique = new HashMap<String, boolean[]>();
			initUniqueCoverageMap();
		}

		return Summary.countTrue(unique.get(MERGE)) / (double) unique.get(MERGE).length;
	}

	/**
	 * Gets the names of members in a String.
	 * @return member names
	 */
	public String getGeneNamesInString()
	{
		String s = "";
		for (GeneAlt gene : members)
		{
			s += " " + gene.getId();
		}
		return s.trim();
	}

	/**
	 * Gets the names of members.
	 * @return member names
	 */
	public List<String> getGeneNames()
	{
		List<String> names = new ArrayList<String>(members.size());
		for (GeneAlt gene : members)
		{
			names.add(gene.getId());
		}
		return names;
	}

	/**
	 * Uses member gene names as id
	 * @return member gene names as id
	 */
	public String getID()
	{
		return getGeneNamesInString();
	}

	/**
	 * Gets the member gene alteration with a name match.
	 * @param id name of the member gene
	 * @return member gene that matches
	 */
	public GeneAlt getGene(String id)
	{
		for (GeneAlt gene : members)
		{
			if (gene.getId().equals(id)) return gene;
		}
		return null;
	}

	/**
	 * Gets a copy of the group.
	 * @return a copy
	 */
	public Group copy()
	{
		Group g = new Group();
		g.members.addAll(members);
		g.candidates.addAll(candidates);
		g.alternatives.putAll(alternatives);
		g.black.addAll(black);
		g.unique.putAll(unique);
		return g;
	}

	// Section oncoprint

	/**
	 * Gets an ordering for the samples to make the oncoprint look nicer.
	 * @return sample ordering for printing oncoprint
	 */
	private List<Integer> getPrintOrdering()
	{
		List<Integer> order = new ArrayList<Integer>();

		for (GeneAlt gene : members)
		{
			Change[] ch = gene.getChanges();

			for (int i = 0; i < ch.length; i++)
			{
				if (ch[i].isAltered() && !order.contains(i)) order.add(i);
			}
		}
		return order;
	}

	/**
	 * Gets the oncoprint of the members in a String.
	 * @return oncoprint
	 */
	public String getPrint()
	{
		List<Integer> order = getPrintOrdering();
		double[] p = calcPvalArray();
		StringBuilder s = new StringBuilder();
		for (GeneAlt gene : members)
		{
			if (s.length() > 0) s.append("\n");
			s.append(gene.gene.getPrint(gene.alt, order)).
				append((gene.getId().length() < 4) ? "  \t" : "\t").
				append((gene.alt == Alteration.ACTIVATING) ? "+" : "-").append("\tp-val: ").
				append(fmt.format(p[members.indexOf(gene)]));
		}
		return s.toString();
	}


	// Section: static methods

	/**
	 * Removes the groups that are already covered by other groups in the given collection.
	 * @param groups groups to filter
	 */
	public static void removeSubsets(Collection<Group> groups)
	{
		for (Group group : new HashSet<Group>(groups))
		{
			for (Group other : groups)
			{
				if (group == other) continue;

				if (group.isSubsetOf(other))
				{
					groups.remove(group);
					break;
				}
			}
		}
	}

	/**
	 * Sorts the given groups to coverage.
	 * @param groups groups to sort
	 */
	public static void sortToCoverage(List<Group> groups)
	{
		Collections.sort(groups, new Comparator<Group>()
		{
			@Override
			public int compare(Group b1, Group b2)
			{
				return new Double(b2.calcCoverage()).compareTo(b1.calcCoverage());
			}
		});
	}
}
