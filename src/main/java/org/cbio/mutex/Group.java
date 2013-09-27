package org.cbio.mutex;

import org.cbio.causality.analysis.Traverse;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.util.ArrayUtil;
import org.cbio.causality.util.FishersCombinedProbability;
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
	 * Genes in the group.
	 */
	List<GeneAlt> members;

	/**
	 * These are the genes that we do not want to consider for expanding.
	 */
	Set<GeneAlt> black;

	/**
	 * These are the genes that we want to remember that we considered those as other candidates.
	 * These genes may not be counted towards multiple hypothesis testing for different reasons.
	 */
	Map<GeneAlt, Set<GeneAlt>> alternatives;

	/**
	 * These are the genes that we considered as alternative for the current members.
	 */
	Map<String, Set<GeneAlt>> candidateMap;

	/**
	 * Boolean arrays showing unique coverage of each  member.
	 */
	Map<String, boolean[]> unique;

	/**
	 * Gistic sets;
	 */
	private Map<String, Set<String>> gistic;

	/**
	 * Common targets of the genes.
	 */
	List<String> targets;

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
	public Group(GeneAlt seed, Map<String, Set<String>> gistic)
	{
		this();
		this.gistic = gistic;
		addGene(seed, new HashSet<GeneAlt>(), false);
		initUniqueCoverageMap();
	}

	/**
	 * Empty constructor that initializes sets and maps.
	 */
	public Group()
	{
		members = new ArrayList<GeneAlt>();
		black = new HashSet<GeneAlt>();
		candidateMap = new HashMap<String, Set<GeneAlt>>();
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
//		adjustToMultipleHypothesisTesting(pvals);
		return pvals;
	}

	/**
	 * Calculates p-values for each gene in the group, assuming that the member at the specified
	 * index is replaced with the given gene.
	 * @return p-values
	 */
	public Map<String, Double> calcPValsAsIf(GeneAlt alternative, int index)
	{
		Map<String, Double> pvals = new HashMap<String, Double>();

		if (members.size() == 2)
		{
			GeneAlt g0 = index != 0 ? members.get(0) : alternative;
			GeneAlt g1 = index != 1 ? members.get(1) : alternative;
			double p = Overlap.calcMutexPval(g0.getBooleanChanges(), g1.getBooleanChanges());

			pvals.put(g0.getId(), p);
			pvals.put(g1.getId(), p);
		}
		else
		{
			for (int i = 0; i < members.size(); i++)
			{
				boolean[] others = i == index? getMergedAlterations(i) :
					getMergedAlterations(i, alternative, index);

				GeneAlt gene = i == index ? alternative : members.get(i);
				pvals.put(gene.getId(), Overlap.calcMutexPval(gene.getBooleanChanges(), others));
			}
		}
		return pvals;
	}

	public double calcOverallPVal()
	{
		return calcOverallPVal(true);
	}

	public double calcOverallPVal(boolean correctForMultipleHypothesisTesting)
	{
		if (size() == 1) return 1;

//		return calcPairsOverallPVal(-1);

//		boolean[][][] b = prepareSetsToTest();
//		double[] pvals = new double[b.length];
//		for (int i = 0; i < pvals.length; i++)
//		{
//			pvals[i] = Overlap.calcCoocPval(b[i]);
//		}
//		double p = Summary.max(pvals);
//		if (size() > 2)
//			p = Math.pow(p, size()-1);
//		return p;

//		double[] pvals = calcPvalArray();
//
//		double pval = 1;
//
//		for (double v : pvals)
//		{
//			pval *= (1 - v);
//		}
//
//		pval = 1 - pval;

		double pval = Summary.max(calcPvalArray());
		if (correctForMultipleHypothesisTesting) pval = adjustToMultipleHypothesisTesting(pval);
		return pval;
	}

	public double calcOverallPVal(GeneAlt alternative, int index)
	{
		assert index < members.size();

		double pval = Summary.max(getPvalArray(calcPValsAsIf(alternative, index)));
		pval = adjustToMultipleHypothesisTesting(pval);
		return pval;
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

	public Set<GeneAlt> getMemberAlternatives(GeneAlt mem, double thr, Set<GeneAlt> cands)
	{
		Set<GeneAlt> set = new HashSet<GeneAlt>();

		for (GeneAlt cand : cands)
		{
			if (Overlap.calcCoocPval(mem.getBooleanChanges(), cand.getBooleanChanges()) < 0.01)
			{
				if (calcOverallPVal(cand, members.indexOf(mem)) <= thr) set.add(cand);
			}
		}
		return set;
	}

	/**
	 * Assuming the given gene is added to the group, calculates the new geometric mean of the
	 * p-values. Does not modify the group.
	 * @param gene gene alteration to consider
	 * @return geometric mean of the new p-values
	 */
	public double calcFuturePVal(GeneAlt gene, Set<GeneAlt> cands)
	{
		List<Set<GeneAlt>> sets = addGene(gene, cands, false);
		double pval = calcOverallPVal();
		removeGene(gene, sets);
		return pval;
	}

	public double calcPairsOverallPVal(int exclude)
	{
		double[] pvals = new double[members.size() * (members.size() - 1) / 2];

		if (members.size() == 2)
		{
			double p = Overlap.calcMutexPval(
				members.get(0).getBooleanChanges(),
				members.get(1).getBooleanChanges());

			pvals[0] = p;
		}
		if (members.size() == 3 && exclude >= 0)
		{
			double p = Overlap.calcMutexPval(
				members.get(exclude == 0 ? 1 : 0).getBooleanChanges(),
				members.get(exclude == 2 ? 1 : 2).getBooleanChanges());

			pvals[0] = p;
		}
		else
		{
			int k = 0;
			for (int i = 0; i < members.size() - 1; i++)
			{
				if (i == exclude) continue;

				for (int j = i + 1; j < members.size(); j++)
				{
					if (j == exclude) continue;

					Set<Integer> skip = new HashSet<Integer>(Arrays.asList(i, j, exclude));

					boolean[] ignore = getMergedAlterations(skip);
					boolean[] use = ArrayUtil.negate(ignore);
//					addInMultiOverlaps(use,
//						members.get(i).getBooleanChanges(),
//						members.get(j).getBooleanChanges());

					pvals[k++] = Overlap.calcMutexPval(
						members.get(i).getBooleanChanges(),
						members.get(j).getBooleanChanges(),
						use);
				}
			}
		}

		double pval = FishersCombinedProbability.pValue(pvals);
		return pval;
//		return adjustToMultipleHypothesisTesting(pval);
	}

//	private void addInMultiOverlaps(boolean[] use, boolean[] i, boolean[] j)
//	{
//		for (int k = 0; k < use.length; k++)
//		{
//			if (!use[k] && i[k] && j[k]) use[k] = true;
//		}
//	}

	/**
	 * These are candidates for alternatives.
	 * @param genes
	 */
	public void remember(Set<String> names, Map<String, GeneAlt> genes, double thr)
	{
		Set<GeneAlt> cands = new HashSet<GeneAlt>();
		for (String name : names)
		{
			if (genes.containsKey(name))
			{
				cands.add(genes.get(name));
			}
		}

		alternatives = new HashMap<GeneAlt, Set<GeneAlt>>();

		for (GeneAlt member : members)
		{
			alternatives.put(member, getMemberAlternatives(member, thr, cands));
		}
	}

	/**
	 * Adds the given gene alteration to the group. Updates the unique coverage map only if the
	 * addition is permanent.
	 * @param gene gene alteration to add
	 * @param permanent whether the addition is final
	 */
	public List<Set<GeneAlt>> addGene(GeneAlt gene, Set<GeneAlt> cands, boolean permanent)
	{
		if (!cands.isEmpty()) cands = filterCandsToQualify(cands, gene);

		members.add(gene);
		candidateMap.put(gene.getId(), cands);

		if (permanent) updateUniqueCoverageMap(gene);

		List<Set<GeneAlt>> candSets = new ArrayList<Set<GeneAlt>>();
		for (Set<GeneAlt> set : candidateMap.values())
		{
			if (set.contains(gene))
			{
				set.remove(gene);
				candSets.add(set);
			}
		}
		return candSets;
	}

	private Set<GeneAlt> filterCandsToQualify(Set<GeneAlt> cands, GeneAlt placed)
	{
		if (cands.isEmpty()) Collections.emptySet();

		// Remove the ones in the same gistic set with a member or placed gene.
		Set<GeneAlt> rem = new HashSet<GeneAlt>();

		// Filter using gistic overlaps

		if (gistic != null)
		{
			for (GeneAlt cand : cands)
			{
				Set<String> set = gistic.get(cand.getId());
				if (set != null)
				{
					if (set.contains(placed.getId())) rem.add(cand);

					for (GeneAlt member : members)
					{
						if (set.contains(member.getId())) rem.add(cand);
					}
				}
			}
		}

		cands = new HashSet<GeneAlt>(cands);
		cands.removeAll(rem);

		List<Set<GeneAlt>> sets = addGene(placed, Collections.<GeneAlt>emptySet(), false);
		double thr = calcOverallPVal(false);
		removeGene(placed, sets);

		boolean[] merged = getMergedAlterations(-1);
		int mCnt = ArrayUtil.countValue(merged, true);
		int n = placed.size();

		Set<GeneAlt> filtered = new HashSet<GeneAlt>();

		for (GeneAlt cand : cands)
		{
			int altCnt = cand.countAltered();

			int o = Math.max(mCnt + altCnt - n, 0);
			if (o > mCnt || o > altCnt)
			{
				System.out.println();
			}
			double p = Overlap.calcMutexPval(n, o, mCnt, altCnt);

			if (p <= thr) filtered.add(cand);
		}

		return filtered;
	}

	/**
	 * Removes a temporary gene from the group.
	 * @param gene gene alteration to remove
	 */
	public void removeGene(GeneAlt gene, List<Set<GeneAlt>> candsToAddTo)
	{
		assert unique == null || !unique.containsKey(gene.getId());
		members.remove(gene);
		candidateMap.remove(gene.getId());

		if (candsToAddTo != null)
		{
			for (Set<GeneAlt> cands : candsToAddTo)
			{
				cands.add(gene);
			}
		}
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
			boolean[] b = gene.getBooleanChangesCopy();
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

//	/**
//	 * Converts the change array to a boolean array where the value is true if altered.
//	 * @param gene gene alteration to convert
//	 * @return coverage array
//	 */
//	private boolean[] getCoverage(GeneAlt gene)
//	{
//		Change[] ch = gene.getChanges();
//		boolean[] b = new boolean[ch.length];
//		for (int i = 0; i < ch.length; i++)
//		{
//			if (ch[i].isAltered()) b[i] = true;
//		}
//		return b;
//	}

	/**
	 * Updates the unique coverage map with the given gene alteration.
	 * @param gene gene alteration to add
	 */
	public void updateUniqueCoverageMap(GeneAlt gene)
	{
		assert !unique.containsKey(gene.getId());

		boolean[] x = gene.getBooleanChangesCopy();
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

		boolean[] x = gene.getBooleanChanges();
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

				if (members.get(j).getBooleanChanges()[k])
				{
					others[k] = true;
					break;
				}
			}
		}
		return others;
	}

	/**
	 * Gets a merged change array for the genes in the group. Skips the gene with the given index.
	 * Assumes the member at the specified <code>index</code> is replaced with the given
	 * <code>alternative</code>.
	 * @param skipIndex index of the gene to skip. Use negative value if no skipping is required
	 * @return merged changes
	 */
	public boolean[] getMergedAlterations(int skipIndex, GeneAlt alternative, int altIndex)
	{
		assert skipIndex != altIndex;

		boolean[] others = new boolean[members.get(0).size()];

		for (int k = 0; k < others.length; k++)
		{
			others[k] = false;

			for (int j = 0; j < members.size(); j++)
			{
				if (j == skipIndex) continue;

				GeneAlt gene = j == altIndex ? alternative : members.get(j);

				if (gene.getBooleanChanges()[k])
				{
					others[k] = true;
					break;
				}
			}
		}
		return others;
	}

	/**
	 * Gets a merged change array for the genes in the group. Skips the gene with the given index.
	 * @param skipIndices indices of the genes to skip.
	 * @return merged changes
	 */
	public boolean[] getMergedAlterations(Set<Integer> skipIndices)
	{
		boolean[] others = new boolean[members.get(0).size()];

		for (int k = 0; k < others.length; k++)
		{
			others[k] = false;

			for (int j = 0; j < members.size(); j++)
			{
				if (skipIndices.contains(j)) continue;

				if (members.get(j).getBooleanChanges()[k])
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
		for (String gene : candidateMap.keySet())
		{
			assert candidateMap.containsKey(gene);

			int candidateSize = candidateMap.get(gene).size() + 1;
			pvals.put(gene, 1 - pow(1 - pvals.get(gene), candidateSize));
		}
	}

	/**
	 * Adjusts the p-values to the multiple hypothesis testing performed.
	 * @param pval overall p-value
	 */
	private double adjustToMultipleHypothesisTesting(double pval)
	{
		int candSize = 0;
		for (Set<GeneAlt> set : candidateMap.values())
		{
			for (GeneAlt member : members)
			{
				if (set.contains(member))
				{
					System.out.println();
				}
				assert !set.contains(member) : "seed: " + members.get(0);
			}
			candSize += set.size();
		}

		candSize++;

//		return pval * candSize;
		return 1 - pow(1 - pval, candSize);
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

		while(calcOverallPVal() > pvalThr)
		{
			removeGene(members.get(members.size() - 1), null);

			if (size() == 1) return;
		}
	}

	/**
	 * Gets the member size of the group.
	 * @return the size
	 */
	public int size()
	{
		int s = members.size();
		if (alternatives != null)
		{
			for (Set<GeneAlt> altSet : alternatives.values())
			{
				s += altSet.size();
			}
		}
		return s;
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
			if (!g.members.contains(gene) && !g.alternativesContain(gene)) return false;
		}
		for (Set<GeneAlt> altSet : alternatives.values())
		{
			for (GeneAlt gene : altSet)
			{
				if (!g.members.contains(gene) && !g.alternativesContain(gene)) return false;
			}
		}

		if (members.size() == g.members.size() && alternatives.size() == g.alternatives.size())
			return calcOverallPVal() >= g.calcOverallPVal();

		return true;
	}

	private boolean alternativesContain(GeneAlt gene)
	{
		for (Set<GeneAlt> set : alternatives.values())
		{
			if (set.contains(gene)) return true;
		}
		return false;
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
			for (GeneAlt al : alternatives.get(gene))
			{
				s += ":" + al.getId();
			}
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

			if (alternatives != null)
				for (GeneAlt al : alternatives.get(gene))
				{
					names.add(al.getId());
				}
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
		for (Set<GeneAlt> altSet : alternatives.values())
		{
			for (GeneAlt gene : altSet)
			{
				if (gene.getId().equals(id)) return gene;
			}
		}
		return null;
	}

	public List<GeneAlt> getAllGenes()
	{
		List<GeneAlt> genes= new ArrayList<GeneAlt>(members);
		for (Set<GeneAlt> altSet : alternatives.values())
		{
			genes.addAll(altSet);
		}
		return genes;
	}

	/**
	 * Gets a copy of the group.
	 * @return a copy
	 */
	public Group copy()
	{
		Group g = new Group();
		g.members.addAll(members);
		g.candidateMap.putAll(candidateMap);
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
			boolean[] ch = gene.getBooleanChanges();

			for (int i = 0; i < ch.length; i++)
			{
				if (ch[i] && !order.contains(i)) order.add(i);
			}
		}
		return order;
	}

	/**
	 * Gets the oncoprint of the members in a String.
	 * @return oncoprint
	 */
	public String getPrint(double thr)
	{
		List<Integer> order = getPrintOrdering();
		Map<String, Double> p = calcPVals();
		StringBuilder s = new StringBuilder();
		for (GeneAlt gene : members)
		{
			if (s.length() > 0) s.append("\n");
			s.append(gene.getPrint(order)).
				append((gene.getId().length() < 4) ? "  \t" : "\t").
				append("\tp-val: ").
				append(fmt.format(p.get(gene.getId())));

			for (GeneAlt cand : alternatives.get(gene))
			{
				s.append("\n");
				s.append(cand.getPrint(order));
			}
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

	public void fetchTragets(Traverse traverse, final Map<String, GeneAlt> genesMap)
	{
		targets = new ArrayList<String>(traverse.getLinkedCommonDownstream(
			new HashSet<String>(getGeneNames())));

		Collections.sort(targets, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				GeneAlt pack1 = genesMap.get(o1);
				GeneAlt pack2 = genesMap.get(o2);

				Integer cnt1 = pack1 == null ? 0 : pack1.countAltered();
				Integer cnt2 = pack2 == null ? 0 : pack2.countAltered();

				if (cnt1.equals(cnt2)) return o1.compareTo(o2);
				else return cnt2.compareTo(cnt1);
			}
		});
	}

	public List<String> getTargets()
	{
		return targets;
	}
}
