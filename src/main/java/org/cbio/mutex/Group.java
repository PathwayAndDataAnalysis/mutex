package org.cbio.mutex;

import org.panda.utility.ArrayUtil;
import org.panda.utility.FormatUtil;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;
import org.panda.utility.statistics.Overlap;

import java.io.*;
import java.util.*;

/**
 * A mutex group of genes.
 * @author Ozgun Babur
 */
public class Group implements Serializable
{
	private static final long serialVersionUID = 374555033512864697L;

	/**
	 * Genes in the group.
	 */
	List<GeneAlt> members;

	GeneAlt candidate;

	/**
	 * These are the genes that we do not want to consider for expanding.
	 */
	Set<GeneAlt> black;

	/**
	 * Common targets of the genes.
	 */
	List<String> targets;

	/**
	 * Seeds to remember. This set is useful after doing group merging.
	 */
	Set<GeneAlt> seedGenes;

	Map<GeneAlt, boolean[]> overlaps;

	boolean[] merge;

	/**
	 * Constructor with the seed gene.
	 * @param seed initial gene alteration
	 */
	public Group(GeneAlt seed)
	{
		this();
		addGene(seed);
	}

	/**
	 * Empty constructor that initializes sets and maps.
	 */
	public Group()
	{
		members = new ArrayList<GeneAlt>();
		black = new HashSet<GeneAlt>();
		overlaps = new HashMap<GeneAlt, boolean[]>();
	}

	/**
	 * Calculates p-values for each gene in the group.
	 * @return p-values
	 */
	public Map<String, Double> calcPVals1()
	{
		Map<String, Double> pvals = new HashMap<String, Double>();

		if (candidate == null)
		{
			if (members.size() == 1)
			{
				pvals.put(members.get(0).getId(), 1D);
			}
			else
			{
				int mergeCnt = ArrayUtil.countValue(merge, true);

				for (GeneAlt member : members)
				{
					int ov = ArrayUtil.countValue(overlaps.get(member), true);
					int a1 = member.getAltCnt();
					int a2 = mergeCnt - a1 + ov;

					double pval = Overlap.calcMutexPval(merge.length, ov, a1, a2);
					pvals.put(member.getId(), pval);
				}
			}
		}
		else
		{
			boolean[] cch = candidate.getBooleanChanges();
			int a2_pre = countMergeWithCandidate(cch);

			for (GeneAlt member : members)
			{
				int ov = countOverlapWithCandidate(overlaps.get(member), member.getBooleanChanges(), cch);
				int a1 = member.getAltCnt();
				int a2 = a2_pre - a1 + ov;

				double pval = Overlap.calcMutexPval(merge.length, ov, a1, a2);
				pvals.put(member.getId(), pval);
			}

			pvals.put(candidate.getId(), Overlap.calcMutexPval(cch, merge));
		}

		return pvals;
	}

	/**
	 * Calculates multiple-hypothesis-corrected p-values for each gene in the group.
	 */
	public Map<String, Double> calcPVals2(Map<String, Double> pvals1)
	{
		if (pvals1 == null) pvals1 = calcPVals1();

		Map<String, Double> pvals2 = new HashMap<String, Double>();

		for (GeneAlt member : members)
		{
			pvals2.put(member.getId(), member.getPvalOfScore(pvals1.get(member.getId())));
		}
		if (candidate != null)
			pvals2.put(candidate.getId(), candidate.getPvalOfScore(pvals1.get(candidate.getId())));

		for (String gene : pvals2.keySet())
		{
			if (pvals1.get(gene) > pvals2.get(gene)) pvals2.put(gene, pvals1.get(gene));
		}

		return pvals2;
	}

	private int countOverlapWithCandidate(boolean[] mov, boolean[] mch, boolean[] cch)
	{
		int cnt = 0;
		for (int i = 0; i < mov.length; i++)
		{
			if (mov[i] || (mch[i] && cch[i])) cnt++;
		}
		return cnt;
	}

	private int countMergeWithCandidate(boolean[] cch)
	{
		int cnt = 0;
		for (int i = 0; i < cch.length; i++)
		{
			if (merge[i] || cch[i]) cnt++;
		}
		return cnt;
	}

	public double calcScore()
	{
		if (size() == 1 && candidate == null) return 1;

		return getMaxValue(calcPVals1());
	}

	/**
	 * Gets the p-values in an array. This array does not contain repeated values, so it is only
	 * good for finding max or min.
	 * @return p-values
	 */
	public double getMaxValue(Map<String, Double> pvalMap)
	{
		double pv = 0;
		for (Double d : pvalMap.values())
		{
			if (d > pv) pv = d;
		}
		return pv;
	}


	/**
	 * Assuming the given gene is added to the group, calculates the new score for the group.
	 * Does not modify the group.
	 * @param gene gene alteration to consider
	 * @return geometric mean of the new p-values
	 */
	public double calcFutureScore(GeneAlt gene)
	{
		this.candidate = gene;
		double score = calcScore();
		this.candidate = null;
		return score;
	}

	public Map<String, Double> calcFuturePvals1(GeneAlt gene)
	{
		this.candidate = gene;
		Map<String, Double> vals = calcPVals1();
		this.candidate = null;
		return vals;
	}

	/**
	 * Assuming the given gene is added to the group, calculates the new pval for the group.
	 * Does not modify the group.
	 * @param gene gene alteration to consider
	 * @return geometric mean of the new p-values
	 */
	public double calcFutureFinalScore(GeneAlt gene)
	{
		this.candidate = gene;
		double pval = calcFinalScore();
		this.candidate = null;
		return pval;
	}

	public double calcFinalScore()
	{
		return getMaxValue(calcPVals2(null));
	}

	/**
	 * Adds the given gene alteration to the group. Updates the unique coverage map only if the
	 * addition is permanent.
	 * @param gene gene alteration to add
	 */
	public void addGene(GeneAlt gene)
	{
		assert !members.contains(gene);

		updateOverlaps(gene);
		members.add(gene);
	}

	public void updateOverlaps(GeneAlt gene)
	{
		boolean[] gch = gene.getBooleanChanges();
		boolean[] gov = new boolean[gch.length];

		for (GeneAlt member : members)
		{
			boolean[] ov = overlaps.get(member);
			boolean[] mch = member.getBooleanChanges();

			for (int i = 0; i < ov.length; i++)
			{
				if (gch[i] && mch[i] && !ov[i]) ov[i] = true;

				if (gch[i] && merge[i]) gov[i] = true;
			}
		}

		overlaps.put(gene, gov);

		if (merge == null) merge = gene.getBooleanChangesCopy();
		else
		{
			boolean[] c = gene.getBooleanChanges();
			for (int i = 0; i < merge.length; i++)
			{
				if (c[i]) merge[i] = true;
			}
		}
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
		// not ok if already a member and not ok if black-listed
		if (black.contains(gene) || members.contains(gene)) return false;

		boolean[] ch = gene.getBooleanChanges();
		for (int i = 0; i < ch.length; i++)
		{
			if (ch[i] && !merge[i])
			{
				return true;
			}
		}

		black.add(gene);
		// does not increase coverage
		return false;
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

		return members.size() < g.members.size() || calcScore() >= g.calcScore();
	}

	/**
	 * Calculates the coverage of the group.
	 * @return the coverage value between 0 and 1
	 */
	public double calcCoverage()
	{
		boolean[] merged = getMergedAlterations(-1);
		return ArrayUtil.countValue(merged, true) / (double) merged.length;
//		return ArrayUtil.countValue(merged, true);
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
		return getGeneNamesInString().replaceAll(" ", "").replaceAll(":", "");
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
		g.black.addAll(black);
		return g;
	}

	/**
	 * Gets an ordering for the samples to make the oncoprint look nicer.
	 * @return sample ordering for printing oncoprint
	 */
	private List<Integer> getPrintOrdering()
	{
		List<Integer> order = new ArrayList<Integer>();

		for (int i = 0; i < members.get(0).getBooleanChanges().length; i++)
		{
			order.add(i);
		}

		final boolean[][] marks = new boolean[members.get(0).getBooleanChanges().length][];

		for (int i = 0; i < marks.length; i++)
		{
			marks[i] = alterationMarks(i);
		}

		final boolean[][] mut = new boolean[members.size()][];
		final boolean[][] cna = new boolean[members.size()][];

		for (int i = 0; i < members.size(); i++)
		{
			mut[i] = members.get(i).getMutated();
			cna[i] = members.get(i).getCNAltered();
			if (cna[i] == null) cna[i] = new boolean[mut[i].length];
		}

		Collections.sort(order, (o1, o2) -> {
			boolean[] m1 = marks[o1];
			boolean[] m2 = marks[o2];

			int c = 0;
			for (int i = 0; i < members.size(); i++)
			{
				if (m1[i] && !m2[i]) c = -1;
				if (!m1[i] && m2[i]) c = 1;
				if (c != 0) break;
			}

			if (c != 0)
			{
				if (getNumberOfInitialPositiveAltOverlap(m1, m2) % 2 == 1) return -c;
				else return c;
			}

			for (int i = 0; i < members.size(); i++)
			{
				if (mut[i][o1] && !mut[i][o2]) return -1;
				if (!mut[i][o1] && mut[i][o2]) return 1;
				if (cna[i][o1] && !cna[i][o2]) return 1;
				if (!cna[i][o1] && cna[i][o2]) return -1;
			}

			return 0;
		});

		return order;
	}

	private boolean[] alterationMarks(int sample)
	{
		boolean[] b = new boolean[members.size()];
		for (int i = 0; i < b.length; i++)
		{
			b[i] = members.get(i).getBooleanChanges()[sample];
		}
		return b;
	}

	private int getNumberOfInitialPositiveAltOverlap(boolean[] m1, boolean[] m2)
	{
		int x = 0;

		for (int i = 0; i < members.size(); i++)
		{
			if (!m1[i] && !m2[i] && x == 0) continue;

			if (m1[i] && m2[i]) x++;
			else break;
		}
		return x;
	}

	/**
	 * Gets the oncoprint of the members in a String.
	 * @return oncoprint
	 */
	public String getPrint()
	{
		return getPrint(false);
	}

	/**
	 * Gets the oncoprint of the members in a String.
	 * @return oncoprint
	 */
	public String getPrint(boolean withTargets)
	{
		return getPrint(null, true, withTargets);
	}

	/**
	 * Gets the oncoprint of the members in a String.
	 * @param withMHT with multiple hypothesis testing
	 * @return oncoprint
	 */
	public String getPrint(Map<String, String> nameConvMap, boolean withMHT,
		boolean withTargets)
	{
		List<Integer> order = getPrintOrdering();
		Map<String, Double> p = calcPVals1();
		Map<String, Double> p2 = null;
		if (withMHT) p2 = calcPVals2(p);
		double score = calcScore();
		StringBuilder s = new StringBuilder();

		String names = getGeneNamesInString();
		if (nameConvMap != null) names = replaceNames(names, nameConvMap);

		s.append("[").append(names).append("]\tcoverage: ").
			append(FormatUtil.roundToSignificantDigits(calcCoverage(), 2)).
			append("\tscore: ").
			append(FormatUtil.roundToSignificantDigits(score, 2));
		if (withMHT) s.append("\tcorrected-score: ").
			append(FormatUtil.roundToSignificantDigits(calcFinalScore(), 2));
		if (withTargets) s.append("\ttargets:").append(getTargetLine(getTargets()));

		for (GeneAlt gene : members)
		{
			s.append("\n").append(gene.getPrint(order)).
				append((gene.getId().length() < 4) ? "  \t" : "\t").
				append("\tp1: ").
				append(FormatUtil.roundToSignificantDigits(p.get(gene.getId()), 2));
			if (withMHT) s.append("\tp2: ").
				append(FormatUtil.roundToSignificantDigits(p2.get(gene.getId()), 2));
		}
		return s.toString();
	}

	private String getTargetLine(List<String> list)
	{
		if (list.size() <= 10) return list.toString();
		else
		{
			String s = "[";

			for (int i = 0; i < 10; i++)
			{
				s += list.get(i) + ", ";
			}

			s += "and more]";
			return s;
		}
	}

	private String replaceNames(String names, Map<String, String> convMap)
	{
		StringBuilder s = new StringBuilder();
		for (String name : names.split(" "))
		{
			s.append(convMap.get(name)).append(" ");
		}
		return s.toString().trim();
	}

	// Section: static methods

	/**
	 * Removes the groups that are already covered by other groups in the given collection.
	 * @param groups groups to filter
	 */
	public static void removeSubsets(Collection<Group> groups)
	{
		for (Group group : new HashSet<>(groups))
		{
			group.initSeeds();
			for (Group other : groups)
			{
				if (group == other) continue;

				if (group.isSubsetOf(other))
				{
					other.mergeSeeds(group);
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
		Collections.sort(groups, (b1, b2) -> new Double(b2.calcCoverage()).compareTo(b1.calcCoverage()));
	}

	public static void serialize(List<Group> groups, String filename) throws IOException
	{
		ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(filename));
		out.writeObject(groups);
		out.close();
	}

	public static List<Group> deserialize(String filename) throws IOException, ClassNotFoundException
	{
		ObjectInputStream in = new ObjectInputStream(new FileInputStream(filename));
		List<Group> groups = (List<Group>) in.readObject();
		in.close();
		return groups;
	}

	public static Set<String> collectGenes(Collection<Group> groups)
	{
		Set<String> genes = new HashSet<>();
		for (Group group : groups)
		{
			genes.addAll(group.getGeneNames());
		}
		return genes;
	}

	public void fetchTragets(DirectedGraph traverse, Map<String, GeneAlt> genesMap)
	{
		Set<String> tars = traverse.getLinkedCommonDownstream(new HashSet<String>(getGeneNames()));

		if (!tars.removeAll(getGeneNames()))
		{
			targets = sortTargetsToFit(getGeneNames(), new ArrayList<String>(tars),
				genesMap);
//			targets = new ArrayList<String>(targets.subList(0, Math.min(targets.size(), 1)));
		}
		else targets = Collections.emptyList();
	}

	private List<String> sortTargetsToFit(Collection<String> mutex,
		Collection<String> comTar, Map<String, GeneAlt> genesMap)
	{
		final Map<String, Double> fit = new HashMap<>();

		for (String tar : comTar)
		{
			if (!genesMap.containsKey(tar)) fit.put(tar, 1D);
			else
			{
				GeneAlt gene = genesMap.get(tar);
				fit.put(tar, calcFutureScore(gene));
			}
		}

		List<String> sorted = new ArrayList<>(comTar);

		Collections.sort(sorted, (o1, o2) -> fit.get(o1).compareTo(fit.get(o2)));

		return sorted;
	}

	public List<String> getTargets()
	{
		return targets;
	}

	public void mergeSeeds(Group other)
	{
		initSeeds();
		other.initSeeds();
		seedGenes.addAll(other.seedGenes);
	}

	private void initSeeds()
	{
		if (this.seedGenes == null)
		{
			seedGenes = new HashSet<>();
			seedGenes.add(this.members.get(0));
		}
	}
}
