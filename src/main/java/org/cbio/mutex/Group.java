package org.cbio.mutex;

import org.cbio.causality.analysis.Graph;
import org.cbio.causality.cocitation.CocitationManager;
import org.cbio.causality.util.ArrayUtil;
import org.cbio.causality.util.FormatUtil;
import org.cbio.causality.util.Overlap;

import java.io.*;
import java.util.*;

/**
 * A mutex group of genes.
 * @author Ozgun Babur
 */
public class Group implements Serializable
{
	private final static CocitationManager coMan = new CocitationManager(
		365, "portal-cache/cocitation-cache");

	private static final long serialVersionUID = 350555030919624296L;

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
		double pval = calcScore();
		this.candidate = null;
		return pval;
	}

	/**
	 * Assuming the given gene is added to the group, calculates the new pval for the group.
	 * Does not modify the group.
	 * @param gene gene alteration to consider
	 * @return geometric mean of the new p-values
	 */
	public double calcFuturePVal(GeneAlt gene)
	{
		this.candidate = gene;
		double pval = calcPvalOfScore(calcScore());
		this.candidate = null;
		return pval;
	}

	public double calcPvalOfScore(double score)
	{
		double worst = 0;

		for (GeneAlt member : members)
		{
			double pval = member.getPvalOfScore(score);
			if (pval > worst) worst = pval;
		}
		if (candidate != null)
		{
			double pval = candidate.getPvalOfScore(score);
			if (pval > worst) worst = pval;
		}
		return worst;
	}

	public double calcPVal()
	{
		return calcPvalOfScore(calcScore());
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
	public String getPrint()
	{
		return getPrint(null);
	}

	/**
	 * Gets the oncoprint of the members in a String.
	 * @return oncoprint
	 */
	public String getPrint(SubtypeAligner sa)
	{
		List<Integer> order = getPrintOrdering();
		Map<String, Double> p = calcPVals1();
		double score = calcScore();
		StringBuilder s = new StringBuilder();

		s.append("[").append(getGeneNamesInString()).append("]\tcover: ").
			append(FormatUtil.roundToSignificantDigits(calcCoverage(), 2)).
			append("\tscore: ").
			append(FormatUtil.roundToSignificantDigits(score, 2)).
			append("\tcorrected-score: ").
			append(FormatUtil.roundToSignificantDigits(calcPVal(), 2)).
			append("\ttargets:").append(getTargets());

		for (GeneAlt gene : members)
		{
			s.append("\n").append(gene.getPrint(order)).
				append((gene.getId().length() < 4) ? "  \t" : "\t").
				append("\tp1: ").
				append(FormatUtil.roundToSignificantDigits(p.get(gene.getId()), 2)).
				append("\tp2: ").
				append(FormatUtil.roundToSignificantDigits(gene.getPvalOfScore(score), 2));
			if (sa != null)
			{
				List<String> subs = sa.getEnrichedSubtypes(gene.getId(), 0.05);
				if (subs != null) s.append("\t").append(subs);
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
		System.out.println("groups before removing subsets = " + groups.size());

		for (Group group : new HashSet<Group>(groups))
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

		System.out.println("groups after removing subsets = " + groups.size());
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
		Set<String> genes = new HashSet<String>();
		for (Group group : groups)
		{
			genes.addAll(group.getGeneNames());
		}
		return genes;
	}

	public void fetchTragets(Graph traverse, Map<String, GeneAlt> genesMap)
	{
		Set<String> tars = traverse.getLinkedCommonDownstream(new HashSet<String>(getGeneNames()));
		if (!tars.removeAll(getGeneNames()))
		{
			targets = sortTargetsToFitAndCocitation(getGeneNames(), new ArrayList<String>(tars),
				genesMap);
			targets = new ArrayList<String>(targets.subList(0, Math.min(targets.size(), 1)));
		}
		else targets = Collections.emptyList();
	}

	private List<String> sortTargetsToFitAndCocitation(Collection<String> mutex,
		Collection<String> comTar, Map<String, GeneAlt> genesMap)
	{
		final Map<String, Double> fit = new HashMap<String, Double>();

		for (String tar : comTar)
		{
			if (!genesMap.containsKey(tar)) fit.put(tar, 1D);
			else
			{
				GeneAlt gene = genesMap.get(tar);
				fit.put(tar, calcFutureScore(gene));
			}
		}

		Map<String, Map<String, Integer>> citMap = new HashMap<String, Map<String, Integer>>();

		for (String gene : mutex)
		{
			citMap.put(gene, coMan.getCocitations(gene));
		}

		List<String> sorted = new ArrayList<String>(comTar);

		final Map<String, Integer> citScore = new HashMap<String, Integer>();

		for (String tar : sorted)
		{
			citScore.put(tar, 0);
			for (Map<String, Integer> cits : citMap.values())
			{
				if (cits != null && cits.containsKey(tar))
				{
					citScore.put(tar, citScore.get(tar) + cits.get(tar));
				}
			}
		}

		Collections.sort(sorted, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				if (fit.get(o1).equals(fit.get(o2)))
				{
					if (citScore.get(o2).equals(citScore.get(o1)))
					{
						return o1.compareTo(o2);
					}
					else return citScore.get(o2).compareTo(citScore.get(o1));
				}
				else return fit.get(o1).compareTo(fit.get(o2));
			}
		});

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
			seedGenes = new HashSet<GeneAlt>();
			seedGenes.add(this.members.get(0));
		}
	}
}
