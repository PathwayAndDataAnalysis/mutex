package org.cbio.mutex;

import org.cbio.causality.data.portal.CBioPortalAccessor;
import org.cbio.causality.data.portal.CaseList;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.util.ArrayUtil;
import org.cbio.causality.util.FDR;
import org.cbio.causality.util.FormatUtil;
import org.cbio.causality.util.Overlap;

import java.io.IOException;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class SubtypeAligner
{
	private PortalDataset dataset;
	private Map<String, GeneAlt> genes;
	List<Group> groups;
	Map<String, List<String>> gene2type;
	Map<String, Map<String, Double>> pvals;
	List<String> names;

	private CBioPortalAccessor portal;
	private boolean[] hypermuts;

	public SubtypeAligner(PortalDataset dataset, List<Group> groups, boolean[] hypers) throws IOException
	{
		this.dataset = dataset;
		this.genes = new HashMap<String, GeneAlt>();

		for (Group group : groups)
		{
			for (GeneAlt member : group.members)
			{
				genes.put(member.getId(), member);
			}
		}

		this.groups = groups;
		this.portal = PortalReader.getPortalAccessor(dataset);
		this.hypermuts = hypers;
		calcEnrichments();
	}


	public SubtypeAligner(PortalDataset dataset, Collection<String> genes) throws IOException
	{
		this.dataset = dataset;
		this.genes = new HashMap<String, GeneAlt>();

		PortalReader pr = new PortalReader();
		Map<String,AlterationPack> packs = pr.readAlterations(dataset, genes);
		hypermuts = AltDistr.getOutlierAltered(packs.values());

		for (AlterationPack pack : packs.values())
		{
			this.genes.put(pack.getId(), new GeneAlt(pack, Alteration.GENOMIC, hypermuts));
		}

		this.portal = PortalReader.getPortalAccessor(dataset);
		calcEnrichments();
	}



	public void calcEnrichments() throws IOException
	{
		this.pvals = new HashMap<String, Map<String, Double>>();
		gene2type = new HashMap<String, List<String>>();
		names = new ArrayList<String>();

		for (int i = 0; i < dataset.subtypeCases.length; i++)
		{
			String subtype = dataset.subtypeCases[i].substring(dataset.subtypeCases[i].lastIndexOf("_") + 1);
			names.add(subtype);

			boolean[] sub = getSubtypeLoc(dataset.subtypeCases[i]);
			sub = trimToNonHyper(sub);

			Map<String, Double> pvals = new HashMap<String, Double>();
			for (String geneName : genes.keySet())
			{
				GeneAlt gene = genes.get(geneName);
				boolean[] alts = gene.getBooleanChanges();
				double pv = Overlap.calcCoocPval(alts, sub);
				pvals.put(geneName, pv);

				if (!this.pvals.containsKey(gene.getId())) this.pvals.put(gene.getId(), new HashMap<String, Double>());
				this.pvals.get(gene.getId()).put(subtype, pv);
			}

			List<String> names = FDR.selectBH(pvals, 0.05);

			for (String name : names)
			{
				if (!gene2type.containsKey(name)) gene2type.put(name, new ArrayList<String>());
				gene2type.get(name).add(subtype);
			}
		}
	}

	public List<String> getAllSubtypes()
	{
		return names;
	}

	public void printOverlapEffect() throws IOException
	{
		int cntP = 0;
		int cntN = 0;

		for (int i = 0; i < dataset.subtypeCases.length; i++)
		{
			boolean[] sub = getSubtypeLoc(dataset.subtypeCases[i]);
			sub = trimToNonHyper(sub);

			for (Group group : groups)
			{
				for (GeneAlt member : group.members)
				{
					boolean[] alts = member.getBooleanChanges();
					double pv1 = Overlap.calcCoocPval(alts, sub);

					alts = getNonOverlappingAlts(member, group);
					double pv2 = Overlap.calcCoocPval(alts, sub);

					if (pv1 < 0.05 || pv2 < 0.05)
					{
						if (pv1 < pv2) cntN++;
						else if (pv1 > pv2) cntP++;

						System.out.println(member + "\t" + group.getID() + "\t" +
							dataset.subtypeCases[i] + "\t" +
							FormatUtil.roundToSignificantDigits(pv1, 2) + "\t" +
							FormatUtil.roundToSignificantDigits(pv2, 2));
					}
				}
			}
		}
		System.out.println("cntP = " + cntP);
		System.out.println("cntN = " + cntN);
	}

	private boolean[] getNonOverlappingAlts(GeneAlt gene, Group group)
	{
		boolean[] ch = gene.getBooleanChangesCopy();
		for (GeneAlt member : group.members)
		{
			if (member == gene) continue;

			boolean[] oth = member.getBooleanChanges();

			for (int i = 0; i < ch.length; i++)
			{
				if (oth[i]) ch[i] = false;
			}
		}
		return ch;
	}

	private boolean[] getSubtypeLoc(String caseID) throws IOException
	{
		CaseList current = portal.getCurrentCaseList();
		for (CaseList cl : portal.getCaseListsForCurrentStudy())
		{
			if (cl.getId().equals(caseID))
			{
				Set<String> samples = new HashSet<String>(Arrays.asList(cl.getCases()));

				boolean[] b = new boolean[current.getCases().length];

				for (int i = 0; i < b.length; i++)
				{
					b[i] = samples.contains(current.getCases()[i]);
				}
				return b;
			}
		}
		return null;
	}

	public List<String> getEnrichedSubtypes(final String gene, double thr)
	{
		if (!pvals.containsKey(gene)) return Collections.emptyList();

		List<String> types = new ArrayList<String>();

		for (String type : pvals.get(gene).keySet())
		{
			Double p = pvals.get(gene).get(type);
			if (p < thr)
			{
				types.add(type);
			}
		}

		Collections.sort(types, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return pvals.get(gene).get(o1).compareTo(pvals.get(gene).get(o2));
			}
		});

		return types;
	}

	private boolean[] trimToNonHyper(boolean[] sub)
	{
		assert sub.length == hypermuts.length;

		int size = ArrayUtil.countValue(hypermuts, false);
		boolean[] b = new boolean[size];
		int i = 0;
		for (int j = 0; j < sub.length; j++)
		{
			if (!hypermuts[j]) b[i++] = sub[j];
		}
		return b;
	}

	public String getMostEnrichedSubtype(String gene, double thr)
	{
		if (!pvals.containsKey(gene)) return null;

		double pv = 1;
		String select = null;
		for (String type : pvals.get(gene).keySet())
		{
			Double p = pvals.get(gene).get(type);
			if (p < thr && pv > p)
			{
				select = type;
				pv = p;
			}
		}
		return select;
	}

	public static void main(String[] args) throws IOException
	{
		String s = "PIK3CA, TP53, MYC, PTK2, CDH1, PIP5K1A, RPS6KB1, PPM1D, EIF4EBP1, TLK2, BTG2, BCAS3, PPP2R5A, GATA3, GRB2, IL2RA, WDR37";
		Set<String> genes = new HashSet<String>(Arrays.asList(s.split(", ")));
		SubtypeAligner sa = new SubtypeAligner(PortalDataset.BRCA_PUB, genes);
		for (String gene : genes)
		{
			System.out.println(gene + "\t" + sa.getEnrichedSubtypes(gene, 0.05));
		}
	}
}
