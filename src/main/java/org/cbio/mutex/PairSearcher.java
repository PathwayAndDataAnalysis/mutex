package org.cbio.mutex;

import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.util.Overlap;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class PairSearcher
{
	/**
	 * Gene alterations.
	 */
	private Map<String, GeneAlt> genes;

	/**
	 * Constructor with network and alterations.
	 * @param packs alterations
	 */
	public PairSearcher(Map<String, AlterationPack> packs)
	{
		this.genes = new HashMap<String, GeneAlt>();

		for (String s : new HashSet<String>(packs.keySet()))
		{
			AlterationPack pack = packs.get(s);

			GeneAlt gene = new GeneAlt(pack, Alteration.GENOMIC, null);
			this.genes.put(gene.getId(), gene);
		}

		System.out.println("filtered to min alt = " + this.genes.size());
	}


	public void search(String filename, Map<String, String> labelMap) throws IOException
	{
		List<GeneAlt[]> pairs = new ArrayList<GeneAlt[]>();
		for (GeneAlt g1 : genes.values())
		{
			for (GeneAlt g2 : genes.values())
			{
				if (g1.getId().compareTo(g2.getId()) < 0)
				{
					pairs.add(new GeneAlt[]{g1, g2});
				}
			}
		}

		final Map<GeneAlt[], Double> pvals = new HashMap<GeneAlt[], Double>();
		for (GeneAlt[] pair : pairs)
		{
			pvals.put(pair, Overlap.calcMutexPval(
				pair[0].getBooleanChanges(), pair[1].getBooleanChanges()));
		}

		Collections.sort(pairs, new Comparator<GeneAlt[]>()
		{
			@Override
			public int compare(GeneAlt[] o1, GeneAlt[] o2)
			{
				return pvals.get(o1).compareTo(pvals.get(o2));
			}
		});

		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

		for (GeneAlt[] pair : pairs)
		{
			writer.write(labelMap.get(pair[0].getId()) + "\t" + labelMap.get(pair[1].getId()) + "\n");
		}

		writer.close();
	}
}
