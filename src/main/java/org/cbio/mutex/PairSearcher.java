package org.cbio.mutex;

import org.panda.utility.Kronometre;
import org.panda.utility.statistics.Overlap;

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
	 */
	public PairSearcher(Map<String, GeneAlt> genes)
	{
		this.genes = genes;
	}

	public void search(String filename) throws IOException
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

		Collections.sort(pairs, (o1, o2) -> pvals.get(o1).compareTo(pvals.get(o2)));

		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

		for (GeneAlt[] pair : pairs)
		{
			writer.write(pair[0].getId() + "\t" + pair[1].getId() + "\n");
		}

		writer.close();
	}

	public static void main(String[] args) throws IOException
	{
		Kronometre kron = new Kronometre();
		Main.dataFileName = "data-simulation/large-dataset/DataMatrix.txt";
		Map<String, GeneAlt> genes = Main.loadAlterations();
		PairSearcher ps = new PairSearcher(genes);
		ps.search("pair-results.txt");
		kron.stop();
		kron.print();
	}
}
