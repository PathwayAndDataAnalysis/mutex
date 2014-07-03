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
public class MutexSimpleRanker
{
	/**
	 * Gene alterations.
	 */
	private Map<String, GeneAlt> genes;

	/**
	 * Constructor with network and alterations.
	 * @param packs alterations
	 */
	public MutexSimpleRanker(Map<String, AlterationPack> packs)
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

	public void search(String filename, int maxGroupSize, Map<String, String> labelMap) throws IOException
	{
		List<Group> groups = new ArrayList<Group>();
		for (GeneAlt gene : genes.values())
		{
			Group group = getBestGroupOfGene(gene, maxGroupSize);
			if (group.size() > 1) groups.add(group);
		}

		Collections.sort(groups, new Comparator<Group>()
		{
			@Override
			public int compare(Group o1, Group o2)
			{
				return new Double(o1.calcScore()).compareTo(o2.calcScore());
			}
		});

		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

		for (Group group : groups)
		{
			List<String> names = group.getGeneNames();
			StringBuilder s = new StringBuilder();
			for (String name : names)
			{
				s.append(labelMap.get(name)).append("\t");
			}
			writer.write(s.toString().trim() + "\n");
		}

		writer.close();
	}

	private Group getBestGroupOfGene(GeneAlt gene, int maxGroupSize)
	{
		Group group = new Group(gene);

		while (group.size() < maxGroupSize && expand(group));
		return group;
	}

	private boolean expand(Group group)
	{
		double current = group.calcScore();
		GeneAlt bestGene = null;
		double bestScore = 1;

		for (GeneAlt gene : genes.values())
		{
			if (group.members.contains(gene)) continue;

			double score = group.calcFutureScore(gene);

			if (score < current && score < bestScore)
			{
				bestScore = score;
				bestGene = gene;
			}
		}

		if (bestGene != null)
		{
			group.addGene(bestGene);
			return true;
		}
		return false;
	}
}
