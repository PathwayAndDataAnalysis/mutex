package org.cbio.mutex;

import org.cbio.causality.analysis.SIFLinker;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.Change;
import org.cbio.causality.util.Overlap;
import org.cbio.causality.util.Summary;
import org.cbio.mutex.mutation.RandomAltGenerator;

import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.util.*;

/**
 * Writes mutex groups as a graph. to visualize in ChiBE.
 * @author Ozgun Babur
 */
public class GraphWriter
{
	/**
	 * Formatter for object values in the graph.
	 */
	static final DecimalFormat fmt = new DecimalFormat("0.00");

	/**
	 * Writes the graph data for the given groups and the given network provider to the given
	 * stream.
	 * @param groups groups of genes
	 * @param linker network helper
	 * @param os stream to write
	 */
	public static void write(List<Group> groups, double coocThr, SIFLinker linker, OutputStream os,
		String graphName, boolean groupCoocCliques) throws IOException
	{
		String s = "type:graph\ttext:" + graphName;
		Set<String> nodes = new HashSet<String>();
		Set<String> pairs = new HashSet<String>();

		for (Group group : groups)
		{
			if (group.size() > 2 || !hasLinkBetween(group, linker))
			{
				// write group id and members

				s += "\ntype:compound\tid:" + group.getID() + "\tmembers:";

				s += group.getGeneNamesInString().replaceAll(" ", ":");

				s += "\ttext:" + "cov: " + fmt.format(group.calcCoverage());
				s += "\ttextcolor:0 0 0\tbgcolor:255 255 255";
			}

			// write member nodes
			for (GeneAlt gene : group.members)
			{
				if (!nodes.contains(gene.getId()))
				{
					s += "\ntype:node\tid:" + gene.getId() + "\ttext:" + gene.getId();
					double rat = gene.getAlteredRatio();
					int v = 255 - (int) (rat * 255);

					boolean activating = isActivated(gene);

					String color = activating ?
						"255 " + v + " " + v : v + " " + v + " 255";
					s += "\tbgcolor:" + color + "\ttooltip:" + fmt.format(rat);
					nodes.add(gene.getId());
				}
			}

			// write relations
			Set<String> genes = new HashSet<String>(group.getGeneNames());
			List<String> relations = linker.linkProgressive(genes, genes, 0);
			for (String rel : relations)
			{
				String pairHash = rel.substring(0, rel.indexOf("\t")) +
					rel.substring(rel.lastIndexOf("\t"));

				if (!pairs.contains(pairHash))
				{
					String[] tok = rel.split("\t");
					GeneAlt src = group.getGene(tok[0]);
					GeneAlt trg = group.getGene(tok[2]);

					s+= "\ntype:edge\tid:" + rel.replaceAll("\t", " ") + "\tsource:" + src.getId() +
						"\ttarget:" + trg.getId() + "\tarrow:Target";

					double pv = Overlap.calcAlterationMutexPval(src.getChanges(), trg.getChanges());

					int v = (int) Math.max(0, 255 - (-Math.log(pv) * 55.4));
					String color = v + " " + v + " " + v;
					s += "\tlinecolor:" + color;

					pairs.add(pairHash);
				}
			}
		}

		// Put the co-occurrence relations
		Set<GeneAlt> set = new HashSet<GeneAlt>();
		Collections.reverse(groups);
		for (Group group : groups)
		{
			set.addAll(group.members);
		}

		List<Set<GeneAlt>> cliques = null;

		if (groupCoocCliques)
		{
			cliques = determineCliques(set, coocThr);
			Map<Set<GeneAlt>, String> cli2id = new HashMap<Set<GeneAlt>, String>();

			int i = 0;
			String idText = "cocNode";
			for (Set<GeneAlt> cli : cliques)
			{
				String id = idText + i++;
				s += "\ntype:node\tid:" + id + "\ttext:";

				s += "\tbgcolor:" + "255 255 255" + "\ttooltip:";
				cli2id.put(cli, id);
			}

			for (Set<GeneAlt> clique : cliques)
			{
				for (GeneAlt gene : clique)
				{
					double pv = calcAveragePVal(clique, gene);

					String cliID = cli2id.get(clique);
					s += "\ntype:edge\tid:" + gene.getId() + " co-occur " + cliID;
					s += "\tsource:" + gene.getId() + "\ttarget:" + cliID;

					int v = (int) Math.max(0, 255 - (-Math.log(pv) * 10));
					String color = v + " " + v + " " + v;
					s += "\tlinecolor:" + color + "\tstyle:Dashed";

					s += "\tarrow:None";
				}
			}
		}

		for (GeneAlt gene1 : set)
		{
			for (GeneAlt gene2 : set)
			{
				if (gene1.getId().compareTo(gene2.getId()) > 0)
				{
					if (groupCoocCliques && isPartOfAClique(gene1, gene2, cliques)) continue;

					double pv = Overlap.calcAlterationCoocPval(
						gene1.getChanges(), gene2.getChanges());

					if (pv < coocThr)
					{
						s += "\ntype:edge\tid:" + gene1.getId() + " co-occur " + gene2.getId();
						s += "\tsource:" + gene1.getId() + "\ttarget:" + gene2.getId();

						int v = (int) Math.max(0, 255 - (-Math.log(pv) * 10));
						String color = v + " " + v + " " + v;
						s += "\tlinecolor:" + color + "\tstyle:Dashed";

						s += "\tarrow:None";
					}
				}
			}
		}
		Collections.reverse(groups);

		// write the graph to the stream
		OutputStreamWriter writer = new OutputStreamWriter(os);
		writer.write(s);
		writer.close();
	}

	private static boolean hasLinkBetween(Group group, SIFLinker linker)
	{
		Set<String> genes = new HashSet<String>(group.getGeneNames());
		return !linker.linkProgressive(genes, genes, 0).isEmpty();
	}

	private static boolean isActivated(GeneAlt gene)
	{
		int cnAc = 0;
		int cnIn = 0;

		for (Change ch : gene.gene.get(Alteration.COPY_NUMBER))
		{
			if (ch == Change.ACTIVATING) cnAc++;
			else if (ch == Change.INHIBITING) cnIn++;
		}

		if (cnAc != cnIn) return cnAc > cnIn;

		cnIn = 0;

		for (Change ch : gene.gene.get(Alteration.MUTATION))
		{
			if (ch == Change.INHIBITING) cnIn++;
		}
		return cnIn * 10 < gene.size();
	}

	private static List<Set<GeneAlt>> determineCliques(Set<GeneAlt> genes, double coocThr)
	{
		List<Set<GeneAlt>> list = new ArrayList<Set<GeneAlt>>();

		for (GeneAlt gene1 : genes)
		{
			for (GeneAlt gene2 : genes)
			{
				if (gene1 == gene2) continue;

				double pv = Overlap.calcAlterationCoocPval(
					gene1.getChanges(), gene2.getChanges());

				if (pv < coocThr)
				{
					Set<GeneAlt> set = new HashSet<GeneAlt>();
					set.add(gene1);
					set.add(gene2);
					list.add(set);
				}
			}
		}

		for (GeneAlt gene : genes)
		{
			for (Set<GeneAlt> set : list)
			{
				boolean allSignif = true;

				for (GeneAlt member : set)
				{
					double pv = Overlap.calcAlterationCoocPval(
						gene.getChanges(), member.getChanges());

					if (pv > coocThr)
					{
						allSignif = false;
						break;
					}
				}

				if (allSignif) set.add(gene);
			}
		}

		List<Set<GeneAlt>> cliques = new ArrayList<Set<GeneAlt>>();

		for (Set<GeneAlt> set : list)
		{
			if (set.size() > 2 && !alreadyContains(cliques, set)) cliques.add(set);
		}
		return cliques;
	}

	private static boolean alreadyContains(List<Set<GeneAlt>> found, Set<GeneAlt> cand)
	{
		for (Set<GeneAlt> set : found)
		{
			if (set.size() >= cand.size() && set.containsAll(cand)) return true;
		}
		return false;
	}

	private static double calcAveragePVal(Set<GeneAlt> set)
	{
		double[] pv = new double[set.size() * (set.size() - 1) / 2];

		int i = 0;
		for (GeneAlt g1 : set)
		{
			for (GeneAlt g2 : set)
			{
				if (g1.getId().compareTo(g2.getId()) < 0)
				{
					pv[i++] = Overlap.calcAlterationCoocPval(g1.getChanges(), g2.getChanges());
				}
			}
		}
		return Summary.geometricMean(pv);
	}

	private static double calcAveragePVal(Set<GeneAlt> set, GeneAlt gene)
	{
		assert set.contains(gene);

		double[] pv = new double[set.size() - 1];

		int i = 0;
		for (GeneAlt g1 : set)
		{
			if (g1.equals(gene)) continue;

			pv[i++] = Overlap.calcAlterationCoocPval(g1.getChanges(), gene.getChanges());
		}
		return Summary.geometricMean(pv);
	}

	private static boolean isPartOfAClique(GeneAlt g1, GeneAlt g2, List<Set<GeneAlt>> cliques)
	{
		for (Set<GeneAlt> clique : cliques)
		{
			if (clique.contains(g1) && clique.contains(g2)) return true;
		}
		return false;
	}
}
