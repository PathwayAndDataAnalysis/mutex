package org.cbio.mutex;

import org.cbio.causality.analysis.SIFLinker;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.util.Overlap;

import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.util.*;

/**
 * Writes mutex groups as a graph.
 * @author Ozgun Babur
 */
public class GraphWriter
{
	/**
	 * Formatter for object values in the graph.
	 */
	static final DecimalFormat fmt = new DecimalFormat("0.0000");

	/**
	 * Writes the graph data for the given groups and the given network provider to the given
	 * stream.
	 * @param groups groups of genes
	 * @param linker network helper
	 * @param os stream to write
	 */
	public static void write(List<Group> groups, SIFLinker linker, OutputStream os,
		String graphName) throws IOException
	{
		String s = "type:graph\ttext:" + graphName;
		Set<String> nodes = new HashSet<String>();
		Set<String> edges = new HashSet<String>();

		for (Group group : groups)
		{
			if (group.size() > 2)
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

					int inh = gene.gene.getAlteredCount(Alteration.INHIBITING);
					int all = gene.gene.getAlteredCount(Alteration.GENOMIC);

					boolean activating =
						gene.alt == Alteration.ACTIVATING || (all - inh) * 10 < all;

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
				if (!edges.contains(rel))
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

					edges.add(rel);
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
		for (GeneAlt gene1 : set)
		{
			for (GeneAlt gene2 : set)
			{
				if (gene1.getId().compareTo(gene2.getId()) > 0)
				{
					double pv = Overlap.calcAlterationCoocPval(
						gene1.getChanges(), gene2.getChanges());

					if (pv < 0.05)
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
}
