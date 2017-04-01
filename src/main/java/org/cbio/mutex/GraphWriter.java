package org.cbio.mutex;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.utility.FormatUtil;
import org.panda.utility.statistics.Overlap;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.*;

/**
 * Writes mutex groups as a graph. to visualize in ChiBE.
 * @author Ozgun Babur
 */
public class GraphWriter
{
	private static final List<String> H_COLS = Arrays.asList("255 255 205", "255 205 255",
		"205 255 255", "205 205 255", "255 205 205", "205 255 205");

	/**
	 * Writes the graph data for the given groups and the given network provider to the given
	 * stream.
	 * @param groups groups of genes
	 * @param network network data
	 */
	public static void write(List<Group> groups, Map<String, GeneAlt> geneMap,
		Network network, String dir, String graphName) throws IOException
	{
		// All genes in mutex groups, ignoring targets
		Set<String> genesInGroups = new HashSet<>();

		Set<String> to = getTargetOnly(groups);
		OutputStreamWriter writer = new OutputStreamWriter(new FileOutputStream(dir + "merged-network.sif"));

		Map<String, String> node2color = new HashMap<>();
		Map<String, String> node2tooltip = new HashMap<>();

		String s = "type:graph\ttext:" + graphName;
//		Set<String> nodes = new HashSet<>();

		Map<String, Set<String>> name2id = new HashMap<>();

		for (Group group : groups)
		{
//			if (group.getGeneNames().contains("PTEN")) continue;

			// write group id and members

			String groupID = group.getID();
			s += "\ntype:compound\tid:" + groupID + "\tmembers:";

			for (String mem : group.getGeneNamesInString().replaceAll(":", " ").split(" "))
			{
				s += mem + groupID + ":";
			}

			s += "\ttext:" + (int) Math.round(group.calcCoverage() * 100) + "%";
			s += "\ttextcolor:0 0 0\tbgcolor:255 255 255";

			// write member nodes
			for (GeneAlt gene : group.members)
			{
				s += writeNodeData(gene, groupID, name2id,
					group.seedGenes.contains(gene), node2color, node2tooltip);
			}

			List<String> memNames = group.getGeneNames();
			genesInGroups.addAll(memNames);
			List<String> tarNames = new ArrayList<>();

			// write target nodes
			for (String tarName : group.getTargets())
			{
				if (!memNames.contains(tarName))
				{
					tarNames.add(tarName);
				}
				// limit targets to one
				if (tarNames.size() == 1) break;
			}

			for (String tarName : tarNames)
			{
				if (geneMap.containsKey(tarName))
					s += writeNodeData(geneMap.get(tarName), groupID, name2id, false, node2color,
						node2tooltip);
				else
				{
					s += writeNodeData(tarName, groupID);
					to.add(tarName);
				}
			}

			// write relations
			Set<String> genes = new HashSet<>(memNames);
			genes.addAll(tarNames);

			List<String> relations = network.linker.linkProgressive(genes, genes, 0);
			for (String rel : relations)
			{
				String[] tok = rel.split("\t");

				String srcID = tok[0] + groupID;
				String trgID = tok[2] + groupID;
				String edgeID = tok[0] + tok[2] + groupID;

				s+= "\ntype:edge\tid:" + edgeID + "\tsource:" + srcID +
					"\ttarget:" + trgID + "\tarrow:Target";

				if (tok[1].equals(SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag()))
				{
					s += "\tlinecolor:50 100 150";
				}
				else if (tok[1].equals(SIFEnum.CONTROLS_EXPRESSION_OF.getTag()))
				{
					s += "\tlinecolor:50 150 50\tstyle:Dashed";
				}
				else
				{
					s += "\tlinecolor:50 50 50";
				}

				writer.write(rel + "\n");
			}
		}

		writer.close();

		// write the graph to the stream
		writer = new OutputStreamWriter(new FileOutputStream(dir + "result-groups.cus"));
		writer.write(s);
		writer.close();

		writer = new OutputStreamWriter(new FileOutputStream(dir + "merged-network.format"));
		writer.write("node\tall-nodes\tcolor\t255 255 255\n");
		for (String id : node2color.keySet())
		{
			writer.write("node\t" + id + "\tcolor\t" + node2color.get(id) + "\n");
			writer.write("node\t" + id + "\ttooltip\t" + node2tooltip.get(id) + "\n");
		}
		for (String target : to)
		{
			writer.write("node\t" + target + "\tbordercolor\t255 255 255\n");
			writer.write("node\t" + target + "\ttextcolor\t155 155 155\n");
		}
		writer.close();

		// Write co-occurrence graph

		writer = new OutputStreamWriter(new FileOutputStream(dir + "co-occurred.sif"));

		for (String sym1 : genesInGroups)
		{
			for (String sym2 : genesInGroups)
			{
				if (sym1.compareTo(sym2) >= 0) continue;

				GeneAlt gene1 = geneMap.get(sym1);
				GeneAlt gene2 = geneMap.get(sym2);

				double pv = Overlap.calcCoocPval(gene1.getBooleanChanges(), gene2.getBooleanChanges());

				if (pv < 0.01)
				{
					writer.write(sym1 + "\tinteracts-with\t" + sym2 + "\n");
				}
			}
		}

		writer.close();
		writer = new OutputStreamWriter(new FileOutputStream(dir + "co-occurred.format"));
		writer.write(  "node\tall-nodes\tcolor\t255 255 255");
		writer.write("\nnode\tall-nodes\tbordercolor\t0 0 0");
		writer.close();
	}

	private static String writeNodeData(GeneAlt gene, String groupID,
		Map<String, Set<String>> name2id, boolean seed, Map<String, String> node2color,
		Map<String, String> node2tooltip)
	{
		String nodeID = gene.getId() + groupID;
		String s = "\ntype:node\tid:" + nodeID + "\ttext:" + gene.getId();
		double rat = gene.getAlteredRatio();
		int v = 255 - (int) (rat * 255);

//		boolean activating = gene.isActivated();
		boolean activating = true;

		String color = activating ?
			"255 " + v + " " + v : v + " " + v + " 255";
		s += "\tbgcolor:" + color + "\ttooltip:" + (int) Math.round(rat * 100) + "%";

		if (!node2color.containsKey(gene.getId()))
		{
			node2color.put(gene.getId(), color);
			node2tooltip.put(gene.getId(), FormatUtil.roundToSignificantDigits(rat, 2) + "");
		}

		if (seed)
		{
			s += "\tbordercolor:150 0 0";
		}

		if (!name2id.containsKey(gene.getId()))
			name2id.put(gene.getId(), new HashSet<>());

		name2id.get(gene.getId()).add(nodeID);
		return s;
	}

	private static String writeNodeData(String name, String groupID)
	{
		String nodeID = name + groupID;
		String s = "\ntype:node\tid:" + nodeID + "\ttext:" + name;

		s += "\tbgcolor:255 255 255";

		return s;
	}

	private static Set<String> getTargetOnly(List<Group> groups)
	{
		Set<String> only = new HashSet<>();
		Set<String> mems = new HashSet<>();
		for (Group group : groups)
		{
			for (GeneAlt member : group.members)
			{
				mems.add(member.getId());
			}
		}
		for (Group group : groups)
		{
			for (String target : group.targets)
			{
				if (!mems.contains(target)) only.add(target);
			}
		}
		return only;
	}
}
