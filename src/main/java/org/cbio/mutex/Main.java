package org.cbio.mutex;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.cbio.causality.analysis.Graph;
import org.cbio.causality.analysis.SIFLinker;
import org.cbio.causality.data.GeneCards;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.network.PathwayCommons;
import org.cbio.causality.network.SPIKE;
import org.cbio.causality.network.SignaLink;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.*;

/**
 * This is the main executor class of the mutex search.
 * @author Ozgun Babur
 */
public class Main
{
	public static final PortalDataset data = PortalDataset.UCEC;
	public static final double MIN_ALTERATION_THR = 0.03;
	public static final double FDR_THR = 0.05;
	public static final int MAX_GROUP_SIZE = 7;
	public static final int RANDOMIZATION_TRIALS1 = 10000;
	public static final int RANDOMIZATION_TRIALS2 = 1000;

	public static final boolean REMOVE_HYPERMUTATED_OUTLIERS = true;
	public static final boolean LIMIT_TO_MUTSIG_AND_GISTIC_GENES = false;
	public static final boolean CROP_GRAPH_TO_MUTSIG_GISTIC_NEIGH = true;

	// For debug purposes
	public static final boolean LOAD_RESULTS_FROM_FILE = false;

	public static void main(String[] args) throws IOException, ClassNotFoundException
	{
		search();
//		checkAGroup();
	}

	public static void search() throws IOException, ClassNotFoundException
	{
		System.out.println("Dataset: " + data.study);

		// load the network

		Graph graphSig = loadPTRGraph();
		Graph graphTR = loadTRGraph();

		Graph graph = new Graph("directed network", "is-upstream-of");
		graph.merge(graphSig);
		graph.merge(graphTR);
		Set<String> symbols = decideSymbols(graph);

		// load the alteration data
		Map<String, AlterationPack> genesMap = readPortal(symbols);

		MutexGreedySearcher searcher = new MutexGreedySearcher(
			graph, genesMap, MIN_ALTERATION_THR, REMOVE_HYPERMUTATED_OUTLIERS);

		List<Group> groups;
		String groupsFilename = "result/" + data.study + ".groups";

		if (!LOAD_RESULTS_FROM_FILE)
		{
			// greedy search for the upstream of each gene
			groups = searcher.search(FDR_THR, MAX_GROUP_SIZE,
				RANDOMIZATION_TRIALS1, RANDOMIZATION_TRIALS2, data.study);

			// sort the list favoring high-coverage
			Group.sortToCoverage(groups);
			Group.serialize(groups, groupsFilename);

			searcher.addOmitted();
		}
		else
		{
			groups = Group.deserialize(groupsFilename);
		}

		SIFLinker linker = new SIFLinker();
		linker.load(graphSig);
		linker.load(graphTR);

		SubtypeAligner sa =
		data == PortalDataset.BRCA ?
			new SubtypeAligner(PortalDataset.BRCA_PUB, Group.collectGenes(groups)) :
			new SubtypeAligner(data, groups, searcher.getHyper());

		// Write the output graph to visualize in ChiBE
		GraphWriter.write(groups, searcher.getGenes(), linker, "result/", data.study, sa);

		// Print textual results
		for (Group group : groups)
		{
			System.out.println(group.getPrint(sa) + "\n");
		}

		System.out.println();
		printAnnotations(getGenes(groups, genesMap));
	}

	private static Set<String> decideSymbols(Graph graph)
	{
		// Filter gens to mutsig and gistic
		Set<String> symbols = graph.getSymbols();
		System.out.println("symbols initial size = " + symbols.size());

		if (LIMIT_TO_MUTSIG_AND_GISTIC_GENES)
		{
			GeneFilterer.filterToMutsigAndGistic(symbols, data, 0.05);
			System.out.println("symbols size filtered by mutsig and gistic = " + symbols.size());
		}

		if (CROP_GRAPH_TO_MUTSIG_GISTIC_NEIGH)
		{
			Set<String> syms = new HashSet<String>(symbols);
			GeneFilterer.filterToMutsigAndGistic(syms, data, 0.05);
			Set<String> neighbors = graph.getNeighbors(syms);
			Set<String> downstream = graph.getDownstream(syms);
			Set<String> upOfDown = graph.getUpstream(downstream);

			syms.addAll(neighbors);
			syms.addAll(upOfDown);

			graph.crop(syms);
			symbols = graph.getSymbols();
			System.out.println("symbols cropped size = " + symbols.size());
		}
		return symbols;
	}

	private static Graph loadTRGraph()
	{
		Graph graphTR = new Graph("transcriptional regulation", SIFEnum.CONTROLS_EXPRESSION_OF.getTag());
		graphTR.merge(PathwayCommons.getGraph(SIFEnum.CONTROLS_EXPRESSION_OF));
		graphTR.merge(SPIKE.getGraphTR());
		graphTR.merge(SignaLink.getGraphTR());
		return graphTR;
	}

	private static Graph loadPTRGraph()
	{
		Graph graphSig = new Graph("signaling", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		graphSig.merge(PathwayCommons.getGraph(SIFEnum.CONTROLS_STATE_CHANGE_OF));
		graphSig.merge(SPIKE.getGraphPostTl());
		graphSig.merge(SignaLink.getGraphPostTl());
		return graphSig;
	}

	private static Map<String, AlterationPack> readPortal(Set<String> symbols) throws IOException
	{
		// load the alteration data
		PortalReader reader = new PortalReader();
		Map<String, AlterationPack> genesMap = reader.readAlterations(data, symbols);
		for (String s : new HashSet<String>(genesMap.keySet()))
		{
			if (!symbols.contains(s)) genesMap.remove(s);
		}
		return genesMap;
	}


	public static void checkAGroup() throws IOException
	{
		Set<String> symbols = new HashSet<String>(Arrays.asList("ERBB2", "IGF1R", "PTEN", "PIK3CA", "AKT1"));

		// load the alteration data
		Map<String, AlterationPack> genesMap = readPortal(symbols);

		List<GeneAlt> genes = new ArrayList<GeneAlt>(genesMap.size());
		for (AlterationPack pack : genesMap.values())
		{
			genes.add(new GeneAlt(pack, Alteration.GENOMIC));
		}

		Group group = new Group(genes.get(0));
		for (int i = 1; i < genes.size(); i++)
		{
			group.addGene(genes.get(i));
		}

		System.out.println("[" + group.getGeneNamesInString() + "]\tcover: " +
			group.calcCoverage() + "\tpval: " + group.calcScore() + "\ttargets:" +
			group.getTargets());

		System.out.println(group.getPrint());
	}

	private static List<String> getGenes(List<Group> groups,
		final Map<String, AlterationPack> genesMap)
	{
		Set<String> genes = new HashSet<String>();
		for (Group group : groups)
		{
			genes.addAll(group.getGeneNames());
		}
		List<String> list = new ArrayList<String>(genes);

		Collections.sort(list, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return new Integer(genesMap.get(o2).getAlteredCount(Alteration.GENOMIC)).compareTo(
					genesMap.get(o1).getAlteredCount(Alteration.GENOMIC));
			}
		});

		return list;
	}

	private static void printAnnotations(List<String> genes)
	{
		System.out.println("Genes in groups: " + genes.size());
		System.out.println(genes + "\n");

		String keyword = data.name;

		List<String> withKw = new ArrayList<String>();
		List<String> withOtherCancer = new ArrayList<String>(genes);
		List<String> notAssoc = new ArrayList<String>();

		for (String gene : genes)
		{
			System.out.print(gene + ": ");
			Set<String> relatedCancers = GeneCards.getRelatedCancers(gene);
			for (String cancer : relatedCancers)
			{
				System.out.print(cancer + ", ");
				if (cancer.contains(keyword) && !withKw.contains(gene)) withKw.add(gene);
			}
			System.out.println();
			if (relatedCancers.isEmpty()) notAssoc.add(gene);
		}
		withOtherCancer.removeAll(withKw);
		withOtherCancer.removeAll(notAssoc);

		System.out.println("\nAlready associated with " + keyword + " cancer: " + withKw.size());
		System.out.println(withKw);
		System.out.println("\nAssociated with other types of cancer: " + withOtherCancer.size());
		System.out.println(withOtherCancer);
		System.out.println("\nNot associated with any cancer: " + notAssoc.size());
		System.out.println(notAssoc);
	}
}
