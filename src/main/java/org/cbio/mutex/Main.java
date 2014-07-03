package org.cbio.mutex;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.cbio.causality.analysis.Graph;
import org.cbio.causality.analysis.SIFLinker;
import org.cbio.causality.data.GeneCards;
import org.cbio.causality.data.portal.BroadAccessor;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;
import org.cbio.causality.network.PathwayCommons;
import org.cbio.causality.network.SPIKE;
import org.cbio.causality.network.SignaLink;
import org.cbio.causality.util.Histogram;
import org.cbio.causality.util.Overlap;
import org.cbio.causality.util.Summary;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

/**
 * This is the main executor class of the mutex search.
 * @author Ozgun Babur
 */
public class Main
{
	public static final PortalDataset data = PortalDataset.GBM;
	public static final double MIN_ALTERATION_THR = data.minAltThr;
	public static final double FDR_THR = 0.05;
	public static final int MAX_GROUP_SIZE = 10;
	public static final int RANDOMIZATION_TRIALS1 = 10000;
	public static final int RANDOMIZATION_TRIALS2 = 1000;

	public static final boolean REMOVE_HYPERMUTATED_OUTLIERS = true;
	public static final boolean LIMIT_TO_MUTSIG_AND_GISTIC_GENES = false;
	public static final boolean CROP_GRAPH_TO_MUTSIG_GISTIC_NEIGH = true;

	// For debug purposes
	public static final boolean LOAD_RESULTS_FROM_FILE = false;
	public static final boolean JUST_WRITE_THE_INPUT_DATA = false;

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

		Graph graph = getGraph(graphSig, graphTR);
		Set<String> symbols = decideSymbols(graph);

		// load the alteration data
		Map<String, AlterationPack> genesMap;
		if (data.name.equals("simulated"))
		{
			genesMap = Simulation.loadSimData(symbols, data);
		}
		else
		{
			genesMap = readPortal(symbols);
		}

		String groupsFilename = "result/" + data.study + ".groups";

		if (data == PortalDataset.SKCM ||
			data == PortalDataset.LUAD ||
//			data == PortalDataset.OV ||
			data == PortalDataset.LUSC)
			applyMutsigThreshold(genesMap, 0.05);

//		printDatasetCharacteristics(genesMap);

		MutexGreedySearcher searcher = LOAD_RESULTS_FROM_FILE ?
			MutexGreedySearcher.deserialize(groupsFilename) :
			new MutexGreedySearcher(graph, genesMap, MIN_ALTERATION_THR, REMOVE_HYPERMUTATED_OUTLIERS);

		if (JUST_WRITE_THE_INPUT_DATA)
		{
			searcher.serialize("cache-" + data.study);
			return;
		}

		if (false)
		{
//			MutexSimpleRanker msr = new MutexSimpleRanker(genesMap);
//			msr.search("ranks-v1.txt", Simulation.extractGroupSize(data), Simulation.getLabelMap(data, genesMap));

			PairSearcher ps = new PairSearcher(genesMap);
			ps.search("pairs-v3.txt", Simulation.getLabelMap(data, genesMap));

//			searcher.writeRankedGroups(MAX_GROUP_SIZE, RANDOMIZATION_TRIALS1, Simulation.getLabelMap(data, genesMap), "mutex-ranked-groups-v3-nograph.txt");
			return;
		}

		List<Group> groups;

		// greedy search for the upstream of each gene
		groups = searcher.search(FDR_THR, MAX_GROUP_SIZE, RANDOMIZATION_TRIALS1,
			RANDOMIZATION_TRIALS2, data.study, LOAD_RESULTS_FROM_FILE ? null : groupsFilename);

		// sort the list favoring high-coverage
		Group.sortToCoverage(groups);
		searcher.addOmitted();

//		groups = Simulation.getTrueGroups(data, searcher);

		SIFLinker linker = new SIFLinker();
		linker.load(graphSig);
		linker.load(graphTR);

		SubtypeAligner sa =
		data == PortalDataset.BRCA ? new SubtypeAligner(PortalDataset.BRCA_PUB, Group.collectGenes(groups)) :
		data == PortalDataset.GBM ? new SubtypeAligner(PortalDataset.GBM_PUB, Group.collectGenes(groups)) :
		data.name.equals("simulated") ? null : new SubtypeAligner(data, groups, searcher.getHyper());

		// Write the output graph to visualize in ChiBE
		GraphWriter.write(groups, searcher.getGenes(), linker, "result/", data.study, sa);
		Map<String, String> nameConvMap = data.name.equals("simulated") ?
			Simulation.getLabelMap(data, genesMap) : null;

		// Print textual results
		for (Group group : groups)
		{
			System.out.println(group.getPrint(sa, nameConvMap, true) + "\n");
		}

		System.out.println();
		printOverlappingCNA(groups, genesMap);
		System.out.println();

		if (data.name.equals("simulated")) Simulation.evaluateSuccess(groups, data, genesMap);
		else printAnnotations(getGenes(groups, genesMap));

//		Oncoprint onco = new Oncoprint(groups, data);
//		onco.write();
	}

	private static Set<String>  decideSymbols(Graph graph) throws FileNotFoundException
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

			if (data == PortalDataset.SIMUL1)
			{
//				syms = SimulatedDataGenerator.getMutsig();
			}
			else GeneFilterer.filterToMutsigAndGistic(syms, data, 0.05);

			Set<String> neighbors = graph.getNeighbors(syms);
			Set<String> downstream = graph.getDownstream(syms);
			Set<String> upOfDown = graph.getUpstream(downstream);

			syms.addAll(neighbors);
			syms.addAll(upOfDown);

			graph.crop(syms);
			symbols = graph.getSymbols();
			System.out.println("symbols cropped size = " + symbols.size());
		}

//		if (data == PortalDataset.SKCM)
//		{
//			GeneFilterer.filterToMutsigAndGistic(symbols, data, 0.05);
//			System.out.println("symbols size filtered by mutsig and gistic = " + symbols.size());
//		}

		return symbols;
	}

	public static Graph getGraph()
	{
		return getGraph(loadPTRGraph(), loadTRGraph());
	}

	public static Graph getGraph(Graph graphSig, Graph graphTR)
	{
		Graph graph = new Graph("directed network", "is-upstream-of");
		graph.merge(graphSig);
		graph.merge(graphTR);
		return graph;
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

	public static void printDatasetCharacteristics(Map<String, AlterationPack> genesMap)
	{
		int sampleSize = genesMap.values().iterator().next().getSize();
		int geneSize = genesMap.size();
		double log2 = Math.log(2);

		int[] cnt = new int[sampleSize];
		for (AlterationPack pack : genesMap.values())
		{
			for (int i = 0; i < sampleSize; i++)
			{
				if (pack.get(Alteration.GENOMIC)[i].isAltered()) cnt[i]++;
			}
		}

		for (int i = 0; i < cnt.length; i++)
		{
			if (cnt[i] == 0) cnt[i]++;
		}

		double[] sampleCov = new double[sampleSize];
		List<Double> geneCov = new ArrayList<Double>();

		double range = 0.5;
		Histogram h1 = new Histogram(range);
		h1.setBorderAtZero(true);
		for (int i = 0; i < sampleSize; i++)
		{
			double rat = Math.log(cnt[i] / (double) geneSize) / log2;
			h1.count(rat);
			sampleCov[i] = rat;
		}
		System.out.println("Sample coverage:");
		h1.print();

		System.out.println("Summary.mean(sampleCov) = " + Summary.mean(sampleCov));
		System.out.println("Summary.stdev(sampleCov) = " + Summary.stdev(sampleCov));

		Histogram h2 = new Histogram(range);
		h2.setBorderAtZero(true);
		for (AlterationPack pack : genesMap.values())
		{
			if (!pack.isAltered()) continue;
			double rat = Math.log(pack.getAlteredRatio(Alteration.GENOMIC)) / log2;
			h2.count(rat);
			geneCov.add(rat);
		}
		System.out.println("Gene coverage:");
		h2.print();

		double[] v = new double[geneCov.size()];
		for (int i = 0; i < v.length; i++)
		{
			v[i] = geneCov.get(i);
		}
		System.out.println("Summary.mean(geneCov) = " + Summary.mean(v));
		System.out.println("Summary.stdev(geneCov) = " + Summary.stdev(v));

		System.exit(0);
	}

	private static void applyMutsigThreshold(Map<String, AlterationPack> genes, double pvalThr)
	{
		String study = data.caseList.substring(0, data.caseList.indexOf("_")).toUpperCase();
		Set<String> mutsig = BroadAccessor.getMutsigGenes(study, pvalThr, false);
		for (AlterationPack pack : genes.values())
		{
			if (!mutsig.contains(pack.getId()))
			{
				for (int i = 0; i < pack.getSize(); i++)
				{
					pack.get(Alteration.MUTATION)[i] = Change.NO_CHANGE;
				}

				pack.complete(Alteration.GENOMIC);
			}
		}
	}

	private static void printOverlappingCNA(List<Group> groups, Map<String, AlterationPack> packMap)
	{
		List<Set<String>> list = new ArrayList<Set<String>>();
		Set<String> genes = new HashSet<String>();
		for (Group group : groups)
		{
			genes.addAll(group.getGeneNames());
		}
		for (String gene : genes)
		{
			list.add(new HashSet<String>(Arrays.asList(gene)));
		}

		for (String gene : genes)
		{
			AlterationPack pack1 = packMap.get(gene);
			if (pack1.get(Alteration.COPY_NUMBER) == null) continue;

			for (Set<String> group : list)
			{
				if (!group.contains(gene))
				{
					boolean allFit = true;
					for (String g2 : group)
					{
						AlterationPack pack2 = packMap.get(g2);
						if (pack2.get(Alteration.COPY_NUMBER) == null ||
							Overlap.calcCoocPval(pack1.get(Alteration.COPY_NUMBER),
								pack2.get(Alteration.COPY_NUMBER)) > 0.001)
						{
							allFit = false;
							break;
						}
					}
					if (allFit) group.add(gene);
				}
			}
		}

		Set<Set<String>> ss = new HashSet<Set<String>>(list);
		System.out.println("Cooccurred CNA groups:");
		for (Set<String> group : ss)
		{
			if (group.size() > 1) System.out.println(group);
		}
	}
}
