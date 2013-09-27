package org.cbio.mutex;

import org.cbio.causality.analysis.SIFLinker;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.*;

/**
 * This is the main executor class of the mutex search.
 * @author Ozgun Babur
 */
public class Main
{
	public static final PortalDataset data = PortalDataset.endometrial;
	public static final double MIN_ALTERATION_THR = 0.07;
	public static final double FDR_THR = 0.05;
	public static final int MAX_GROUP_SIZE = 5;
	public static final int RANDOMIZATION_TRIALS = 10;

	public static final boolean REMOVE_HYPERMUTATED_OUTLIERS = true;
	public static final boolean LIMIT_TO_MUTSIG_AND_GISTIC_GENES = false;
	public static final boolean USE_GISTIC_GROUPING = true;

	public static void main(String[] args) throws IOException
	{
		// load the network

		SIFLinker linker = new SIFLinker();
		linker.load(Main.class.getResourceAsStream("PC.txt"),
			"controls-state-change", "controls-expression", "controls-degradation");
		linker.load(Main.class.getResourceAsStream("SPIKE.txt"), "is-upstream-of");

		// Filter gens to mutsig and gistic
		Set<String> symbols = linker.traverse.getSymbols();
		System.out.println("symbols initial size = " + symbols.size());

		if (LIMIT_TO_MUTSIG_AND_GISTIC_GENES)
		{
			GeneFilterer.filterToMutsigAndGistic(symbols, data);
			System.out.println("symbols size filtered by mutsig and gistic= " + symbols.size());
		}

		Map<String, Set<String>> gistic = USE_GISTIC_GROUPING ?
			GeneFilterer.getGisticSets(data) : null;

		// load the alteration data
		PortalReader reader = new PortalReader();
		final Map<String,AlterationPack> genesMap = reader.readAlterations(data, symbols);
		for (String s : new HashSet<String>(genesMap.keySet()))
		{
			if (!symbols.contains(s)) genesMap.remove(s);
		}

		// greedy search for the upstream of each gene

		MutexGreedySearcher searcher = new MutexGreedySearcher(
			linker.traverse, genesMap, MIN_ALTERATION_THR, gistic, REMOVE_HYPERMUTATED_OUTLIERS);

		List<Group> groups = searcher.search(FDR_THR, MAX_GROUP_SIZE, RANDOMIZATION_TRIALS);
		double thr = getWorstPval(groups);

//		// pair search
//		MutexPairSearcher searcher = new MutexPairSearcher(genesMap, MIN_ALTERATION_THR);
//		List<Group> groups = searcher.search(FDR_THR, RANDOMIZATION_TRIALS);

		// clean subsets in the result
		Group.removeSubsets(groups);

		// sort the list favoring high-coverage
		Group.sortToCoverage(groups);

		// Write the output graph to visualize in ChiBE
		GraphWriter.write(groups, 0.01, linker,
			new FileOutputStream("/home/ozgun/Desktop/mutex/" + data.name + ".cus"),
//			new FileOutputStream("C:\\Users\\ozgun\\Desktop\\mutex\\" + data.name + ".cus"),
			data.name, false);

		// Print textual results
		for (Group group : groups)
		{
			System.out.print("[" + group.getGeneNamesInString() + "]\tcover: " + group.calcCoverage() +
				"\tpval: " + group.calcOverallPVal() + "\ttargets:" + group.getTargets());

			System.out.println();

			System.out.println(group.getPrint(thr));
			System.out.println();
		}
	}

	private static double getWorstPval(List<Group> groups)
	{
		double worst = 0;

		for (Group group : groups)
		{
			double pv = group.calcOverallPVal();
			if (pv > worst)
			{
				worst = pv;
			}
		}
		return worst;
	}
}
