package org.cbio.mutex;

import org.cbio.causality.analysis.SIFLinker;
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
	public static final PortalDataset data = PortalDataset.colon;
	public static final double MIN_ALTERATION_THR = 0.03;
	public static final double PVAL_THR = 0.05;
	public static final double FDR_THR = 0.05;
	public static final int MAX_GROUP_SIZE = 6;
	public static final int RANDOMIZATION_TRIALS = 10;

	public static void main(String[] args) throws IOException
	{
		// load the network

		SIFLinker linker = new SIFLinker();
		linker.load(Main.class.getResourceAsStream("PC.txt"),
			"controls-state-change", "controls-expression", "controls-degradation");
		linker.load(Main.class.getResourceAsStream("SPIKE.txt"), "is-upstream-of");

		// load the alteration data
		PortalReader reader = new PortalReader();
		Map<String,AlterationPack> genesMap = reader.readAlterations(
			data, linker.traverse.getSymbols());

		// greedy search for the upstream of each gene

		MutexGreedySearcher searcher = new MutexGreedySearcher(
			linker.traverse, genesMap, MIN_ALTERATION_THR);

		List<Group> groups = searcher.search(
			PVAL_THR, FDR_THR, MAX_GROUP_SIZE, RANDOMIZATION_TRIALS);

		// clean subsets in the result
		Group.removeSubsets(groups);

		// sort the list favoring high-coverage
		Group.sortToCoverage(groups);

		// Write the output graph to visualize in ChiBE
		GraphWriter.write(groups, linker,
			new FileOutputStream("/home/ozgun/Desktop/mutex/greedy/" + data.name + ".cus"),
			data.name);

		// Print textual results
		for (Group group : groups)
		{
			List<String> dwstr = new ArrayList<String>(linker.traverse.getLinkedCommonDownstream(
				new HashSet<String>(group.getGeneNames())));

			Collections.sort(dwstr);

			System.out.print(group.getGeneNames() + "\tcover: " + group.calcCoverage() +
				"\ttargets:");

			for (String s : dwstr)
			{
				System.out.print(" " + s);
			}
			System.out.println();

			System.out.println(group.getPrint());
			System.out.println();
		}
	}
}
