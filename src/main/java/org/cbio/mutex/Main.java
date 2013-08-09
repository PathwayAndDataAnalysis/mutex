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
	public static final double MIN_ALTERATION_THR = 0.03;
	public static final double FDR_THR = 0.05;
	public static final int MAX_GROUP_SIZE = 10;
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
		final Map<String,AlterationPack> genesMap = reader.readAlterations(
			data, linker.traverse.getSymbols());

		// greedy search for the upstream of each gene

		MutexGreedySearcher searcher = new MutexGreedySearcher(
			linker.traverse, genesMap, MIN_ALTERATION_THR);

		List<Group> groups = searcher.search(FDR_THR, MAX_GROUP_SIZE, RANDOMIZATION_TRIALS);

//		// pair search
//		MutexPairSearcher searcher = new MutexPairSearcher(genesMap, MIN_ALTERATION_THR);
//		List<Group> groups = searcher.search(FDR_THR, RANDOMIZATION_TRIALS);

		// clean subsets in the result
		Group.removeSubsets(groups);

		// sort the list favoring high-coverage
		Group.sortToCoverage(groups);

		// Write the output graph to visualize in ChiBE
		GraphWriter.write(groups, 0.001, linker,
			new FileOutputStream("/home/ozgun/Desktop/mutex/greedy/" + data.name + ".cus"),
			data.name, true);

		// Print textual results
		for (Group group : groups)
		{
			List<String> dwstr = new ArrayList<String>(linker.traverse.getLinkedCommonDownstream(
				new HashSet<String>(group.getGeneNames())));

			Collections.sort(dwstr, new Comparator<String>()
			{
				@Override
				public int compare(String o1, String o2)
				{
					AlterationPack pack1 = genesMap.get(o1);
					AlterationPack pack2 = genesMap.get(o2);

					Integer cnt1 = pack1 == null ? 0 : pack1.getAlteredCount(Alteration.GENOMIC);
					Integer cnt2 = pack2 == null ? 0 : pack2.getAlteredCount(Alteration.GENOMIC);

					if (cnt1.equals(cnt2)) return o1.compareTo(o2);
					else return cnt2.compareTo(cnt1);
				}
			});

			System.out.print(group.getGeneNames() + "\tcover: " + group.calcCoverage() +
				"\tpval: " + group.calcOverallPVal() + "\ttargets:");

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
