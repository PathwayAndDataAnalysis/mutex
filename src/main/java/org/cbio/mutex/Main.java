package org.cbio.mutex;

import org.cbio.causality.analysis.SIFLinker;
import org.cbio.causality.model.AlterationPack;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.List;
import java.util.Map;

/**
 * This is the executor class of the search.
 * @author Ozgun Babur
 */
public class Main
{
	public static void main(String[] args) throws IOException
	{
		// decide parameters of current analysis
		PortalDataset data = PortalDataset.colon;
		String outFile = "/home/ozgun/Desktop/mutex/greedy/" + data.name + ".cus";
		double minAlt = 0.03;
		double pvalThr = 0.05;

		// load the network

		SIFLinker linker = new SIFLinker();
		linker.load(Main.class.getResourceAsStream("network.txt"),
			"STATE_CHANGE", "TRANSCRIPTION", "DEGRADATION");

		// load the alteration data
		PortalReader reader = new PortalReader();
		Map<String,AlterationPack> genesMap = reader.readAlterations(
			data, linker.traverse.getSymbols());

		// greedy search for the upstream of each gene
		MutexSearcher searcher = new MutexSearcher(linker.traverse, genesMap, minAlt);
		List<Group> groups = searcher.search(pvalThr, minAlt, 5);

		// clean subsets in the result
		Group.removeSubsets(groups);

		// sort the list favoring high-coverage
		Group.sortToCoverage(groups);

		GraphWriter.write(groups, linker, new FileOutputStream(outFile), data.name);

		// Print oncoprints for the groups
		for (Group group : groups)
		{
			System.out.println(group.getGeneNames() + "\t" + group.calcCoverage());
			System.out.println(group.getPrint());
			System.out.println();
		}
	}
}
