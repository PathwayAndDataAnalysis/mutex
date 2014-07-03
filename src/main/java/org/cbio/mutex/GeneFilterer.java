package org.cbio.mutex;

import org.cbio.causality.data.portal.BroadAccessor;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class GeneFilterer
{
	public static void filterToMutsigAndGistic(Set<String> symbols, PortalDataset dataset,
		double fdrThr)
	{
		String study = dataset.caseList.substring(0, dataset.caseList.indexOf("_")).toUpperCase();
		Set<String> genes = BroadAccessor.getMutsigGenes(study, fdrThr, true);
		genes.addAll(BroadAccessor.getExpressionVerifiedGistic(study, fdrThr));

//		Set<String> check = new HashSet<String>(Arrays.asList("RPS6KB1 TLK2".split(", ")));
//		check.removeAll(genes);
//		System.out.println(check);
//		if (true) System.exit(0);

		symbols.retainAll(genes);
	}

	public static Map<String, Set<String>> getGisticSets(PortalDataset dataset, double fdrThr)
	{
		String study = dataset.caseList.substring(0, dataset.caseList.indexOf("_")).toUpperCase();
		List<Set<String>> sets = BroadAccessor.getGisticGeneSets(study, 0.5);

		Set<String> mutsig = getMutsig(dataset, fdrThr);

		Map<String, Set<String>> map = new HashMap<String, Set<String>>();

		for (Set<String> set : sets)
		{
			set.removeAll(mutsig);
			for (String s : set)
			{
				map.put(s, set);
			}
		}

		return map;
	}

	public static Set<String> getMutsig(PortalDataset dataset, double fdrThr)
	{
		String study = dataset.caseList.substring(0, dataset.caseList.indexOf("_")).toUpperCase();
		Set<String> genes = BroadAccessor.getMutsigGenes(study, fdrThr, true);
		return genes;
	}

	public static List<String> getGisticExpVerified(PortalDataset dataset, double fdrThr)
	{
		String study = dataset.caseList.substring(0, dataset.caseList.indexOf("_")).toUpperCase();
		List<String> genes = BroadAccessor.getExpressionVerifiedGistic(study, fdrThr);
		return genes;
	}
}
