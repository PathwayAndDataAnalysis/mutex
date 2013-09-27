package org.cbio.mutex;

import org.cbio.causality.data.portal.BroadAccessor;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class GeneFilterer
{
	public static void filterToMutsigAndGistic(Set<String> symbols, PortalDataset dataset)
	{
		String study = dataset.caseList.substring(0, dataset.caseList.indexOf("_")).toUpperCase();
		Set<String> genes = BroadAccessor.getMutsigGenes(study, 0.5);
		genes.addAll(BroadAccessor.getGisticGenes(study, 0.5));
		symbols.retainAll(genes);

	}

	public static Map<String, Set<String>> getGisticSets(PortalDataset dataset)
	{
		String study = dataset.caseList.substring(0, dataset.caseList.indexOf("_")).toUpperCase();
		List<Set<String>> sets = BroadAccessor.getGisticGeneSets(study, 0.5);

		Set<String> mutsig = getMutsig(dataset);

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

	public static Set<String> getMutsig(PortalDataset dataset)
	{
		String study = dataset.caseList.substring(0, dataset.caseList.indexOf("_")).toUpperCase();
		Set<String> genes = BroadAccessor.getMutsigGenes(study, 0.5);
		return genes;
	}
}
