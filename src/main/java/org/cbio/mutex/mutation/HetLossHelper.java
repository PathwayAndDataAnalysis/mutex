package org.cbio.mutex.mutation;

import org.cbio.causality.data.portal.CBioPortalAccessor;
import org.cbio.causality.data.portal.CBioPortalManager;
import org.cbio.causality.data.portal.CaseList;
import org.cbio.causality.data.portal.GeneticProfile;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class HetLossHelper
{
	public static Set<String>[] collectNonExistingGenes(Collection<String> symbols,
		CBioPortalAccessor acc)
	{
		CBioPortalManager cman = new CBioPortalManager();
		Set<String>[] sets = null;

		for (String symbol : symbols)
		{
			String[] data = cman.getDataForGene(symbol,
				acc.getGeneticProfileContainingName("copy-number"), acc.getCurrentCaseList());

			for (int i = 0; i < data.length; i++)
			{
				if (sets == null) sets = new Set[data.length];
				if (sets[i] == null) sets[i] = new HashSet<String>();

				double val = MutCount.isNaN(data[i]) ? 0 : Double.parseDouble(data[i]);
				if (val != 0)
				{
					sets[i].add(symbol);
				}
			}
		}
		return sets;
	}
}
