package org.cbio.mutex.mutation;

import org.cbio.causality.data.portal.CBioPortalAccessor;
import org.cbio.causality.data.portal.CBioPortalManager;
import org.cbio.causality.data.portal.GeneticProfile;

import java.util.Collection;

/**
 * @author Ozgun Babur
 */
public class MutAndCopyNumberRelater
{
	public static void printPlot(Collection<String> symbols, CBioPortalAccessor acc)
	{
		CBioPortalManager cman = new CBioPortalManager();
		GeneticProfile mutGP = acc.getGeneticProfileContainingName("mutation");
		GeneticProfile cnGP = acc.getGeneticProfileContainingName("copy-number");

		int[] sum = new int[5];
		int[] size = new int[5];

		for (String symbol : symbols)
		{
			String[] muts = cman.getDataForGene(symbol, mutGP, acc.getCurrentCaseList());
			String[] cns = cman.getDataForGene(symbol, cnGP, acc.getCurrentCaseList());

			for (int i = 0; i < muts.length; i++)
			{
				if (MutCount.isNaN(cns[i])) continue;

				int j = cns[i].equals("-2") ? 0 :
					cns[i].equals("-1") ? 1 :
						cns[i].equals("0") ? 2 :
							cns[i].equals("1") ? 3 :
								cns[i].equals("2") ? 4 : -1;


				if (!MutCount.isNaN(muts[i]))
				{
					sum[j] += 1;
					if (j == 0) System.out.println("muts[i] = " + muts[i]);
				}
				size[j] += 1;
			}
		}

		for (int i = 0; i < sum.length; i++)
		{
			System.out.println((i - 2) + "\t" + sum[i] + "\t" + size[i] + "\t" + (sum[i] / (double) size[i]));
		}
	}
}
