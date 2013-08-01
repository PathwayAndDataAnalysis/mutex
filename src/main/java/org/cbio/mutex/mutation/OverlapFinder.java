package org.cbio.mutex.mutation;

import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;

import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class OverlapFinder
{
	public static double getOverlapDif(AlterationPack pack1, AlterationPack pack2,
		ProbDistHelper pdh)
	{
		double dif = 0;

		for (int i = 0; i < pack1.getSize(); i++)
		{
			double pv1 = pdh.getProb(pack1.getId(), i);
			double pv2 = pdh.getProb(pack2.getId(), i);

			double e = - (pv1 * pv2);

			if (pack1.getChange(Alteration.MUTATION, i).isAltered() &&
				pack2.getChange(Alteration.MUTATION, i).isAltered())
			{
				e += 1;
			}

			dif += e;
		}

		return dif;
	}

	public static double getAltDif(AlterationPack pack, ProbDistHelper pdh, Set<String>[] lost)
	{

		int alteredCount = 0;
		double dif = 0;

		for (int i = 0; i < pack.getSize(); i++)
		{
			if (lost[i].contains(pack.getId())) continue;

			dif += pdh.getRandProb(pack.getId(), i);
			alteredCount += pack.getChange(Alteration.MUTATION, i).isAltered() ? 1 : 0;
		}

		return alteredCount - dif;
	}
}
