package org.cbio.mutex.mutation;

import org.cbio.causality.idmapping.Length;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.util.Histogram;
import org.cbio.causality.util.Histogram2D;
import org.cbio.causality.util.Progress;
import org.cbio.causality.util.Summary;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class PairwiseStat
{
	public static Histogram getOverlapDistribution(List<AlterationPack> packs, Map<String, String[]> mutMap,
		Set<String>[] lostGenes)
	{
		ProbDistHelper pdh = new ProbDistHelper(packs, mutMap, lostGenes);

		Histogram h = new Histogram(1);

		Progress p = new Progress(packs.size());
		double sum = 0;
		int cnt = 0;

		for (AlterationPack pack1 : packs)
		{
			p.tick();
			for (AlterationPack pack2 : packs)
			{
				if (pack1.getId().compareTo(pack2.getId()) >= 0) continue;

				double dif = OverlapFinder.getOverlapDif(pack1, pack2, pdh);
				h.count(dif);

				sum += dif;
				cnt++;
			}
		}

		System.out.println("average deviation of overlap =  " + (sum / cnt));

		return h;
	}

	public static Histogram getAlterationDistribution(List<AlterationPack> packs,
		Map<String, String[]> mutMap, Set<String>[] lostGenes)
	{
		ProbDistHelper pdh = new ProbDistHelper(packs, mutMap,lostGenes);

		double sum = 0;
		Histogram h = new Histogram(1);

		for (AlterationPack pack : packs)
		{
			double dif = OverlapFinder.getAltDif(pack, pdh, lostGenes);
			h.count(dif);

			sum += dif;

			if (dif < -10) System.out.println(pack.getId() + "\t" +
				Length.of(pack.getId()) + "\t" + pack.getAlteredCount(Alteration.MUTATION) +
				"\t" + Summary.sum(pdh.randMap.get(pack.getId())));
		}

		System.out.println("average = " + sum / packs.size());
		return h;
	}

	public static List<AlterationPack> filterToAlterationCounts(List<AlterationPack> packs,
		Map<String, String[]> mutMap, Set<String>[] lostGenes)
	{
		List<AlterationPack> list = new ArrayList<AlterationPack>();

		ProbDistHelper pdh = new ProbDistHelper(packs, mutMap, lostGenes);

		for (AlterationPack pack : packs)
		{
			double dif = OverlapFinder.getAltDif(pack, pdh, lostGenes);
			if (dif > 5) list.add(pack);
		}

		return list;
	}
}
