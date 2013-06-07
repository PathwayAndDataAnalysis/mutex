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

/**
 * @author Ozgun Babur
 */
public class PairwiseStat
{
	public static void printOverlapStats(List<AlterationPack> packs, Map<String, String[]> mutMap)
	{
		ProbDistHelper pdh = new ProbDistHelper(packs, mutMap);

		int cntCooc = 0;
		int cntMutex = 0;

		Histogram h = new Histogram(1);

		Progress p = new Progress(packs.size());

		for (AlterationPack pack1 : packs)
		{
			p.tick();
			for (AlterationPack pack2 : packs)
			{
				if (pack1.getId().compareTo(pack2.getId()) >= 0) continue;

				double dif = OverlapFinder.getOverlapDif(pack1, pack2, pdh);
				h.count(dif);

				if (Math.abs(dif) > 5)
				{
					if (dif > 0) cntCooc++;
					else cntMutex++;
				}
			}
		}

		h.print();

		System.out.println("cntCooc = " + cntCooc);
		System.out.println("cntMutex = " + cntMutex);
	}

	public static List<AlterationPack> filterToAlterationCounts(List<AlterationPack> packs,
		Map<String, String[]> mutMap)
	{
		List<AlterationPack> list = new ArrayList<AlterationPack>();

		ProbDistHelper pdh = new ProbDistHelper(packs, mutMap);

//		double sum = 0;
//		Histogram h = new Histogram(1);
		Histogram2D h = new Histogram2D(1);

		for (AlterationPack pack : packs)
		{
			double dif = OverlapFinder.getAltDif(pack, pdh);
			double abs = OverlapFinder.getAltAbsDif(pack, pdh);
			h.count(dif, abs);

//			h.count(dif);
//			sum += dif;

			if (dif > 5) list.add(pack);

//			if (dif < -30) System.out.println("pack.getId() = " + pack.getId() + "\t" +
//				Length.of(pack.getId()) + "\t" + pack.getAlteredCount(Alteration.MUTATION) +
//				"\t" + Summary.sum(pdh.randMap.get(pack.getId())));
		}

//		System.out.println("average = " + sum / packs.size());
		h.plot();
		return list;
	}
}
