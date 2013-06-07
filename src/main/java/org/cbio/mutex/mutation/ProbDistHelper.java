package org.cbio.mutex.mutation;

import org.cbio.causality.idmapping.Length;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Ozgun Babur
 */
public class ProbDistHelper
{
	Map<String, double[]> probMap = new HashMap<String, double[]>();
	Map<String, double[]> randMap = new HashMap<String, double[]>();
//	Map<String, double[]> drivMap = new HashMap<String, double[]>();
	Map<String, String[]> mutMap;

	public ProbDistHelper(List<AlterationPack> packs, Map<String, String[]> mutMap)
	{
		this.mutMap = mutMap;

		Map<Integer, Integer> counts = new HashMap<Integer, Integer>();

		for (int i = 0; i < packs.get(0).getSize(); i++)
		{
			int cnt = MutCount.getTotalMutationCount(mutMap, i);
			counts.put(i, cnt);
		}

		int totalSeqLength = Length.getTotalLength(mutMap.keySet());

		for (AlterationPack pack : packs)
		{
			int altCnt = pack.countAltered(Alteration.MUTATION);
			double[] pv = new double[pack.getSize()];

			for (int i = 0; i < pv.length; i++)
			{
				pv[i] = calcProbOfAtLeastOneHit(Length.of(pack.getId()), totalSeqLength, counts.get(i));
			}

			double[] randpv = new double[pv.length];
			System.arraycopy(pv, 0, randpv, 0, pv.length);
			randMap.put(pack.getId(), randpv);

			double total = 0;
			for (double v : pv) total += v;

			double diff = altCnt - total;

			if (diff > 0)
			{
				diff /= pv.length;
				for (int i = 0; i < pv.length; i++)
				{
					pv[i] += diff;
				}
			}
			else if (diff < 0)
			{
				double factor = altCnt / total;
				for (int i = 0; i < pv.length; i++)
				{
					pv[i] *= factor;
				}
			}

			probMap.put(pack.getId(), pv);
		}
	}

	public double calcProbOfAtLeastOneHit(int length, int totalLength, int times)
	{
		return 1 - Math.pow((totalLength - length) / (double) totalLength, times);
	}

	public double getProb(String symbol, int pos)
	{
		if (!probMap.containsKey(symbol))
			throw new IllegalArgumentException("Symbol not known: " + symbol);

		double[] pv = probMap.get(symbol);

		if (pv.length <= pos) throw new IllegalArgumentException("Index too large: " + pos);

		return pv[pos];
	}

	public double getRandProb(String symbol, int pos)
	{
		if (!randMap.containsKey(symbol))
			throw new IllegalArgumentException("Symbol not known: " + symbol);

		double[] pv = randMap.get(symbol);

		if (pv.length <= pos) throw new IllegalArgumentException("Index too large: " + pos);

		return pv[pos];
	}
}
