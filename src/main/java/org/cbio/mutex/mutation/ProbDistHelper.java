package org.cbio.mutex.mutation;

import org.cbio.causality.idmapping.Length;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class ProbDistHelper
{
	Map<String, double[]> probMap = new HashMap<String, double[]>();
	Map<String, double[]> randMap = new HashMap<String, double[]>();
//	Map<String, double[]> drivMap = new HashMap<String, double[]>();
	Map<String, String[]> mutMap;
	Set<String>[] lostGenes;

	public ProbDistHelper(List<AlterationPack> packs, Map<String, String[]> mutMap, Set<String>[] lostGenes)
	{
		this.mutMap = mutMap;
		this.lostGenes = lostGenes;

		Map<Integer, Integer> counts = new HashMap<Integer, Integer>();

		for (int i = 0; i < packs.get(0).getSize(); i++)
		{
			int cnt = MutCount.getTotalMutationCount(mutMap, i, lostGenes[i]);
			counts.put(i, cnt);
		}

		int[] totalSeqLength = new int[packs.get(0).getSize()];
		for (int i = 0; i < totalSeqLength.length; i++)
		{
			Set<String> set = new HashSet<String>(mutMap.keySet());
			set.removeAll(lostGenes[i]);
			totalSeqLength[i] = Length.getTotalLength(set);
		}

		for (AlterationPack pack : packs)
		{
			int altCnt = pack.countAltered(Alteration.MUTATION);
			double[] pv = new double[pack.getSize()];

			for (int i = 0; i < pv.length; i++)
			{
				pv[i] = calcProbOfAtLeastOneHit(Length.of(pack.getId()), totalSeqLength[i], counts.get(i));
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

		if (lostGenes[pos].contains(symbol)) return 0;

		double[] pv = probMap.get(symbol);

		if (pv.length <= pos) throw new IllegalArgumentException("Index too large: " + pos);

		return pv[pos];
	}

	public double getRandProb(String symbol, int pos)
	{
		if (!randMap.containsKey(symbol))
			throw new IllegalArgumentException("Symbol not known: " + symbol);

		if (lostGenes[pos].contains(symbol)) return 0;

		double[] pv = randMap.get(symbol);

		if (pv.length <= pos) throw new IllegalArgumentException("Index too large: " + pos);

		return pv[pos];
	}
}
