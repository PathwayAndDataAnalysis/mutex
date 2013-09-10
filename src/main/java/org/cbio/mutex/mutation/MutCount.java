package org.cbio.mutex.mutation;

import org.cbio.causality.data.portal.CBioPortalAccessor;
import org.cbio.causality.idmapping.Length;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.util.Summary;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class MutCount
{
	public static Map<String, String[]> prepareMutMap(List<AlterationPack> packs,
		CBioPortalAccessor cacc)
	{
		Map<String, String[]> mutMap = new HashMap<String, String[]>();

		for (AlterationPack pack : packs)
		{
			String[] muts = cacc.getManager().getDataForGene(pack.getId(),
				cacc.getGeneticProfileContainingName("mutations"), cacc.getCurrentCaseList());
			mutMap.put(pack.getId(), muts);
		}
		return mutMap;
	}

	public static int getTotalMutationCount(Map<String, String[]> mutMap, int sample)
	{
		int total = 0;
		for (String[] muts : mutMap.values())
		{
			total += getMutationCount(muts[sample]);
		}
		return total;
	}

	public static int getTotalMutationCount(Map<String, String[]> mutMap, int sample, Set<String> ignore)
	{
		int total = 0;
		for (String symbol : mutMap.keySet())
		{
			if (ignore.contains(symbol)) continue;

			String[] muts = mutMap.get(symbol);
			total += getMutationCount(muts[sample]);

		}
		return total;
	}

	public static int[] getSampleMutCounts(Map<String, String[]> mutMap)
	{
		int[] cnts = new int[mutMap.values().iterator().next().length];

		for (int i = 0; i < cnts.length; i++)
		{
			for (String symbol : mutMap.keySet())
			{
				cnts[i] += getMutationCount(symbol, mutMap, i);
			}
		}
		return cnts;
	}
	
	public static int getMutationCount(String symbol, Map<String, String[]> mutMap, int sample)
	{
		return getMutationCount(mutMap.get(symbol)[sample]);
	}

	public static int getMutationCount(String mutStr)
	{
		if (isNaN(mutStr)) return 0;

		return mutStr.split(",").length;
	}

	final static String[] NULLS = new String[]{"", "NaN", "NA", "null"};

	public static boolean isNaN(String s)
	{
		if (s == null) return true;
		s = s.trim();
		for (String val : NULLS)
		{
			if (s.equalsIgnoreCase(val)) return true;
		}
		return false;
	}

	public static void countAvgLengths(Map<String, String[]> mutMap)
	{
		Map<Integer, List<Integer>> map = new HashMap<Integer, List<Integer>>();

		for (String symbol : mutMap.keySet())
		{
			String[] muts = mutMap.get(symbol);

			for (String mut : muts)
			{
				int cnt = getMutationCount(mut);

				if (!map.containsKey(cnt)) map.put(cnt, new ArrayList<Integer>());
				map.get(cnt).add(Length.of(symbol));
			}
		}

		for (Integer len : map.keySet())
		{
			System.out.println(len + "\t" + Summary.mean(map.get(len)) + "\t" + map.get(len).size());
		}
	}
}
