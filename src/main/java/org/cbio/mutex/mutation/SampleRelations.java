package org.cbio.mutex.mutation;

import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class SampleRelations
{
	public static Map<Integer, Set<String>> getMutatedPerSample(List<AlterationPack> packs)
	{
		Map<Integer, Set<String>> map = new HashMap<Integer, Set<String>>();
		for (int i = 0; i < packs.get(0).getSize(); i++)
		{
			Set<String> set = new HashSet<String>();

			for (AlterationPack pack : packs)
			{
				if (pack.getChange(Alteration.MUTATION, i).isAltered())
				{
					set.add(pack.getId());
				}
			}

			map.put(i, set);
		}
		return map;
	}

	public static void printSampleDistAsSif(List<AlterationPack> packs)
	{
		printSampleDistAsSif(getMutatedPerSample(packs));
	}

	public static void printSampleDistAsSif(Map<Integer, Set<String>> map)
	{
		for (int i = 0; i < map.size() - 1; i++)
		{
			Set<String> setI = map.get(i);

			for (int j = i + 1; j < map.size(); j++)
			{
				Set<String> setJ = map.get(j);

				Set<String> com = new HashSet<String>(setI);
				com.retainAll(setJ);

				if (com.size() > 100)
				{
					System.out.println(i + "\tBINDS_TO\t" + j);
				}
			}

		}

	}
}
