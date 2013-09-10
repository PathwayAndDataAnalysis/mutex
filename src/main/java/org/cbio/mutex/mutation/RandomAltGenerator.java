package org.cbio.mutex.mutation;

import org.cbio.causality.idmapping.Length;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;

import java.util.*;

/**
 * @author Ozgun Babur
 */
public class RandomAltGenerator
{
	public static Map<String, String[]> generate(List<String> symbols, int[] mutCnt,
		Set<String>[] lost)
	{
		Map<String, String[]> mutMap = new HashMap<String, String[]>(symbols.size());

		List<Integer> cumLength = new ArrayList<Integer>(symbols.size());

		for (String symbol : symbols)
		{
			mutMap.put(symbol, new String[mutCnt.length]);
		}

		Random r = new Random();

		for (int i = 0; i < mutCnt.length; i++)
		{
			cumLength.clear();
			int sum = 0;
			List<String> current = new ArrayList<String>();

			for (String symbol : symbols)
			{
				if (lost[i].contains(symbol)) continue;
				current.add(symbol);

				Integer length = Length.of(symbol);

				assert length > 0;

				sum += length;
				cumLength.add(sum);
			}

			int total = cumLength.get(cumLength.size() - 1);


			for (int j = 0; j < mutCnt[i]; j++)
			{
				int val = r.nextInt(total) + 1;

				int hit = getHitIndex(cumLength, val);

				String[] muts = mutMap.get(current.get(hit));
				if (muts[i] == null) muts[i] = "x";
				else muts[i] += ",x";
			}
		}

		return mutMap;
	}

	private static int getHitIndex(List<Integer> cumLength, int val)
	{
		int min = 0;
		int max = cumLength.size();
		int index = max / 2;

		boolean forassert = true;

		while (!isHit(cumLength, val, index))
		{
			assert forassert;

			if (cumLength.get(index) < val)
			{
				if (max - index == 1)
				{
					index = max;
					forassert = false;
				}
				else
				{
					min = index;
					index = (max + index) / 2;
				}
			}
			else
			{
				if (index - min == 1)
				{
					index = min;
					forassert = false;
				}
				else
				{
					max = index;
					index = (min + index) / 2;
				}
			}
		}

		return index;
	}

	private static boolean isHit(List<Integer> cumLength, int val, int index)
	{
		return val <= cumLength.get(index) &&
			(index == 0 || val > cumLength.get(index - 1));
	}

	public static List<AlterationPack> prepareAltpAcksFromMutMap(Map<String, String[]> mutMap)
	{
		List<AlterationPack> packs = new ArrayList<AlterationPack>(mutMap.size());
		int sampleSize = mutMap.values().iterator().next().length;

		for (String symbol : mutMap.keySet())
		{
			String[] muts = mutMap.get(symbol);
			AlterationPack pack = new AlterationPack(symbol);
			Change[] ch = new Change[sampleSize];
			for (int i = 0; i < ch.length; i++)
			{
				ch[i] = muts[i] == null || muts[i].isEmpty() ?
					Change.NO_CHANGE : Change.UNKNOWN_CHANGE;
			}
			pack.put(Alteration.MUTATION, ch);
			packs.add(pack);
		}
		return packs;
	}
}
