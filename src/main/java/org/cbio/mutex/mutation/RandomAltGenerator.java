package org.cbio.mutex.mutation;

import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class RandomAltGenerator
{
	public static List<AlterationPack> generate(Set<String> symbols, int sampleSize)
	{
		List<AlterationPack> packs = new ArrayList<AlterationPack>();

		for (String symbol : symbols)
		{
			AlterationPack pack = new AlterationPack(symbol);
			Change[] ch = new Change[sampleSize];
			for (int i = 0; i < ch.length; i++)
			{
				ch[i] = Change.NO_CHANGE;
			}
			pack.put(Alteration.MUTATION, ch);
		}
		return null; // todo implement
	}
}
