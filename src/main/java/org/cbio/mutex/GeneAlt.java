package org.cbio.mutex;

import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;
import org.cbio.causality.util.ArrayUtil;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

/**
 * A gene alteration data.
 * @author Ozgun Babur
 */
public class GeneAlt implements Cloneable, Serializable
{
	/**
	 * Gene ID.
	 */
	String id;

	/**
	 * The array of alterations.
	 */
	int[] alterations;

	/**
	 * Changes array shuffled in a sticky way.
	 */
	public boolean[] shuf;

	/**
	 * Changes array as booleans.
	 */
	boolean[] ch;

	/**
	 * Count of altered samples.
	 */
	private int altCnt;

	/**
	 * This is the estimated null distribution of p-values.
	 */
	List<Double> randScores;

	/**
	 * This is the estimated null distribution of group scores with this gene.
	 */
	private List<Double> randScoresSave;

	private static final long serialVersionUID = 2664760285698573701L;

	/**
	 * Constructor with parameters.
	 */
	public GeneAlt(String[] data)
	{
		this.id = data[0];
		this.alterations = new int[data.length - 1];
		for (int i = 0; i < alterations.length; i++)
		{
			alterations[i] = Integer.parseInt(data[i + 1]);
		}
	}

	/**
	 * Constructor for subtype cases
	 * @param pack
	 */
	public GeneAlt(AlterationPack pack)
	{
		this.id = pack.getId();
		this.alterations = convertAlterations(pack, null);
	}

	/**
	 * Gets the sample values in a boolean array.
	 * @return changes in a boolean array
	 */
	public boolean[] getBooleanChangesCopy()
	{
		boolean[] ch = getBooleanChanges();
		boolean[] cop = new boolean[ch.length];
		System.arraycopy(ch, 0, cop, 0, ch.length);
		return cop;
	}

	/**
	 * Gets the sample values in a boolean array.
	 * @return changes in a boolean array
	 */
	public boolean[] getBooleanChanges()
	{
		if (ch == null)
		{
			if (shuf == null)
			{
				ch = new boolean[alterations.length];

				for (int i = 0; i < ch.length; i++)
				{
					ch[i] = alterations[i] != 0;
				}
			}
			else
			{
				ch = new boolean[shuf.length];
				System.arraycopy(shuf, 0, ch, 0, ch.length);
			}
			altCnt = ArrayUtil.countValue(ch, true);
		}

		return ch;
	}

	public boolean[] getMutated()
	{
		boolean[] b = new boolean[alterations.length];
		for (int i = 0; i < b.length; i++)
		{
			b[i] = alterations[i] == Letter.MUT.code ||
				alterations[i] == Letter.AMP_MUT.code ||
				alterations[i] == Letter.DEL_MUT.code;
		}
		return b;
	}

	public boolean[] getCNAltered()
	{
		boolean[] b = new boolean[alterations.length];
		for (int i = 0; i < b.length; i++)
		{
			b[i] = alterations[i] == Letter.AMP.code ||
				alterations[i] == Letter.DEL.code ||
				alterations[i] == Letter.AMP_MUT.code ||
				alterations[i] == Letter.DEL_MUT.code;
		}
		return b;
	}

	public int getAltCnt()
	{
		getBooleanChanges();
		return altCnt;
	}

	/**
	 * Gets the ID of the gene.
	 * @return id of the gene
	 */
	public String getId()
	{
		return id;
	}

	/**
	 * Two gene alterations are equal only if their ID are equal.
	 * @param obj the object to check
	 * @return true if equal
	 */
	@Override
	public boolean equals(Object obj)
	{
		return obj instanceof GeneAlt && this.id != null && ((GeneAlt) obj).id != null &&
			this.id.equals(((GeneAlt) obj).id);
	}

	/**
	 * Gets the hash code to use in maps.
	 * @return hash code
	 */
	@Override
	public int hashCode()
	{
		return id.hashCode();
	}

	/**
	 * Gets size of the samples.
	 * @return size of the samples
	 */
	public int size()
	{
		return alterations.length;
	}

	/**
	 * Gets the ratio of altered samples.
	 * @return alteration ratio
	 */
	public double getAlteredRatio()
	{
		return Math.round((countAltered() / (double) getBooleanChanges().length) * 1E10) / 1E10;
	}

	/**
	 * Gets the count of altered samples.
	 * @return alteration count
	 */
	public int countAltered()
	{
		return ArrayUtil.countValue(getBooleanChanges(), true);
	}

	public String getPrint(List<Integer> order)
	{
		assert new HashSet<Integer>(order).size() == order.size();

		Letter[] let = Letter.convertToLetter(alterations);

		StringBuilder buf = new StringBuilder();
		for (Integer o : order)
		{
			buf.append(let[o].print);
		}
		for (int i = 0; i < let.length; i++)
		{
			if (!order.contains(i))
				buf.append(let[i].print);
		}

		buf.append("  ").append(getId());
		return buf.toString();
	}

	public void shuffle()
	{
		Boolean[] bool = new Boolean[getBooleanChanges().length];
		for (int i = 0; i < bool.length; i++)
		{
			bool[i] = ch[i];
		}
		Collections.shuffle(Arrays.asList(bool));
		for (int i = 0; i < bool.length; i++)
		{
			ch[i] = bool[i];
		}
	}

	public void unshuffle()
	{
		ch = null;
	}

	public void unshuffleSticky()
	{
		shuf = null;
		ch = null;
		randScores = randScoresSave;
	}

	public void shuffleSticky()
	{
		shuffle();
		shuf = new boolean[ch.length];
		System.arraycopy(ch, 0, shuf, 0, ch.length);
		randScoresSave = randScores;
		randScores = null;
	}

	@Override
	public String toString()
	{
		return getId();
	}

	/**
	 * Random permutation scores sorted ascending. Smaller score is more significant.
	 * @param randScores
	 */
	public void setRandScores(List<Double> randScores)
	{
		this.randScores = randScores;
	}

	public double getPvalOfScore(double score)
	{
		int i = 0;
		for (double randScore : randScores)
		{
			if (randScore <= score) i++;
			else break;
		}

		return i / (double) randScores.size();
	}

	public static int[] convertAlterations(AlterationPack pack, boolean[] hyper)
	{
		int[] alts = new int[hyper == null ? pack.getSize() : ArrayUtil.countValue(hyper, false)];

		int j = 0;
		for (int i = 0; i < pack.getSize(); i++)
		{
			if (hyper != null && hyper[i])
			{
				continue;
			}

			int code;

			if (pack.get(Alteration.MUTATION)[i].isAltered())
			{
				if (pack.get(Alteration.COPY_NUMBER)[i] == Change.ACTIVATING) code = Letter.AMP_MUT.code;
				else if (pack.get(Alteration.COPY_NUMBER)[i] == Change.INHIBITING) code = Letter.DEL_MUT.code;
				else code = Letter.MUT.code;
			}
			else
			{
				if (pack.get(Alteration.COPY_NUMBER)[i] == Change.ACTIVATING) code = Letter.AMP.code;
				else if (pack.get(Alteration.COPY_NUMBER)[i] == Change.INHIBITING) code = Letter.DEL.code;
				else code = 0;
			}

			alts[j++] = code;
		}
		assert j == alts.length;
		return alts;
	}
}
