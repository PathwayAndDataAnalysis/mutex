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
	 * Pack of gene alterations.
	 */
	private AlterationPack gene;

	/**
	 * Related alteration type.
	 */
	private Alteration alt;

	/**
	 * Changes array shuffled in a sticky way.
	 */
	private boolean[] shuf;

	/**
	 * Changes array.
	 */
	private boolean[] ch;

	/**
	 * Count of altered samples.
	 */
	private int altCnt;

	/**
	 * Negative of changes array.
	 */
	private boolean[] neg;

	/**
	 * A marking for hyper mutated samples. The two boolean arrays ch and neg shouldn't contain data
	 * for these samples, so their sizes should be smaller.
	 */
	private boolean[] hyper;

	/**
	 * This is the estimated null distribution of group scores with this gene.
	 */
	double[] randScores;

	/**
	 * This is the estimated null distribution of group scores with this gene.
	 */
	private double[] randScoresSave;

	/**
	 * Constructor with parameters.
	 * @param gene pack of gene alterations
	 * @param alt the type of alteration to use
	 */
	public GeneAlt(AlterationPack gene, Alteration alt)
	{
		this(gene, alt, null);
	}

	/**
	 * Constructor with parameters.
	 * @param gene pack of gene alterations
	 * @param alt the type of alteration to use
	 */
	public GeneAlt(AlterationPack gene, Alteration alt, boolean[] hyper)
	{
		this.gene = gene;
		this.alt = alt;
		this.hyper = hyper;
	}

	/**
	 * Gets the related change array in the pack.
	 * @return the related change array
	 */
	private Change[] getChanges()
	{
		return gene.get(alt);
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
				Change[] c = getChanges();
				ch = new boolean[c.length];

				for (int i = 0; i < c.length; i++)
				{
					ch[i] = c[i].isAltered();
				}

				if (hyper != null) ch = crop(ch);
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

	public int getAltCnt()
	{
		getBooleanChanges();
		return altCnt;
	}

	private boolean[] crop(boolean[] ch)
	{
		assert ch.length == hyper.length;

		boolean[] c = new boolean[ArrayUtil.countValue(hyper, false)];

		int k = 0;
		for (int i = 0; i < hyper.length; i++)
		{
			if (!hyper[i]) c[k++] = ch[i];
		}
		return c;
	}

	private String[] crop(String[] let)
	{
		assert let.length == hyper.length;

		String[] c = new String[ArrayUtil.countValue(hyper, false)];

		int k = 0;
		for (int i = 0; i < hyper.length; i++)
		{
			if (!hyper[i]) c[k++] = let[i];
		}
		return c;
	}

	/**
	 * Gets the sample values in a boolean array.
	 * @return changes in a boolean array
	 */
	public boolean[] getNegativeChanges()
	{
		if (neg == null) neg = ArrayUtil.negate(getBooleanChanges());
		return neg;
	}

	/**
	 * Gets the ID of the gene.
	 * @return id of the gene
	 */
	public String getId()
	{
		return gene.getId();
	}

	/**
	 * Two gene alterations are equal only if the pack and the alteration are equal.
	 * @param obj the object to check
	 * @return true if equal
	 */
	@Override
	public boolean equals(Object obj)
	{
		if (obj instanceof GeneAlt)
		{
			GeneAlt other = (GeneAlt) obj;
			return gene.equals(other.gene) && alt.equals(other.alt);
		}
		return false;
	}

	/**
	 * Gets the hash code to use in maps.
	 * @return hash code
	 */
	@Override
	public int hashCode()
	{
		return gene.hashCode() + alt.hashCode();
	}

	/**
	 * Gets size of the samples.
	 * @return size of the samples
	 */
	public int size()
	{
		return getBooleanChanges().length;
//		return gene.getSize();
	}

	/**
	 * Gets the ratio of altered samples.
	 * @return alteration ratio
	 */
	public double getAlteredRatio()
	{
		return countAltered() / (double) getBooleanChanges().length;
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

		String[] let = getLetterChanges();

		assert getBooleanChanges().length == let.length;

		StringBuilder buf = new StringBuilder();
		for (Integer o : order)
		{
			buf.append(let[o]);
		}
		for (int i = 0; i < let.length; i++)
		{
			if (!order.contains(i))
				buf.append(let[i]);
		}

		buf.append("  ").append(getId());
		return buf.toString();
	}

	private String[] getLetterChanges()
	{
		Change[] c = getChanges();
		Change[] mut = gene.get(Alteration.MUTATION);
		Change[] cna = gene.get(Alteration.COPY_NUMBER);

		String[] let = new String[c.length];

		for (int i = 0; i < c.length; i++)
		{
			let[i] = c[i].isAltered() ? mut[i].isAltered() ? "M" :
				cna[i] == Change.ACTIVATING ? "A" : cna[i] == Change.INHIBITING ? "D" :
					"." : ".";
		}

		return crop(let);
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
			if (neg != null) neg[i] = !ch[i];
		}
	}

	public void unshuffle()
	{
		ch = null;
		neg = null;
	}

	public void unshuffleSticky()
	{
		shuf = null;
		ch = null;
		neg = null;
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
	public void setRandScores(double[] randScores)
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

		return i / (double) randScores.length;
	}

	public double getMinPval()
	{
		return getPvalOfScore(randScores[0]);
	}

	public GeneAlt copy()
	{
		return new GeneAlt(gene, alt, hyper);
	}
}
