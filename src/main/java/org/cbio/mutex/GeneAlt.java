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
	public boolean[] shuf;

	/**
	 * Changes array.
	 */
	boolean[] ch;

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
	List<Double> randScores;

	/**
	 * This is the estimated null distribution of group scores with this gene.
	 */
	private List<Double> randScoresSave;

	private static final long serialVersionUID = 466072025998513732L;

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
	public Change[] getChanges(Alteration alt)
	{
		return crop(gene.get(alt));
	}

	public AlterationPack getGene()
	{
		return gene;
	}

	public boolean[] getHyper()
	{
		return hyper;
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
	public boolean[] getBooleanChanges(Alteration alt)
	{
		Change[] ch = getChanges(alt);
		if (ch == null) return null;

		boolean[] b = new boolean[ch.length];
		for (int i = 0; i < ch.length; i++)
		{
			b[i] = ch[i].isAltered();
		}
		return b;
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
				Change[] c = getChanges(alt);
				ch = new boolean[c.length];

				for (int i = 0; i < c.length; i++)
				{
					ch[i] = c[i].isAltered();
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

	public int getAltCnt()
	{
		getBooleanChanges();
		return altCnt;
	}

	private Change[] crop(Change[] ch)
	{
		if (hyper == null || ch == null) return ch;
		assert ch.length == hyper.length;

		Change[] c = new Change[ArrayUtil.countValue(hyper, false)];

		int k = 0;
		for (int i = 0; i < hyper.length; i++)
		{
			if (!hyper[i]) c[k++] = ch[i];
		}
		return c;
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
		Change[] c = getChanges(alt);
		Change[] mut = getChanges(Alteration.MUTATION);
		Change[] cna = getChanges(Alteration.COPY_NUMBER);
		if (cna == null)
		{
			cna = new Change[mut.length];
			for (int i = 0; i < cna.length; i++)
			{
				cna[i] = Change.NO_CHANGE;
			}
		}

		String[] let = new String[c.length];

		for (int i = 0; i < c.length; i++)
		{
			let[i] = c[i].isAltered() ? (mut[i].isAltered() ?
				(cna[i] == Change.ACTIVATING ? "B" : cna[i] == Change.INHIBITING ? "E" : "M") :
				(cna[i] == Change.ACTIVATING ? "A" : cna[i] == Change.INHIBITING ? "D" : "." ))
				: ".";
		}

		return let;
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

	public double getMinPval()
	{
		return getPvalOfScore(randScores.get(0));
	}

	public GeneAlt copy()
	{
		return new GeneAlt(gene, alt, hyper);
	}
}
