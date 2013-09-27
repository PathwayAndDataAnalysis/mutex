package org.cbio.mutex;

import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;
import org.cbio.causality.util.ArrayUtil;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

/**
 * A gene alteration data.
 * @author Ozgun Babur
 */
public class GeneAlt implements Cloneable
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
	 * Changes array.
	 */
	private boolean[] ch;

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
			Change[] c = getChanges();
			ch = new boolean[c.length];

			for (int i = 0; i < c.length; i++)
			{
				ch[i] = c[i].isAltered();
			}

			if (hyper != null) ch = crop(ch);
		}

		return ch;
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

	public void removeMinorCopyNumberAlts()
	{
		assert alt == Alteration.GENOMIC : "Use this method only with genomic alteration. " +
			"Remember that this method will alter genomic changes array.";

		int up = 0;
		int dw = 0;
		for (Change c : gene.get(Alteration.COPY_NUMBER))
		{
			if (c == Change.ACTIVATING) up++;
			else if (c == Change.INHIBITING) dw++;
		}

//		if (up == dw)
//		{
//			System.out.println("Gene " + gene.getId() + " is equally altered (up: " + up +
//				", dw: " + dw + "). Choosing down");
//		}

		Change c = dw < up ? Change.ACTIVATING : Change.INHIBITING;

		for (int i = 0; i < gene.getSize(); i++)
		{
			gene.get(Alteration.GENOMIC)[i] = gene.get(Alteration.MUTATION)[i].isAltered() ||
				gene.get(Alteration.COPY_NUMBER)[i] == c ? c : Change.NO_CHANGE;
		}

		this.ch = null;
	}

	public boolean isActivated()
	{
		int cnAc = 0;
		int cnIn = 0;

		for (Change ch : gene.get(Alteration.COPY_NUMBER))
		{
			if (ch == Change.ACTIVATING) cnAc++;
			else if (ch == Change.INHIBITING) cnIn++;
		}

		if (cnAc != cnIn) return cnAc > cnIn;

		cnIn = 0;

		for (Change ch : gene.get(Alteration.MUTATION))
		{
			if (ch == Change.INHIBITING) cnIn++;
		}
		return cnIn * 10 < size();
	}

	public String getPrint(List<Integer> order)
	{
		assert new HashSet<Integer>(order).size() == order.size();

		boolean[] ch = getBooleanChanges();
		StringBuilder buf = new StringBuilder();
		for (Integer o : order)
		{
			buf.append(ch[o] ? "x" : ".");
		}
		for (int i = 0; i < ch.length; i++)
		{
			if (!order.contains(i))
				buf.append(ch[i] ? "x" : ".");
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
			if (neg != null) neg[i] = !ch[i];
		}
	}

	public void unshuffle()
	{
		ch = null;
		neg = null;
	}

	@Override
	public String toString()
	{
		return getId();
	}
}
