package org.cbio.mutex;

import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;

import java.util.Arrays;
import java.util.Collections;

/**
 * A gene alteration data.
 * @author Ozgun Babur
 */
public class GeneAlt implements Cloneable
{
	/**
	 * Pack of gene alterations.
	 */
	AlterationPack gene;

	/**
	 * Related alteration type.
	 */
	Alteration alt;

	/**
	 * Changes array.
	 */
	boolean[] ch;

	/**
	 * Constructor with parameters.
	 * @param gene pack of gene alterations
	 * @param alt the type of alteration to use
	 */
	public GeneAlt(AlterationPack gene, Alteration alt)
	{
		this.gene = gene;
		this.alt = alt;
	}

	/**
	 * Gets the related change array in the pack.
	 * @return the related change array
	 */
	public Change[] getChanges()
	{
		return gene.get(alt);
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
		}

		return ch;
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
		return gene.getSize();
	}

	/**
	 * Gets the ratio of altered samples.
	 * @return alteration ratio
	 */
	public double getAlteredRatio()
	{
		return gene.getAlteredRatio(alt);
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
}
