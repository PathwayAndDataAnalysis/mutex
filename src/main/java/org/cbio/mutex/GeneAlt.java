package org.cbio.mutex;

import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;

/**
 * A gene alteration data.
 * @author Ozgun Babur
 */
public class GeneAlt
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
}
