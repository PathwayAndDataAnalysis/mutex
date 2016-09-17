package org.cbio.mutex;

import org.panda.utility.ArrayUtil;
import org.panda.utility.statistics.Summary;

import java.io.Serializable;
import java.util.*;

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
	 * If present, mapping of types to the indices in the alteration array.
	 */
	private Map<String, int[]> typeMap;

	/**
	 * If present, subsets of the alteration array divided to types.
	 */
	private Map<String, boolean[]> typeAlts;

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
			try
			{
				alterations[i] = Integer.parseInt(data[i + 1]);
			} catch (NumberFormatException e){
				System.out.println("Matrix value is not integer. Gene = " + id + ", column = " + i);
				throw e;
			}
		}
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

	public double getMutatedRatio()
	{
		boolean[] b = getMutated();
		return Math.round((Summary.countTrue(b) / (double) b.length) * 1E10) / 1E10;
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
	 * Sets the tissue to sample association map, if ever exists.
	 */
	public void setTypeMap(Map<String, int[]> map)
	{
		this.typeMap = map;
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

	public void shuffle(Random r)
	{
		if (typeMap == null)
		{
			boolean[] ch = getBooleanChanges();
			ArrayUtil.shuffle(ch, r);
		}
		else
		{
			if (typeAlts == null)
			{
				boolean[] ch = getBooleanChanges();
				typeAlts = new HashMap<>();
				typeMap.keySet().forEach(t ->
				{
					int[] ind = typeMap.get(t);
					boolean[] b = new boolean[ind.length];
					for (int i = 0; i < ind.length; i++)
					{
						b[i] = ch[ind[i]];
					}

					if (!ArrayUtil.isUniform(b)) typeAlts.put(t, b);
				});
			}

			typeAlts.keySet().forEach(t ->
			{
				boolean[] val = typeAlts.get(t);
				ArrayUtil.shuffle(val, r);
				int[] ind = typeMap.get(t);
				for (int i = 0; i < val.length; i++)
				{
					ch[ind[i]] = val[i];
				}
			});
		}
	}

	/**
	 * Beware. There is no un-shuffle of this operation.
	 */
	public void shufflePermanent()
	{
		Integer[] arr = new Integer[alterations.length];
		for (int i = 0; i < arr.length; i++)
		{
			arr[i] = alterations[i];
		}
		Collections.shuffle(Arrays.asList(arr));
		for (int i = 0; i < arr.length; i++)
		{
			alterations[i] = arr[i];
		}
		ch = null;
	}

	public void unshuffle()
	{
		ch = null;
		typeAlts = null;
	}

	public void unshuffleSticky()
	{
		shuf = null;
		ch = null;
		typeAlts = null;
		randScores = randScoresSave;
		randScoresSave = null;
	}

	public void shuffleSticky()
	{
		if (randScoresSave != null) throw new RuntimeException("The method shuffleSticky cannot " +
			"be called while a previous call exists. Please call unshuffleSticky before calling" +
			"this method second time.");

		shuffle(new Random());
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
}
