package org.cbio.mutex;

import org.cbio.causality.util.FishersCombinedProbability;
import org.cbio.causality.util.Overlap;
import org.cbio.causality.util.RandomnessChecker;
import org.cbio.mutex.Group;
import org.junit.Test;

import java.util.Random;

/**
 * @author Ozgun Babur
 */
public class IndependentPValsTest
{
	@Test
	public void testRandomness()
	{
		boolean[][] alts = new boolean[2][100];

		Group g = new Group();

		RandomnessChecker rc = new RandomnessChecker();

		for (int i = 0; i < 10000; i++)
		{
			for (int j = 0; j < alts.length; j++)
			{
				for (int k = 0; k < alts[j].length; k++)
				{
					alts[j][k] = Math.random() <= 1. / alts.length;
				}
			}

			double[] pvals = g.calcIndependentPVals(alts);

			double pval = FishersCombinedProbability.pValue(pvals);
			rc.add(pval);
		}

		System.out.println(rc.getStatusForThreshold(0.05));
	}

	@Test
	public void testMutexPVal()
	{
		boolean[][] alts = new boolean[2][100];

		for (int j = 0; j < alts.length; j++)
		{
			for (int k = 0; k < alts[j].length; k++)
			{
				alts[j][k] = Math.random() <= 0.5;
			}
		}

		double orig_pval = Overlap.calcMutexPval(alts[0], alts[1], null);
		System.out.println("orig_pval = " + orig_pval);

		RandomnessChecker rc = new RandomnessChecker();

		for (int i = 0; i < 10000; i++)
		{
			randomize(alts);
			double pval = Overlap.calcMutexPval(alts[0], alts[1], null);
			rc.add(pval);
		}

		System.out.println(rc.getStatusForThreshold(orig_pval));

	}

	@Test
	public void testIndependentMutexPVal()
	{
		boolean[][] alts = new boolean[3][200];

		for (int j = 0; j < alts.length; j++)
		{
			for (int k = 0; k < alts[j].length; k++)
			{
				alts[j][k] = Math.random() <= 1D / alts.length;
			}
		}

		Group g = new Group();

		double orig_pval = FishersCombinedProbability.pValue(g.calcIndependentPVals(alts));
		System.out.println("orig_pval = " + orig_pval);

		RandomnessChecker rc = new RandomnessChecker();

		for (int i = 0; i < 10000; i++)
		{
			randomize(alts);
			double pval = FishersCombinedProbability.pValue(g.calcIndependentPVals(alts));
			rc.add(pval);
		}

		System.out.println(rc.getStatusForThreshold(orig_pval));
	}

	static final Random r = new Random();
	public void randomize(boolean[][] alts)
	{
		int m = alts[0].length;

		for (boolean[] alt : alts)
		{
			for (int i = 0; i < m; i++)
			{
				int a = r.nextInt(m);
				boolean t = alt[a];
				alt[a] = alt[i];
				alt[i] = t;
			}
		}
	}
}
