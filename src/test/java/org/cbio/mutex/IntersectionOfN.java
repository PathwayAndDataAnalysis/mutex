package org.cbio.mutex;

import org.panda.utility.ArrayUtil;
import org.panda.utility.Progress;
import org.panda.utility.statistics.Overlap;
import org.panda.utility.statistics.Summary;
import org.panda.utility.statistics.UniformityChecker;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author Ozgun Babur
 */
public class IntersectionOfN
{
	List<Boolean>[] g;

	public IntersectionOfN(int size, double... alt)
	{
		g = new List[alt.length];

		for (int i = 0; i < alt.length; i++)
		{
			g[i] = generateGene(size, (int) (alt[i] * size));
		}
	}

	private ArrayList<Boolean> generateGene(int size, int altered)
	{
		ArrayList<Boolean> list = new ArrayList<Boolean>(size);

		for (int i = 0; i < size; i++)
		{
			list.add(i < altered);
		}
		return list;
	}

	public boolean[] toPrim(List<Boolean> gene)
	{
		boolean[] b = new boolean[gene.size()];

		for (int i = 0; i < b.length; i++)
		{
			b[i] = gene.get(i);
		}

		return b;
	}

	public boolean[] uniteButOne(int exclude)
	{
		boolean[] b = new boolean[g[0].size()];

		for (int i = 0; i < b.length; i++)
		{
			for (int j = 0; j < g.length; j++)
			{
				if (j == exclude) continue;
				if (g[j].get(i))
				{
					b[i] = true;
					break;
				}
			}
		}
		return b;
	}

	private void shuffle()
	{
		for (List<Boolean> gene : g)
		{
			Collections.shuffle(gene);
		}
	}

	public List<Double> spitMulti(int cnt)
	{
		List<Double> list = new ArrayList<>(cnt);

		Progress prg = new Progress(cnt);
		for (int i = 0; i < cnt; i++)
		{
			list.add(spit2());
			prg.tick();
		}
		return list;
	}

	public double spit()
	{
//		return Overlap.calcCoocPval(toBoolArray());

		boolean[][][] b = new boolean[g.length][g.length][];
		for (int i = 0; i < g.length; i++)
		{
			shuffle();
			for (int j = 0; j < g.length; j++)
			{
				b[i][j] = toPrim(g[j]);
				if (i != j) b[i][j] = ArrayUtil.negate(b[i][j]);
			}
		}
		double[] pv = new double[g.length];
		for (int i = 0; i < g.length; i++)
		{
			pv[i] = Overlap.calcCoocPval(b[i]);
		}
		double max = Math.round(1E10 * Summary.max(pv)) / 1E10;
		double p = Math.pow(max, g.length);
		return p;
	}

	public double spit2()
	{
		shuffle();
		double[] pv = new double[g.length];
		for (int i = 0; i < g.length; i++)
		{
			pv[i] = Overlap.calcMutexPval(toPrim(g[i]), uniteButOne(i));
		}
		double max = Math.round(1E10 * Summary.max(pv)) / 1E10;
		double p = Math.pow(max, g.length - 1);
//		p = Math.pow(p, 0.5-1) * Math.pow(1-p, 0.5-1);
		return p;
	}

	public static void main(String[] args)
	{
		IntersectionOfN inFi = new IntersectionOfN(1000, 0.3, 0.2,  0.3);
		List<Double> pvals = inFi.spitMulti(10000);
		UniformityChecker.plot(pvals);

//		System.out.println();
//		Frequency f = new Frequency();
//		f.count(p);
//		f.print();
	}
}
