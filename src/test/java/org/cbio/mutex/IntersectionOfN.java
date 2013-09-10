package org.cbio.mutex;

import org.cbio.causality.util.*;

import java.util.ArrayList;
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

	private void shuffle()
	{
		for (List<Boolean> gene : g)
		{
			Collections.shuffle(gene);
		}
	}

	public double[] spitMulti(int cnt)
	{
		double[] p = new double[cnt];

		Progress prg = new Progress(p.length);
		for (int i = 0; i < cnt; i++)
		{
//			p[i] = maxOfN(g.length);
			p[i] = spit();
			prg.tick();
		}
		return p;
	}

	private static void takePower(double[] v, double expo)
	{
		for (int i = 0; i < v.length; i++)
		{
			v[i] = Math.pow(v[i], expo);
		}
	}

	private boolean[][] toBoolArray()
	{
		boolean[][] b = new boolean[g.length][];
		for (int i = 0; i < b.length; i++)
		{
			b[i] = toPrim(g[i]);
//			if (i > 0) b[i] = ArrayUtil.negate(b[i]);
		}
		return b;
	}

	public double maxOfN(int n)
	{
		double[] pval = new double[n];
		for (int i = 0; i < n; i++)
		{
			pval[i] = spit();
		}
		return Math.pow(Summary.max(pval), n);
	}

	public double spit()
	{
		shuffle();
//		return Overlap.calcCoocPval(toBoolArray());

		boolean[][][] b = new boolean[g.length][g.length][];
		for (int i = 0; i < g.length; i++)
		{
			for (int j = 0; j < g.length; j++)
			{
				b[i][j] = toPrim(g[j]);
				if (j!=i) b[i][j] = ArrayUtil.negate(b[i][j]);
			}
		}
		double[] pv = new double[g.length];
		for (int i = 0; i < g.length; i++)
		{
			pv[i] = Overlap.calcCoocPval(b[i]);
		}
		return Math.pow(Summary.max(pv), g.length);
	}

	public static void main(String[] args)
	{
		IntersectionOfN inFi = new IntersectionOfN(100, 0.3, 0.5, 0.2);
		double[] p = inFi.spitMulti(10000);
//		takePower(p, 3);
		DiscretePvalHisto h = new DiscretePvalHisto(p, 0.05);
		h.plot();

		System.out.println();
		Frequency f = new Frequency();
		f.count(p);
		f.print();
	}
}
