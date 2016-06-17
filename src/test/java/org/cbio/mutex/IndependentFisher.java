package org.cbio.mutex;

import org.panda.utility.Progress;
import org.panda.utility.statistics.Overlap;
import org.panda.utility.statistics.Summary;
import org.panda.utility.statistics.UniformityChecker;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author Ozgun Babur
 */
public class IndependentFisher
{
	boolean[] g1;
	List<Boolean> g2;

	public IndependentFisher(int size, int alt1, int alt2)
	{
		List<Boolean> g11 = generateGene(size, alt1);
		g1 = toPrim(g11);
		g2 = generateGene(size, alt2);
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

	public double spit()
	{
		Collections.shuffle(g2);
		return Overlap.calcMutexPval(g1, toPrim(g2));
	}

	public double[] spitMulti(int cnt)
	{
		double[] p = new double[cnt];
		for (int i = 0; i < cnt; i++)
		{
			p[i] = spit();
		}
		return p;
	}

	public double spitMax(int cnt)
	{
		return Summary.max(spitMulti(cnt));
	}

	public double[] generateMaxPvals(int size, int cnt)
	{
		double[] p = new double[size];

		Progress prg = new Progress(size);
		for (int i = 0; i < size; i++)
		{
			p[i] = spitMax(cnt);
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

	public static void main(String[] args)
	{
		IndependentFisher inFi = new IndependentFisher(1000, 400, 400);
		int dim = 2;
		double[] p = inFi.generateMaxPvals(10000, dim);
		takePower(p, dim);

		List<Double> list = new ArrayList<>(p.length);
		for (double v : p) list.add(v);
		UniformityChecker.plot(list);
	}
}
