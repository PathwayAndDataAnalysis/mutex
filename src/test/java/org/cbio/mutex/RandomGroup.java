package org.cbio.mutex;

import org.panda.utility.ArrayUtil;
import org.panda.utility.Progress;
import org.panda.utility.statistics.Histogram;
import org.panda.utility.statistics.Overlap;
import org.panda.utility.statistics.Summary;
import org.panda.utility.statistics.UniformityChecker;

import java.text.DecimalFormat;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class RandomGroup
{
	List<List<Boolean>> group;

	public RandomGroup(int size, int... altCnt)
	{
		group = new ArrayList<List<Boolean>>();
		for (int cnt : altCnt)
		{
			group.add(generateGene(size, cnt));
		}
	}

	public void shuffle()
	{
		for (List<Boolean> gene : group)
		{
			Collections.shuffle(gene);
		}
	}

	public boolean[][] getOneAgainstOthers(int one)
	{
		assert one < group.size();

		boolean[] single = toPrim(one);
		boolean[] others = mergeOthers(one);

		return new boolean[][]{single, others};
	}

	public boolean[] toPrim(int index)
	{
		assert index < group.size();

		List<Boolean> gene = group.get(index);

		boolean[] b = new boolean[gene.size()];

		for (int i = 0; i < b.length; i++)
		{
			b[i] = gene.get(i);
		}

		return b;
	}

	public boolean[][] negate()
	{
		boolean[][] b = new boolean[group.size()][];

		for (int i = 0; i < group.size(); i++)
		{
			b[i] = ArrayUtil.negate(toPrim(i));
		}
		return b;
	}

	public boolean[] mergeOthers(int exclude)
	{
		assert exclude < group.size();

		boolean[] b = new boolean[group.get(exclude).size()];

		for (int i = 0; i < b.length; i++)
		{
			b[i] = false;

			for (int j = 0; j < group.size(); j++)
			{
				if (j == exclude) continue;

				if (group.get(j).get(i))
				{
					b[i] = true;
					break;
				}
			}
		}

		return b;
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

	private double[] calcPvals()
	{
		double[] pvals = new double[group.size()];

		for (int i = 0; i < group.size(); i++)
		{
			boolean[][] b = getOneAgainstOthers(i);

			pvals[i] = Overlap.calcMutexPval(b[0], b[1]);
		}
		return pvals;
	}

	private double[] calcPvalsOf3()
	{
		boolean[][] neg = negate();

		return new double[]{
			Overlap.calcCoocPval(toPrim(0), neg[1], neg[2]),
			Overlap.calcCoocPval(toPrim(1), neg[0], neg[2]),
			Overlap.calcCoocPval(toPrim(2), neg[1], neg[0]),
		};
	}

	public List<Double>[] generatePvalsDist(int trials)
	{
		List<Double>[] list = new List[group.size()];
		for (int i = 0; i < group.size(); i++)
		{
			list[i] = new ArrayList<Double>();
		}

		Progress prog = new Progress(trials);

		for (int i = 0; i < trials; i++)
		{
			shuffle();
			double[] p = calcPvalsOf3();

			for (int j = 0; j < p.length; j++)
			{
				list[j].add(p[j]);
			}
			prog.tick();
		}
		return list;
	}

	DecimalFormat f = new DecimalFormat("0.00");
	private void printCorMatrix(List<Double>[] list)
	{
		for (int i = 0; i < list.length; i++)
		{
			System.out.println();
			for (int j = 0; j < list.length; j++)
			{
				if (i == j) System.out.print("\t    ");
				else
					System.out.print("\t" + f.format(cor(list[i], list[j])));
			}
		}
	}

	private void printWorstDistr(List<Double>[] lists)
	{
		List<Double> v = getWorstPvals(lists);
		UniformityChecker.plot(v);
	}

	private List<Double> getWorstPvals(List<Double>[] lists)
	{
		int size = lists[0].size();
		List<Double> v = new ArrayList<>(size);

		for (int i = 0; i < size; i++)
		{
			double x = 0;
			for (List<Double> list : lists)
			{
				if (list.get(i) > x) x = list.get(i);
			}
			v.add(0D);//todo
//			v.add(x);
		}
		return v;
	}

	private void makeUniform(double[] v)
	{
		for (int i = 0; i < v.length; i++)
		{
			v[i] = Math.pow(v[i], (group.size() - 1));
		}
	}

	private double cor(List<Double> list1, List<Double> list2)
	{
		double m1 = Summary.meanOfDoubles(list1);
		double m2 = Summary.meanOfDoubles(list2);

		double nom = 0;
		double s1 = 0;
		double s2 = 0;

		for (int i = 0; i < list1.size(); i++)
		{
			double v1 = list1.get(i) - m1;
			double v2 = list2.get(i) - m2;
			nom += v1 * v2;
			s1 += v1 * v1;
			s2 += v2 * v2;
		}

		s1 = Math.sqrt(s1);
		s2 = Math.sqrt(s2);

		return nom / (s1 * s2);
	}

	private void printFrequency(double[] v)
	{
		Map<Double, Integer> map = new HashMap<Double, Integer>();

		for (double v1 : v)
		{
			if (!map.containsKey(v1)) map.put(v1, 0);
			map.put(v1, map.get(v1) + 1);
		}

		List<Double> list = new ArrayList<Double>(map.keySet());
		Collections.sort(list);

		for (Double d : list)
		{
			System.out.println(d + "\t" + map.get(d));
		}
	}

//	private List<Double> selectLeastSignificant(List<Double>[] lists)
//	{
//		List<Double> select = new ArrayList<Double>();
//
//		for (int i = 0; i < lists[0].size(); i++)
//		{
//			double v =
//		}
//	}

	public static void main(String[] args)
	{
		RandomGroup group = new RandomGroup(100, 40, 40, 40);
		List<Double>[] lists = group.generatePvalsDist(10000);
		List<Double> pv = group.getWorstPvals(lists);
//		group.makeUniform(pv);
		UniformityChecker.plot(pv);
	}

//	public static void main(String[] args)
//	{
//		int x = 10;
//
//	}
}
