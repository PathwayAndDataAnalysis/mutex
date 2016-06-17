package org.cbio.mutex;

import org.panda.utility.graph.Graph;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

/**
 * @deprecated
 */
public class PrintPerformances
{
	private static final int VERSION = 1; // or 2 or 3

	public static void main(String[] args) throws IOException, InterruptedException
	{
		PrintPerformances t = new PrintPerformances();
//		t.printTrueCounts();
		t.printDiscoveryGraphs();
	}

	private void printDiscoveryGraphs() throws IOException
	{
		Graph graph = new Network();
		final Map<String, String> labelMap = null;
//		final Map<String, String> labelMap = Simulation.getLabelMap(dataset,
//			Simulation.loadSimData(graph.getSymbols(), dataset).keySet());

		final String base = "/home/ozgun/Documents/mutex-comparison/";
		CustomResult[] cus = new CustomResult[]{
			new CustomResult()
			{
				@Override
				public String getName()
				{
					return "Mutex";
				}

				@Override
				public String getFilename1()
				{
					return "data-simulation/large-dataset/with-network/ranked-groups.txt";
				}

				@Override
				public String getFilename2()
				{
					return "data-simulation/small-dataset/with-network/ranked-groups.txt";
				}

				@Override
				public String[] getGroupMembers(String line)
				{
					if (line.startsWith("S")) return null;
					line = line.substring(line.indexOf("\t") + 1);
					line = line.substring(line.indexOf("\t") + 1);
					return line.split("\t");
				}
			},
			new CustomResult()
			{
				@Override
				public String getName()
				{
					return "Pair search";
				}

				@Override
				public String getFilename1()
				{
					return "pairs-v1.txt";
				}

				@Override
				public String getFilename2()
				{
					return "pairs-v2.txt";
				}

				@Override
				public String[] getGroupMembers(String line)
				{
					return line.split("\t");
				}
			},
			new CustomResult()
			{
				@Override
				public String getName()
				{
					return "RME";
				}

				@Override
				public String getFilename1()
				{
					return base + "rmeMod/simdata1/topModules";
				}

				@Override
				public String getFilename2()
				{
					return base + "rmeMod/simdata2/topModules";
				}

				@Override
				public String[] getGroupMembers(String line)
				{
					return line.substring(line.lastIndexOf("\t") + 1).split(",");
				}
			},
			new CustomResult()
			{
				@Override
				public String getName()
				{
					return "Dendrix";
				}

				@Override
				public String getFilename1()
				{
					return base + "Dendrix/results-v1.txt";
				}

				@Override
				public String getFilename2()
				{
					return base + "Dendrix/results-v2.txt";
				}

				@Override
				public String[] getGroupMembers(String line)
				{
					return line.substring(line.indexOf("\t") + 1, line.lastIndexOf("\t")).split("\t");
				}
			},
			new CustomResult()
			{
				@Override
				public String getName()
				{
					return "MDPFinder";
				}

				@Override
				public String getFilename1()
				{
					return base + "MDPfinder/results-v1.txt";
				}

				@Override
				public String getFilename2()
				{
					return base + "MDPfinder/results-v2.txt";
				}

				@Override
				public String[] getGroupMembers(String line)
				{
					if (line.startsWith("gene")) return null;
					line = line.split("\t")[0];
					return line.split(",");
				}
			},
			new CustomResult()
			{
				@Override
				public String getName()
				{
					return "Multi-dendrix";
				}

				@Override
				public String getFilename1()
				{
					return base + "multi-dendrix/results-v1.txt";
				}

				@Override
				public String getFilename2()
				{
					return base + "multi-dendrix/results-v2.txt";
				}

				@Override
				public String[] getGroupMembers(String line)
				{
					return line.substring(0, line.lastIndexOf("\t")).split("\t");
				}
			},
			new CustomResult()
			{
				@Override
				public String getName()
				{
					return "MEMo";
				}

				@Override
				public String getFilename1()
				{
					return null;
				}

				@Override
				public String getFilename2()
				{
					return base + "MEMo/cancer_data/simulation2/MemoReport.txt";
				}

				@Override
				public String[] getGroupMembers(String line)
				{
					if (line.startsWith("Module ID")) return null;

					line = line.split("\t")[1];
					String[] tok = line.split(",");
					String[] result = new String[tok.length - 1];

					for (int i = 0; i < result.length; i++)
					{
						result[i] = tok[i].trim().substring(0, tok[i].trim().indexOf(" "));
						result[i] = labelMap.get(result[i]);
					}
					return result;
				}
			},
			new CustomResult()
			{
				@Override
				public String getName()
				{
					return "MEMo-with-supplement";
				}

				@Override
				public String getFilename1()
				{
					return null;
				}

				@Override
				public String getFilename2()
				{
					return base + "MEMo/cancer_data/simulation2-suppl/MemoReport.txt";
				}

				@Override
				public String[] getGroupMembers(String line)
				{
					if (line.startsWith("Module ID")) return null;

					line = line.split("\t")[1];
					String[] tok = line.split(",");
					String[] result = new String[tok.length - 1];

					for (int i = 0; i < result.length; i++)
					{
						result[i] = tok[i].trim().substring(0, tok[i].trim().indexOf(" "));
						result[i] = labelMap.get(result[i]);
					}
					return result;
				}
			},
			new CustomResult()
			{
				@Override
				public String getName()
				{
					return "ME";
				}

				@Override
				public String getFilename1()
				{
					return null;
				}

				@Override
				public String getFilename2()
				{
					return base + "muex/results.txt";
				}

				@Override
				public String[] getGroupMembers(String line)
				{
					if (line.startsWith("G")) return null;
					line = line.substring(line.indexOf("\t") + 1);
					String[] tok = line.split("\t");
					String[] result = new String[3];
					System.arraycopy(tok, 0, result, 0, 3);
					return result;
				}
			},
			new CustomResult()
			{
				@Override
				public String getName()
				{
					return "Mutex (no reduction)";
				}

				@Override
				public String getFilename1()
				{
					return "data-simulation/large-dataset/without-network/ranked-groups.txt";
				}

				@Override
				public String getFilename2()
				{
					return "data-simulation/small-dataset/without-network/ranked-groups.txt";
				}

				@Override
				public String[] getGroupMembers(String line)
				{
					if (line.startsWith("S")) return null;
					line = line.substring(line.indexOf("\t") + 1);
					line = line.substring(line.indexOf("\t") + 1);
					return line.split("\t");
				}
			},
		};

		Map<String, List<double[]>> map = new HashMap<String, List<double[]>>();

		for (CustomResult cu : cus)
		{
			List<double[]> points = getROCPoints(cu);
//			List<double[]> points = getDiscoveryGraph(cu);
			if (points == null) continue;

			map.put(cu.getName(), points);

//			System.out.println("\t" + cu.getName());
//			for (double[] point : points)
//			{
//				System.out.println(point[0] + "\t" + point[1]);
//			}
//			System.out.println();
		}

		int size = 0;
		for (List<double[]> list : map.values())
		{
			if (list.size() > size) size = list.size();
		}

		for (CustomResult cu : cus)
		{
			if (map.containsKey(cu.getName())) System.out.print("\t" + cu.getName() + "\t");
		}
		for (int i = 0; i < size; i++)
		{
			System.out.println();

			for (CustomResult cu : cus)
			{
				if (map.containsKey(cu.getName()))
				{
					List<double[]> list = map.get(cu.getName());
					if (list.size() > i)
						System.out.print(list.get(i)[0] + "\t" + list.get(i)[1] + "\t");
					else
						System.out.print("\t\t");
				}
			}
		}
		System.out.println();
	}


	private List<double[]> getROCPoints(CustomResult cr) throws FileNotFoundException
	{
		double CP = VERSION == 1 ? 150 : 60;
		double CN = VERSION == 1 ? 680 : 96;

		Set<String> ts = new HashSet<String>();
		Set<String> fs = new HashSet<String>();

		List<double[]> plot = new ArrayList<double[]>();

		String filename = VERSION == 1 ? cr.getFilename1() :
			VERSION == 2 ? cr.getFilename2() : null;

		if (filename == null) return null;
		if (!new File(filename).exists())
		{
			System.out.println("File not found: " + filename);
			return null;
		}

		Scanner sc = new Scanner(new File(filename));
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			String[] members = cr.getGroupMembers(line);
			if (members == null) continue;
			for (String gene : members)
			{
				if (gene.startsWith("T")) ts.add(gene);
				else fs.add(gene);
			}

			double x = fs.size() / CN;
			double y = ts.size() / CP;

			if (plot.isEmpty() ||
				x > plot.get(plot.size() - 1)[0] ||
				y > plot.get(plot.size() - 1)[1])
			{
				plot.add(new double[]{x, y});
			}
		}
		return plot;
	}

	private List<double[]> getDiscoveryGraph(CustomResult cr) throws FileNotFoundException
	{
		double CP = VERSION == 1 ? 150 : 60;
		double CN = VERSION == 1 ? 680 : 96;

		Set<String> ts = new HashSet<String>();
		Set<String> fs = new HashSet<String>();

		String filename = VERSION == 1 ? cr.getFilename1() :
			VERSION == 2 ? cr.getFilename2() : null;

		if (filename == null) return null;
		if (!new File(filename).exists())
		{
			System.out.println("File not found: " + filename);
			return null;
		}

		Map<Double, Integer> counts = new HashMap<Double, Integer>();
		for (int i = 0; i <= 100; i++)
		{
			double d =  i / 100D;
			counts.put(d, 0);
		}

		Scanner sc = new Scanner(new File(filename));
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			String[] members = cr.getGroupMembers(line);
			if (members == null) continue;
			for (String gene : members)
			{
				if (gene.startsWith("T")) ts.add(gene);
				else fs.add(gene);
			}

			double fdr = fs.size() / (double) (fs.size() + ts.size());

			for (Double v : counts.keySet())
			{
				if (v >= fdr && counts.get(v) < ts.size()) counts.put(v, ts.size());
			}
		}

		List<double[]> res = new ArrayList<double[]>();

		for (int i = 0; i <= 100; i++)
		{
			double d =  i / 100D;
//			res.add(new double[]{d, counts.get(d) / CP});
			res.add(new double[]{d, counts.get(d)});
		}
		return res;
	}

	private void printTrueCounts() throws FileNotFoundException
	{
		int t = 0;
		int f = 0;

		Scanner sc = new Scanner(new File("/home/ozgun/Documents/mutex-comparison/muex/data/SimData.txt"));
		sc.nextLine();
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			if (line.startsWith("T")) t++;
			else f++;
		}
		System.out.println("t = " + t);
		System.out.println("f = " + f);
	}

	private interface CustomResult
	{
		String getName();
		String getFilename1();
		String getFilename2();
		String[] getGroupMembers(String line);
	}
}
