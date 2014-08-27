package org.cbio.mutex;

import org.cbio.causality.analysis.Graph;
import org.cbio.causality.data.GeneCards;

import java.io.*;
import java.util.*;

/**
 * Processes the results of many analysis and prepares an integrated result graph.
 *
 * @author Ozgun Babur
 */
public class Integrator
{
	/**
	 * Parent directory of analysis results.
	 */
	static final String dir = "data/";

	/**
	 * Names of directories containing the studies to include.
	 */
	static final String[] study = new String[]{"adrenocortical", "breast", "colorectal",
		"endometrial-cna", "endometrial-mut", "glioblastoma", "head-and-neck", "kidney-papillary", "kidney-renal",
		"leukemia", "lower-grade-glioma", "lung-adeno", "lung-squamous", "ovarian", "prostate",
		"skin", "stomach", "thyroid"};

	public static void main(String[] args) throws IOException
	{
		run();
	}

	static void run() throws IOException
	{
		Map<String, Integer> geneScore = new HashMap<String, Integer>();
		Map<String, Integer> pairScore = new HashMap<String, Integer>();
		Map<String, Set<String>> genesMap = new HashMap<String, Set<String>>();
		Map<String, List<List<String>>> groupsMap = new HashMap<String, List<List<String>>>();

		for (String s : study)
		{
			System.out.println("\nStudy = " + s);
			File f = new File(dir + s);
			assert f.exists() && f.isDirectory();

			if (!new File(f.getPath() + File.separator + "fdr-guide.txt").exists())
			{
				System.out.println("Not complete!");
				continue;
			}

			double cutoff = getBestCutoff(f.getPath());
			List<List<String>> groups = readResults(f.getPath(), cutoff);
			groupsMap.put(s, groups);
			System.out.println("Groups = " + groups.size());

			Set<String> genes = new HashSet<String>();
			Set<String> pairs = new HashSet<String>();

			for (List<String> group : groups)
			{
				for (String g1 : group)
				{
					genes.add(g1);
					for (String g2 : group)
					{
						if (g1.compareTo(g2) >= 0) continue;

						String key = key(g1, g2);
						pairs.add(key);
					}
				}
			}

			for (String gene : genes)
			{
				if (!geneScore.containsKey(gene)) geneScore.put(gene, 1);
				else geneScore.put(gene, geneScore.get(gene) + 1);
			}
			for (String pair : pairs)
			{
				if (!pairScore.containsKey(pair)) pairScore.put(pair, 1);
				else pairScore.put(pair, pairScore.get(pair) + 1);
			}

			genesMap.put(s, genes);
		}

		int brightestTone = 100;
		int extraTone = 200;
		List<String> recurrent = prepareSIF(geneScore, pairScore, "integrated", 2, 1,
			brightestTone, extraTone);

		printCancerAssociations(genesMap, recurrent);
		printTargets(groupsMap);
	}

	private static void printTargets(Map<String, List<List<String>>> groupsMap)
	{
		Graph graph = new Network();

		final Map<String, Map<String, List<List<String>>>> targetMap =
			new HashMap<String, Map<String, List<List<String>>>>();

		for (String s : groupsMap.keySet())
		{
			for (List<String> group : groupsMap.get(s))
			{
				Set<String> targets = graph.getLinkedCommonDownstream(new HashSet<String>(group));

				cropToMembersIfContains(targets, group);

				for (String target : targets)
				{
					if (!targetMap.containsKey(target))
						targetMap.put(target, new HashMap<String, List<List<String>>>());

					if (!targetMap.get(target).containsKey(s))
						targetMap.get(target).put(s, new ArrayList<List<String>>());

					targetMap.get(target).get(s).add(group);
				}
			}
		}

		List<String> targets = new ArrayList<String>(targetMap.keySet());

		Collections.sort(targets, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return new Integer(targetMap.get(o2).size()).compareTo(targetMap.get(o1).size());
			}
		});

		for (int i = 0; i < 100; i++)
		{
			String target = targets.get(i);

			System.out.println("\ntarget = " + target);
			System.out.println("Upstream neighbors = " + graph.getUpstream(Collections.singleton(target)).size());
			System.out.println("Studies = " + targetMap.get(target).size());
			for (String s : targetMap.get(target).keySet())
			{
				System.out.print(s + ":");
				for (List<String> group : targetMap.get(target).get(s))
				{
					System.out.print(" " + group);
				}
				System.out.println();
			}
		}
	}

	private static void cropToMembersIfContains(Set<String> targets, List<String> group)
	{
		Set<String> s = new HashSet<String>(group);
		s.retainAll(targets);
		if (!s.isEmpty()) targets.retainAll(s);
	}

	private static void printCancerAssociations(Map<String, Set<String>> genesMap, List<String> recurrent)
	{
		for (String s : study)
		{
			if (!genesMap.containsKey(s)) continue;

			System.out.print(s + ": ");
			List<String> list = new ArrayList<String>(recurrent);
			list.retainAll(genesMap.get(s));
			System.out.println(list.toString());
		}
		System.out.println();
	}

	/**
	 * Generates a unique key for the pair, using the given two gene names.
	 */
	static String key(String s1, String s2)
	{
		if (s1.compareTo(s2) < 0) return s1 + " " + s2;
		else return s2 + " " + s1;
	}

	/**
	 * Reads the best cutoff score for the analysis results.
	 */
	static double getBestCutoff(String dir) throws FileNotFoundException
	{
		Scanner sc = new Scanner(new File(dir + "/fdr-guide.txt"));
		double bestCutoff = -Double.MAX_VALUE;
		double bestScore = -Double.MAX_VALUE;
		double bestFDR = 1;

		sc.nextLine();
		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			double cutoff = Double.parseDouble(token[0]);
			int size = Integer.parseInt(token[1]);
			double fdr = Double.parseDouble(token[2]);
			double tp = Double.parseDouble(token[3]);
			double fp = Double.parseDouble(token[4]);
			double score = Double.parseDouble(token[5]);

			if (bestScore < score)
			{
				bestScore = score;
				bestCutoff = cutoff;
				bestFDR = fdr;
			}
		}

		System.out.println("bestFDR = " + bestFDR);
		return bestCutoff;
	}

	/**
	 * Reads result groups.
	 */
	static List<List<String>> readResults(String dir, double cutoff) throws FileNotFoundException
	{
		List<List<String>> groups = new ArrayList<List<String>>();

		Scanner sc = new Scanner(new File(dir + "/ranked-groups.txt"));
		sc.nextLine();

		while (sc.hasNextLine())
		{
			String[] row = sc.nextLine().split("\t");
			double score = Double.parseDouble(row[0]);
			if (score > cutoff) break;

			List<String> group = new ArrayList<String>(Arrays.asList(row).subList(1, row.length));
			groups.add(group);
		}
		return groups;
	}

	/**
	 * Writes the sif file to visualize using ChiBE.
	 * @return recurrent genes
	 */
	private static List<String> prepareSIF(final Map<String, Integer> geneScore,
		Map<String, Integer> pairScore, String file, int recurrenceThr, int bindExistingThr,
		int minCol, int extraCol) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(file + ".sif"));

		Set<String> wroteGenes = new HashSet<String>();

		for (String pair : pairScore.keySet())
		{
			int score = pairScore.get(pair);
			if (score < recurrenceThr) continue;
			String[] g = pair.split(" ");
			writer.write(g[0] + "\tinteracts-with\t" + g[1] + "\n");
			wroteGenes.add(g[0]);
			wroteGenes.add(g[1]);
		}

		Set<String> strongGenes = new HashSet<String>(wroteGenes);

		for (String gene : geneScore.keySet())
		{
			if (!wroteGenes.contains(gene))
			{
				int score = geneScore.get(gene);
				if (score < recurrenceThr) continue;
				writer.write(gene + '\n');
				wroteGenes.add(gene);
			}
		}

		Set<String> extras = new HashSet<String>();
		for (String g1 : wroteGenes)
		{
			for (String g2 : wroteGenes)
			{
				if (g1.compareTo(g2) >= 0) continue;
				if (strongGenes.contains(g1) && strongGenes.contains(g2)) continue;

				String key = key(g1, g2);

				if (pairScore.containsKey(key) && pairScore.get(key) < recurrenceThr &&
					pairScore.get(key) >= bindExistingThr)
				{
					extras.add(key);
					key = key.replaceAll(" ", "\tinteracts-with\t");
					writer.write(key + "\n");
				}
			}
		}

		writer.close();

		int max = Math.max(getMax(geneScore), getMax(pairScore));
		System.out.println("\nmax value = " + max);

		writer = new BufferedWriter(new FileWriter(file + ".format"));
		writer.write("node\tall-nodes\tcolor\t255 255 255\n");
		writer.write("node\tall-nodes\tborderwidth\t2\n");
		writer.write("edge\tall-edges\twidth\t2\n");

		for (String gene : wroteGenes)
		{
			String col = getColor(geneScore.get(gene), recurrenceThr, max, minCol);
			writer.write("node\t" + gene + "\ttextcolor\t" + col + "\n");
			writer.write("node\t" + gene + "\tbordercolor\t" + col + "\n");
		}
		for (String pair : pairScore.keySet())
		{
			int score = pairScore.get(pair);
			if (score >= recurrenceThr)
			{
				String col = getColor(score, recurrenceThr, max, minCol);
				String edge = pair.replaceAll(" ", " interacts-with ");
				writer.write("edge\t" + edge + "\tcolor\t" + col + "\n");
			}
			else if (extras.contains(pair))
			{
				String col = extraCol + " " + extraCol + " " + extraCol;
				String edge = pair.replaceAll(" ", " interacts-with ");
				writer.write("edge\t" + edge + "\tcolor\t" + col + "\n");
				writer.write("edge\t" + edge + "\twidth\t1\n");
			}
		}
		writer.close();

		List<String> list = new ArrayList<String>(wroteGenes);
		Collections.sort(list);
		System.out.println("list.size() = " + list.size());
		for (String gene : list)
		{
			System.out.print(gene + ": ");
			Set<String> relatedCancers = GeneCards.getRelatedCancers(gene);
			for (String cancer : relatedCancers)
			{
				System.out.print(cancer + ", ");
			}
			System.out.println();
		}

		Collections.sort(list, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return geneScore.get(o2).compareTo(geneScore.get(o1));
			}
		});

		for (String gene : list)
		{
			System.out.println(gene + "\t" + geneScore.get(gene));
		}

		return list;
	}

	private static int getMax(Map<String, Integer> map)
	{
		int max = 0;
		for (Integer i : map.values())
		{
			if (max < i) max = i;
		}
		return max;
	}

	private static String getColor(int val, int minVal, int maxVal, int minCol)
	{
		assert val >= minVal;

		int dif = maxVal - minVal;
		double rat = (maxVal - val) / (double) dif;

		int col = (int) Math.round(minCol * rat);
		return col + " " + col + " " + col;
	}
}
