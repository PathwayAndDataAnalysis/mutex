package org.cbio.mutex;

import org.panda.utility.ArrayUtil;
import org.panda.utility.statistics.Summary;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

/**
 * @author Ozgun Babur
 */
public class HyperAlteredRemover
{
	public static void cleanHyper(String inFile, String outFile) throws IOException
	{
		Map<String, int[]> geneMap = new HashMap<>();

		Scanner sc = new Scanner(new File(inFile));
		String line = sc.nextLine();
		String[] sample = line.substring(line.indexOf("\t") + 1).split("\t");

		while (sc.hasNextLine())
		{
			line = sc.nextLine();
			String gene = line.substring(0, line.indexOf("\t"));
			int[] vals = convert(line.substring(line.indexOf("\t") + 1));
			geneMap.put(gene, vals);
		}

		double[] cnt = new double[sample.length];
		Arrays.fill(cnt, 0);

		for (String gene : geneMap.keySet())
		{
			int[] vals = geneMap.get(gene);
			for (int i = 0; i < vals.length; i++)
			{
				if (vals[i] > 0) cnt[i]++;
			}
		}

		boolean[] hyper = Summary.markOutliers(cnt, true);

		System.out.println("Hyper altered samples = " + ArrayUtil.countValue(hyper, true));

		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));

		for (int i = 0; i < sample.length; i++)
		{
			if (!hyper[i]) writer.write("\t" + sample[i]);
		}

		for (String gene : geneMap.keySet())
		{
			writer.write("\n" + gene);
			int[] vals = geneMap.get(gene);

			for (int i = 0; i < sample.length; i++)
			{
				if (!hyper[i]) writer.write("\t" + vals[i]);
			}
		}

		writer.close();


	}

	private static int[] convert(String valString)
	{
		String[] s = valString.split("\t");
		int[] vals = new int[s.length];
		for (int i = 0; i < vals.length; i++)
		{
			vals[i] = Integer.parseInt(s[i]);
		}
		return vals;
	}

	public static void main(String[] args) throws IOException
	{
		String dir = "/home/ozgun/Documents/TCGA/SARC/mutex/UPS/";
		cleanHyper(dir + "DataMatrix.txt", dir + "DataMatrix-clean.txt");
	}
}
