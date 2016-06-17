package org.cbio.mutex;

import java.io.*;
import java.util.*;

/**
 * This class is for converting one of the alteration matrix formats to Mutex format.
 * @author Ozgun Babur
 */
public class GAMConverter
{
	public static void convert(String inputFile, String outDir) throws IOException
	{
		Map<String, boolean[]> mutMap = new HashMap<String, boolean[]>();
		Map<String, boolean[]> ampMap = new HashMap<String, boolean[]>();
		Map<String, boolean[]> delMap = new HashMap<String, boolean[]>();

		Scanner sc = new Scanner(new File(inputFile));

		String line = sc.nextLine();
		String samples = line.substring(line.indexOf("\tTCGA"));
		int size = samples.substring(1).split("\t").length;
		System.out.println("size = " + size);

		while (sc.hasNextLine())
		{
			line = sc.nextLine();
			String[] token = line.split("\t");

			boolean[] b = convertToBoolean(line.substring(line.indexOf("\"\t") + 2));
			String[] genes = token[2].replaceAll("\"", "").split(",");

			Map<String, boolean[]> map;

			if (token[1].equals("HYPER_NODE_AMPLIFICATION")) map = ampMap;
			else if (token[1].equals("HYPER_NODE_DELETION")) map = delMap;
			else if (token[1].equals("HYPER_NODE_MUTATION")) map = mutMap;
			else throw new RuntimeException("Illegal type: " + token[1]);

			for (String gene : genes)
			{
				if (gene.isEmpty()) continue;

				map.put(gene, b);
			}
		}

		Set<String> genes = new HashSet<String>(mutMap.keySet());
		genes.addAll(ampMap.keySet());
		genes.addAll(delMap.keySet());

		new File(outDir).mkdirs();

		BufferedWriter writer = new BufferedWriter(new FileWriter(outDir + "DataMatrix.txt"));

		writer.write(samples);

		for (String gene : genes)
		{
			writer.write("\n" + gene);

			for (int i = 0; i < size; i++)
			{
				writer.write("\t");
				if (mutMap.containsKey(gene) && mutMap.get(gene)[i])
				{
					if (ampMap.containsKey(gene) && ampMap.get(gene)[i]) writer.write("4");
					else if (delMap.containsKey(gene) && delMap.get(gene)[i]) writer.write("5");
					else writer.write("1");
				}
				else
				{
					if (ampMap.containsKey(gene) && ampMap.get(gene)[i]) writer.write("2");
					else if (delMap.containsKey(gene) && delMap.get(gene)[i]) writer.write("3");
					else writer.write("0");
				}
			}
		}

		writer.close();
	}

	private static boolean[] convertToBoolean(String line)
	{
		String[] s = line.split("\t");
		boolean[] b = new boolean[s.length];
		for (int i = 0; i < b.length; i++)
		{
			b[i] = s[i].equals("1");
		}
		return b;
	}

	public static void main(String[] args) throws IOException
	{
		String dir = "/home/ozgun/Documents/TCGA/SARC/GAM/UPS/";
		convert("/home/ozgun/Documents/TCGA/SARC/GAM_gistic_specific/GAM_UPS.txt", dir);
	}
}
