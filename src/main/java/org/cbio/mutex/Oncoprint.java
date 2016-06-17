package org.cbio.mutex;

import java.io.*;
import java.util.*;

/**
 * This class generates SVG code for a given oncoprint.
 *
 * @author Ozgun Babur
 */
public class Oncoprint
{
	static int id = 0;
	public static void main(String[] args) throws IOException
	{
//		String s =  "MM................................................................................................\n" +
//					"..MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMD................................................................\n" +
//					".......................DDDDDDDDDDDMMDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD.......\n" +
//					"......................M....................................................................AAAAAAA";
//		write(s, "temp.svg");

		String dir = "/home/ozgun/Documents/TCGA/SARC/mutex/All/";
		generate((
			"SPTBN4\n" +
			"FRS2\n" +
			"EFNB3\n" +
			"TP53\n" +
			"MDM2\n" +
			"CTDSP2").split("\n"), dir + "DataMatrix.txt", dir + "oncoprint/story3.svg", false);
//		generateAllOncoprints(dir);
	}

	public static void generateAllOncoprints(String dir, boolean includeUnaltered) throws IOException
	{
		new File(dir + "/oncoprint").mkdirs();
		Scanner sc = new Scanner(new File(dir + "/oncoprint.txt"));
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			if (line.startsWith("["))
			{
				String s = line.substring(0, line.indexOf("\t"));
				s = s.substring(1, s.length() - 1);
				String[] genes = s.split(" ");
				generate(genes, dir + "DataMatrix.txt", dir + "oncoprint/" +
					Arrays.toString(genes) + ".svg", includeUnaltered);
			}
		}
	}

	public static void generate(String[] genes, String matrixFile, String outFile,
		boolean includeUnaltered) throws IOException
	{
		int[][] matrix = readFromMatrixFile(genes, matrixFile);
		List<Integer> order = getPrintOrdering(matrix);
		String oncoprint = convertMatrixToString(matrix, order, includeUnaltered);
		write(oncoprint, outFile);
	}

	private static void write(String oncoString, String filename) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));
		writer.write(header);

		int row = 0;
		String[] lines = oncoString.split("\n");
		for (String letters : lines)
		{
			for (int i = 0; i < letters.length(); i++)
			{
				String let = letters.substring(i, i+1);

				placeRect(i, row, let, writer);
			}
			row++;
		}

		writer.write("</svg>");
		writer.close();
	}

	private static void placeRect(int col, int row, String type, Writer writer) throws IOException
	{
		double height = 23;
		double width = 2;

		double x = col * (width);
		double y = row * (height + 4);

		String style =
			type.equals(".") || type.equals("M") ? "fill:#d3d3d3" :
			type.equals("A") || type.equals("B") ? "fill:#ff0000":
			type.equals("D") || type.equals("E") ? "fill:#0000ff": null;

		String s =
			"    <rect\n" +
			"       style=\"" + style + "\"\n" +
			"       x=\"" + x + "\"\n" +
			"       height=\"" + height + "\"\n" +
			"       width=\"" + width + "\"\n" +
			"       y=\"" + y + "\"\n" +
			"       id=\"rect" + (++id) + "\" />\n";

		writer.write(s);

		if (type.equals("M") || type.equals("B") || type.equals("E"))
		{
			s = "    <rect\n" +
				"       style=\"fill:#008000\"\n" +
				"       x=\"" + x + "\"\n" +
				"       height=\"7.66\"\n" +
				"       width=\"" + width + "\"\n" +
				"       y=\"" + (y + 7.67) + "\"\n" +
				"       id=\"rect" + (++id) + "\" />\n";
			writer.write(s);
		}

	}

	private static String convertMatrixToString(int[][] matrix, List<Integer> order,
		boolean includeUnaltered)
	{
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < matrix.length; i++)
		{
			if (i > 0) sb.append("\n");

			for (int j = 0; j < order.size(); j++)
			{
				int sample = order.get(j);
				if (!includeUnaltered)
				{
					boolean b = false;
					for (int k = 0; k < matrix.length; k++)
					{
						if (matrix[k][sample] != 0)
						{
							b = true;
							break;
						}
					}
					if (!b) break;
				}

				sb.append(letterCode(matrix[i][sample]));
			}
		}
		return sb.toString();
	}

	private static String letterCode(int alt)
	{
		switch (alt)
		{
			case 0: return ".";
			case 1: return "M";
			case 2: return "A";
			case 3: return "D";
			case 4: return "B";
			case 5: return "E";
			default: throw new IllegalArgumentException("Illegal alteration value: " + alt);
		}
	}

	private static int[][] readFromMatrixFile(String[] genes, String file) throws FileNotFoundException
	{
		Set<String> set = new HashSet<>(Arrays.asList(genes));
		int[][] x = new int[genes.length][];
		Scanner sc = new Scanner(new File(file));
		sc.nextLine();

		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			String sym = line.substring(0, line.indexOf("\t"));
			if (set.contains(sym))
			{
				int[] alts = convert(line.substring(line.indexOf("\t") + 1));

				for (int i = 0; i < genes.length; i++)
				{
					if (genes[i].equals(sym)) x[i] = alts;
				}
			}
		}
		return x;
	}

	/**
	 * Gets an ordering for the samples to make the oncoprint look nicer.
	 * @return sample ordering for printing oncoprint
	 */
	private static List<Integer> getPrintOrdering(int[][] matrix)
	{
		List<Integer> order = new ArrayList<>();

		int sampleSize = matrix[0].length;
		final int geneSize = matrix.length;

		for (int i = 0; i < sampleSize; i++)
		{
			order.add(i);
		}

		final boolean[][] marks = new boolean[sampleSize][];

		for (int i = 0; i < sampleSize; i++)
		{
			marks[i] = alterationMarks(matrix, i);
		}

		final boolean[][] mut = new boolean[geneSize][];
		final boolean[][] cna = new boolean[geneSize][];

		for (int i = 0; i < geneSize; i++)
		{
			mut[i] = getMutated(matrix[i]);
			cna[i] = getCNAltered(matrix[i]);
		}

		Collections.sort(order, (o1, o2) ->
		{
			boolean[] m1 = marks[o1];
			boolean[] m2 = marks[o2];

			int c = 0;
			for (int i = 0; i < geneSize; i++)
			{
				if (m1[i] && !m2[i]) c = -1;
				if (!m1[i] && m2[i]) c = 1;
				if (c != 0) break;
			}

			if (c != 0)
			{
				if (getNumberOfInitialPositiveAltOverlap(m1, m2) % 2 == 1) return -c;
				else return c;
			}

			for (int i = 0; i < geneSize; i++)
			{
				if (mut[i][o1] && !mut[i][o2]) return -1;
				if (!mut[i][o1] && mut[i][o2]) return 1;
				if (cna[i][o1] && !cna[i][o2]) return 1;
				if (!cna[i][o1] && cna[i][o2]) return -1;
			}

			return 0;
		});

		return order;
	}

	private static boolean[] alterationMarks(int[][] matrix, int sample)
	{
		boolean[] b = new boolean[matrix.length];
		for (int i = 0; i < b.length; i++)
		{
			b[i] = matrix[i][sample] != 0;
		}
		return b;
	}

	private static int getNumberOfInitialPositiveAltOverlap(boolean[] m1, boolean[] m2)
	{
		int x = 0;

		for (int i = 0; i < m1.length; i++)
		{
			if (!m1[i] && !m2[i] && x == 0) continue;

			if (m1[i] && m2[i]) x++;
			else break;
		}
		return x;
	}

	private static boolean[] getMutated(int[] alts)
	{
		boolean[] m = new boolean[alts.length];
		for (int i = 0; i < alts.length; i++)
		{
			m[i] = alts[i] == 1 || alts[i] == 4 || alts[i] == 5;
		}
		return m;
	}

	private static boolean[] getCNAltered(int[] alts)
	{
		boolean[] c = new boolean[alts.length];
		for (int i = 0; i < alts.length; i++)
		{
			c[i] = alts[i] == 2 || alts[i] == 3 || alts[i] == 4 || alts[i] == 5;
		}
		return c;
	}

	private static int[] convert(String altStr)
	{
		String[] token = altStr.split("\t");
		int[] alts = new int[token.length];
		for (int i = 0; i < alts.length; i++)
		{
			alts[i] = Integer.parseInt(token[i]);
		}
		return alts;
	}

	public static final String header = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n" +
		"<svg\n" +
		"   xmlns:dc=\"http://purl.org/dc/elements/1.1/\"\n" +
		"   xmlns:cc=\"http://creativecommons.org/ns#\"\n" +
		"   xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"\n" +
		"   xmlns:svg=\"http://www.w3.org/2000/svg\"\n" +
		"   xmlns=\"http://www.w3.org/2000/svg\"\n" +
		"   xmlns:sodipodi=\"http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd\"\n" +
		"   xmlns:inkscape=\"http://www.inkscape.org/namespaces/inkscape\"\n" +
		"   version=\"1.1\"\n" +
		"   height=\"81\"\n" +
		"   width=\"570.16998\"\n" +
		"   id=\"svg2\"\n" +
		"   inkscape:version=\"0.48.3.1 r9886\"\n" +
		"   sodipodi:docname=\"temp.svg\">\n" +
		"  <metadata\n" +
		"     id=\"metadata1946\">\n" +
		"    <rdf:RDF>\n" +
		"      <cc:Work\n" +
		"         rdf:about=\"\">\n" +
		"        <dc:format>image/svg+xml</dc:format>\n" +
		"        <dc:type\n" +
		"           rdf:resource=\"http://purl.org/dc/dcmitype/StillImage\" />\n" +
		"        <dc:title></dc:title>\n" +
		"      </cc:Work>\n" +
		"    </rdf:RDF>\n" +
		"  </metadata>\n" +
		"  <defs\n" +
		"     id=\"defs1944\">\n" +
		"    <marker\n" +
		"       inkscape:stockid=\"StopL\"\n" +
		"       orient=\"auto\"\n" +
		"       refY=\"0\"\n" +
		"       refX=\"0\"\n" +
		"       id=\"StopL\"\n" +
		"       style=\"overflow:visible\">\n" +
		"      <path\n" +
		"         id=\"path4366\"\n" +
		"         d=\"M 0,5.65 0,-5.65\"\n" +
		"         style=\"fill:none;stroke:#000000;stroke-width:1pt\"\n" +
		"         transform=\"scale(0.8,0.8)\"\n" +
		"         inkscape:connector-curvature=\"0\" />\n" +
		"    </marker>\n" +
		"    <marker\n" +
		"       style=\"overflow:visible\"\n" +
		"       id=\"DistanceEnd\"\n" +
		"       refX=\"0\"\n" +
		"       refY=\"0\"\n" +
		"       orient=\"auto\"\n" +
		"       inkscape:stockid=\"DistanceEnd\">\n" +
		"      <g\n" +
		"         id=\"g2301\">\n" +
		"        <path\n" +
		"           style=\"fill:none;stroke:#ffffff;stroke-width:1.14999998;stroke-linecap:square\"\n" +
		"           d=\"M 0,0 -2,0\"\n" +
		"           id=\"path2316\"\n" +
		"           inkscape:connector-curvature=\"0\" />\n" +
		"        <path\n" +
		"           style=\"fill:#000000;fill-rule:evenodd;stroke:none\"\n" +
		"           d=\"M 0,0 -13,4 -9,0 -13,-4 0,0 z\"\n" +
		"           id=\"path2312\"\n" +
		"           inkscape:connector-curvature=\"0\" />\n" +
		"        <path\n" +
		"           style=\"fill:none;stroke:#000000;stroke-width:1;stroke-linecap:square\"\n" +
		"           d=\"M 0,-4 0,40\"\n" +
		"           id=\"path2314\"\n" +
		"           inkscape:connector-curvature=\"0\" />\n" +
		"      </g>\n" +
		"    </marker>\n" +
		"    <inkscape:path-effect\n" +
		"       effect=\"spiro\"\n" +
		"       id=\"path-effect3423\"\n" +
		"       is_visible=\"true\" />\n" +
		"    <marker\n" +
		"       inkscape:stockid=\"StopL\"\n" +
		"       orient=\"auto\"\n" +
		"       refY=\"0\"\n" +
		"       refX=\"0\"\n" +
		"       id=\"StopL-6\"\n" +
		"       style=\"overflow:visible\">\n" +
		"      <path\n" +
		"         inkscape:connector-curvature=\"0\"\n" +
		"         id=\"path4366-9\"\n" +
		"         d=\"M 0,5.65 0,-5.65\"\n" +
		"         style=\"fill:none;stroke:#000000;stroke-width:1pt\"\n" +
		"         transform=\"scale(0.8,0.8)\" />\n" +
		"    </marker>\n" +
		"    <inkscape:path-effect\n" +
		"       effect=\"spiro\"\n" +
		"       id=\"path-effect3423-0\"\n" +
		"       is_visible=\"true\" />\n" +
		"  </defs>\n" +
		"  <sodipodi:namedview\n" +
		"     pagecolor=\"#ffffff\"\n" +
		"     bordercolor=\"#666666\"\n" +
		"     borderopacity=\"1\"\n" +
		"     objecttolerance=\"10\"\n" +
		"     gridtolerance=\"10\"\n" +
		"     guidetolerance=\"10\"\n" +
		"     inkscape:pageopacity=\"0\"\n" +
		"     inkscape:pageshadow=\"2\"\n" +
		"     inkscape:window-width=\"1855\"\n" +
		"     inkscape:window-height=\"1176\"\n" +
		"     id=\"namedview1942\"\n" +
		"     showgrid=\"false\"\n" +
		"     fit-margin-top=\"0\"\n" +
		"     fit-margin-left=\"0\"\n" +
		"     fit-margin-right=\"0\"\n" +
		"     fit-margin-bottom=\"0\"\n" +
		"     inkscape:snap-global=\"false\"\n" +
		"     inkscape:zoom=\"1.56829\"\n" +
		"     inkscape:cx=\"199.74182\"\n" +
		"     inkscape:cy=\"62.318951\"\n" +
		"     inkscape:window-x=\"65\"\n" +
		"     inkscape:window-y=\"24\"\n" +
		"     inkscape:window-maximized=\"1\"\n" +
		"     inkscape:current-layer=\"svg2\" />\n";
}
