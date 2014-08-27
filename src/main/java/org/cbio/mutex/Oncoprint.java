package org.cbio.mutex;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;

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
		String s =
					"MM................................................................................................\n" +
					"..MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMD................................................................\n" +
					".......................DDDDDDDDDDDMMDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD.......\n" +
					"......................M....................................................................AAAAAAA";

		BufferedWriter writer = new BufferedWriter(new FileWriter("temp.svg"));
		writer.write(header);

		int row = 0;
		String[] lines = s.split("\n");
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
