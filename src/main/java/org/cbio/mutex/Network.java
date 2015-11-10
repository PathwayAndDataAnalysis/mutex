package org.cbio.mutex;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.cbio.causality.analysis.Graph;
import org.cbio.causality.analysis.SIFLinker;
import org.cbio.causality.network.PathwayCommons;
import org.cbio.causality.network.SPIKE;
import org.cbio.causality.network.SignaLink;

import java.io.*;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class Network extends Graph
{
	SIFLinker linker;

	public Network(String filename) throws FileNotFoundException
	{
		super("directed network", "is-upstream-of");
		Graph graphSig = loadPTRGraph(filename);
		Graph graphTR = loadTRGraph(filename);
		merge(graphSig);
		merge(graphTR);
		load(new FileInputStream(filename), Collections.<String>emptySet(),
			new HashSet<String>(Arrays.asList("is-upstream-of")));

		linker = new SIFLinker();
		linker.load(this);
		linker.load(graphSig);
		linker.load(graphTR);
	}

	public Network() throws FileNotFoundException
	{
		super("directed network", "is-upstream-of");
		Graph graphSig = loadPTRGraph(null);
		Graph graphTR = loadTRGraph(null);
		merge(graphSig);
		merge(graphTR);
		linker = new SIFLinker();
		linker.load(graphSig);
		linker.load(graphTR);
	}

	private static Graph loadTRGraph(String filename) throws FileNotFoundException
	{
		Graph graphTR = new Graph("transcriptional regulation", SIFEnum.CONTROLS_EXPRESSION_OF.getTag());

		if (filename == null)
		{
			graphTR.merge(PathwayCommons.getGraph(SIFEnum.CONTROLS_EXPRESSION_OF));
			graphTR.merge(SPIKE.getGraphTR());
			graphTR.merge(SignaLink.getGraphTR());
		}
		else
		{
			graphTR.load(new FileInputStream(filename), Collections.<String>emptySet(),
				new HashSet<String>(Arrays.asList(SIFEnum.CONTROLS_EXPRESSION_OF.getTag())));
		}
		return graphTR;
	}

	private static Graph loadPTRGraph(String filename) throws FileNotFoundException
	{
		Graph graphSig = new Graph("signaling", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		if (filename == null)
		{
			graphSig.merge(PathwayCommons.getGraph(SIFEnum.CONTROLS_STATE_CHANGE_OF));
			graphSig.merge(SPIKE.getGraphPostTl());
			graphSig.merge(SignaLink.getGraphPostTl());
		}
		else
		{
			graphSig.load(new FileInputStream(filename), Collections.<String>emptySet(),
				new HashSet<String>(Arrays.asList(SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag())));
		}
		return graphSig;
	}

	public static void main(String[] args) throws FileNotFoundException
	{
		try
		{
			mergeNetworkFiles();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		System.exit(0);
		Network n = new Network();
		final Map<String, Integer> degreeMap = new HashMap<String, Integer>();
		for (String sym : n.getSymbols())
		{
			degreeMap.put(sym, n.getUpstream(sym).size());
		}
		List<String> genes = new ArrayList<String>(degreeMap.keySet());
		Collections.sort(genes, new Comparator<String>()
		{
			@Override
			public int compare(String o1, String o2)
			{
				return degreeMap.get(o2).compareTo(degreeMap.get(o1));
			}
		});

		for (int i = 0; i < 100; i++)
		{
			String gene = genes.get(i);
			System.out.println((i+ 1) + "\t" + degreeMap.get(gene) + "\t" + gene);
		}
	}

	private static void mergeNetworkFiles() throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter("Network.sif"));
		Scanner sc = new Scanner(new File("Network-PTM.sif"));
		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			writer.write(token[0] + "\t" + token[1] + "\t" + token[2] + "\n");
		}
		sc = new Scanner(new File("Network-TR.sif"));
		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			writer.write(token[0] + "\t" + token[1] + "\t" + token[2] + "\n");
		}
		writer.close();
	}
}
