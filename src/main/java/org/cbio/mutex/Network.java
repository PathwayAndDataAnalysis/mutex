package org.cbio.mutex;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.resource.network.PathwayCommons;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;
import org.panda.utility.graph.SIFLinker;

import java.io.*;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class Network extends DirectedGraph
{
	SIFLinker linker;

	public Network(String filename) throws FileNotFoundException
	{
		super("directed network", "is-upstream-of");
		linker = new SIFLinker();
		addResource(filename);
	}

	public Network() throws FileNotFoundException
	{
		this(null);
	}

	public void addResource(String filename) throws FileNotFoundException
	{
		linker.load(loadPTRGraph(filename));
		linker.load(loadTRGraph(filename));
		linker.load(loadIsUpstreamOfGraph(filename));
		merge(linker.graph);
	}

	private static DirectedGraph loadTRGraph(String filename) throws FileNotFoundException
	{
		DirectedGraph graphTR = new DirectedGraph("transcriptional regulation", SIFEnum.CONTROLS_EXPRESSION_OF.getTag());

		if (filename == null)
		{
			graphTR.merge(PathwayCommons.get().getGraph(SIFEnum.CONTROLS_EXPRESSION_OF));
//			graphTR.merge(SPIKE.getGraphTR());
//			graphTR.merge(SignaLink.getGraphTR());
		}
		else
		{
			graphTR.load(new FileInputStream(filename), Collections.singleton(SIFEnum.CONTROLS_EXPRESSION_OF.getTag()));
		}
		return graphTR;
	}

	private static DirectedGraph loadPTRGraph(String filename) throws FileNotFoundException
	{
		DirectedGraph graphSig = new DirectedGraph("signaling", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		if (filename == null)
		{
			graphSig.merge(PathwayCommons.get().getGraph(SIFEnum.CONTROLS_STATE_CHANGE_OF));
//			graphSig.merge(SPIKE.getGraphPostTl());
//			graphSig.merge(SignaLink.getGraphPostTl());
		}
		else
		{
			graphSig.load(new FileInputStream(filename), Collections.singleton(SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag()));
		}
		return graphSig;
	}

	private static DirectedGraph loadIsUpstreamOfGraph(String filename) throws FileNotFoundException
	{
		DirectedGraph graph = new DirectedGraph("signaling relations", "is-upstream-of");

		if (filename != null) graph.load(new FileInputStream(filename), Collections.singleton("is-upstream-of"));

		return graph;
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
		final Map<String, Integer> degreeMap = new HashMap<>();
		for (String sym : n.getSymbols())
		{
			degreeMap.put(sym, n.getUpstream(sym).size());
		}
		List<String> genes = new ArrayList<String>(degreeMap.keySet());
		Collections.sort(genes, (o1, o2) -> degreeMap.get(o2).compareTo(degreeMap.get(o1)));

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
