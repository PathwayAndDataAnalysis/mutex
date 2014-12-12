package org.cbio.mutex;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.cbio.causality.analysis.Graph;
import org.cbio.causality.analysis.SIFLinker;
import org.cbio.causality.network.PathwayCommons;
import org.cbio.causality.network.SPIKE;
import org.cbio.causality.network.SignaLink;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;

/**
 * @author Ozgun Babur
 */
public class Network extends Graph
{
	SIFLinker linker;

	public Network(String filename) throws FileNotFoundException
	{
		super("directed network", "is-upstream-of");
		load(new FileInputStream(filename), Collections.<String>emptySet(),
			new HashSet<String>(Arrays.asList(
				SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag(),
				SIFEnum.CONTROLS_EXPRESSION_OF.getTag())));

		linker = new SIFLinker();
		linker.load(this);
	}

	public Network()
	{
		super("directed network", "is-upstream-of");
		Graph graphSig = loadPTRGraph();
		Graph graphTR = loadTRGraph();
		merge(graphSig);
		merge(graphTR);
		linker = new SIFLinker();
		linker.load(graphSig);
		linker.load(graphTR);
	}

	private static Graph loadTRGraph()
	{
		Graph graphTR = new Graph("transcriptional regulation", SIFEnum.CONTROLS_EXPRESSION_OF.getTag());
		graphTR.merge(PathwayCommons.getGraph(SIFEnum.CONTROLS_EXPRESSION_OF));
		graphTR.merge(SPIKE.getGraphTR());
		graphTR.merge(SignaLink.getGraphTR());
		return graphTR;
	}

	private static Graph loadPTRGraph()
	{
		Graph graphSig = new Graph("signaling", SIFEnum.CONTROLS_STATE_CHANGE_OF.getTag());
		graphSig.merge(PathwayCommons.getGraph(SIFEnum.CONTROLS_STATE_CHANGE_OF));
		graphSig.merge(SPIKE.getGraphPostTl());
		graphSig.merge(SignaLink.getGraphPostTl());
		return graphSig;
	}

}
