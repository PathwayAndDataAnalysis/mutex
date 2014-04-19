package org.cbio.mutex.mutation;

import org.cbio.causality.data.portal.CBioPortalAccessor;
import org.cbio.causality.data.portal.PortalDataset;
import org.cbio.causality.idmapping.Length;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.network.HPRD;
import org.cbio.causality.util.Histogram;

import java.io.IOException;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class Main
{
	public static void main(String[] args) throws IOException
	{
		PortalDataset dataset = PortalDataset.PROSTATE_TCGA_MUT_CN;

		CBioPortalAccessor cacc = new CBioPortalAccessor(dataset);
		CBioPortalAccessor.setCacheDir("/home/ozgun/Projects/causality/portal-cache");

		List<AlterationPack> packs = new ArrayList<AlterationPack>();

		Set<String> symbols = HPRD.getAllSymbols();
		symbols.retainAll(Length.getSymbols());
		System.out.println("Length.getSymbols() = " + Length.getSymbols().size());

		for (String symbol : symbols)
		{
			AlterationPack pack = cacc.getAlterations(symbol);
			if (pack != null && pack.isAltered(Alteration.MUTATION))
			{
				packs.add(pack);
			}
		}
		Map<String, String[]> mutMap = MutCount.prepareMutMap(packs, cacc);

		symbols = new HashSet<String>();
		for (AlterationPack pack : packs)
		{
			symbols.add(pack.getId());
		}

		Set<String>[] lost = HetLossHelper.collectNonExistingGenes(symbols, cacc);

		MutAndCopyNumberRelater.printPlot(symbols, cacc);
		if (true) return;

//		MutCount.countAvgLengths(mutMap);
//		SampleRelations.printSampleDistAsSif(packs);


		Histogram his = new Histogram(1);
		int times = 0;
		for (int i = 0; i < times; i++)
		{
			int[] mutCnts = MutCount.getSampleMutCounts(mutMap);
			Map<String, String[]> mutMap2 = RandomAltGenerator.generate(new ArrayList<String>(symbols), mutCnts, lost);
			List<AlterationPack> packs2 = RandomAltGenerator.prepareAltpAcksFromMutMap(mutMap2);
//			Histogram h = PairwiseStat.getAlterationDistribution(packs2, mutMap2, lost);
			Histogram h = PairwiseStat.getOverlapDistribution(packs2, mutMap2, lost);
			his.add(h);
		}

//		Histogram h = PairwiseStat.getAlterationDistribution(packs, mutMap, lost);
		Histogram h = PairwiseStat.getOverlapDistribution(packs, mutMap, lost);
		h.print();
//		h.printTogether(his, 1D / times);
	}
}
