package org.cbio.mutex.mutation;

import org.cbio.causality.data.portal.CBioPortalAccessor;
import org.cbio.causality.data.portal.PortalDataset;
import org.cbio.causality.hprd.HPRD;
import org.cbio.causality.idmapping.Length;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;

import java.io.IOException;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class Main
{
	public static void main(String[] args) throws IOException
	{
		PortalDataset dataset = PortalDataset.ENDOMETRIAL_MUT;

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

//		symbols = new HashSet<String>();
//		for (AlterationPack pack : packs)
//		{
//			symbols.add(pack.getId());
//		}

//		MutCount.countAvgLengths(mutMap);
//		SampleRelations.printSampleDistAsSif(packs);

//		PairwiseStat.printOverlapStats(packs, mutMap);
		packs = PairwiseStat.filterToAlterationCounts(packs, mutMap);
		mutMap = MutCount.prepareMutMap(packs, cacc);
		PairwiseStat.printOverlapStats(packs, mutMap);
	}
}
