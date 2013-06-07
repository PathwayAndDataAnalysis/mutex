package org.cbio.mutex;

import org.cbio.causality.data.portal.CBioPortalAccessor;
import org.cbio.causality.data.portal.CancerStudy;
import org.cbio.causality.data.portal.CaseList;
import org.cbio.causality.data.portal.GeneticProfile;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;

/**
 * Reads a dataset from cBio Portal. Uses local files as cache.
 * @author Ozgun Babur
 */
public class PortalReader
{
	public Map<String, AlterationPack> readAlterations(PortalDataset data, Set<String> syms) throws IOException
	{
		// cBio portal configuration

		System.out.println();

		if (new File(data.filename).exists())
		{
			return AlterationPack.readFromFile(data.filename);
		}
		else if (getClass().getResourceAsStream(data.filename) != null)
		{
			return AlterationPack.readFromFile(new InputStreamReader(getClass().getResourceAsStream(
				data.filename)));
		}

		CBioPortalAccessor cBioPortalAccessor = new CBioPortalAccessor();
		CancerStudy cancerStudy = cBioPortalAccessor.getCancerStudies().get(data.study); // GBM
		cBioPortalAccessor.setCurrentCancerStudy(cancerStudy);

		List<GeneticProfile> geneticProfilesForCurrentStudy =
			cBioPortalAccessor.getGeneticProfilesForCurrentStudy();
		List<GeneticProfile> gp = new ArrayList<GeneticProfile>();
		for (int prof : data.profile)
		{
			gp.add(geneticProfilesForCurrentStudy.get(prof));
		}
		cBioPortalAccessor.setCurrentGeneticProfiles(gp);

		List<CaseList> caseLists = cBioPortalAccessor.getCaseListsForCurrentStudy();
		cBioPortalAccessor.setCurrentCaseList(caseLists.get(data.caseList));

		System.out.println("syms.size() = " + syms.size());
		long time = System.currentTimeMillis();

		Map<String, AlterationPack> map = new HashMap<String, AlterationPack>();

		for (String sym : syms)
		{
			AlterationPack alt = cBioPortalAccessor.getAlterations(sym);
			if (alt.isAltered(Alteration.GENOMIC))
			{
				map.put(sym, alt);
			}
		}
		System.out.println("map.size() = " + map.size());
		time = System.currentTimeMillis() - time;
		System.out.println("read in " + (time / 1000D) + " seconds");

		AlterationPack.writeToFile(map, data.filename);

		return map;
	}
}
