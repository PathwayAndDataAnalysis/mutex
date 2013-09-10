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

		if (new File(data.filename).exists())
		{
			return AlterationPack.readFromFile(data.filename,
				Alteration.MUTATION, Alteration.COPY_NUMBER);
		}
		else if (getClass().getResourceAsStream(data.filename) != null)
		{
			return AlterationPack.readFromFile(new InputStreamReader(getClass().getResourceAsStream(
				data.filename)), Alteration.MUTATION, Alteration.COPY_NUMBER);
		}

		CBioPortalAccessor cBioPortalAccessor = new CBioPortalAccessor();
		CancerStudy cancerStudy = findCancerStudy(cBioPortalAccessor.getCancerStudies(), data.study); // GBM
		cBioPortalAccessor.setCurrentCancerStudy(cancerStudy);

		List<GeneticProfile> geneticProfilesForCurrentStudy =
			cBioPortalAccessor.getGeneticProfilesForCurrentStudy();
		List<GeneticProfile> gp = new ArrayList<GeneticProfile>();
		for (String prof : data.profile)
		{
			gp.add(findProfile(geneticProfilesForCurrentStudy, prof));
		}
		cBioPortalAccessor.setCurrentGeneticProfiles(gp);

		List<CaseList> caseLists = cBioPortalAccessor.getCaseListsForCurrentStudy();
		cBioPortalAccessor.setCurrentCaseList(findCaseList(caseLists, data.caseList));

		System.out.println("syms.size() = " + syms.size());
		long time = System.currentTimeMillis();

		Map<String, AlterationPack> map = new HashMap<String, AlterationPack>();

		for (String sym : syms)
		{
			AlterationPack alt = cBioPortalAccessor.getAlterations(sym);
			if (alt != null && alt.isAltered(Alteration.GENOMIC))
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

	private CancerStudy findCancerStudy(List<CancerStudy> list, String id)
	{
		for (CancerStudy study : list)
		{
			if (study.getStudyId().equals(id)) return study;
		}
		return null;
	}

	private CaseList findCaseList(List<CaseList> list, String id)
	{
		for (CaseList cl : list)
		{
			if (cl.getId().equals(id)) return cl;
		}
		return null;
	}

	private GeneticProfile findProfile(List<GeneticProfile> list, String id)
	{
		for (GeneticProfile profile : list)
		{
			if (profile.getId().equals(id)) return profile;
		}
		return null;
	}
}
