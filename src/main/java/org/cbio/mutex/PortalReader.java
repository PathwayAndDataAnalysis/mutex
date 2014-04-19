package org.cbio.mutex;

import org.cbio.causality.analysis.CNVerifier;
import org.cbio.causality.data.portal.*;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;
import org.cbio.causality.util.ArrayUtil;
import org.cbio.causality.util.Histogram;
import org.cbio.causality.util.StudentsT;
import org.cbio.causality.util.Summary;

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
	CBioPortalAccessor accessor;
	ExpDataManager expMan;
	CNVerifier cnVerifier;

	public Map<String, AlterationPack> readAlterations(PortalDataset data, Collection<String> syms)
		throws IOException
	{
		accessor = getPortalAccessor(data);
		expMan = getExpDataMan(data, accessor);
		cnVerifier = new CNVerifier(expMan, 0.05);

		// cBio portal configuration

//		String filename = data.study + ".txt";
//
//		if (new File(filename).exists())
//		{
//			return AlterationPack.readFromFile(filename);
//		}
//		else if (getClass().getResourceAsStream(filename) != null)
//		{
//			return AlterationPack.readFromFile(new InputStreamReader(getClass().getResourceAsStream(
//				filename)));
//		}

		System.out.println("syms.size() = " + syms.size());
		long time = System.currentTimeMillis();

		Map<String, AlterationPack> map = new HashMap<String, AlterationPack>();

//		int[] cnts = null;
		List<String> sorted = new ArrayList<String>(syms);
		Collections.sort(sorted);
		for (String sym : sorted)
		{
			AlterationPack alt = accessor.getAlterations(sym);

			cnVerifier.verify(alt);

			if (alt != null && alt.get(Alteration.COPY_NUMBER) != null &&
				alt.get(Alteration.MUTATION) != null)
			{
				removeMinorCopyNumberAlts(alt);
				alt.complete(Alteration.GENOMIC);

				if (alt.isAltered(Alteration.GENOMIC)) map.put(sym, alt);

//				if (cnts == null) cnts = new int[alt.getSize()];
//				for (int i = 0; i < alt.getSize(); i++)
//				{
//					if (alt.getChange(Alteration.GENOMIC, i).isAltered()) cnts[i]++;
//				}
			}
		}
		System.out.println("map.size() = " + map.size());
		time = System.currentTimeMillis() - time;
		System.out.println("read in " + (time / 1000D) + " seconds");

//		AlterationPack.writeToFile(map, filename);

//		System.out.println("Median altered gene count = " + Summary.median(cnts));
//		Histogram h = new Histogram(50);
//		h.setBorderAtZero(true);
//		for (int i = 0; i < cnts.length; i++)
//		{
//			h.count(cnts[i]);
//		}
//		h.print();

		Histogram h = new Histogram(1);
		h.setBorderAtZero(true);
		for (AlterationPack pack : map.values())
		{
			h.count(Math.log(pack.getAlteredRatio(Alteration.GENOMIC)) / Math.log(2));
		}
		h.print();

		return map;
	}

	public void removeMinorCopyNumberAlts(AlterationPack gene)
	{
		if (gene.get(Alteration.COPY_NUMBER) == null) return;

		int up = 0;
		int dw = 0;
		for (Change c : gene.get(Alteration.COPY_NUMBER))
		{
			if (c == Change.ACTIVATING) up++;
			else if (c == Change.INHIBITING) dw++;
		}

		if (up == 0 && dw == 0) return;

		if (up == dw)
		{
			System.out.println("Gene " + gene.getId() + " is equally altered (up: " + up +
				", dw: " + dw + "). Choosing none");
		}

		if (up > 0 && dw > 0)
		{
			Change c = dw < up ? Change.ACTIVATING : dw > up ? Change.INHIBITING : null;

			Change[] changes = gene.get(Alteration.COPY_NUMBER);
			for (int i = 0; i < gene.getSize(); i++)
			{
				changes[i] = changes[i] == c ? c : Change.NO_CHANGE;
			}
		}
	}


	public static CBioPortalAccessor getPortalAccessor(PortalDataset data) throws IOException
	{
		CBioPortalAccessor cBioPortalAccessor = new CBioPortalAccessor();

		CancerStudy cancerStudy = findCancerStudy(cBioPortalAccessor.getCancerStudies(), data.study); // GBM
		cBioPortalAccessor.setCurrentCancerStudy(cancerStudy);

		List<GeneticProfile> geneticProfilesForCurrentStudy =
			cBioPortalAccessor.getGeneticProfilesForCurrentStudy();
		List<GeneticProfile> gp = new ArrayList<GeneticProfile>();
		for (String prof : data.profile)
		{
			if (prof.endsWith("mrna")) continue;
			gp.add(findProfile(geneticProfilesForCurrentStudy, prof));
		}
		cBioPortalAccessor.setCurrentGeneticProfiles(gp);

		List<CaseList> caseLists = cBioPortalAccessor.getCaseListsForCurrentStudy();
		cBioPortalAccessor.setCurrentCaseList(findCaseList(caseLists, data.caseList));
		return cBioPortalAccessor;
	}

	private ExpDataManager getExpDataMan(PortalDataset data, CBioPortalAccessor acc) throws IOException
	{
		for (String prof : data.profile)
		{
			if (prof.endsWith("mrna"))
			{
				return new ExpDataManager(
					findProfile(acc.getGeneticProfilesForCurrentStudy(), prof),
					acc.getCurrentCaseList());
			}
		}
		return null;
	}

	private static CancerStudy findCancerStudy(List<CancerStudy> list, String id)
	{
		for (CancerStudy study : list)
		{
			if (study.getStudyId().equals(id)) return study;
		}
		return null;
	}

	private static CaseList findCaseList(List<CaseList> list, String id)
	{
		for (CaseList cl : list)
		{
			if (cl.getId().equals(id)) return cl;
		}
		return null;
	}

	private static GeneticProfile findProfile(List<GeneticProfile> list, String id)
	{
		for (GeneticProfile profile : list)
		{
			if (profile.getId().equals(id)) return profile;
		}
		return null;
	}

	public CBioPortalAccessor getAccessor()
	{
		return accessor;
	}
}
