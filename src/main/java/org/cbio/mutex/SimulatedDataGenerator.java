package org.cbio.mutex;

import org.cbio.causality.analysis.Graph;
import org.cbio.causality.data.portal.*;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;
import org.cbio.causality.util.CollectionUtil;

import java.io.IOException;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class SimulatedDataGenerator
{
	private PortalDataset sampleData;
	private Graph graph;
	private static final PortalDataset simData = PortalDataset.SIMUL;
	private static final int SAMPLE_SIZE = 100;
	private static final int TOTAL_GROUPS = 10;
	private static final int GROUP_SIZE = 5;
	private static final int NOISE_ALTS_PER_SAMPLE = 300;
	private static final Alteration ALT = Alteration.MUTATION;
	private static final double[] PROBS = new double[]{35, 25, 20, 15, 5};
	private List<List<String>> groups;
	private Map<List<String>, List<String>> pollMap;
	private Map<String, AlterationPack> packMap;

	private void initAlterationData()
	{
		packMap = new HashMap<String, AlterationPack>();
		for (String gene : graph.getSymbols())
		{
			AlterationPack pack = new AlterationPack(gene);
			Change[] changes = new Change[SAMPLE_SIZE];
			for (int i = 0; i < changes.length; i++)
			{
				changes[i] = Change.NO_CHANGE;
			}
			pack.put(ALT, changes);
		}
	}

	private void alterEachGroupForEachSample()
	{
		for (int i = 0; i < SAMPLE_SIZE; i++)
		{
			Set<String> altered = new HashSet<String>();

			for (List<String> group : groups)
			{
				if (CollectionUtil.intersects(group, altered)) continue;

				String gene = selectOne(altered);
				AlterationPack pack = packMap.get(gene);

				assert pack.getChange(ALT, i) == Change.NO_CHANGE;
				pack.get(ALT)[i] = Change.ACTIVATING;
				altered.add(gene);
			}
		}
	}

	private void addNoise()
	{
		for (int i = 0; i < SAMPLE_SIZE; i++)
		{

		}
	}

	private void decideGroups()
	{
		groups = new ArrayList<List<String>>();

		for (String gene : graph.getSymbols())
		{
			Set<String> upstream = graph.getUpstream(gene);
			if (upstream.size() > GROUP_SIZE)
			{
				List<String> group = new ArrayList<String>();
				List<String> list = new ArrayList<String>(upstream);
				Collections.shuffle(list);
				group.addAll(list.subList(0, GROUP_SIZE));
				groups.add(group);
			}

			if (groups.size() == TOTAL_GROUPS) break;
		}

		pollMap = new HashMap<List<String>, List<String>>();

		for (List<String> group : groups)
		{
			List<String> poll = new ArrayList<String>();

			for (int i = 0; i < group.size(); i++)
			{
				String gene = group.get(i);

				for (int j = 0; j < PROBS[i]; j++)
				{
					poll.add(gene);
				}
			}

			pollMap.put(group, poll);
		}
	}

	public CBioPortalAccessor getPortalAccessor()
	{
		try
		{
			CBioPortalAccessor acc = new CBioPortalAccessor();
			acc.setWithoutCheck(new CancerStudy(simData.study, simData.study, simData.study),
				new CaseList(simData.caseList, simData.caseList, prepareCaseArray()),
				getProfiles());

			return acc;
		}
		catch (IOException e)
		{
			e.printStackTrace();
			return null;
		}
	}

	private String[] prepareCaseArray()
	{
		String[] s = new String[SAMPLE_SIZE];
		for (int i = 0; i < s.length; i++)
		{
			s[i] = "S" + i;
		}
		return s;
	}

	private GeneticProfile[] getProfiles()
	{
		GeneticProfile[] gp = new GeneticProfile[simData.profile.length];
		for (int i = 0; i < gp.length; i++)
		{
			gp[i] = new GeneticProfile(simData.profile[i], simData.profile[i], simData.profile[i],
				ProfileType.MUTATION_EXTENDED);
		}
		return gp;
	}

	private String selectOne(Set<String> group)
	{
		Collections.shuffle(pollMap.get(group));
		return pollMap.get(group).get(0);
	}

}
