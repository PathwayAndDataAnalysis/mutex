package org.cbio.mutex;

import org.cbio.causality.analysis.Graph;
import org.cbio.causality.data.portal.*;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;
import org.cbio.causality.util.CollectionUtil;

import java.io.*;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class SimulatedDataGenerator
{
	private Graph graph;
	private static final PortalDataset simData = PortalDataset.SIMUL;
	private static final int SAMPLE_SIZE = 100;
	private static final int TOTAL_GROUPS = 20;
	private static final int GROUP_SIZE = 5;
	private static final Alteration ALT = Alteration.MUTATION;
	private static final double[] PROBS = new double[]{50, 20, 15, 10, 5};
	private List<List<String>> groups;
	private Map<List<String>, List<String>> pollMap;
	private Map<String, AlterationPack> packMap;
	private List<String> genes;
	private static final String GROUP_FILE = "sim-groups.txt";

	public SimulatedDataGenerator(Graph graph)
	{
		this.graph = graph;
	}

	public void generate() throws IOException
	{
		decideGroups();
		writeGroups();
		initAlterationData();
		alterEachGroupForEachSample();
		addNoise();
		for (AlterationPack pack : packMap.values())
		{
			pack.complete(Alteration.GENOMIC);
		}
		saveDataToCache();
		for (List<String> group : groups)
		{
			for (String gene : group)
			{
				packMap.remove(gene);
			}
		}
		Main.printDatasetCharacteristics(packMap);
	}

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
			packMap.put(gene, pack);
		}
		genes = new ArrayList<String>(packMap.keySet());
		Collections.sort(genes);
	}

	private void alterEachGroupForEachSample()
	{
		for (int i = 0; i < SAMPLE_SIZE; i++)
		{
			Set<String> altered = new HashSet<String>();

			for (List<String> group : groups)
			{
				if (CollectionUtil.intersects(group, altered)) continue;

				String gene = selectOne(group);
				AlterationPack pack = packMap.get(gene);

				assert pack.getChange(ALT, i) == Change.NO_CHANGE;
				pack.get(ALT)[i] = Change.ACTIVATING;
				altered.add(gene);
			}
		}
	}

	private void addNoise()
	{
		List<List<String>> groups = getAlteredGenes();
		for (int i = 0; i < SAMPLE_SIZE; i++)
		{
			for (String gene : groups.get(i))
			{
				packMap.get(gene).get(ALT)[i] = Change.ACTIVATING;
			}
		}
	}

	private List<List<String>> getAlteredGenes()
	{
		Random r = new Random();
		List<List<String>> groups = new ArrayList<List<String>>();
		List<String> mutated = new LinkedList<String>();

		double[] geneCov = new double[genes.size()];
		for (int i = 0; i < geneCov.length; i++)
		{
			geneCov[i] = getCoverageOfAGene(r);
			for (int j = 0; j < Math.round(geneCov[i] * SAMPLE_SIZE); j++)
			{
				mutated.add(genes.get(i));
			}
		}

		Collections.shuffle(mutated);

		double[] sampleCov = new double[SAMPLE_SIZE];
		for (int i = 0; i <SAMPLE_SIZE; i++)
		{
			sampleCov[i] = getCoverageOfASample(r);
		}

		for (int i = 0; i <SAMPLE_SIZE; i++)
		{
			double total = 0;
			for (int j = i; j < SAMPLE_SIZE; j++)
			{
				total += sampleCov[j];
			}

			int size = (int) Math.round((sampleCov[i] / total) * mutated.size());
			List<String> group = new LinkedList<String>();

			for (int j = 0; j < size; j++)
			{
				String select = null;
				for (int k = 0; k < mutated.size(); k++)
				{
					if (!group.contains(mutated.get(0)))
					{
						select = mutated.remove(0);
						break;
					}
					else
					{
						mutated.add(mutated.remove(0));
					}
				}

				if (select != null)
				{
					group.add(select);
				}
				else
				{
					break;
				}
			}

			groups.add(group);
		}

		return groups;
	}

	private double getCoverageOfAGene(Random r)
	{
		return Math.pow(2, (r.nextGaussian() * 1.6635857655985917) - 7.019409115929619);
	}

	private double getCoverageOfASample(Random r)
	{
		return Math.pow(2, (r.nextGaussian() * 1.7877393339701853) - 6.89664352676453);
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

			assert group.size() == PROBS.length;

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

	private static GeneticProfile[] getProfiles()
	{
		GeneticProfile[] gp = new GeneticProfile[simData.profile.length];
		for (int i = 0; i < gp.length; i++)
		{
			gp[i] = new GeneticProfile(simData.profile[i], simData.profile[i], simData.profile[i],
				ProfileType.MUTATION_EXTENDED);
		}
		return gp;
	}

	private String selectOne(List<String> group)
	{
		Collections.shuffle(pollMap.get(group));
		return pollMap.get(group).get(0);
	}

	private void saveDataToCache() throws IOException
	{
		String dir = "portal-cache/" + simData.profile[0] + "/" + simData.caseList + "/";
		File f = new File(dir);

		if (!f.exists()) f.mkdirs();

		for (AlterationPack pack : packMap.values())
		{
			BufferedWriter writer = new BufferedWriter(new FileWriter(dir + pack.getId()));

			Change[] ch = pack.get(ALT);

			for (int i = 0; i < ch.length; i++)
			{
				if (i > 0) writer.write("\t");
				writer.write(ch[i].isAltered() ? "A1B" : "NaN");
			}

			writer.close();
		}

		CBioPortalAccessor acc = getPortalAccessor();

		for (AlterationPack pack : packMap.values())
		{
			AlterationPack p2 = acc.getAlterations(pack.getId());
			p2.complete(ALT);

			assert pack.countAltered(ALT) == p2.countAltered(ALT);
		}
	}

	public static CBioPortalAccessor getPortalAccessor()
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

	private static String[] prepareCaseArray()
	{
		String[] s = new String[SAMPLE_SIZE];
		for (int i = 0; i < s.length; i++)
		{
			s[i] = "S" + i;
		}
		return s;
	}

	public static List<List<String>> readGroups() throws FileNotFoundException
	{
		List<List<String>> groups = new ArrayList<List<String>>();
		Scanner sc = new Scanner(new File(GROUP_FILE));
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			if (!line.isEmpty())
			{
				groups.add(new ArrayList<String>(Arrays.asList(line.split("\t"))));
			}
		}
		return groups;
	}

	public static Set<String> getMutsig() throws FileNotFoundException
	{
		Set<String> set = new HashSet<String>();
		for (List<String> group : readGroups())
		{
			set.add(group.get(0));
		}
		return set;
	}

	private void writeGroups() throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(GROUP_FILE));

		for (List<String> group : groups)
		{
			for (int i = 0; i < group.size(); i++)
			{
				if (i > 0) writer.write("\t");
				writer.write(group.get(i));
			}
			writer.write("\n");
		}

		writer.close();
	}

	public static Map<String, AlterationPack> loadSimData(Set<String> syms)
	{
		Map<String, AlterationPack> packs = new HashMap<String, AlterationPack>();
		CBioPortalAccessor acc = getPortalAccessor();
		for (String sym : syms)
		{
			AlterationPack pack = acc.getAlterations(sym);
			if (pack != null)
			{
				pack.complete(Alteration.GENOMIC);
				packs.put(sym, pack);
			}
		}
		return packs;
	}

	public static Map<String, String> anonymizeGroupMembers() throws FileNotFoundException
	{
		List<List<String>> groups = readGroups();
		Map<String, String> nameMap = new HashMap<String, String>();
		int A = (int) 'A';
		for (int i = 0; i < groups.size(); i++)
		{
			List<String> group = groups.get(i);
			String pre = (char)(A + i) + "";
			for (int j = 0; j < group.size(); j++)
			{
				String gene = group.get(j);
				if (nameMap.containsKey(gene)) nameMap.put(gene, nameMap.get(gene) + "-" + pre + j);
				else nameMap.put(gene, pre + j);
			}
		}

		for (List<String> group : groups)
		{
			for (String gene : group)
			{
				System.out.print(nameMap.get(gene) + "\t");
			}
			System.out.println();
		}
		System.out.println("nameMap = " + nameMap.size());
		return nameMap;
	}

	public static void evaluateSuccess(List<Group> inferred) throws FileNotFoundException
	{
		List<List<String>> groups = readGroups();

		Set<String> allSim = new HashSet<String>();
		for (List<String> group : groups)
		{
			allSim.addAll(group);
		}

		Set<String> allInf = new HashSet<String>();
		for (Group group : inferred)
		{
			allInf.addAll(group.getGeneNames());
		}

		int allInfSize = allInf.size();
		System.out.println("allInfSize = " + allInfSize);
		Set<String> temp = new HashSet<String>(allInf);
		temp.removeAll(allSim);
		int falseSize = temp.size();
		System.out.println("falseSize = " + falseSize);

		int trueSize = allInfSize - falseSize;
		System.out.println("trueSize = " + trueSize);

		double recovery = trueSize / (double) allSim.size();
		System.out.println("recovery = " + recovery);

		double fdr = falseSize / (double) allInfSize;
		System.out.println("fdr = " + fdr);

		List<List<String>> tiers = new ArrayList<List<String>>();
		for (int i = 0; i < groups.get(0).size(); i++)
		{
			tiers.add(new ArrayList<String>());
		}
		for (List<String> group : groups)
		{
			for (int i = 0; i < group.size(); i++)
			{
				tiers.get(i).add(group.get(i));
			}
		}
		for (int i = 0; i < tiers.size(); i++)
		{
			System.out.println("Tier " + i + ": " +
				CollectionUtil.countOverlap(tiers.get(i), allInf) + "/" + tiers.get(i).size());
		}
	}

	public static void main(String[] args) throws FileNotFoundException
	{
		anonymizeGroupMembers();
	}
}
