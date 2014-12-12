package org.cbio.mutex;

import org.cbio.causality.analysis.Graph;
import org.cbio.causality.data.portal.*;
import org.cbio.causality.idmapping.EntrezGene;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;
import org.cbio.causality.util.Histogram;
import org.cbio.causality.util.TermCounter;

import java.io.*;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class Simulation
{
	private static final Alteration ALT = Alteration.MUTATION;

	private List<List<String>> decideGroups(Graph graph, Set<String> genes, int groupSize,
		int totalGroups)
	{
		List<List<String>> groups = new ArrayList<List<String>>();
		Set<String> selected = new HashSet<String>();

		List<String> syms = new ArrayList<String>(graph.getSymbols());
		Collections.shuffle(syms);

		for (String gene : syms)
		{
			Set<String> upstream = new HashSet<String>(graph.getUpstream(gene));
			upstream.retainAll(genes);
			upstream.removeAll(selected);
			if (upstream.size() >= groupSize)
			{
				List<String> group = new ArrayList<String>();
				List<String> list = new ArrayList<String>(upstream);
				Collections.shuffle(list);
				group.addAll(list.subList(0, groupSize));
				selected.addAll(group);
				groups.add(group);
			}

			if (groups.size() == totalGroups) break;
		}

		if (groups.size() < totalGroups)
		{
			System.err.println("Cannot create all groups. Desired: " + totalGroups + ", created: " + groups.size());
			System.exit(0);
		}

		return groups;
	}

	private List<Set<String>> generateSamples(Map<String, GeneAlt> geneMap, List<List<String>> groups)
	{
		for (GeneAlt alt : geneMap.values())
		{
			alt.shuffleSticky();
		}

		Map<GeneAlt, List<List<GeneAlt>>> gene2groups = new HashMap<GeneAlt, List<List<GeneAlt>>>();
		for (List<String> group : groups)
		{
			for (String id : group)
			{
				GeneAlt gene = geneMap.get(id);

				if (!gene2groups.containsKey(gene)) gene2groups.put(gene, new ArrayList<List<GeneAlt>>());

				List<GeneAlt> g = new ArrayList<GeneAlt>();
				for (String id2 : group)
				{
					g.add(geneMap.get(id2));
				}
				gene2groups.get(gene).add(g);
			}
		}

		System.out.println("-- before mutex ----\n");
		printGroups(geneMap, groups);

		Random r = new Random();

		for (int j = 0; j < 2; j++)
		{
			for (GeneAlt gene : gene2groups.keySet())
			{
				for (int i = 0; i < gene.shuf.length; i++)
				{
					if (overlaps(gene, i, gene2groups))
					{
						int newInd;
						do
						{
							newInd = r.nextInt(gene.shuf.length);
						}
						while (gene.shuf[newInd]);

						gene.shuf[newInd] = true;
						gene.shuf[i] = false;
					}
				}
			}
		}

		for (GeneAlt gene : geneMap.values())
		{
			gene.ch = null;
		}

		System.out.println("\n\n-- after mutex ----\n");
		printGroups(geneMap, groups);

		int sampleSize = geneMap.values().iterator().next().shuf.length;
		List<Set<String>> samples = new ArrayList<Set<String>>(sampleSize);

		for (int i = 0; i < sampleSize; i++)
		{
			Set<String> set = new HashSet<String>();
			for (GeneAlt gene : geneMap.values())
			{
				if (gene.shuf[i]) set.add(gene.getId());
			}
			samples.add(set);
		}

		return samples;
	}

	private boolean overlaps(GeneAlt gene, int index, Map<GeneAlt, List<List<GeneAlt>>> gene2groups)
	{
		if (!gene.shuf[index]) return false;

		for (List<GeneAlt> group : gene2groups.get(gene))
		{
			for (GeneAlt g2 : group)
			{
				if (gene == g2) continue;

				if (g2.shuf[index]) return true;
			}
		}
		return false;
	}

	private void printGroups(Map<String, GeneAlt> geneMap, List<List<String>> groups)
	{
		for (List<String> list : groups)
		{
			Group group = new Group();
			for (String id : list)
			{
				group.addGene(geneMap.get(id));
			}
			String s = group.getPrint(null, null, false);
			System.out.println(s.substring(0, s.indexOf("\n")));
		}
	}

	private void printCoverages(Collection<String> genes, List<Set<String>> samples,
		List<List<String>> groups)
	{
		System.out.println("genes.size() = " + genes.size());

		Set<String> groupGenes = new HashSet<String>();
		for (List<String> group : groups)
		{
			groupGenes.addAll(group);
		}

		double range = 0.01;
		Histogram h = new Histogram(range);
		Histogram h_group = new Histogram(range);
		Histogram h_nongr = new Histogram(range);
		h.setBorderAtZero(true);
		h_group.setBorderAtZero(true);
		h_nongr.setBorderAtZero(true);
		for (String gene : genes)
		{
			int c = 0;
			for (Set<String> sample : samples)
			{
				if (sample.contains(gene))
				{
					c++;
				}
			}
			double cov = c / (double) samples.size();
			h.count(cov);
			if (groupGenes.contains(gene)) h_group.count(cov);
			else h_nongr.count(cov);
		}
		System.out.println("Gene coverage:");
		h.print();
		System.out.println("Gene coverage - in and out of groups:");
		h_group.printTogether(h_nongr);

		h = new Histogram(range);
		h.setBorderAtZero(true);
		for (Set<String> sample : samples)
		{
			h.count(sample.size() / (double) genes.size());
		}
		System.out.println("Sample coverage:");
		h.print();
	}

	private void printCoverages(Map<String, AlterationPack> packMap, PortalDataset simData) throws FileNotFoundException
	{
		Collection<String> genes = new HashSet<String>();
		List<Set<String>> samples = new ArrayList<Set<String>>();
		int sampleSize = packMap.values().iterator().next().getSize();
		for (int i = 0; i < sampleSize; i++)
		{
			samples.add(new HashSet<String>());
		}
		List<List<String>> groups = readGroups(simData);

		for (String gene : packMap.keySet())
		{
			genes.add(gene);
			AlterationPack pack = packMap.get(gene);
			for (int i = 0; i < sampleSize; i++)
			{
				if (pack.get(ALT)[i].isAltered()) samples.get(i).add(gene);
			}
		}

		printCoverages(genes, samples, groups);
	}

	private void printGroupsStats(List<List<String>> groups)
	{
		Map<String, Integer> cnt = new HashMap<String, Integer>();
		for (List<String> group : groups)
		{
			for (String gene : group)
			{
				if (cnt.containsKey(gene)) cnt.put(gene, cnt.get(gene) + 1);
				else cnt.put(gene, 1);
			}
		}
		TermCounter tc = new TermCounter();
		for (Integer c : cnt.values())
		{
			tc.addTerm(c + "");
		}
		tc.print();
	}

	private void generate(Graph graph, Map<String, GeneAlt> geneMap, PortalDataset simData) throws IOException
	{
		Set<String> genes = geneMap.keySet();
		List<List<String>> groups = decideGroups(graph, genes, extractGroupSize(simData), extractNumberOfGroups(simData));
		writeGroups(groups, simData);
		printGroupsStats(groups);
		List<Set<String>> samples = generateSamples(geneMap, groups);
		saveDataToCache(transferToAlterationData(samples, genes), simData);
		printCoverages(genes, samples, groups);
	}

	private void generateNoWrite(Graph graph, Map<String, GeneAlt> geneMap, PortalDataset simData)
	{
		Set<String> genes = geneMap.keySet();
		List<List<String>> groups = decideGroups(graph, genes, extractGroupSize(simData), extractNumberOfGroups(simData));
		printGroupsStats(groups);
		List<Set<String>> samples = generateSamples(geneMap, groups);
		printCoverages(genes, samples, groups);
	}

	private static int extractSampleSize(PortalDataset simData)
	{
		int i = simData.study.indexOf("_ss");
		return Integer.parseInt(simData.study.substring(i + 3));
	}

	private int extractNumberOfGroups(PortalDataset simData)
	{
		int i = simData.study.indexOf("_ng");
		return Integer.parseInt(simData.study.substring(i + 3, simData.study.indexOf("_", i + 1)));
	}

	public static int extractGroupSize(PortalDataset simData)
	{
		int i = simData.study.indexOf("_gs");
		return Integer.parseInt(simData.study.substring(i + 3, simData.study.indexOf("_", i + 1)));
	}

	private Map<String, AlterationPack> transferToAlterationData(List<Set<String>> samples, Set<String> genes)
	{
		Map<String, AlterationPack> packMap = new HashMap<String, AlterationPack>();
		for (String gene : genes)
		{
			AlterationPack pack = new AlterationPack(gene);
			Change[] changes = new Change[samples.size()];
			for (int i = 0; i < changes.length; i++)
			{
				changes[i] = samples.get(i).contains(gene) ? Change.ACTIVATING : Change.NO_CHANGE;
			}
			pack.put(ALT, changes);
			packMap.put(gene, pack);
		}

		return packMap;
	}


	private void writeGroups(List<List<String>> groups, PortalDataset simData) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(getGroupFilename(simData)));

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

	public static List<List<String>> readGroups(PortalDataset simData) throws FileNotFoundException
	{
		List<List<String>> groups = new ArrayList<List<String>>();
		Scanner sc = new Scanner(new File(getGroupFilename(simData)));
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

	public static CBioPortalAccessor getPortalAccessor(PortalDataset simData) throws IOException
	{
		CBioPortalAccessor acc = new CBioPortalAccessor();
		acc.setWithoutCheck(new CancerStudy(simData.study, simData.study, simData.study),
			new CaseList(simData.caseList, simData.caseList, prepareCaseArray(extractSampleSize(simData))),
			getProfiles(simData));

		return acc;
	}

	private static GeneticProfile[] getProfiles(PortalDataset simData)
	{
		GeneticProfile[] gp = new GeneticProfile[simData.profile.length];
		for (int i = 0; i < gp.length; i++)
		{
			gp[i] = new GeneticProfile(simData.profile[i], simData.profile[i], simData.profile[i],
				ProfileType.MUTATION_EXTENDED);
		}
		return gp;
	}

	private static String[] prepareCaseArray(int sampleSize)
	{
		String[] s = new String[sampleSize];
		for (int i = 0; i < s.length; i++)
		{
			s[i] = "S" + i;
		}
		return s;
	}

	public static Map<String, AlterationPack> loadSimData(Set<String> syms, PortalDataset simData) throws IOException
	{
		Map<String, AlterationPack> packs = new HashMap<String, AlterationPack>();
		CBioPortalAccessor acc = getPortalAccessor(simData);
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

	private void saveDataToCache(Map<String, AlterationPack> packMap, PortalDataset simData) throws IOException
	{
		String dir = "portal-cache/" + simData.profile[0] + "/" + simData.caseList + "/";
		File f = new File(dir);


		if (f.exists()) System.err.println("Simulated data cache already exists!");
		else f.mkdirs();

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

		CBioPortalAccessor acc = getPortalAccessor(simData);

		for (AlterationPack pack : packMap.values())
		{
			AlterationPack p2 = acc.getAlterations(pack.getId());
			p2.complete(ALT);

			assert pack.countAltered(ALT) == p2.countAltered(ALT);
		}
	}

	private static String getGroupFilename(PortalDataset simData)
	{
		return "groups-" + simData.study + ".txt";
	}

	public static Map<String, String> getLabelMap(PortalDataset simData,
		Set<String> ids) throws FileNotFoundException
	{
		List<List<String>> groups = readGroups(simData);
		Map<String, String> nameMap = new HashMap<String, String>();

		for (int i = 0; i < groups.size(); i++)
		{
			List<String> group = groups.get(i);
			String pre = "T" + (i + 1);
			for (int j = 0; j < group.size(); j++)
			{
				String gene = group.get(j);
				String nm = pre + "." + (j + 1);
				if (nameMap.containsKey(gene)) nameMap.put(gene, nameMap.get(gene) + "-" + nm.substring(1));
				else nameMap.put(gene, nm);
			}
		}

		int i = 0;
		for (String gene : ids)
		{
			if (!nameMap.containsKey(gene))
			{
				nameMap.put(gene, "F" + (++i));
			}
		}

		return nameMap;
	}

	protected Graph convertGraph(Graph g1, PortalDataset dataset) throws FileNotFoundException
	{
		Map<String, String> map = getLabelMap(dataset, g1.getSymbols());

		Graph g2 = new Graph("Simulated graph", "is-upstream-of");

		for (String source : g1.getOneSideSymbols(true))
		{
			for (String target : g1.getDownstream(source))
			{
				assert map.containsKey(source);
				assert map.containsKey(target);

				g2.putRelation(map.get(source), map.get(target), true);
			}
		}
		return g2;
	}

	public static void evaluateSuccess(List<Group> inferred, PortalDataset simData) throws IOException
	{
		System.out.println("inferred groups size = " + inferred.size());
		Set<String> allSim = getGenesInGroups(simData);

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

		CBioPortalAccessor acc = getPortalAccessor(simData);
		Map<Integer, Set<String>> tierMap = new HashMap<Integer, Set<String>>();
		for (String gene : allSim)
		{
			AlterationPack pack = acc.getAlterations(gene);
			int t = pack.getAlteredCount(ALT);
			if (!tierMap.containsKey(t)) tierMap.put(t, new HashSet<String>());
			tierMap.get(t).add(gene);
		}
		int sampleSize = extractSampleSize(simData);
		System.out.println("Tier recovery rates:\nA#\texists\trecovered\tratio");
		for (int i = 0; i < sampleSize; i++)
		{
			if (tierMap.containsKey(i))
			{
				Set<String> found = new HashSet<String>(tierMap.get(i));
				found.retainAll(allInf);
				int exists = tierMap.get(i).size();
				int foundSize = found.size();
				System.out.println(i + "\t" + exists + "\t" + foundSize + "\t" + (foundSize / (double) exists));
			}
		}
	}

	private static Set<String> getGenesInGroups(PortalDataset simData) throws FileNotFoundException
	{
		List<List<String>> groups = readGroups(simData);
		Set<String> trueGenes = new HashSet<String>();
		for (List<String> group : groups)
		{
			trueGenes.addAll(group);
		}
		return trueGenes;
	}

	public static int[] getTrueAndFalseCounts(Set<String> genes, PortalDataset simData) throws FileNotFoundException
	{
		Set<String> trues = getGenesInGroups(simData);
		int tc = 0;
		int fc = 0;
		for (String gene : genes)
		{
			if (trues.contains(gene)) tc++;
			else fc++;
		}
		return new int[]{tc, fc};
	}

	public static Set<String>[] getTrueAndFalses(Set<String> genes, PortalDataset simData) throws FileNotFoundException
	{
		Set<String> trues = getGenesInGroups(simData);
		Set<String> ts = new HashSet<String>();
		Set<String> fs = new HashSet<String>();
		for (String gene : genes)
		{
			if (trues.contains(gene)) ts.add(gene);
			else fs.add(gene);
		}
		return new Set[]{ts, fs};
	}

	public void writeDataToMatrixForRME(PortalDataset simData, Graph graph, String filename) throws IOException
	{
		Map<String, AlterationPack> packs = loadSimData(graph.getSymbols(), simData);
		String[] cases = getPortalAccessor(simData).getCurrentCaseList().getCases();
		Map<String, String> labelMap = getLabelMap(simData, packs.keySet());
		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

		for (String aCase : cases) writer.write("\t" + aCase);

		for (String gene : packs.keySet())
		{
			AlterationPack pack = packs.get(gene);
			writer.write("\n" + labelMap.get(gene));

			for (Change ch : pack.get(ALT))
			{
				writer.write("\t" + (ch.isAltered() ? 1 : 0));
			}
		}

		writer.close();
	}

	public void writeDataToMatrixForMDPFinder(PortalDataset simData, Graph graph, String filename) throws IOException
	{
		Map<String, AlterationPack> packs = loadSimData(graph.getSymbols(), simData);
		String[] cases = getPortalAccessor(simData).getCurrentCaseList().getCases();
		Map<String, String> labelMap = getLabelMap(simData, packs.keySet());
		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

		List<AlterationPack> packList = new ArrayList<AlterationPack>(packs.values());

		for (AlterationPack pack : packList) writer.write("\t" + labelMap.get(pack.getId()));

		int i = 0;
		for (String aCase : cases)
		{
			writer.write("\n" + aCase);

			for (AlterationPack pack : packList)
			{
				writer.write("\t" + (pack.get(ALT)[i].isAltered() ? 1 : 0));
			}
			i++;
		}

		writer.close();
	}

	public void writeDataToMatrixForDendrix(PortalDataset simData, Graph graph, String matrixFile, String geneFile) throws IOException
	{
		Map<String, AlterationPack> packs = loadSimData(graph.getSymbols(), simData);
		String[] cases = getPortalAccessor(simData).getCurrentCaseList().getCases();
		Map<String, String> labelMap = getLabelMap(simData, packs.keySet());
		BufferedWriter writer = new BufferedWriter(new FileWriter(matrixFile));

		List<String> newNames = new ArrayList<String>(labelMap.values());
		Collections.sort(newNames);

		for (int i = 0; i < cases.length; i++)
		{
			if (i > 0) writer.write("\n");
			writer.write(cases[i]);

			for (AlterationPack pack : packs.values())
			{
				if (pack.get(ALT)[i].isAltered()) writer.write("\t" + labelMap.get(pack.getId()));
//				if (pack.get(ALT)[i].isAltered()) writer.write("\t" + pack.getId());
			}
		}
		writer.close();
		writer = new BufferedWriter(new FileWriter(geneFile));
		for (AlterationPack pack : packs.values())
		{
			writer.write(labelMap.get(pack.getId()) + "\n");
//			writer.write(pack.getId() + "\n");
		}
		writer.close();
	}

	public void writeDataToMatrixForMEMo(PortalDataset simData, Graph graph, String mafFile,
		String mutsigFile, String casesFile, String cnaFile) throws IOException
	{
		Map<String, AlterationPack> packs = loadSimData(graph.getSymbols(), simData);
		String[] cases = getPortalAccessor(simData).getCurrentCaseList().getCases();
		Map<String, String> labelMap = getLabelMap(simData, packs.keySet());
		BufferedWriter writer = new BufferedWriter(new FileWriter(mafFile));
		writer.write("Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_position\tEnd_position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tdbSNP_RS\tdbSNP_Val_Status\tTumor_Sample_Barcode\tMatched_Norm_Sample_Barcode\tMatch_Norm_Seq_Allele1\tMatch_Norm_Seq_Allele2\tTumor_Validation_Allele1\tTumor_Validation_Allele2\tMatch_Norm_Validation_Allele1\tMatch_Norm_Validation_Allele2\tVerification_Status\tValidation_Status\tMutation_Status\tSequencing_Phase\tSequence_Source\tValidation_Method\tScore\tBAM_file\tSequencer\tchromosome_name_WU\tstart_WU\tstop_WU\treference_WU\tvariant_WU\ttype_WU\tgene_name_WU\ttranscript_name_WU\ttranscript_species_WU\ttranscript_source_WU\ttranscript_version_WU\tstrand_WU\ttranscript_status_WU\ttrv_type_WU\tc_position_WU\tamino_acid_change_WU\tucsc_cons_WU\tdomain_WU\tall_domains_WU\tdeletion_substructures_WU\tGenome_Change_Broad\tAnnotation_Transcript_Broad\tcDNA_Change_Broad\tCodon_Change_Broad\tProtein_Change_Broad\tFAM:variant\tFAM:GE.rank\tFAM:CNA\tFAM:OV.variant.samples\tFAM:OV.gene.samples\tFAM:mapping.issue\tFAM:FImpact\tFAM:FI.score\tFAM:Func.region\tFAM:bindsite.protein\tFAM:bindsite.DNA/RNA\tFAM:bindsite.sm.mol\tFAM:CancerGenes\tFAM:TS\tFAM:OG\tFAM:COSMIC.mutations\tFAM:COSMIC.cancers\tFAM:Uniprot.regions\tFAM:TS.interacts\tFAM:OG.interacts\tFAM:Pfam.domain\tFAM:link.var\tFAM:link.MSA\tFAM:link.PDB");

		List<String> newNames = new ArrayList<String>(labelMap.values());
		Collections.sort(newNames);

		for (int i = 0; i < cases.length; i++)
		{
			for (AlterationPack pack : packs.values())
			{
				if (pack.get(ALT)[i].isAltered())
				{
					writer.write("\n" + pack.getId() + "\t" + EntrezGene.getID(pack.getId()));
					writer.write("\tgenome.wustl.edu\t36\t9\t23691590\t23691590\t+\tmissense\tSNP\tC\tC\tT\tUnknown\tUnknown\t");
					writer.write(cases[i]);
					writer.write("\tTCGA-02-0083-10A-01W\tC\tC\t\t\t\t\t\tValid\tSomatic\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tmissense\t\tG167D\t\t\t\t\t\t\t\t\t\tG167D\t\t\t\t\t\tH\t3.32\t\t\t\t\t\t\t\t\t\t\t\t\t\thttp://mutationassessor.org?cm=var&var=9,23691590,C,T&fts=all\thttp://mutationassessor.org/?cm=msa&ty=f&p=ELAV2_HUMAN&rb=127&re=195&var=G167D\thttp://mutationassessor.org/pdb.php?prot=ELAV2_HUMAN&from=127&to=195&var=G167D");
				}
			}
		}
		writer.close();
		writer = new BufferedWriter(new FileWriter(mutsigFile));
		writer.write("rank\tgene\tN\tn\tnVal\tnVer\tCpG\tC+G\tA+T\tIndel\tp\tq");
		int k = 0;
		for (AlterationPack pack : packs.values())
		{
			writer.write("\n" + (++k) + "\t" + pack.getId() + "\t145177\t48\t48\t0\t18\t17\t10\t3\t<1E-11\t<1E-8");
		}
		writer.close();
		writer = new BufferedWriter(new FileWriter(casesFile));
		writer.write("#cancer_type_id: simulated\n" +
			"#stable_id: simulated\n" +
			"#case_list_name: Simulated\n" +
			"#case_list_description: Simulated");
		for (String aCase : cases)
		{
			writer.write("\n" + aCase);
		}
		writer.close();
		writer = new BufferedWriter(new FileWriter(cnaFile));
		writer.write("Entrez_Gene_ID");
		for (String aCase : cases) writer.write("\t" + aCase);
		boolean printed = false;
		for (AlterationPack pack : packs.values())
		{
			if (!printed && labelMap.get(pack.getId()).startsWith("F")) {System.out.println(EntrezGene.getID(pack.getId())); printed = true;}
			writer.write("\n" + EntrezGene.getID(pack.getId()));
			for (int i = 0; i < cases.length; i++)
			{
				writer.write("\t0");
			}
		}
		writer.close();
	}

	private void writeSupplementInteractionGraphForMEMo(PortalDataset simData, String filename) throws IOException
	{
		List<List<String>> groups = readGroups(simData);
		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

		for (List<String> group : groups)
		{
			for (int i = 0; i < group.size() - 1; i++)
			{
				for (int j = i+1; j < group.size(); j++)
				{
					writer.write(group.get(i) + "\tinteracts-with\t" + group.get(j) + "\n");
				}
			}
		}

		writer.close();
	}

	private void printMEMoGroups(String filename, PortalDataset simData, Graph graph) throws IOException
	{
		Map<String, AlterationPack> packs = loadSimData(graph.getSymbols(), simData);
		Map<String, String> labelMap = getLabelMap(simData, packs.keySet());

		Scanner sc = new Scanner(new File(filename));
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();

			if (line.startsWith("Mo")) System.out.println(line);
			else
			{
				String[] tok = line.split("\t");

				String gline = "";
				for (String s : tok[1].split(","))
				{
					s = s.trim();
					if (s.isEmpty()) continue;
					s = s.split(" ")[0];
					gline += labelMap.get(s) + " ";
				}

				tok[1] = gline;
				for (String s : tok)
				{
					System.out.print(s + "\t");
				}
				System.out.println();
			}

		}
	}

	public static void main(String[] args) throws IOException, ClassNotFoundException
	{
		Simulation sim = new Simulation();
		PortalDataset data = PortalDatasetEnum.SIMUL2.data;

//		String dir = "/home/ozgun/Documents/mutex-comparison/MEMo/cancer_data/simulation2-suppl/";
//		sim.printMEMoGroups(dir + "MemoReport.txt", data, Main.getGraph());
//		sim.writeDataToMatrixForRME(data, Main.getGraph(), "/home/ozgun/Documents/mutex-comparison/rmeMod/simdata3/SimData-v3.txt");
//		sim.writeDataToMatrixForMDPFinder(data, Main.getGraph(), "/home/ozgun/Documents/mutex-comparison/MDPfinder/data/simulated3/SimData-v3.txt");
//		sim.writeDataToMatrixForDendrix(data, Main.getGraph(), "/home/ozgun/Documents/mutex-comparison/Dendrix/temp/SimData-v3.txt", "/home/ozgun/Documents/mutex-comparison/Dendrix/temp/analyzed_genes-v3.txt");
//		String dir = "/home/ozgun/Documents/mutex-comparison/MEMo/cancer_data/simulation3/";
//		sim.writeDataToMatrixForMEMo(data, Main.getGraph(), dir + "data_mutations_MAF.txt", dir + "sig_genes_phase_1_2.txt", dir + "cases_all_three.txt", dir + "data_CNA_RAE.txt");
//		sim.writeSupplementInteractionGraphForMEMo(data, "/home/ozgun/Documents/mutex-comparison/MEMo/data/simnetwork.sif");

//		MutexGreedySearcher searcher = MutexGreedySearcher.deserialize("cache-brca_tcga_pub");
//		sim.generate(Main.getGraph(), searcher.getGenes(), data);
//		sim.generateNoWrite(Main.getGraph(), searcher.getGenes(), data);

//		Map<String, AlterationPack> packMap = loadSimData(Main.getGraph().getSymbols(), data);
//		sim.printCoverages(packMap, data);

		Graph graph = sim.convertGraph(new Network(), data);
		graph.write("data/simulation2/network.txt");
	}
}