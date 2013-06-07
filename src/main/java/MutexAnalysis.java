import org.cbio.causality.analysis.SIFLinker;
import org.cbio.causality.analysis.Traverse;
import org.cbio.causality.data.portal.CBioPortalAccessor;
import org.cbio.causality.data.portal.CancerStudy;
import org.cbio.causality.data.portal.CaseList;
import org.cbio.causality.data.portal.GeneticProfile;
import org.cbio.causality.idmapping.HGNC;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;
import org.cbio.causality.util.Overlap;
import org.cbio.causality.util.Progress;
import org.cbio.causality.util.Summary;
import org.biopax.paxtools.model.level3.ProteinReference;
import org.biopax.paxtools.model.level3.Xref;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

/**
 * @author Ozgun Babur
 * @deprecated
 */
public class MutexAnalysis
{
	static final DecimalFormat fmt = new DecimalFormat("0.0000");
	static final double THR = 0.05;
	static final int MAX_GROUP_SIZE = 4;

	public static void main(String[] args) throws IOException, CloneNotSupportedException
	{
		MutexAnalysis ma = new MutexAnalysis();
		ma.run();
	}

	public void run() throws IOException, CloneNotSupportedException
	{
		// Load network

		SIFLinker linker = new SIFLinker();
		linker.load(getClass().getResourceAsStream("network.txt"),
			"STATE_CHANGE", "TRANSCRIPTION", "DEGRADATION");

		Traverse trav = linker.traverse;

		// Load alteration data

		Map<String, AlterationPack> map = readAlterations(trav.getSymbols());
		for (String s : new HashSet<String>(map.keySet()))
		{
			AlterationPack pack = map.get(s);

			// Filter out genes with les than 3% alteration
			if (pack.getAlteredRatio(Alteration.GENOMIC) < 0.03) map.remove(s);
		}

//		if (true)
//		{
//			AltBundle bun = new AltBundle(map.get("RB1"), Alteration.INHIBITING, 1);
//			bun.add(map.get("CDK4"), Alteration.ACTIVATING);
//			bun.add(map.get("CDKN2A"), Alteration.INHIBITING);
//
//			bun.findBestUseMapConf(0.05);
//
//			System.out.println(bun.getPrint());
//			return;
//		}

		System.out.println("altered genes = " + map.size());

		List<AltBundle> modules = new ArrayList<AltBundle>();
		Progress p = new Progress(map.size());

		for (String seedName : map.keySet())
		{
			p.tick();

			AlterationPack seed = map.get(seedName);
			if (!seed.isAltered(Alteration.GENOMIC)) continue;

			for (int i = MAX_GROUP_SIZE - 1; i > 0; i--)
			{
				List<AltBundle> groups = getMutexGroups(seed, map, trav, THR, i);

				for (AltBundle group : groups) if (!isSubset(group, modules)) modules.add(group);
			}
		}
		System.out.println();


		Collections.sort(modules, new Comparator<AltBundle>()
		{
			@Override
			public int compare(AltBundle b1, AltBundle b2)
			{
				return new Double(b2.calcCoverage()).compareTo(b1.calcCoverage());
			}
		});

		for (AltBundle bun : modules)
		{
			Set<String> genes = new HashSet<String>(bun.getAllGenes());
			List<String> rels = linker.linkProgressive(genes, genes, 0);
			if (rels.isEmpty()) continue;

			bun.sortToMostAltered();
			System.out.println(bun);
			System.out.println(bun.getPrint());

			for (String rel : rels) System.out.println(rel);
			System.out.println();
		}

//		System.out.println(getGraph(modules, linker));
	}

	private Map<String, AlterationPack> readAlterations(Set<String> syms) throws IOException
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

	private String getSymbol(ProteinReference pr)
	{
		for (Xref xref : pr.getXref())
		{
			if (xref.getDb().equals("HGNC"))
			{
				String id = xref.getId().substring(xref.getId().indexOf(":") + 1);
				return HGNC.getSymbol(id);
			}
		}
		return null;
	}

	public void printAlterationCountsPerSample(Map<String, AlterationPack> map)
	{
		AlterationPack aPack = map.values().iterator().next();
		int size = aPack.getSize();
		Set<Alteration> types = aPack.getAlterationTypes();
		
		int[][] cnt = new int[types.size()][size];

		int i = 0;
		for (Alteration type : types)
		{
			for (int j = 0; j < cnt[i].length; j++)
			{
				for (String sm : map.keySet())
				{
					AlterationPack pack = map.get(sm);
					if (pack.getChange(type, j).isAltered()) cnt[i][j]++;
				}
			}

			i++;
		}

		System.out.print("Sample");
		for (Alteration type : types)
		{
			System.out.print("\t" + type);
		}
		for (int j = 0; j < cnt[0].length; j++)
		{
			System.out.print("\n" + j);
			for (i = 0; i < cnt.length; i++)
			{
				System.out.print("\t" + cnt[i][j]);
			}
			System.out.println();
		}
	}

	private boolean isSubset(AltBundle b, Collection<AltBundle> col)
	{
		for (AltBundle b2 : col)
		{
			if (isSubset(b, b2)) return true;
		}
		return false;
	}
	
	private boolean isSubset(AltBundle b1, AltBundle b2)
	{
		return b2.alts.containsAll(b1.alts);
	}
	
	class AltBundle implements Comparable, Cloneable
	{
		String id;
		List<AlterationPack> alts;
		double pval;
		Map<AlterationPack, Alteration> useMap;
		int candidateSize;

		AltBundle(AlterationPack seed, Alteration addKey, int candidateSize)
		{
			this.candidateSize = candidateSize;
			useMap = new HashMap<AlterationPack, Alteration>();
			this.alts = new ArrayList<AlterationPack>();
			add(seed, addKey);
		}

		public void add(AlterationPack pack, Alteration addKey)
		{
			alts.add(pack);
			useMap.put(pack, addKey);
		}

		private void updateID()
		{
			id = "";
			for (int i = 0; i < alts.size(); i++)
			{
				id += " " + alts.get(i).getId();
			}
			id = id.trim();
		}

		private double calcWorstPVal()
		{
			return getWorstPval(calcPVals());
		}
		
		public double calcAveragePValScore()
		{
			return Summary.geometricMean(calcPVals());
		}
		
		private double[] calcPVals()
		{
			double[] pval = new double[alts.size()];

			if (alts.size() == 2)
			{
				double p = Overlap.calcAlterationMutexPval(
					alts.get(0).get(useMap.get(alts.get(0))),
					alts.get(1).get(useMap.get(alts.get(1))));

				pval[0] = p;
				pval[1] = p;
			}
			else
			{
				int x = 0;

				for (int i = 0; i < alts.size(); i++)
				{
					Change[] others = new Change[alts.get(0).getSize()];

					for (int k = 0; k < others.length; k++)
					{
						others[k] = Change.NO_CHANGE;

						for (int j = 0; j < alts.size(); j++)
						{
							if (j == i) continue;

							if (alts.get(j).get(useMap.get(alts.get(j)))[k].isAltered())
							{
								others[k] = Change.ACTIVATING;
								break;
							}
						}
					}

					pval[x++] = Overlap.calcAlterationMutexPval(
						alts.get(i).get(useMap.get(alts.get(i))), others);
				}

				assert x == pval.length;
			}
			adjustToMultipleHypothesisTesting(pval);
			return pval;
		}

		private void adjustToMultipleHypothesisTesting(double[] pval)
		{
			// First pval belongs to the seed, don't touch it
			for (int i = 1; i < pval.length; i++)
			{
				pval[i] = 1 - pow(1 - pval[i], candidateSize);
			}
		}

		private double pow(double v, int exp)
		{
			double e = 1;
			for (int i = 0; i < exp; i++)
			{
				e *= v;
			}
			return e;
		}
		
		private int pow(int v, int exp)
		{
			int e = 1;
			for (int i = 0; i < exp; i++)
			{
				e *= v;
			}
			return e;
		}
		
		public double getWorstPval(double[] pval)
		{
			if (pval.length == 1) return pval[0];

			double w = 0;

			for (double pv : pval)
			{
				if (w < pv) w = pv;
			}
			return w;
		}

		public double calcCoverage()
		{
			int total = alts.get(0).get(useMap.get(alts.get(0))).length;
			int altered = 0;
			
			for (int i = 0; i < total; i++)
			{
				for (AlterationPack alt : alts)
				{
					if (alt.get(useMap.get(alt))[i].isAltered())
					{
						altered++;
						break;
					}
				}
			}
			return altered / (double) total;
		}
		
		@Override
		public int compareTo(Object o)
		{
			if (o instanceof AltBundle)
			{
				AltBundle bundle = (AltBundle) o;
				if (alts.size() == bundle.alts.size())
					return new Double(Math.abs(pval)).compareTo(Math.abs(bundle.pval));
				return new Integer(bundle.alts.size()).compareTo(alts.size());
			}
			else return 0;
		}
		
		public double absPVal()
		{
			return Math.abs(pval);
		}

		public boolean findBestUseMapConf(double thr)
		{
			double bestAvg = 1;
			int bestIndex = -1;
			
			int k = pow(2, alts.size());
			for (int i = 0; i < k; i++)
			{
				for (int j = 0; j < alts.size(); j++)
				{
					if ((i / pow(2, j)) % 2 == 0) useMap.put(alts.get(j), Alteration.ACTIVATING);
					else useMap.put(alts.get(j), Alteration.INHIBITING);
				}

				double[] pv = calcPVals();
				
				if (getWorstPval(pv) < thr)
				{
					double avg = Summary.geometricMean(pv);
					if (avg < bestAvg) bestIndex = i;
				}
			}

			if (bestIndex > -1)
			{
				for (int j = 0; j < alts.size(); j++)
				{
					if ((bestIndex / pow(2, j)) % 2 == 0) useMap.put(alts.get(j), Alteration.ACTIVATING);
					else useMap.put(alts.get(j), Alteration.INHIBITING);
				}
				return true;
			}
			else return false;
		}
		
		@Override
		public String toString()
		{
			return "Genes: " + id + "\tAvg p-val: " + fmt.format(calcAveragePValScore()) + "\tseed: " + alts.get(0).getId() +
				"\t coverage: " + fmt.format(calcCoverage());
		}

		public String getPrint()
		{
			return getPrint(getPrintOrdering(alts));
		}
		
		public AlterationPack getSeed()
		{
			return alts.get(0);
		}
		
		public AlterationPack getAltPack(String id)
		{
			for (AlterationPack alt : alts)
			{
				if (alt.getId().equals(id)) return alt;
			}
			return null;
		}
		
		public String getPrint(List<Integer> order)
		{
			double[] p = calcPVals();
			int i = 0;
			StringBuilder s = new StringBuilder();
			for (AlterationPack alt : alts)
			{
				if (s.length() > 0) s.append("\n");
				s.append(alt.getPrint(useMap.get(alt), order)).
					append((alt.getId().length() < 4) ? "  \t" : "\t").
					append((useMap.get(alt) == Alteration.ACTIVATING) ? "+" : "-").append("\tp-val: ").
					append(fmt.format(p[alts.indexOf(alt)]));

				i++;
			}
			return s.toString();
		}
		
		protected void sortToMostAltered()
		{
			Collections.sort(alts, new Comparator<AlterationPack>()
			{
				@Override
				public int compare(AlterationPack alt1, AlterationPack alt2)
				{
					return new Integer(alt2.countAltered(useMap.get(alt2))).compareTo(alt1.countAltered(useMap.get(alt1)));
				}
			});

			updateID();
		}
		
		private List<Integer> getPrintOrdering(List<AlterationPack> alts)
		{
			List<Integer> order = new ArrayList<Integer>();

			for (AlterationPack alt : alts)
			{
				Change[] ch = alt.get(useMap.get(alt));

				for (int i = 0; i < ch.length; i++)
				{
					if (ch[i].isAltered() && !order.contains(i)) order.add(i);
				}
			}
			return order;
		}

		public List<String> getAllGenes()
		{
			List<String> list = new ArrayList<String>(alts.size());
			int i = 0;
			for (AlterationPack alt : alts)
			{
				list.add(alt.getId());
				i++;
			}
			return list;
		}
		
		public String getGeneNamesInString()
		{
			String s = "";
			for (String g : getAllGenes())
			{
				s += " " + g;
			}
			return s;
		}
		
		@Override
		public int hashCode()
		{
			int h = 0;
			for (AlterationPack alt : alts)
			{
				h += alt.hashCode();
			}
			return h;
		}

		@Override
		public boolean equals(Object o)
		{
			if (o instanceof AltBundle)
			{
				AltBundle bundle = (AltBundle) o;
				
				if (alts.size() == bundle.alts.size())
				{
					if (alts.containsAll(bundle.alts)) return true;
				}
			}
			return false;
		}

		@Override
		protected Object clone() throws CloneNotSupportedException
		{
			AltBundle b = (AltBundle) super.clone();
			b.alts = new ArrayList<AlterationPack>(alts);
			b.useMap = new HashMap<AlterationPack, Alteration>(useMap);
			return b;
		}
		
		public boolean contains(AlterationPack pack)
		{
			return this.alts.contains(pack);
		}
	}

	private List<AltBundle> getMutexGroups(AlterationPack seed, Map<String, AlterationPack> altMap,
		Traverse trav, double thr, int size)
	{
		List<AltBundle> groups = new ArrayList<AltBundle>();

		Set<String> neigh = trav.goBFS(Collections.singleton(seed.getId()), 
			Collections.EMPTY_SET, false);

		neigh.retainAll(altMap.keySet());

		List<String> nList = new ArrayList<String>(neigh);

		for (String n : neigh)
		{
			if (!altMap.get(n).isAltered(Alteration.GENOMIC)) nList.remove(n);
		}

		if (nList.size() < size) return Collections.emptyList();


		int[] i = new int[size];
		init(i);
		
		do
		{
			increment(i, nList.size());

			AltBundle bun = new AltBundle(seed, Alteration.GENOMIC, nList.size());

			for (int index : i)
			{
				bun.add(altMap.get(nList.get(index)), Alteration.GENOMIC);
			}

			if (bun.findBestUseMapConf(thr)) groups.add(bun);

		} while (!iterationComplete(i, nList.size()));

		return groups;
	}

	private void init(int[] i)
	{
		for (int j = 0; j < i.length; j++)
		{
			i[j] = j;
		}
		i[i.length-1]--;
	}
	
	private boolean iterationComplete(int[] i, int size)
	{
		for (int j = 0; j < i.length; j++)
		{
			if (i[j] < size - i.length + j) return false;
		}
		return true;
	}

	private void increment(int[] i, int size)
	{
		boolean reachedEnd = false;

		for (int j = i.length-1; j >= 0; j--)
		{
			if (i[j] < size - i.length + j)
			{
				i[j]++;

				if (reachedEnd)
				{
					for (int k = j+1; k < i.length; k++)
					{
						i[k] = i[j] + (k-j);
					}
				}

				break;
			}
			else
			{
				reachedEnd = true;
			}
		}
	}

	private String getGraph(List<AltBundle> bundles, SIFLinker linker)
	{
		String s = "type:graph\ttext:mutex";
		Set<String> nodes = new HashSet<String>();
		Set<String> edges = new HashSet<String>();

		for (AltBundle bundle : bundles)
		{
			if (bundle.alts.size() > 2)
			{
				s += "\ntype:compound\tid:" + bundle.getGeneNamesInString() + "\tmembers";

				for (AlterationPack pack : bundle.alts)
				{
					s += ":" + pack.getId();
				}

				s += "\ttext:" + "cov: " + fmt.format(bundle.calcCoverage());
				s += "\ttextcolor:0 0 0\tbgcolor:255 255 255";
			}
			
			for (AlterationPack pack : bundle.alts)
			{
				if (!nodes.contains(pack.getId()))
				{
					s += "\ntype:node\tid:" + pack.getId() + "\ttext:" + pack.getId();
					Alteration alt = bundle.useMap.get(pack);
					double rat = pack.getAlteredRatio(alt);
					int v = 255 - (int) (rat * 255);
					
					int inh = pack.getAlteredCount(Alteration.INHIBITING);
					int all = pack.getAlteredCount(Alteration.GENOMIC);
					
					boolean activating = alt == Alteration.ACTIVATING || (all - inh) * 10 < all;

					String color = activating ? 
						"255 " + v + " " + v : v + " " + v + " 255";
					s += "\tbgcolor:" + color + "\ttooltip:" + fmt.format(rat);
					nodes.add(pack.getId());
				}
			}

			Set<String> genes = new HashSet<String>(bundle.getAllGenes());
			List<String> relations = linker.linkProgressive(genes, genes, 0);
			for (String rel : relations)
			{
				if (!edges.contains(rel))
				{
					String[] tok = rel.split("\t");
					AlterationPack src = bundle.getAltPack(tok[0]);
					AlterationPack trg = bundle.getAltPack(tok[2]);
					
					s+= "\ntype:edge\tid:" + rel.replaceAll("\t", " ") + "\tsource:" + src.getId() +
						"\ttarget:" + trg.getId() + "\tarrow:Target";
					
					double pv = Overlap.calcAlterationMutexPval(
						src.get(bundle.useMap.get(src)), trg.get(bundle.useMap.get(trg)));
					
					int v = (int) Math.max(0, 255 - (-Math.log(pv) * 55.4));
					String color = v + " " + v + " " + v;
					s += "\tlinecolor:" + color;
					
					edges.add(rel);
				}				
			}
		}
		
		// Put the co-occurrence relations
		
		Set<AlterationPack> set = new HashSet<AlterationPack>();
		Map<AlterationPack, Alteration> useMap = new HashMap<AlterationPack, Alteration>();
		Collections.reverse(bundles);
		for (AltBundle bundle : bundles)
		{
			set.addAll(bundle.alts);
			useMap.putAll(bundle.useMap);
		}
		for (AlterationPack pack1 : set)
		{
			for (AlterationPack pack2 : set)
			{
				if (pack1.getId().compareTo(pack2.getId()) > 0)
				{
					double pv = Overlap.calcAlterationCoocPval(
						pack1.get(useMap.get(pack1)), pack2.get(useMap.get(pack2)));

					if (pv < 0.05)
					{
						s += "\ntype:edge\tid:" + pack1.getId() + " co-occur " + pack2.getId();
						s += "\tsource:" + pack1.getId() + "\ttarget:" + pack2.getId();

						int v = (int) Math.max(0, 255 - (-Math.log(pv) * 10));
						String color = v + " " + v + " " + v;
						s += "\tlinecolor:" + color + "\tstyle:Dashed";

						s += "\tarrow:None";
					}
				}
			}
		}
		Collections.reverse(bundles);
		return s;
	}

	static class Dataset
	{
		public Dataset(String filename, int study, int caseList, int[] profile)
		{
			this.filename = filename;
			this.study = study;
			this.caseList = caseList;
			this.profile = profile;
		}

		String filename;
		int study;
		int caseList;
		int[] profile;
	}
	
	public static final Dataset glioblastoma = new Dataset(
		"Glioblastoma.txt", 7, 0, new int[]{1, 3});
	public static final Dataset ovarian = new Dataset(
		"Ovarian.txt", 16, 0, new int[]{0, 12});
	public static final Dataset breast = new Dataset(
		"Breast.txt", 3, 1, new int[]{0, 8});
	public static final Dataset colon = new Dataset(
		"Colon.txt", 5, 0, new int[]{0, 10});

	static final Dataset data = glioblastoma;
}

