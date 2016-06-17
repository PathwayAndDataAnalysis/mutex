package org.cbio.mutex;

import org.panda.resource.tcga.CNAReader;
import org.panda.resource.tcga.ExpressionReader;
import org.panda.resource.tcga.MutSigReader;
import org.panda.resource.tcga.MutationReader;
import org.panda.utility.ArrayUtil;
import org.panda.utility.statistics.Summary;
import org.panda.utility.statistics.TTest;

import java.io.*;
import java.util.*;

/**
 * This class is used to prepare the alteration data matrix from muation (maf), CNA (gistic), and
 * expression (RNAseq) data files.
 *
 * @author Ozgun Babur
 */
public class MatrixFromTCGAFiles
{
	public static void prepareMulti(String dir, String outDir) throws IOException
	{
		for (File f : new File(dir).listFiles())
		{
			if (f.isDirectory())
			{
				prepare(f.getPath(), outDir + "/" + f.getName());
			}
		}
	}

	public static void prepare(String dir, String outDir) throws IOException
	{
		prepare(dir + "/mutation.maf", dir + "/copynumber.txt", dir + "/expression.txt",
			dir, dir + "/scores-amplified.txt", dir + "/scores-deleted.txt", outDir);
	}

	public static void prepare(String mutFile, String CNAFile, String expFile,
		String mutScoreFile, String ampScoreFile, String delScoreFile, String outDir) throws IOException
	{
		MutationReader mutR = new MutationReader(mutFile);
		CNAReader cnaR = new CNAReader(CNAFile);
		ExpressionReader expR = new ExpressionReader(expFile);

		Set<String> sampleSet = mutR.getSamples();
		sampleSet.retainAll(cnaR.getSamples());
		sampleSet.retainAll(expR.getSamples());

		String[] samples = sampleSet.toArray(new String[sampleSet.size()]);
		Arrays.sort(samples);

		System.out.println("samples.size() = " + samples.length);

		Set<String> genes = new HashSet<>(mutR.getGenes());
		Set<String> cnaGenes = new HashSet<>(cnaR.getGenes());
		cnaGenes.retainAll(expR.getGenes());
		genes.addAll(cnaGenes);

		System.out.println("genes = " + genes.size());

		Map<String, Double> mutsigScores = MutSigReader.readPValues(mutScoreFile);
		Map<String, Double> ampScores = readCNAScores(ampScoreFile);
		Map<String, Double> delScores = readCNAScores(delScoreFile);
		List<Gene> rankedGenes = getScoredGenes(genes, samples, cnaR, expR, mutsigScores, ampScores, delScores);

		if (!new File(outDir + "/whole").exists()) new File(outDir + "/whole").mkdirs();
		BufferedWriter wM = new BufferedWriter(new FileWriter(outDir + "/whole/DataMatrix.txt"));
		BufferedWriter wR = new BufferedWriter(new FileWriter(outDir + "/RankedGenes.txt"));

		for (String sample : samples)
		{
			wM.write("\t" + sample);
		}

		wR.write(Gene.getHeader());

		for (Gene gene : rankedGenes)
		{
			String name = gene.name;
			boolean[] mut = mutR.getGeneAlterationArray(name, samples);
			double[] exp = expR.getGeneAlterationArray(name, samples);
			int[] cna = cnaR.getExpVerifiedCNA(name, samples, exp, 0.05);

			if ((mut == null || ArrayUtil.countValue(mut, true) == 0) &&
				(cna == null || ArrayUtil.countValue(cna, 0) == cna.length)) continue;

			if (!considerExpressed(exp, getAltered(mut, cna), 6, 0.05))
			{
//				System.out.println(name + " is considered not expressed.");
				continue;
			}

			wM.write("\n" + name);
			for (int i = 0; i < samples.length; i++)
			{
				boolean m = mut != null && mut[i];
				boolean a = cna != null && cna[i] > 0;
				boolean d = cna != null && cna[i] < 0;

				int v = m ? a ? 4 : d ? 5 : 1 : a ? 2 : d ? 3 : 0;

				wM.write("\t" + v);
			}

			wR.write("\n" + gene.toString());
		}

		wM.close();
		wR.close();
	}

	private static List<Gene> getScoredGenes(Set<String> geneSet, String[] samples,
		CNAReader cnaR, ExpressionReader expR,
		Map<String, Double> mutsigScores, Map<String, Double> ampScores, Map<String, Double> delScores)
	{
		List<Gene> genes = new ArrayList<>();
		for (String name : geneSet)
		{
			Gene gene = new Gene(name);
			gene.mutPval = mutsigScores.containsKey(name) ? mutsigScores.get(name) : 1;
			gene.ampGisticPval = ampScores.containsKey(name) ? ampScores.get(name) : 1;
			gene.delGisticPval = delScores.containsKey(name) ? delScores.get(name) : 1;
			double[] exp = expR.getGeneAlterationArray(name, samples);
			gene.ampExpPval = cnaR.getCNAPval(name, samples, true, exp);
			gene.delExpPval = cnaR.getCNAPval(name, samples, false, exp);
			gene.complete();
			genes.add(gene);
		}
		Collections.sort(genes);
		return genes;
	}

	private static Map<String, Double> readCNAScores(String filename) throws FileNotFoundException
	{
		Map<String, Double> map = new HashMap<>();
		Scanner sc = new Scanner(new File(filename));
		sc.nextLine();
		String[] token = sc.nextLine().split("\t");
		double[] pvals = new double[token.length - 1];
		for (int i = 1; i < token.length; i++)
		{
			pvals[i - 1] = Double.parseDouble(token[i]);
		}
		sc.nextLine();
		sc.nextLine();
		while (sc.hasNextLine())
		{
			token = sc.nextLine().split("\t");
			for (int i = 1; i < token.length; i++)
			{
				map.put(token[i], pvals[i - 1]);
			}
		}
		return map;
	}

	/**
	 * This methods checks if the mean expression of the gene is above a threshold, OR
	 * the altered samples are differently expressed. If one of those holds, then the gene is
	 * considered expressed in the dataset.
	 * @param exp
	 * @param altered
	 * @param meanExpThr
	 * @param diffExpPvalThr
	 * @return
	 */
	private static boolean considerExpressed(double[] exp, boolean[] altered, double meanExpThr,
		double diffExpPvalThr)
	{
		double mean = Summary.mean(exp);
		if (mean > meanExpThr) return true;

		double[] e1 = ArrayUtil.subset(exp, altered);
		if (e1.length == 0) return false;
		double[] e0 = ArrayUtil.subset(exp, ArrayUtil.negate(altered));
		if (e0.length == 0) return false;

		double pval = TTest.getPValOfMeanDifference(e1, e0);
		return pval <= diffExpPvalThr;
	}

	private static boolean[] getAltered(boolean[] mut, int[] cna)
	{
		boolean[] alt = new boolean[mut.length];
		System.arraycopy(mut, 0, alt, 0, mut.length);

		if (cna == null) return alt;

		for (int i = 0; i < cna.length; i++)
		{
			if (cna[i] != 0) alt[i] = true;
		}
		return alt;
	}

	static class Gene implements Comparable
	{
		String name;
		double mutPval = 1;
		double ampGisticPval = 1;
		double ampExpPval = 1;
		double delGisticPval = 1;
		double delExpPval = 1;
		double ampPval = 1;
		double delPval = 1;

		public Gene(String name)
		{
			this.name = name;
		}

		@Override
		public int compareTo(Object o)
		{
			if (o instanceof Gene)
			{
				Gene g = (Gene) o;

				return new Double(score()).compareTo(g.score());
			}
			else return 0;
		}

		public double score()
		{
			return Math.min(mutPval, Math.min(ampPval, delPval));
		}

		@Override
		public String toString()
		{
			return name + "\t" + mutPval + "\t" +
				ampPval + "\t" + ampGisticPval + "\t" + ampExpPval + "\t" +
				delPval + "\t" + delGisticPval + "\t" + delExpPval ;
		}

		public void complete()
		{
			ampPval = Math.max(ampExpPval, ampGisticPval);
			delPval = Math.max(delExpPval, delGisticPval);
		}

		static String getHeader()
		{
			return "Gene\tMutSig\tAmp\tAmp Gistic\tAmp exp\tDel\tDel Gistic\tDel exp";
		}
	}

	public static void main(String[] args) throws IOException
	{
//		String dir = "/home/ozgun/Documents/TCGA/SARC/mutex/";
//		String mutFile = dir + "genome.wustl.edu_SARC.IlluminaHiSeq_DNASeq_automated.1.1.0.somatic.maf";
//		String cnaFile = dir + "Gistic/UPS/all_thresholded.by_genes.txt";
//		String expFile = dir + "SARC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt";
//
//		prepare(mutFile, cnaFile, expFile, dir + "UPS/");

//		String study = "STES";
//		String dir = "/home/babur/Documents/TCGA/" + study;
//		String out = "/home/babur/Documents/mutex/TCGA/" + study;
//		prepare(dir, out);

		String dir = "/home/babur/Documents/TCGA";
		String out = "/home/babur/Documents/mutex/TCGA";
		prepareMulti(dir, out);
	}
}
