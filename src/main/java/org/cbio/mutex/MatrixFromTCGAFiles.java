package org.cbio.mutex;

import org.cbio.causality.data.tcgafile.CNAReader;
import org.cbio.causality.data.tcgafile.ExpressionReader;
import org.cbio.causality.data.tcgafile.MutationReader;
import org.cbio.causality.util.ArrayUtil;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * This class is used to prepare the alteration data matrix from muation (maf), CNA (gistic), and
 * expression (RNAseq) data files.
 *
 * @author Ozgun Babur
 */
public class MatrixFromTCGAFiles
{
	public static void prepare(String mutFile, String CNAFile, String expFile, String outDir) throws IOException
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

		Set<String> genes = new HashSet<String>(mutR.getGenes());
		Set<String> cnaGenes = new HashSet<String>(cnaR.getGenes());
		cnaGenes.retainAll(expR.getGenes());
		genes.addAll(cnaGenes);

		System.out.println("genes = " + genes.size());

		BufferedWriter writer = new BufferedWriter(new FileWriter(outDir + "DataMatrix.txt"));

		for (String sample : samples)
		{
			writer.write("\t" + sample);
		}

		for (String gene : genes)
		{
			boolean[] mut = mutR.getGeneAlterationArray(gene, samples);
			double[] exp = expR.getGeneAlterationArray(gene, samples);
			int[] cna = cnaR.getExpVerifiedCNA(gene, samples, exp, 0.05);

			if ((mut == null || ArrayUtil.countValue(mut, true) == 0) &&
				(cna == null || ArrayUtil.countValue(cna, 0) == cna.length)) continue;


			writer.write("\n" + gene);
			for (int i = 0; i < samples.length; i++)
			{
				boolean m = mut != null && mut[i];
				boolean a = cna != null && cna[i] > 0;
				boolean d = cna != null && cna[i] < 0;

				int v = m ? a ? 4 : d ? 5 : 1 : a ? 2 : d ? 3 : 0;

				writer.write("\t" + v);
			}
		}

		writer.close();

	}

	public static void main(String[] args) throws IOException
	{
		String dir = "/home/ozgun/Documents/TCGA/SARC/mutex/";
		String mutFile = dir + "genome.wustl.edu_SARC.IlluminaHiSeq_DNASeq_automated.1.1.0.somatic.maf";
		String cnaFile = dir + "Gistic/UPS/all_thresholded.by_genes.txt";
		String expFile = dir + "SARC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt";

		prepare(mutFile, cnaFile, expFile, dir + "UPS/");
	}
}
