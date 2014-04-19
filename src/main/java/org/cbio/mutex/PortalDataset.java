package org.cbio.mutex;

/**
 * A dataset from the portal.
 * @author Ozgun Babur
 */
public class PortalDataset
{
	/**
	 * Constructor with parameters.
	 * @param study index of cancer study in portal
	 * @param caseList index of case-list in portal
	 * @param profile indices of the desired genomic profiles
	 */
	public PortalDataset(String name, String study, String caseList, String[] profile, String[] subtypes)
	{
		this.name = name;
		this.study = study;
		this.caseList = caseList;
		this.profile = profile;
		this.subtypeCases = subtypes;
	}

	String name;
	String study;
	String caseList;
	String[] profile;
	String[] subtypeCases;

	@Override
	public String toString()
	{
		String s = "Study: \"" + study + "\" ";
		return null;
	}

	public static final PortalDataset GBM = new PortalDataset("glioblastoma", "gbm_tcga", "gbm_tcga_3way_complete", new String[]{"gbm_tcga_mutations", "gbm_tcga_gistic", "gbm_tcga_rna_seq_v2_mrna"}, new String[]{});
	public static final PortalDataset OV = new PortalDataset("ovarian", "ov_tcga", "ov_tcga_cnaseq", new String[]{"ov_tcga_mutations", "ov_tcga_gistic", "ov_tcga_mrna_merged_median_Zscores"}, new String[]{});
	public static final PortalDataset BRCA = new PortalDataset("breast", "brca_tcga", "brca_tcga_3way_complete", new String[]{"brca_tcga_mutations", "brca_tcga_gistic", "brca_tcga_rna_seq_v2_mrna"}, new String[]{});
	public static final PortalDataset BRCA_PUB = new PortalDataset("breast-pub", "brca_tcga_pub", "brca_tcga_pub_complete", new String[]{"brca_tcga_pub_mutations", "brca_tcga_pub_gistic", "brca_tcga_pub_mrna"}, new String[]{"brca_tcga_pub_basal", "brca_tcga_pub_her2", "brca_tcga_pub_luma", "brca_tcga_pub_lumb", "brca_tcga_pub_claudin"});
	public static final PortalDataset COADREAD = new PortalDataset("colorectal", "coadread_tcga_pub", "coadread_tcga_pub_3way_complete", new String[]{"coadread_tcga_pub_mutations", "coadread_tcga_pub_gistic", "coadread_tcga_pub_rna_seq_mrna"}, new String[]{});
	public static final PortalDataset UCEC = new PortalDataset("endometrial", "ucec_tcga_pub", "ucec_tcga_pub_3way_complete", new String[]{"ucec_tcga_pub_mutations", "ucec_tcga_pub_gistic", "ucec_tcga_pub_rna_seq_v2_mrna"}, new String[]{"ucec_tcga_pub_cnhigh", "ucec_tcga_pub_cnlow", "ucec_tcga_pub_msi", "ucec_tcga_pub_pole"});
	public static final PortalDataset THCA = new PortalDataset("thyroid", "thca_tcga", "thca_tcga_3way_complete", new String[]{"thca_tcga_mutations", "thca_tcga_gistic", "thca_tcga_rna_seq_v2_mrna"}, new String[]{});
	public static final PortalDataset SIMUL = new PortalDataset("simulated", "simulated", "simulated_caseset", new String[]{"simulated_profile"}, new String[]{});
}
