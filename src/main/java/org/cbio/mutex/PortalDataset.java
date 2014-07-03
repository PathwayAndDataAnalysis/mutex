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
	public PortalDataset(String name, double minAltThr, String study, String caseList, String[] profile, String[] subtypes)
	{
		this.name = name;
		this.minAltThr = minAltThr;
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
	double minAltThr;

	private static final String S1 = "ng50_gs3_ss958";
	private static final String S2 = "ng20_gs3_ss463";
	private static final String S3 = "ng10_gs6_ss463";

	public static final PortalDataset GBM = new PortalDataset("glioblastoma", 0.02, "gbm_tcga", "gbm_tcga_3way_complete", new String[]{"gbm_tcga_mutations", "gbm_tcga_gistic", "gbm_tcga_rna_seq_v2_mrna"}, new String[]{});
	public static final PortalDataset GBM_PUB = new PortalDataset("glioblastoma-pub", 0.01, "gbm_tcga_pub", "gbm_tcga_pub_3way_complete", new String[]{"gbm_tcga_pub_mutations", "gbm_tcga_pub_cna_consensus", "gbm_tcga_pub_mrna"}, new String[]{"gbm_tcga_pub_expr_classical", "gbm_tcga_pub_expr_mesenchymal", "gbm_tcga_pub_expr_neural", "gbm_tcga_pub_expr_proneural"});
	public static final PortalDataset OV = new PortalDataset("ovarian", 0.05, "ov_tcga", "ov_tcga_3way_complete", new String[]{"ov_tcga_mutations", "ov_tcga_gistic", "ov_tcga_rna_seq_v2_mrna"}, new String[]{});
	public static final PortalDataset BRCA = new PortalDataset("breast", 0.03, "brca_tcga", "brca_tcga_3way_complete", new String[]{"brca_tcga_mutations", "brca_tcga_gistic", "brca_tcga_rna_seq_v2_mrna"}, new String[]{});
	public static final PortalDataset BRCA_PUB = new PortalDataset("breast-pub", 0.03, "brca_tcga_pub", "brca_tcga_pub_complete", new String[]{"brca_tcga_pub_mutations", "brca_tcga_pub_gistic", "brca_tcga_pub_mrna"}, new String[]{"brca_tcga_pub_basal", "brca_tcga_pub_her2", "brca_tcga_pub_luma", "brca_tcga_pub_lumb", "brca_tcga_pub_claudin"});
	public static final PortalDataset COADREAD = new PortalDataset("colorectal", 0.02,  "coadread_tcga_pub", "coadread_tcga_pub_3way_complete", new String[]{"coadread_tcga_pub_mutations", "coadread_tcga_pub_gistic", "coadread_tcga_pub_rna_seq_mrna"}, new String[]{});
	public static final PortalDataset UCEC = new PortalDataset("endometrial", 0.03, "ucec_tcga_pub", "ucec_tcga_pub_3way_complete", new String[]{"ucec_tcga_pub_mutations", "ucec_tcga_pub_gistic", "ucec_tcga_pub_rna_seq_v2_mrna"}, new String[]{"ucec_tcga_pub_cnhigh", "ucec_tcga_pub_cnlow", "ucec_tcga_pub_msi", "ucec_tcga_pub_pole"});
	public static final PortalDataset THCA = new PortalDataset("thyroid", 0.01, "thca_tcga", "thca_tcga_3way_complete", new String[]{"thca_tcga_mutations", "thca_tcga_gistic", "thca_tcga_rna_seq_v2_mrna"}, new String[]{});
	public static final PortalDataset STAD = new PortalDataset("stomach", 0.03, "stad_tcga", "stad_tcga_3way_complete", new String[]{"stad_tcga_mutations", "stad_tcga_gistic", "stad_tcga_rna_seq_v2_mrna"}, new String[]{});
	public static final PortalDataset SKCM = new PortalDataset("skin", 0.03, "skcm_tcga", "skcm_tcga_3way_complete", new String[]{"skcm_tcga_mutations", "skcm_tcga_gistic", "skcm_tcga_rna_seq_v2_mrna"}, new String[]{});
	public static final PortalDataset LAML = new PortalDataset("leukemia", 0.01, "laml_tcga", "laml_tcga_3way_complete", new String[]{"laml_tcga_mutations", "laml_tcga_gistic", "laml_tcga_rna_seq_v2_mrna"}, new String[]{});
	public static final PortalDataset ACC = new PortalDataset("adrenocortical", 0.01, "acc_tcga", "acc_tcga_3way_complete", new String[]{"acc_tcga_mutations", "acc_tcga_gistic", "acc_tcga_rna_seq_v2_mrna"}, new String[]{});
	public static final PortalDataset LGG = new PortalDataset("lower-grade-glioma", 0.01, "lgg_tcga", "lgg_tcga_3way_complete", new String[]{"lgg_tcga_mutations", "lgg_tcga_gistic", "lgg_tcga_rna_seq_v2_mrna"}, new String[]{});
	public static final PortalDataset HNSC = new PortalDataset("head-and-neck", 0.03, "hnsc_tcga", "hnsc_tcga_3way_complete", new String[]{"hnsc_tcga_mutations", "hnsc_tcga_gistic", "hnsc_tcga_rna_seq_v2_mrna"}, new String[]{});
	public static final PortalDataset LUAD = new PortalDataset("lung-adeno", 0.05, "luad_tcga_pub", "luad_tcga_pub_3way_complete", new String[]{"luad_tcga_pub_mutations", "luad_tcga_pub_gistic", "luad_tcga_pub_rna_seq_v2_mrna"}, new String[]{});
	public static final PortalDataset LUSC = new PortalDataset("lung-squamous", 0.03, "lusc_tcga_pub", "lusc_tcga_pub_3way_complete", new String[]{"lusc_tcga_pub_mutations", "lusc_tcga_pub_gistic", "lusc_tcga_pub_rna_seq_mrna"}, new String[]{});
	public static final PortalDataset KIRC = new PortalDataset("kidney-renal", 0.01, "kirc_tcga", "kirc_tcga_3way_complete", new String[]{"kirc_tcga_mutations", "kirc_tcga_gistic", "kirc_tcga_rna_seq_v2_mrna"}, new String[]{});
	public static final PortalDataset KIRP = new PortalDataset("kidney-papillary", 0.02, "kirp_tcga", "kirp_tcga_3way_complete", new String[]{"kirp_tcga_mutations", "kirp_tcga_gistic", "kirp_tcga_rna_seq_v2_mrna"}, new String[]{});
	public static final PortalDataset PRAD = new PortalDataset("prostate", 0.03, "prad_tcga", "prad_tcga_3way_complete", new String[]{"prad_tcga_mutations", "prad_tcga_gistic", "prad_tcga_rna_seq_v2_mrna"}, new String[]{});
	public static final PortalDataset SARC = new PortalDataset("sarcoma", 0.01, "sarc_mskcc", "sarc_mskcc_complete", new String[]{"sarc_mskcc_mutations", "sarc_mskcc_cna", "sarc_mskcc_mrna"}, new String[]{"sarc_mskcc_ddlps", "sarc_mskcc_gist", "sarc_mskcc_lms", "sarc_mskcc_mfh", "sarc_mskcc_myxoid", "sarc_mskcc_pleo", "sarc_mskcc_synovial"});
	public static final PortalDataset SIMUL1 = new PortalDataset("simulated", 0.03, "simulated_" + S1, "simulated_caseset", new String[]{"simulated_profile_" + S1}, new String[]{});
	public static final PortalDataset SIMUL2 = new PortalDataset("simulated", 0.03, "simulated_" + S2, "simulated_caseset", new String[]{"simulated_profile_" + S2}, new String[]{});
	public static final PortalDataset SIMUL3 = new PortalDataset("simulated", 0.03, "simulated_" + S3, "simulated_caseset", new String[]{"simulated_profile_" + S3}, new String[]{});
}
