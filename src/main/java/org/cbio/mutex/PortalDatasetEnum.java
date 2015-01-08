package org.cbio.mutex;

/**
 * A dataset from the portal.
 * @author Ozgun Babur
 */
public enum PortalDatasetEnum
{
	// TCGA datasets

	GBM_PUB("GBM-pub", 0.01,
	"gbm_tcga_pub", "gbm_tcga_pub_3way_complete", new String[]{"gbm_tcga_pub_mutations", "gbm_tcga_pub_cna_consensus", "gbm_tcga_pub_mrna"},
	new String[]{"gbm_tcga_pub_expr_classical", "gbm_tcga_pub_expr_mesenchymal", "gbm_tcga_pub_expr_neural", "gbm_tcga_pub_expr_proneural"}, null),

	GBM("GBM", 0.01,
		"gbm_tcga", "gbm_tcga_3way_complete", new String[]{"gbm_tcga_mutations", "gbm_tcga_gistic", "gbm_tcga_rna_seq_v2_mrna"},
		new String[]{}, GBM_PUB),

	OV("OV", 0.06,
		"ov_tcga", "ov_tcga_3way_complete", new String[]{"ov_tcga_mutations", "ov_tcga_gistic", "ov_tcga_rna_seq_v2_mrna"},
		new String[]{}, null),

	BRCA_PUB("BRCA-pub", 0.03,
		"brca_tcga_pub", "brca_tcga_pub_complete", new String[]{"brca_tcga_pub_mutations", "brca_tcga_pub_gistic", "brca_tcga_pub_mrna"},
		new String[]{"brca_tcga_pub_basal", "brca_tcga_pub_her2", "brca_tcga_pub_luma", "brca_tcga_pub_lumb", "brca_tcga_pub_claudin"}, null),

	BRCA("BRCA", 0.04,
		"brca_tcga", "brca_tcga_3way_complete", new String[]{"brca_tcga_mutations", "brca_tcga_gistic", "brca_tcga_rna_seq_v2_mrna"},
		new String[]{}, BRCA_PUB),

	COADREAD("COADREAD", 0.02,
		"coadread_tcga_pub", "coadread_tcga_pub_3way_complete", new String[]{"coadread_tcga_pub_mutations", "coadread_tcga_pub_gistic", "coadread_tcga_pub_rna_seq_mrna"},
		new String[]{}, null),

	UCEC("UCEC", 0.04,
		"ucec_tcga_pub", "ucec_tcga_pub_3way_complete", new String[]{"ucec_tcga_pub_mutations", "ucec_tcga_pub_gistic", "ucec_tcga_pub_rna_seq_v2_mrna"},
		new String[]{"ucec_tcga_pub_cnhigh", "ucec_tcga_pub_cnlow", "ucec_tcga_pub_msi", "ucec_tcga_pub_pole"}, null),

	THCA("THCA", 0.01,
		"thca_tcga", "thca_tcga_3way_complete", new String[]{"thca_tcga_mutations", "thca_tcga_gistic", "thca_tcga_rna_seq_v2_mrna"},
		new String[]{}, null),

	STAD("STAD", 0.03,
		"stad_tcga_pub", "stad_tcga_pub_3way_complete", new String[]{"stad_tcga_pub_mutations", "stad_tcga_pub_gistic", "stad_tcga_pub_rna_seq_v2_mrna"},
		new String[]{}, null),

	SKCM("SKCM", 0.06,
		"skcm_tcga", "skcm_tcga_3way_complete", new String[]{"skcm_tcga_mutations", "skcm_tcga_gistic", "skcm_tcga_rna_seq_v2_mrna"},
		new String[]{}, null),

	LAML("LAML", 0.01,
		"laml_tcga", "laml_tcga_3way_complete", new String[]{"laml_tcga_mutations", "laml_tcga_gistic", "laml_tcga_rna_seq_v2_mrna"},
		new String[]{}, null),

	ACC("ACC", 0.02,
		"acc_tcga", "acc_tcga_3way_complete", new String[]{"acc_tcga_mutations", "acc_tcga_gistic", "acc_tcga_rna_seq_v2_mrna"},
		new String[]{}, null),

	LGG("LGG", 0.01,
		"lgg_tcga", "lgg_tcga_3way_complete", new String[]{"lgg_tcga_mutations", "lgg_tcga_gistic", "lgg_tcga_rna_seq_v2_mrna"},
		new String[]{}, null),

	HNSC("HNSC", 0.04,
		"hnsc_tcga", "hnsc_tcga_3way_complete", new String[]{"hnsc_tcga_mutations", "hnsc_tcga_gistic", "hnsc_tcga_rna_seq_v2_mrna"},
		new String[]{}, null),

	LUAD("LUAD", 0.05,
		"luad_tcga", "luad_tcga_3way_complete", new String[]{"luad_tcga_mutations", "luad_tcga_gistic", "luad_tcga_rna_seq_v2_mrna"},
		new String[]{}, null),

	LUSC("LUSC", 0.04,
		"lusc_tcga_pub", "lusc_tcga_pub_3way_complete", new String[]{"lusc_tcga_pub_mutations", "lusc_tcga_pub_gistic", "lusc_tcga_pub_rna_seq_mrna"},
		new String[]{}, null),

	KIRC("KIRC", 0.01,
		"kirc_tcga", "kirc_tcga_3way_complete", new String[]{"kirc_tcga_mutations", "kirc_tcga_gistic", "kirc_tcga_rna_seq_v2_mrna"},
		new String[]{}, null),

	KIRP("KIRP", 0.01,
		"kirp_tcga", "kirp_tcga_3way_complete", new String[]{"kirp_tcga_mutations", "kirp_tcga_gistic", "kirp_tcga_rna_seq_v2_mrna"},
		new String[]{}, null),

	KICH("KICH", 0.01,
		"kich_tcga", "kich_tcga_3way_complete", new String[]{"kich_tcga_mutations", "kich_tcga_gistic", "kich_tcga_rna_seq_v2_mrna"},
		new String[]{}, null),

	PRAD("PRAD", 0.01,
		"prad_tcga", "prad_tcga_3way_complete", new String[]{"prad_tcga_mutations", "prad_tcga_gistic", "prad_tcga_rna_seq_v2_mrna"},
		new String[]{}, null),

	CESC("CESC", 0.01,
		"cesc_tcga", "cesc_tcga_3way_complete", new String[]{"cesc_tcga_mutations", "cesc_tcga_gistic", "cesc_tcga_rna_seq_v2_mrna"},
		new String[]{}, null),

	LIHC("LIHC", 0.01,
		"lihc_tcga", "lihc_tcga_3way_complete", new String[]{"lihc_tcga_mutations", "lihc_tcga_gistic", "lihc_tcga_rna_seq_v2_mrna"},
		new String[]{}, null),

	BLCA("BLCA", 0.01,
		"blca_tcga", "blca_tcga_3way_complete", new String[]{"blca_tcga_mutations", "blca_tcga_gistic", "blca_tcga_rna_seq_v2_mrna"},
		new String[]{}, null),

	UCS("UCS", 0.01,
		"ucs_tcga", "ucs_tcga_3way_complete", new String[]{"ucs_tcga_mutations", "ucs_tcga_gistic", "ucs_tcga_rna_seq_v2_mrna"},
		new String[]{}, null),

//	SARC("sarcoma", 0.01,
//		"sarc_mskcc", "sarc_mskcc_complete", new String[]{"sarc_mskcc_mutations", "sarc_mskcc_cna", "sarc_mskcc_mrna"},
//		new String[]{"sarc_mskcc_ddlps", "sarc_mskcc_gist", "sarc_mskcc_lms", "sarc_mskcc_mfh", "sarc_mskcc_myxoid", "sarc_mskcc_pleo", "sarc_mskcc_synovial"}, null),

	// Simulated datasets

	SIMUL1("simulated", 0.03, "simulated_" + "ng50_gs3_ss958", "simulated_caseset", new String[]{"simulated_profile_" + "ng50_gs3_ss958"}, new String[]{}, null),
	SIMUL2("simulated", 0.03, "simulated_" + "ng20_gs3_ss463", "simulated_caseset", new String[]{"simulated_profile_" + "ng20_gs3_ss463"}, new String[]{}, null),
	SIMUL3("simulated", 0.03, "simulated_" + "ng10_gs6_ss463", "simulated_caseset", new String[]{"simulated_profile_" + "ng10_gs6_ss463"}, new String[]{}, null);

	PortalDataset data;

	/**
	 * Constructor with parameters.
	 * @param study index of cancer study in portal
	 * @param caseList index of case-list in portal
	 * @param profile indices of the desired genomic profiles
	 */
	PortalDatasetEnum(String name, double minAltThr, String study, String caseList,
		String[] profile, String[] subtypes, PortalDatasetEnum subtypeProvider)
	{
		this.data = new PortalDataset(name, minAltThr, study, caseList, profile, subtypes, subtypeProvider == null ? null : subtypeProvider.data);
	}

	public static PortalDataset find(String name)
	{
		for (PortalDatasetEnum data : values())
		{
			if (data.data.name.equals(name)) return data.data;
		}
		return null;
	}

	public static PortalDataset findByStudyID(String studyID)
	{
		for (PortalDatasetEnum data : values())
		{
			if (data.data.study.equals(studyID)) return data.data;
		}
		return null;
	}
}
