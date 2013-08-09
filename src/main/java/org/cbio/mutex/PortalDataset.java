package org.cbio.mutex;

/**
 * A dataset from the portal.
 * @author Ozgun Babur
 */
public class PortalDataset
{
	/**
	 * Constructor with parameters.
	 * @param name an identifier name of the dataset
	 * @param filename filename for local cash
	 * @param study index of cancer study in portal
	 * @param caseList index of case-list in portal
	 * @param profile indices of the desired genomic profiles
	 */
	public PortalDataset(String name, String filename,
		String study, String caseList, String... profile)
	{
		this.name = name;
		this.filename = filename;
		this.study = study;
		this.caseList = caseList;
		this.profile = profile;
	}

	String name;
	String filename;

	String study;
	String caseList;
	String[] profile;

	public static final PortalDataset glioblastoma = new PortalDataset(
		"glioblastoma", "Glioblastoma.txt", "gbm_tcga_pub", "gbm_tcga_pub_cnaseq", "gbm_tcga_pub_mutations", "gbm_tcga_pub_cna_rae");
	public static final PortalDataset ovarian = new PortalDataset(
		"ovarian", "Ovarian.txt", "ov_tcga_pub", "ov_tcga_pub_cna_seq", "ov_tcga_pub_mutations", "ov_tcga_pub_gistic");
	public static final PortalDataset breast = new PortalDataset(
		"breast", "Breast.txt", "brca_tcga_pub", "brca_tcga_pub_cnaseq", "brca_tcga_pub_mutations", "brca_tcga_pub_gistic");
	public static final PortalDataset colon = new PortalDataset(
		"colon", "Colon.txt", "coadread_tcga_pub", "coadread_tcga_pub_cna_seq", "coadread_tcga_pub_mutations", "coadread_tcga_pub_gistic");
	public static final PortalDataset endometrial = new PortalDataset(
		"endometrial", "Endometrial.txt", "ucec_tcga_pub", "ucec_tcga_pub_cnaseq", "ucec_tcga_pub_mutations", "ucec_tcga_pub_gistic");
}
