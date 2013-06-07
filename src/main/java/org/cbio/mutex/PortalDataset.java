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
	public PortalDataset(String name, String filename, int study, int caseList, int[] profile)
	{
		this.name = name;
		this.filename = filename;
		this.study = study;
		this.caseList = caseList;
		this.profile = profile;
	}

	String name;
	String filename;
	int study;
	int caseList;
	int[] profile;

	public static final PortalDataset glioblastoma = new PortalDataset(
		"glioblastoma", "Glioblastoma.txt", 7, 0, new int[]{1, 3});
	public static final PortalDataset ovarian = new PortalDataset(
		"ovarian", "Ovarian.txt", 16, 0, new int[]{0, 12});
	public static final PortalDataset breast = new PortalDataset(
		"breast", "Breast.txt", 3, 1, new int[]{0, 8});
	public static final PortalDataset colon = new PortalDataset(
		"colon", "Colon.txt", 5, 0, new int[]{0, 10});
}
