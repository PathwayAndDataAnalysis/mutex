package org.cbio.mutex;

/**
 * A dataset from the portal.
 * @author Ozgun Babur
 */
public class PortalDataset implements Cloneable
{
	String name;
	String study;
	String caseList;
	String[] profile;
	String[] subtypeCases;
	double minAltThr;
	boolean[] hyper;
	PortalDataset subtypeProvider;

	/**
	 * Constructor with parameters.
	 * @param study index of cancer study in portal
	 * @param caseList index of case-list in portal
	 * @param profile indices of the desired genomic profiles
	 */
	PortalDataset(String name, double minAltThr, String study, String caseList,
		String[] profile, String[] subtypes, PortalDataset subtypeProvider)
	{
		this.name = name;
		this.minAltThr = minAltThr;
		this.study = study;
		this.caseList = caseList;
		this.profile = profile;
		this.subtypeCases = subtypes;
		this.subtypeProvider = subtypeProvider;
	}

	@Override
	protected Object clone() throws CloneNotSupportedException
	{
		PortalDataset cln = (PortalDataset) super.clone();
		if (profile != null) cln.profile = profile.clone();
		if (subtypeCases != null) cln.subtypeCases = subtypeCases.clone();
		if (hyper != null) cln.hyper = hyper.clone();
		return cln;
	}
}
