package org.cbio.mutex;

import org.cbio.causality.data.portal.CBioPortalAccessor;
import org.cbio.causality.data.portal.CaseList;
import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.Change;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * This class prepares text files to visualize at http://www.cbioportal.org/public-portal/tools.jsp.
 *
 * @author Ozgun Babur
 */
public class Oncoprint
{
	private List<Group> groups;
	PortalDataset dataset;
	public static final String DIR = "result/oncoprint/";
	public static String dir;
	CBioPortalAccessor portal;

	public Oncoprint(List<Group> groups, PortalDataset dataset) throws IOException
	{
		this.groups = groups;
		this.dataset = dataset;
		dir = DIR + dataset.study + "/";
		if (!new File(dir).exists()) new File(dir).mkdirs();
		portal = PortalReader.getPortalAccessor(dataset);
	}

	public void write() throws IOException
	{
		CaseList cases = portal.getCurrentCaseList();
		for (Group group : groups)
		{
			writeData(group, cases);
		}
	}

	private List<String> getCaseIDs(CaseList cases, boolean[] hyper)
	{
		assert cases.getCases().length == hyper.length;

		List<String> list = new ArrayList<String>();
		for (int i = 0; i < hyper.length; i++)
		{
			if (!hyper[i]) list.add(cases.getCases()[i]);
		}
		return list;
	}

	private void writeData(Group group, CaseList cases) throws IOException
	{
		String gID = group.getID();

		BufferedWriter writer = new BufferedWriter(new FileWriter(dir + gID + "-mut.txt"));

		writer.write("Hugo_Symbol\tEntrez_Gene_Id\tsample_id\tProtein_Change");
		boolean[] hyper = group.members.get(0).getHyper();

		for (GeneAlt gene : group.members)
		{
			Change[] muts = gene.getGene().get(Alteration.MUTATION);

			for (int i = 0; i < muts.length; i++)
			{
				if (hyper[i]) continue;

				if (muts[i].isAltered()) writer.write("\n" + gene.getId() + "\t0\t" + cases.getCases()[i] + "\tK1E");
			}
		}
		writer.close();

		writer = new BufferedWriter(new FileWriter(dir + gID + "-cna.txt"));

		writer.write("Hugo_Symbol\tEntrez_Gene_Id");

		for (int i = 0; i < hyper.length; i++)
		{
			if (!hyper[i]) writer.write("\t" + cases.getCases()[i]);
		}

		for (GeneAlt gene : group.members)
		{
			writer.write("\n" + gene.getId() + "\t0");
			Change[] cnas = gene.getGene().get(Alteration.COPY_NUMBER);

			for (int i = 0; i < hyper.length; i++)
			{
				if (!hyper[i]) writer.write("\t" + (cnas[i] == Change.ACTIVATING ? "2" :
					cnas[i] == Change.INHIBITING ? "-2" : "0"));
			}
		}

		writer.close();
	}
}
