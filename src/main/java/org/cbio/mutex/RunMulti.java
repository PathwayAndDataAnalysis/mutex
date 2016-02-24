package org.cbio.mutex;

import java.io.File;
import java.io.IOException;

/**
 * Created by babur on 2/22/16.
 */
public class RunMulti
{
	public static void main(String[] args) throws IOException, ClassNotFoundException
	{
		if (args.length < 1) return;
		String dir = args[0];


		for (File file : new File(dir).listFiles())
		{
			if (file.isDirectory())
			{
				String[] a = new String[args.length];
				a[0] = file.getPath();

				System.arraycopy(args, 1, a, 1, args.length - 1);

				if (new File(file.getPath() + "/parameters.txt").exists())
				{
					Main.main(a);
				}
				else
				{
					main(a);
				}
			}
		}
	}
}
