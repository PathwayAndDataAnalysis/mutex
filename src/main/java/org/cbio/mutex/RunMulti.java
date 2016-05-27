package org.cbio.mutex;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

/**
 * Created by babur on 2/22/16.
 */
public class RunMulti
{
	public static void main(String[] args) throws IOException, ClassNotFoundException
	{
		if (args.length < 1) return;
		String dir = args[0];

		Set<String> avoid = args.length > 1 ? new HashSet<>(Arrays.asList(args).subList(1, args.length)) :
			Collections.<String>emptySet();

		for (File file : new File(dir).listFiles())
		{
			if (file.isDirectory() && !avoid.contains(file.getName()))
			{
				String[] a = new String[args.length];
				a[0] = file.getPath();

				System.arraycopy(args, 1, a, 1, args.length - 1);

				if (Files.exists(Paths.get(file.getPath() + "/parameters.txt")) &&
					!Files.exists(Paths.get(file.getPath() + "/ranked-groups.txt")))
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
