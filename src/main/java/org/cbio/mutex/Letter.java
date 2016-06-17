package org.cbio.mutex;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Ozgun Babur
 */
public enum Letter
{
	NO_ALT(0, ".", "No alteration"),
	MUT(1, "M", "Mutated"),
	AMP(2, "A", "Amplified"),
	DEL(3, "D", "Deleted"),
	AMP_MUT(4, "B", "Amplified and mutated"),
	DEL_MUT(5, "E", "Deleted and mutated");

	int code;
	String description;
	String print;

	private Letter(int code, String print, String description)
	{
		this.code = code;
		this.print = print;
		this.description = description;
	}

	private static Map<Integer, Letter> map;

	static
	{
		map = new HashMap<>();
		for (Letter let : Letter.values())
		{
			map.put(let.code, let);
		}
	}

	static Letter letterOf(int code)
	{
		return map.get(code);
	}

	static Letter[] convertToLetter(int[] alts)
	{
		Letter[] letters = new Letter[alts.length];
		for (int i = 0; i < alts.length; i++)
		{
			letters[i] = letterOf(alts[i]);
		}
		return letters;
	}
}
