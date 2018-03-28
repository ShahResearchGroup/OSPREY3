package edu.duke.cs.osprey.tools;

import java.util.HashMap;
import java.util.Map;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.RCTuple;

public class ResidueIndexMap {

	private Map<Integer, Integer> designIndexToPDBIndex = new HashMap<>();
	private Map<Integer, Integer> PDBIndexToDesignIndex = new HashMap<>();
	private Map<String, String> RCToAAName = new HashMap<>();
	
	public static ResidueIndexMap createResidueIndexMap(ConfSpace conformationSpace)
	{
		return new ResidueIndexMap(conformationSpace);
	}
	
	public int designIndexToPDBIndex(int designIndex)
	{
		if(!designIndexToPDBIndex.containsKey(designIndex))
		{
			return -1;
		}
		return designIndexToPDBIndex.get(designIndex);
	}
	
	public int PDBIndexToDesignIndex(int PDBIndex)
	{
		if(!PDBIndexToDesignIndex.containsKey(PDBIndex))
		{
			System.err.println("unrecognized PDBIndex.");
		}
		return PDBIndexToDesignIndex.get(PDBIndex);
	}
	
	private ResidueIndexMap (ConfSpace conformationSpace)
	{
		for(PositionConfSpace positionSpace : conformationSpace.posFlex)
		{
			int PDBIndex = positionSpace.res.getPDBIndex();
			int designIndex = positionSpace.designIndex;
			designIndexToPDBIndex.put(designIndex, PDBIndex);
			PDBIndexToDesignIndex.put(PDBIndex, designIndex);
			for(RC rotamer: positionSpace.RCs)
			{
				RCToAAName.put(positionSpace.designIndex+":"+rotamer.RCIndex, rotamer.AAType);
			}
		}
	}
	
	public String getSequenceOfRCTuple(RCTuple conf)
	{
		String output = "";
		for(int i = 0; i < conf.size(); i++)
		{
			output+=conf.pos.get(i)+":"+RCToAAName.get(conf.pos.get(i)+":"+conf.RCs.get(i))+", ";
		}
		return output;
	}
}
