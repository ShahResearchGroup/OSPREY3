package edu.duke.cs.osprey.tools;

import java.util.HashMap;
import java.util.Map;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.sparse.sequence.Sequence;

public class ResidueIndexMap {

	private Map<Integer, Integer> designIndexToPDBIndex = new HashMap<>();
	private Map<Integer, Integer> PDBIndexToDesignIndex = new HashMap<>();
	private Map<String, String> RCToAAName = new HashMap<>();
	private Map<Integer, Map<String, Integer>> AAAssignmentToRC = new HashMap<>();
	private Map<Integer, Map<String, Integer>> numConfsAtAAAtPos = new HashMap<>();
	
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
				if(!AAAssignmentToRC.containsKey(designIndex))
					AAAssignmentToRC.put(designIndex, new HashMap<>());
				if(!AAAssignmentToRC.get(designIndex).containsKey(rotamer.AAType))
					AAAssignmentToRC.get(designIndex).put(rotamer.AAType, rotamer.RCIndex);
				if(!numConfsAtAAAtPos.containsKey(designIndex))
					numConfsAtAAAtPos.put(designIndex, new HashMap<>());
				if(!numConfsAtAAAtPos.get(designIndex).containsKey(rotamer.AAType))
					numConfsAtAAAtPos.get(designIndex).put(rotamer.AAType,0);
				numConfsAtAAAtPos.get(designIndex).
					put(rotamer.AAType,numConfsAtAAAtPos.get(designIndex).get(rotamer.AAType) + 1);
			}
		}
	}
	
	public Sequence getSequenceOfRCTuple(RCTuple conf)
	{
		Sequence output = new Sequence();
		for(int i = 0; i < conf.size(); i++)
		{
			output.addAssignment(conf.pos.get(i), RCToAAName.get(conf.pos.get(i)+":"+conf.RCs.get(i)));
		}
		return output;
	}
	
	public int getConfsForAAAtPos(int designIndex, String AA)
	{
		return numConfsAtAAAtPos.get(designIndex).get(AA);
	}
	
	public RCTuple createTemplateConfForSequence(Sequence seq)
	{
		RCTuple output = new RCTuple();
		for(Integer residueIndex : seq.residues())
		{
			output.addRC(residueIndex, AAAssignmentToRC.get(residueIndex).get(seq.getAAAt(residueIndex)));
		}
		return output;
	}
}
