package edu.duke.cs.osprey.sparse;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Map;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.sparse.sequence.Sequence;
import edu.duke.cs.osprey.tools.ResidueIndexMap;

public class Subproblem {
	private Set<Integer> LSet;
	private Set<Integer> lambdaSet;
	private Set<Integer> MSet;
	private Set<Integer> MULambdaSet;
	List<ConformationProcessor> processors = new ArrayList<>();
	public Subproblem leftSubproblem;
	public Subproblem rightSubproblem;
	private ResidueIndexMap residueIndexMap;
	protected RCs localConfSpace;
	private Map<Integer, Map<Integer, Integer>> PDBRCToSubspaceRCMap;
	private Map<Integer, Map<Integer, Integer>> SubspaceRCToPDBRCMap;
	
	public Subproblem (RCs superSpace, TreeEdge sparseTree, ResidueIndexMap resMap) {
		this(superSpace, sparseTree, resMap, null);
	}
	
	public Subproblem (RCs superSpace, TreeEdge sparseTree, ResidueIndexMap resMap, RCTuple initialConf) {
		localConfSpace = superSpace;
		residueIndexMap = resMap;
		
		TreeEdge curEdge = sparseTree;
		lambdaSet = curEdge.getLambda();
		MSet = curEdge.getM();
		LSet = curEdge.getL();
		MULambdaSet = new HashSet<Integer>();
		MULambdaSet.addAll(lambdaSet);
		MULambdaSet.addAll(MSet);
		initSubSpaceRCMap();
		if(initialConf != null)
			localConfSpace = superSpace.returnSubspace(initialConf);

		generateSubproblems(sparseTree);
		
	}

	private void initSubSpaceRCMap () {
		PDBRCToSubspaceRCMap = new HashMap<>();
		SubspaceRCToPDBRCMap = new HashMap<>();
		for(int PDBIndex : MULambdaSet)
		{
			int RCSpaceIndex = residueIndexMap.PDBIndexToDesignIndex(PDBIndex);
			if(localConfSpace == null)
			{
				System.out.println("Null.");
			}
			for(int RCSpaceRCIndex = 0; RCSpaceRCIndex < localConfSpace.getNum(RCSpaceIndex); RCSpaceRCIndex++)
			{
				int PDBSpaceRCIndex = localConfSpace.get(RCSpaceIndex, RCSpaceRCIndex);
				mapPDBRCToRCSpace(PDBIndex, PDBSpaceRCIndex, RCSpaceIndex, RCSpaceRCIndex);
			}
		}
	}

	private void mapPDBRCToRCSpace (int PDBIndex, int PDBSpaceRCIndex, int RCSpaceIndex, int RCSpaceRCIndex) {
		if(!PDBRCToSubspaceRCMap.containsKey(PDBIndex))
			PDBRCToSubspaceRCMap.put(PDBIndex, new HashMap<>());
		PDBRCToSubspaceRCMap.get(PDBIndex).put(PDBSpaceRCIndex, RCSpaceRCIndex);		
		if(!SubspaceRCToPDBRCMap.containsKey(RCSpaceIndex))
			SubspaceRCToPDBRCMap.put(RCSpaceIndex, new HashMap<>());
		SubspaceRCToPDBRCMap.get(RCSpaceIndex).put(RCSpaceRCIndex, PDBSpaceRCIndex);
		
	}
	
	public int mapSubproblemConfToIndex(RCTuple queryAssignment) {
		if(queryAssignment.size() < 1 || MSet.size() < 1)
			return 0;
		int largestDesignIndex = Collections.max(queryAssignment.pos)+1;
		int[] RCSpaceRCIndices = new int[largestDesignIndex];

		Arrays.fill(RCSpaceRCIndices,-1);
		int[] RCSpaceSizes = new int[largestDesignIndex];
		Arrays.fill(RCSpaceSizes,1);
		for(int tupleIndex = 0; tupleIndex < queryAssignment.size(); tupleIndex++) {
			int pos = queryAssignment.pos.get(tupleIndex);
			int PDBIndex = residueIndexMap.designIndexToPDBIndex(pos);
			int RC = queryAssignment.RCs.get(tupleIndex);
			if(MSet.contains(PDBIndex)) {
				
				int subspaceRCIndex = getRCSpaceRCIndex(PDBIndex, RC);
				RCSpaceRCIndices[pos] = subspaceRCIndex;
				RCSpaceSizes[pos] = localConfSpace.getNum(pos);
			}
		}

		
		int newSubproblemIndex = 0;
		int newMultiplier = 1;
		for(int i = 0; i < RCSpaceRCIndices.length; i++) {
			if(RCSpaceRCIndices[i] > -1)
			{
				newSubproblemIndex += RCSpaceRCIndices[i]*newMultiplier;
				newMultiplier *= RCSpaceSizes[i];
			}
		}
		return newSubproblemIndex;
	}

	public int mapSubproblemConfToIndex(int[] localConf)
	{
		int subproblemIndex = 0;
		int multiplier = 1;
		for(int i = 0; i < localConf.length; i++)
		{
			int PDBIndex = residueIndexMap.designIndexToPDBIndex(i);
			int subspaceRCIndex = getRCSpaceRCIndex(PDBIndex, localConf[i]);
			subproblemIndex += subspaceRCIndex*multiplier;
			multiplier *= localConfSpace.getNum(i);
		}
		
		return subproblemIndex;
	}
	
	public int mapSubproblemSeqToIndex(Sequence seq)
	{
		return mapSubproblemConfToIndex(residueIndexMap.createTemplateConfForSequence(seq));
	}
	
	private int getRCSpaceRCIndex(int PDBIndex, int PDBSpaceRCIndex)
	{
		return PDBRCToSubspaceRCMap.get(PDBIndex).get(PDBSpaceRCIndex);
	}
	
	private void generateSubproblems (TreeEdge sparseTree) {


		TreeEdge curEdge = sparseTree;
		//ordered array for recursion
		if(curEdge.leftChild != null)
		{
			leftSubproblem = new Subproblem(localConfSpace, sparseTree.leftChild.getCofEdge(), residueIndexMap);
		}
		if(curEdge.rightChild != null)
		{	
			rightSubproblem = new Subproblem(localConfSpace, sparseTree.rightChild.getCofEdge(), residueIndexMap);
		}
	}

	public Sequence getSequenceForConf(RCTuple conf)
	{
		return residueIndexMap.getSequenceOfRCTuple(conf);
	}
	
	public void addConformationProcessor(ConformationProcessor processor)
	{
		if(processor.recurse())
		{
			if(leftSubproblem != null)
				leftSubproblem.addConformationProcessor(processor);
			if(rightSubproblem != null)
				rightSubproblem.addConformationProcessor(processor);
		}
		processors.add(processor);
	}

	public BigInteger getTotalConformations()
	{
		return localConfSpace.unprunedConfsFromRCs();
	}
	
	public BigInteger getTotalLocalConformations()
	{
		BigInteger numConformations = BigInteger.ONE;
		for(int PDBIndex : MULambdaSet)
		{
			int designIndex = residueIndexMap.PDBIndexToDesignIndex(PDBIndex);
			numConformations = numConformations.multiply(BigInteger.valueOf(localConfSpace.get(designIndex).length));
		}
		return numConformations;
	}
	
	public BigInteger getTotalMConformations()
	{
		BigInteger numConformations = BigInteger.ONE;
		for(int PDBIndex : MSet)
		{
			int designIndex = residueIndexMap.PDBIndexToDesignIndex(PDBIndex);
			numConformations = numConformations.multiply(BigInteger.valueOf(localConfSpace.get(designIndex).length));
		}
		return numConformations;
	}
	
	public BigInteger getSubtreeTESS()
	{
		BigInteger numConformations = getTotalLocalConformations();

		if(leftSubproblem != null)
			numConformations = numConformations.add(leftSubproblem.getSubtreeTESS());
		if(rightSubproblem != null)
			numConformations = numConformations.add(rightSubproblem.getSubtreeTESS());
		return numConformations;
	}

	public void preprocess () {
		if(leftSubproblem != null)
			leftSubproblem.preprocess();
		if(rightSubproblem != null)
			rightSubproblem.preprocess();
		int[] currentConf = new int[localConfSpace.getNumPos()];
		recursivelyProcessTuples(0,currentConf);
	}



	protected void recursivelyProcessTuples (int position, int[] currentConf) {
		if(position >= localConfSpace.getNumPos())
		{
			RCTuple confTuple = new RCTuple(currentConf);
			for(ConformationProcessor proc : processors)
			{
				proc.processConformation(confTuple);
			}
			return;
		}
		
		int PDBIndex = residueIndexMap.designIndexToPDBIndex(position);
		if(!MULambdaSet.contains(PDBIndex))
		{
			currentConf[position] = -1;
			recursivelyProcessTuples(position+1, currentConf);
			return;
		}
		
		for(int i = 0; i < localConfSpace.getNum(position); i++)
		{
			currentConf[position] = localConfSpace.get(position, i);
			recursivelyProcessTuples(position+1, currentConf);
		}
	}

	private String printSubspaceConf (int[] currentConf) {
		String output = "(";
		for(int i = 0; i < currentConf.length-1; i++)
		{
			output+=i+":"+currentConf[i]+", ";
		}
		output = output+(currentConf.length-1)+":"+currentConf[currentConf.length-1]+")";
		return output;
	}

	private String printConf (int[] currentConf) {
		String output = "(";
		for(int i = 0; i < currentConf.length-1; i++)
		{
			int PDBIndex = residueIndexMap.designIndexToPDBIndex(i);
			if(MULambdaSet.contains(PDBIndex))
				output+=PDBIndex+":"+currentConf[i]+", ";
		}
		int finalPDBIndex = residueIndexMap.designIndexToPDBIndex(currentConf.length-1);
		if(MULambdaSet.contains(finalPDBIndex))
			output += (finalPDBIndex)+":"+currentConf[currentConf.length-1];
		output += ")";
		return output;
	}

	public boolean isMULambdaAssignment (RCTuple queryAssignment) {
		return setEquals(queryAssignment.pos, MULambdaSet);
	}
	
	public boolean isMAssignment (RCTuple queryAssignment) {
		return setEquals(queryAssignment.pos, MSet);
	}
	
	private boolean setEquals(Collection<Integer> a, Collection<Integer> b)
	{
		if(a.size() != b.size())
			return false;
		for(int index : a)
		{
			if(!b.contains(residueIndexMap.designIndexToPDBIndex(index)))
				return false;
		}
		return true;
	}
	
	public boolean isMULambdaAssignment (Sequence queryAssignment) {
		return setEquals(queryAssignment.residues(), MULambdaSet);
	}
	
	public boolean isMAssignment (Sequence queryAssignment) {
		return setEquals(queryAssignment.residues(), MSet);
	}

	public RCTuple extractSubproblemAssignment (RCTuple queryAssignment) {
		return filterTuple(MULambdaSet, queryAssignment);
	}
	
	public RCTuple extractSubproblemLambdaAssignment (RCTuple queryAssignment) {
		return filterTuple(lambdaSet, queryAssignment);
	}


	public RCTuple extractSubproblemMAssignment (RCTuple queryAssignment) {
		return filterTuple(MSet, queryAssignment);
	}
	
	
	public RCTuple extractSubproblemLAssignment (RCTuple queryAssignment) {
		return filterTuple(LSet, queryAssignment);
	}
	
	private RCTuple filterTuple (Set<Integer> residueSet, RCTuple queryAssignment) {
		RCTuple filteredTuple = new RCTuple();

		for(int i = 0; i < queryAssignment.size(); i++)
		{
			int residueIndex = queryAssignment.pos.get(i);
			if(residueSet.contains(residueIndexMap.designIndexToPDBIndex(residueIndex)))
				filteredTuple = filteredTuple.addRC(residueIndex, queryAssignment.RCs.get(i));
		}		
		if(residueSet.size() != filteredTuple.size())
		{
			System.err.println("Tuple filter failed. Size is wrong.");
			System.out.println("lambdaSet: "+residueSet);
			filteredTuple = new RCTuple();

			for(int i = 0; i < queryAssignment.size(); i++)
			{
				int residueIndex = queryAssignment.pos.get(i);
				int PDBIndex = residueIndexMap.designIndexToPDBIndex(residueIndex);
				if(residueSet.contains(PDBIndex))
				{
					System.out.println("Adding RC at design index "+residueIndex+", whose PDB Index is "+PDBIndex);
					filteredTuple = filteredTuple.addRC(residueIndex, queryAssignment.RCs.get(i));
				}
			}	
		}
		assert(residueSet.size() == filteredTuple.size());

		return filteredTuple;
	}
	
	
	
	public Sequence extractSubproblemAssignment (Sequence queryAssignment) {
		return filterTuple(MULambdaSet, queryAssignment);
	}
	
	public Sequence extractSubproblemLambdaAssignment (Sequence queryAssignment) {
		return filterTuple(lambdaSet, queryAssignment);
	}
	
	public Sequence extractSubproblemMAssignment (Sequence queryAssignment) {
		return filterTuple(MSet, queryAssignment);
	}
	
	public Sequence extractSubproblemLAssignment (Sequence queryAssignment) {
		return filterTuple(LSet, queryAssignment);
	}
	
	private Sequence filterTuple(Set<Integer> positions, Sequence queryAssignment)
	{
		Sequence filteredTuple = new Sequence();

		for(Integer residueIndex : queryAssignment.residues())
		{
			if(positions.contains(residueIndexMap.designIndexToPDBIndex(residueIndex)))
				filteredTuple.addAssignment(residueIndex, queryAssignment.getAAAt(residueIndex));
		}		
		if(positions.size() != filteredTuple.size())
		{
			System.err.println("Tuple filter failed. Size is wrong.");
			System.out.println("input set: "+positions);
			filteredTuple = new Sequence();

			for(Integer residueIndex : queryAssignment.residues())
			{
				int PDBIndex = residueIndexMap.designIndexToPDBIndex(residueIndex);
				if(positions.contains(PDBIndex))
				{
					System.out.println("Adding RC at design index "+residueIndex+", whose PDB Index is "+PDBIndex);
					filteredTuple.addAssignment(residueIndex, queryAssignment.getAAAt(residueIndex));
				}
			}	
		}
		assert(positions.size() == filteredTuple.size());

		return filteredTuple;
	}
	
	public String printTreeMol()
	{
		return printTreeMol("");
	}
	
	public String printTreeDesign()
	{
		return printTreeDesign("");
	}
	
    public String printTreeDesign(String prefix)
    {
        String out = prefix + "[";
        for(int i : lambdaSet)
            out+=residueIndexMap.PDBIndexToDesignIndex(i)+", ";
        out+="]";

        out += " - L Set:[";

        for(int i : LSet)
            out+=residueIndexMap.PDBIndexToDesignIndex(i)+", ";
        out+="]";

        out += " - M Set:[";

        for(int i : MSet)
            out+=residueIndexMap.PDBIndexToDesignIndex(i)+", ";
        out+="]";

//        boolean showHeaps = false;
//        if(showHeaps)
//        {
//            out+="heaps: \n";
//            for(PriorityQueue<Conf> heap : A2)
//            {
//                if(heap.size() > 0)
//                    out+=""+heap+"\n";
//            }
//        }
        if(leftSubproblem != null)
            out += "\n"+leftSubproblem.printTreeDesign(prefix+"+L--");
        if(rightSubproblem != null)
            out += "\n"+rightSubproblem.printTreeDesign(prefix+"+R--");
        return out;
    }
	
    public String printTreeMol(String prefix)
    {
        String out = prefix + "[";
        for(int i : lambdaSet)
            out+=i+", ";
        out+="]";

        out += " - L Set:[";

        for(int i : LSet)
            out+=i+", ";
        out+="]";

        out += " - M Set:[";

        for(int i : MSet)
            out+=i+", ";
        out+="]";
//        boolean showHeaps = false;
//        if(showHeaps)
//        {
//            out+="heaps: \n";
//            for(PriorityQueue<Conf> heap : A2)
//            {
//                if(heap.size() > 0)
//                    out+=""+heap+"\n";
//            }
//        }
        if(leftSubproblem != null)
            out += leftSubproblem.printTreeMol(prefix+"\n+L--");
        if(rightSubproblem != null)
            out += rightSubproblem.printTreeMol(prefix+"\n+R--");
        return out;
    }
    
    public String toString()
    {
    	String out = "[";
        for(int i : lambdaSet)
            out+=residueIndexMap.PDBIndexToDesignIndex(i)+", ";
        out+="]";

        out += " - L Set:[";

        for(int i : LSet)
            out+=residueIndexMap.PDBIndexToDesignIndex(i)+", ";
        out+="]";

        out += " - M Set:[";

        for(int i : MSet)
            out+=residueIndexMap.PDBIndexToDesignIndex(i)+", ";
        out+="]";
        return out;
    }

	public BigInteger getLambdaConformations () {
		BigInteger numConformations = BigInteger.ONE;
		for(int PDBIndex : lambdaSet)
		{
			int designIndex = residueIndexMap.PDBIndexToDesignIndex(PDBIndex);
			numConformations = numConformations.multiply(BigInteger.valueOf(localConfSpace.get(designIndex).length));
		}
		return numConformations;
	}
	
	public int getNumConfsForSequence (Sequence seq) {
		int numConfs = 1;
		for(int residueIndex : seq.residues())
		{
			numConfs *= residueIndexMap.getConfsForAAAtPos(residueIndex, seq.getAAAt(residueIndex));
		}
		// TODO Auto-generated method stub
		return numConfs;
	}

	public boolean isRoot () {
		return MSet.size() < 1;
	}

	public boolean isInternalNode () {
		return lambdaSet.size() < 1;
	}

	public boolean isValidConf (RCTuple assignment) {
		RCTuple lambdaAssignment =  extractSubproblemLambdaAssignment(assignment);
		return lambdaAssignment.size() == lambdaSet.size()
				|| (isInternalNode() && leftSubproblem.isValidConf(assignment));
	}
	
	public boolean isValidSequence (Sequence assignment) {
		Sequence lambdaAssignment =  extractSubproblemLambdaAssignment(assignment);
		return lambdaAssignment.size() == lambdaSet.size()
				|| (isInternalNode() && leftSubproblem.isValidSequence(assignment));
	}





}
