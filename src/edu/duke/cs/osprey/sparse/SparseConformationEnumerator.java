package edu.duke.cs.osprey.sparse;

import java.util.Iterator;
import edu.duke.cs.osprey.confspace.RCTuple;

/***
 * This class enumerates conformations in order. 
 * @author Jon
 *
 */

public class SparseConformationEnumerator implements ConformationProcessor, Iterator<RCTuple>{
	
	private SubproblemConfEnumerator root;
	
	public SparseConformationEnumerator(BranchDecomposedProblem searchProblem, PartialConformationEnergyFunction eFunc)
	{
		root = new SubproblemConfEnumerator(searchProblem.getRoot(), eFunc);
	}


	@Override
	public boolean hasNext () {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public RCTuple next () {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public void processConformation (RCTuple conformation) {
		// TODO Auto-generated method stub
		
	}


	@Override
	public boolean recurse () {
		// TODO Auto-generated method stub
		return false;
	}

}
