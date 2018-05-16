package edu.duke.cs.osprey.sparse;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.Minimizer.Result;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;

public class PartialConformationEnergyFunction {
	
	private EnergyFunction fullEnergyFunction;
	private ConfSpace conformations;
	private EnergyMatrix eMat;
	
	public PartialConformationEnergyFunction(SearchProblem problem, EnergyFunction termE, ConfSpace conformationSpace)
	{
		eMat = problem.emat;
		fullEnergyFunction = termE;
		conformations = conformationSpace;
	}
	
	public double computePartialEnergy(RCTuple partialAssignment)
	{
		return computePartialEnergyGivenPriorConformation(null, partialAssignment);
	}
	
	public double computePartialEnergyGivenPriorConformation(RCTuple priorConformation, RCTuple partialAssignment)
	{
		RCTuple combinedAssignment = new RCTuple();
		double MEnergy = 0;
		if(priorConformation!=null)
		{
			MEnergy = eMat.getInternalEnergy(priorConformation);
			combinedAssignment = combinedAssignment.combineRC(priorConformation);
		}
		combinedAssignment = combinedAssignment.combineRC(partialAssignment);
		
		return eMat.getInternalEnergy(combinedAssignment) - MEnergy;
	}
	

	
	public double computePartialEnergyOld(RCTuple priorConformation, RCTuple partialAssignment)
	{
		RCTuple combinedAssignment = new RCTuple();
		if(priorConformation!=null)
			combinedAssignment = combinedAssignment.combineRC(priorConformation);
		combinedAssignment = combinedAssignment.combineRC(partialAssignment);
		System.out.println("Energy:"+fullEnergyFunction.getEnergy());

		MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(fullEnergyFunction,conformations,combinedAssignment);
		System.out.println("Energy:"+fullEnergyFunction.getEnergy());

        DoubleMatrix1D bestDOFVals;

        if(mof.getNumDOFs()>0){//there are continuously flexible DOFs to minimize
            CCDMinimizer ccdMin = new CCDMinimizer(mof,true);
            bestDOFVals = ccdMin.minimize().dofValues;
        }
        else//molecule is already in the right, rigid conformation
            bestDOFVals = DoubleFactory1D.dense.make(0);

        System.out.println("Energy:"+fullEnergyFunction.getEnergy());

        return mof.getEnergyAndReset(bestDOFVals);
	}

}
