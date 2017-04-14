/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.energy.forcefield;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.Serializable;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.Map;
import java.util.Set;

import edu.duke.cs.osprey.restypes.DAminoAcidHandler;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.HashCalculator;
import edu.duke.cs.osprey.tools.StringParsing;

/**
 * Manages the EEF1 solvation parameters;
 * 		EEF1 handles only natural amino acids, so solvation energies are computed for proteins and ligands (natural amino acids) only, and not cofactors;
 * 		Some important notes about assumptions of the EEF1 model are given in the comments to getSolvationParameters()
 * 
*/






public class EEF1 implements Serializable {
	
	private static final long serialVersionUID = -4783417295676415124L;
	
	public static final double trigConst = 2.0/(4.0*Math.PI*Math.sqrt(Math.PI));
	
	//Variables to store the EEF1 solvation paramaters;
	//These values are specific to eef1parm.dat, and may have to be modified for different
	//		versions of the parameter file
	final int numAtTypesEEF1 = 17;
	String groupEEF1names[] = new String[numAtTypesEEF1];
	String atTypeEEF1names[] = new String[numAtTypesEEF1];
	double atEEF1Vol[] = new double[numAtTypesEEF1];
	double dGiRef[] = new double[numAtTypesEEF1];
	double dGiFree[] = new double[numAtTypesEEF1];
	double dHiRef[] = new double[numAtTypesEEF1];
	double dCpiRef[] = new double[numAtTypesEEF1];
	double lambdai[] = new double[numAtTypesEEF1];
	double vdWri[] = new double[numAtTypesEEF1];	

	
	//constructor
	public EEF1(){
            
	}
	
	//Reads the solvation parameters for EEF1 from the file eef1parm.dat;
	//If a different versio of the paramater file is used, some changes to the
	//		current function may be necessary
	public void readEEF1parm() throws Exception {
		
		try (BufferedReader bufread = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("/config/eef1parm.dat")))) {
			
			String curLine = null;
			int tmpInt = 0;

			//Skip over the first line of header info
			curLine = bufread.readLine();
			
			//Read each parameter line, which is in the format:
			//	Group type   Atom type   Volume   DGi(ref)   DGi(free)   DHi(ref)   DCpi(ref)   lambda_i   vdWri
			curLine = bufread.readLine();		
			tmpInt = 0; // temporary integer
			// Until we're at a blank line (or until we've read numAtomTypes)
			while (!(curLine==null)) {
				groupEEF1names[tmpInt] = StringParsing.getToken(curLine,1);  // snag group name
				atTypeEEF1names[tmpInt] = StringParsing.getToken(curLine,2);  // snag atom type
				atEEF1Vol[tmpInt] = (new Double(StringParsing.getToken(curLine,3))).doubleValue();
				dGiRef[tmpInt] = (new Double(StringParsing.getToken(curLine,4))).doubleValue();
				dGiFree[tmpInt] = (new Double(StringParsing.getToken(curLine,5))).doubleValue();
				dHiRef[tmpInt] = (new Double(StringParsing.getToken(curLine,6))).doubleValue();
				dCpiRef[tmpInt] = (new Double(StringParsing.getToken(curLine,7))).doubleValue();
				lambdai[tmpInt] = (new Double(StringParsing.getToken(curLine,8))).doubleValue();
				vdWri[tmpInt] = (new Double(StringParsing.getToken(curLine,9))).doubleValue();
				tmpInt++;
				curLine = bufread.readLine();
			}
		}
	}
	
	//Returns the dG(ref), dG(free), atom volume, and lambda parameters for 
	//		the atom specified by moleculeAtomNumber atom1;
	//The mapping is based on the group types in eef1parm.dat: different versions
	//		of the parameter file may require modifications in this function;
	//The assumption is made that the temperature in Kelvin is 298.15, since this is
	//		the temperature for which dGi(free) and dGi(ref) are derived in eef1parm.dat
	//Note that some of the C and S groups are *extended* heavy atoms (heavy atoms which include
	//		the bonded hydrogen atoms), so in the original EEF1 and CHARMM implementation, the
	//		corresponding vdW radii are increased accordingly. We thus use the vdW radii from the
	//		extended Param 19 (see Appendix B, Neria & Karplus, 1996) parameters for the solvation energy computation only;
	//The groups in eef1 are based on the natural amino acids. Additional groups may have to be included, 
	//		and the current groups may have to be modified if new non-natural amino acids or other
	//		atom types or molecules (such as AMP or ATP) are included
	
	public static class SolvParams {
		public double dGref;
		public double dGfree;
		public double volume;
		public double lambda;
		public double radius;
	}
	
	public boolean getSolvationParameters(Atom at1, SolvParams solvparams) {
				
		int solvGroupIndex = getSolvGroupIndex(at1);
		
		//System.out.print(m.residue[at1.moleculeResidueNumber].fullName+" "+at1.name);
		
		if (solvGroupIndex==-1) //group not found
			return false;
		else {
			solvparams.dGref = dGiRef[solvGroupIndex];
			solvparams.dGfree = dGiFree[solvGroupIndex];
			solvparams.volume = atEEF1Vol[solvGroupIndex];
			solvparams.lambda = lambdai[solvGroupIndex];
			solvparams.radius = vdWri[solvGroupIndex];
			
			//System.out.println("           "+groupEEF1names[solvGroupIndex]);
			
			return true;
		}		
	}
	
	private static Set<String> warnedAtomTypes = new HashSet<>();
	private static void warnAtomType(String type) {
		if (warnedAtomTypes.add(type)) {
			System.err.println("WARNING: couldn't find solvation parameters for atom type: " + type + ", using default values");
		}
	}
	
	public void getSolvationParametersOrDefaults(Atom atom, SolvParams solvparams) {
		boolean success = getSolvationParameters(atom, solvparams);
		if (!success) {
			
			// if there's no params, don't crash, use defaults instead
			warnAtomType(atom.forceFieldType);
			
			solvparams.dGref = 0;
			solvparams.dGfree = 0;
			solvparams.volume = 0;
			solvparams.lambda = 1;
			solvparams.radius = 0;
		}
	}
	
	public static class SolvPairParams {
		public double radius1;
		public double lambda1;
		public double alpha1;
		public double radius2;
		public double lambda2;
		public double alpha2;
	}
	
	public void getSolvationPairParams(Atom atom1, Atom atom2, double scale, SolvPairParams out) {
		
		// start with defaults
		double dGfree1 = 0;
		double volume1 = 0;
		out.lambda1 = 1;
		out.radius1 = 0;
		
		double dGfree2 = 0;
		double volume2 = 0;
		out.lambda2 = 1;
		out.radius2 = 0;
		
		// lookup atom 1
		int i1 = getSolvGroupIndex(atom1);
		if (i1 == -1) {
			warnAtomType(atom1.forceFieldType);
		} else {
			dGfree1 = dGiFree[i1];
			volume1 = atEEF1Vol[i1];
			out.lambda1 = lambdai[i1];
			out.radius1 = vdWri[i1];
		}
		
		// lookup atom 2
		int i2 = getSolvGroupIndex(atom2);
		if (i2 == -1) {
			warnAtomType(atom2.forceFieldType);
		} else {
			dGfree2 = dGiFree[i2];
			volume2 = atEEF1Vol[i2];
			out.lambda2 = lambdai[i2];
			out.radius2 = vdWri[i2];
		}
		
		// calc mixed params
		scale *= trigConst;
		out.alpha1 = scale*dGfree1*volume2/out.lambda1;
		out.alpha2 = scale*dGfree2*volume1/out.lambda2;
	}
	
	private class SolvGroupKey {
		
		public final String atomName;
		public final String elementType;
		public final String AAname;
		public final int numBoundH;
		public final boolean isCarboxylO;
		
		private final int hashCode;
		
		public SolvGroupKey(Atom atom) {
			
			atomName = atom.name;
			elementType = atom.elementType;
			AAname = atom.res.template.name;
			numBoundH = getNumBoundH(atom);
			isCarboxylO = isCarboxylO(atom);
			
			hashCode = HashCalculator.combineHashes(
				atomName.hashCode(),
				elementType.hashCode(),
				AAname.hashCode(),
				numBoundH,
				isCarboxylO ? 2 : 1
			);
		}
		
		@Override
		public int hashCode() {
			return hashCode;
		}
		
		@Override
		public boolean equals(Object obj) {
			SolvGroupKey other = (SolvGroupKey)obj;
			return this.isCarboxylO == other.isCarboxylO
				&& this.numBoundH == other.numBoundH
				&& this.AAname.equals(other.AAname)
				&& this.elementType.equals(other.elementType)
				&& this.atomName.equals(other.atomName);
		}
	}
	
	private static HashMap<SolvGroupKey,Integer> solvGroupIndices = new HashMap<>();
	
	//Determines the solvation group index for the given atom;
	//This function assumes that the predfefined groups are as given in eef1parm.dat;
	//		Modifcations to this functions will be required if different group types are used
	private int getSolvGroupIndex(Atom atom) {
		SolvGroupKey key = new SolvGroupKey(atom);
		Integer index = solvGroupIndices.get(key);
		if (index == null) {
			index = getSolvGroupIndex(key);
			solvGroupIndices.put(key, index);
		}
		return index;
		//return getSolvGroupIndex(atom.name, atom.elementType, atom.res.template.name, getNumBoundH(atom), isCarboxylO(atom));
	}
	
	private int getSolvGroupIndex(SolvGroupKey key) {
		return getSolvGroupIndex(key.atomName, key.elementType, key.AAname, key.numBoundH, key.isCarboxylO);
	}
	
	private int getSolvGroupIndex(String atomName, String elementType, String AAname, int numBoundH, boolean isCarboxylO) {
		
		if(DAminoAcidHandler.isDAminoAcidName(AAname))
			AAname = DAminoAcidHandler.getLEquivalent(AAname);//look up parameters by L-amino acid name (has same energetics)
		
		boolean aromatic = isAromatic(atomName, elementType, AAname);
		
		if (elementType.equalsIgnoreCase("C")) {
			
			if ( atomName.equalsIgnoreCase("CG") && 
					( AAname.equalsIgnoreCase("TYR") || AAname.equalsIgnoreCase("PHE") || AAname.equalsIgnoreCase("HIS") 
							|| AAname.equalsIgnoreCase("HIP") || AAname.equalsIgnoreCase("HID") || AAname.equalsIgnoreCase("HIE")) )
				return findSolvGroup("CR"); //CG of (TYR or PHE or HIS)
			
			else if ( ( atomName.equalsIgnoreCase("CG") || atomName.equalsIgnoreCase("CD2") || atomName.equalsIgnoreCase("CE2") )
					&& AAname.equalsIgnoreCase("TRP") )
				return findSolvGroup("CR"); // (CG or CD2 or CE2) of TRP
			
			else if ( atomName.equalsIgnoreCase("CZ") && AAname.equalsIgnoreCase("ARG") )
				return findSolvGroup("CR"); //CZ of PHE
			
			else if ( (aromatic) && (numBoundH==1) ) // extended aromatic C with 1 H, so CR1E group
				return findSolvGroup("CR1E");
			
			else if ( (!aromatic) && (numBoundH==1) ) // extended aliphatic C with 1 H, so CH1E group
				return findSolvGroup("CH1E");
			
			else if ( (!aromatic) && (numBoundH==2) ) // extended aliphatic C with 2 H, so CH2E group
				return findSolvGroup("CH2E");
			
			else if ( (!aromatic) && (numBoundH==3) ) // extended aliphatic C with 3 H, so CH3E group
				return findSolvGroup("CH3E");
			
			else // default C group
				return findSolvGroup("C");
		}
		else if (elementType.equalsIgnoreCase("N")) {
			
			if (AAname.equalsIgnoreCase("PRO")) // backbone N of PRO, so N group
				return findSolvGroup("N");
			
			else if (numBoundH==1) // any N with 1 H
				return findSolvGroup("NH1");
			
			else if (AAname.equalsIgnoreCase("ARG") && (numBoundH==2)) // guanidinium N in (side-chain of) ARG, so NC2 group
				return findSolvGroup("NC2");
			
			else if ( (aromatic) && (numBoundH==0) ) // aromatic (side-chain) N with 0 H, so NR group
				return findSolvGroup("NR");
			
			else if (numBoundH==2) // any (side-chain) N with 2 H, so NH2 group
				return findSolvGroup("NH2");
			
			else if (numBoundH==3) // any N with 3 H, so NH3 group
				return findSolvGroup("NH3");
			
			else //unknown nitrogen group
				return -1;
		}
		else if (elementType.equalsIgnoreCase("O")) {
			
			if (numBoundH==1) // hydroxyl oxygen, so OH1 group
				return findSolvGroup("OH1");
			
			else if (isCarboxylO) // carboxyl oxygen, so OC group
				return findSolvGroup("OC");
			
			else //all other oxygens (carbonyl oxygen)
				return findSolvGroup("O");
		}
		else if (elementType.equalsIgnoreCase("S")) {
			
			if (numBoundH==1) // extended S with 1 H, so SH1E group
				return findSolvGroup("SH1E");
			
			else // any other S, so S group
				return findSolvGroup("S");
		}
		else //group not found
			return -1;
	}
	
	//Find the solvation group index for the group name grName
	private int findSolvGroup(String grName){
		for (int i=0; i<groupEEF1names.length; i++){
			if (groupEEF1names[i].equalsIgnoreCase(grName))
				return i;
		}
		return -1; //group name not found
	}
	
	//Determines if ths given atom is a carboxyl oxygen (there should be only one carbon bound to the
	//		given oxygen and exactly two oxygens (including the given one) should be bound to that carbon)
	private boolean isCarboxylO(Atom at1){
		
		if (!at1.elementType.equalsIgnoreCase("O")) //not O
			return false;
		else {
			
			int numBoundC = 0;
			boolean isCarboxylO = false;
			if (at1.bonds!=null){ //check bonds
				for (int i=0; i<at1.bonds.size(); i++){		
					
					Atom at2 = at1.bonds.get(i);
					
					if (at2.elementType.equalsIgnoreCase("C")){ //found a bound C
						
						numBoundC++;
						if (numBoundC>1) //more than 1 C bound
							return false;
						else
							isCarboxylO = isCOhelper(at2);
					}
				}
			}
			return isCarboxylO;			
		}
	}
	
	//Used by isCarboxylO(); see comments above
	private boolean isCOhelper(Atom at2){		
		
		int numBoundO = 0;
		if (at2.bonds!=null){ //check bonds
			for (int j=0; j<at2.bonds.size(); j++){
				if (at2.bonds.get(j).elementType.equalsIgnoreCase("O"))
					numBoundO++;
			}
		}
		
		if (numBoundO==2)
			return true;
		else
			return false;
	}
	
	//Determines if the given heavy atom is aromatic
	private boolean isAromatic(String atomName, String elementType, String AAname) {
		
		if (DAminoAcidHandler.isDAminoAcidName(AAname)) {
			AAname = DAminoAcidHandler.getLEquivalent(AAname);//look up parameters by L-amino acid name (has same energetics)
		}
		
		boolean isHis = (AAname.equalsIgnoreCase("HIS") || AAname.equalsIgnoreCase("HIP") || AAname.equalsIgnoreCase("HID") || AAname.equalsIgnoreCase("HIE"));
		
		if ( !(AAname.equalsIgnoreCase("PHE") || AAname.equalsIgnoreCase("TYR") 
					|| AAname.equalsIgnoreCase("TRP") || isHis) ) //not PHE and not TYR and not TRP or HIS, so not aromatic
				return false;
		else {//PHE or TYR or TRP or HIS, so check if aromatic atom
			
			if (!elementType.equalsIgnoreCase("C")){ //not a carbon
				if ( elementType.equalsIgnoreCase("N") && 
						( !atomName.equalsIgnoreCase("N") && (AAname.equalsIgnoreCase("TRP") || isHis) ) ) //N in the ring of TRP or HIS
					return true;
				else
					return false;
			}
			else {//a carbon, so check if CA or CB or C
				if (atomName.equalsIgnoreCase("CA") || atomName.equalsIgnoreCase("CB") 
						|| atomName.equalsIgnoreCase("C")) //CA or CB or C, so cannot be aromatic
					return false;
				else //aromatic
					return true;
			}
		}
	}
	
	//Get the number of H bound to the given atom
	private int getNumBoundH(Atom at1){
		
		int numBoundH = 0;
		for (int i=0; i<at1.bonds.size(); i++){
			if (at1.bonds.get(i).elementType.equalsIgnoreCase("H"))
				numBoundH++;
		}
		
		return numBoundH;
	}
	
	private Map<ResidueTemplate,Double> internalSolvEnergies = new IdentityHashMap<>();
	
	public double getInternalEnergy(Residue res) {
		
		// check the cache first
		Double energy = internalSolvEnergies.get(res.template);
		if (energy != null) {
			return energy;
		}
		
		// add up all the dGref terms for all the atoms
		energy = 0.0;
		for (Atom atom : res.atoms) {
			if (!atom.isHydrogen()) {
				int index = getSolvGroupIndex(atom);
				if (index == -1) {
					warnAtomType(atom.forceFieldType);
				} else {
					energy += dGiRef[index];
				}
			}
		}
			
		internalSolvEnergies.put(res.template, energy);
		
		return energy;
	}
}
