/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import edu.duke.cs.osprey.sparse.Subproblem;
import edu.duke.cs.osprey.sparse.sequence.Sequence;

/**
 *
 * @author mhall44
 */
public class RCTuple implements Serializable {
    
	private static final long serialVersionUID = 3773470316855163174L;
	
	//a tuple of RCs
    public ArrayList<Integer> pos;//which flexible positions
    public ArrayList<Integer> RCs;//the RCs themselves (residue-specific numbering, as in the TupleMatrices)
    private boolean immutable = false;
    
    // A non-painfully annoying parallel array implemtation of an ordered set of ordered pairs.
    private ArrayList<RCPair> RCPairs = new ArrayList<>();;

    
    public RCTuple(){
        //empty pos, RCs (basically tuple of nothing
        pos = new ArrayList<>();
        RCs = new ArrayList<>();
    }
    
    public RCTuple(boolean immutable){
        //empty pos, RCs (basically tuple of nothing
    	this.immutable = immutable;
        pos = new ArrayList<>();
        RCs = new ArrayList<>();
    }
    
    
    public RCTuple(ArrayList<Integer> pos, ArrayList<Integer> RCs) {
        this.pos = pos;
        this.RCs = RCs;
    }
    
    
    //one-RC tuple
    public RCTuple(int pos1, int RC1){
    	this();
    	set(pos1, RC1);
    }
    
    
    //a pair
    public RCTuple(int pos1, int RC1, int pos2, int RC2){
    	this();
    	set(pos1, RC1, pos2, RC2);
    }
    
    //Sometimes we'll want to generate an RC tuple from a conformation, specified as RCs for all positions
    //in order.  
    //In this case, negative values are not (fully) defined, so the tuple contains all positions
    //with positive values in conf
    public RCTuple(int[] conf){
    	this();
    	set(conf);
    }
    
    public RCTuple copy() {
    	RCTuple copy = new RCTuple();
    	copy.set(this);
    	return copy;
	}
    
    /**
     * Combine two RCTuples. Should throw an error and crash if an RC is overwritten with 
     * different RC.
     * @param source
     */
    public RCTuple combineRC(RCTuple source) {
    	RCTuple output = copy();
    	for(int RCIndex = 0; RCIndex < source.size(); RCIndex++) {
    		int sourcePos = source.pos.get(RCIndex);
    		int sourceRC = source.RCs.get(RCIndex);
    		assert(sourceRC != -1); 
    		
    		int correspondingRC = lookupRC(sourcePos);
    		if(correspondingRC != -1 && correspondingRC != sourceRC)
    		{
    			System.err.println("Error combining RCs: would overwrite "
    					+this+" with "+source+", which differs at position "+sourcePos);
    		}
    		assert(correspondingRC == -1 || correspondingRC == sourceRC);
    		if(correspondingRC == -1)
    		{
    			output.pos.add(sourcePos);
    			output.RCs.add(sourceRC);
    		}
    	}

    	return output;
    }
    
    public boolean consistentWith(RCTuple otherTuple) {
    	for(int RCIndex = 0; RCIndex < otherTuple.size(); RCIndex++) {
    		int sourcePos = otherTuple.pos.get(RCIndex);
    		int sourceRC = otherTuple.RCs.get(RCIndex);
    		assert(sourceRC != -1); 
    		
    		int correspondingRC = lookupRC(sourcePos);
    		if(correspondingRC != -1 && correspondingRC != sourceRC)
    		{
    			System.out.println("Inconsistent RCs: "
    					+this+" != "+otherTuple+". They differ at position "+sourcePos);
    			return false;
    		}
    	}

    	return true;
    }
    
    public boolean consistentWith(Sequence otherSequence, Subproblem sourceProblem) {
    	Sequence selfSequence = sourceProblem.getSequenceForConf(this);
    	return selfSequence.consistentWith(otherSequence);
    }
    
    private int lookupRC(int queryPos)
    {
    	for (int i = 0; i < pos.size(); i++){
    		if(pos.get(i) == queryPos)
    			return RCs.get(i);
    	}
    	return -1;
    }


	public void set(int pos, int rc) {
    	this.pos.clear();
    	this.RCs.clear();
    	
        this.pos.add(pos);
        this.RCs.add(rc);
    }
    
    public void set(int pos1, int rc1, int pos2, int rc2) {
    	this.pos.clear();
    	this.RCs.clear();
    	
        this.pos.add(pos1);
        this.RCs.add(rc1);
        
        this.pos.add(pos2);
        this.RCs.add(rc2);
    }
    
    public void set(int[] conf) {
    	pos.clear();
    	RCs.clear();
        for(int posNum=0; posNum<conf.length; posNum++){
            if(conf[posNum]>=0){//RC fully defined
                pos.add(posNum);
                RCs.add(conf[posNum]);
            }
        }
    }
    
    public void set(RCTuple other) {
    	pos.clear();
    	RCs.clear();
    	pos.addAll(other.pos);
    	RCs.addAll(other.RCs);
    }
    
    public int size() {
    	return pos.size();
    }
    
    
    public boolean isSameTuple(RCTuple tuple2){
    	
    	// short circuit: same instance must have same value
    	if (this == tuple2) {
    		return true;
    	}
    	
        //do the two tuple objects specify the same tuple of RCs?
        if( (pos.size()!=RCs.size()) || (tuple2.pos.size()!=tuple2.RCs.size()) )
            throw new RuntimeException("ERROR: Ill-defined RC tuple");
        
        
        if(pos.size()!=tuple2.pos.size())
            return false;
        
        //tuples are well-defined and same size...check position by position
        for(int index=0; index<pos.size(); index++){
            if(pos.get(index)!=tuple2.pos.get(index))
                return false;
            if(RCs.get(index)!=tuple2.RCs.get(index))
                return false;
        }
        
        //if we get here they're the same
        return true;
    }
    
    
    public String stringListing(){
        //Listing the RCs in a string
        String ans = "";
        
        for(int posNum=0; posNum<pos.size(); posNum++){
            ans = ans + "Res " + pos.get(posNum) + " RC " + RCs.get(posNum) + " ";
        }
        
        return ans;
    }
    
    public String stringListingShort(){
    	if(RCPairs.size() < 1)
    	{
    		for(int RCIndex = 0; RCIndex < size(); RCIndex++)
    		{
    			RCPairs.add(new RCPair(pos.get(RCIndex),RCs.get(RCIndex)));
    		}
    	}
    	Collections.sort(RCPairs);
        //Listing the RCs in a string
        String ans = "(";
        
        for(int RCIndex=0; RCIndex<RCPairs.size(); RCIndex++){
        	RCPair pair = RCPairs.get(RCIndex);
            ans = ans + pair.pos + ":" + pair.RC + ", ";
        }
        
        return ans+")";
    }
    
    public String toString(){
    	return stringListingShort();
    }
    
    
    public RCTuple subtractMember(int index){
        //Make a copy of this RCTuple with the given member removed
        //index is an index in pos and RCs
        ArrayList<Integer> newPos = new ArrayList<>();
        ArrayList<Integer> newRCs = new ArrayList<>();
        
        for(int ind=0; ind<pos.size(); ind++){
            if(ind!=index){
                newPos.add(pos.get(ind));
                newRCs.add(RCs.get(ind));
            }
        }
        
        return new RCTuple(newPos,newRCs);
    }
    
    @SuppressWarnings("unchecked")
	public RCTuple addRC(int addedPos, int addedRC){
        //Make a copy of this RCTuple with (addPos,addRC) added
        ArrayList<Integer> newPos = (ArrayList<Integer>) pos.clone();
        ArrayList<Integer> newRCs = (ArrayList<Integer>) RCs.clone();
        
        newPos.add(addedPos);
        newRCs.add(addedRC);
        
        return new RCTuple(newPos,newRCs);
    }
    
    public int hashCode()
    {
    	return toString().hashCode();
    }
    
    private class RCPair implements Comparable<RCPair>{

    	int pos;
    	int RC;
    	public RCPair(int position, int RCAssignment) {
    		pos = position;
    		RC = RCAssignment;
    	}
    	

		@Override
		public int compareTo (RCPair arg0) {
			return pos - arg0.pos;
		}
    }
}
