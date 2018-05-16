package edu.duke.cs.osprey.sparse;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.TreeSet;

import edu.duke.cs.osprey.confspace.RCTuple;

public class SubproblemConfEnumerator implements ConformationProcessor {

	private PartialConformationEnergyFunction energyFunction;
	private ArrayList<Map<String, PriorityQueue<ScoredAssignment>>> lambdaHeaps; 
	private ArrayList<PriorityQueue<ScoredAssignment>> templateHeaps;
	private Subproblem sourceProblem;
	private SubproblemConfEnumerator leftSubproblem;
	private SubproblemConfEnumerator rightSubproblem;
	private ChildConfManager childConfs;
	private static RCTuple emptyConf = new RCTuple();
	public static boolean debugOutput = true;
	
	public SubproblemConfEnumerator (Subproblem subproblem, PartialConformationEnergyFunction eFunc) {
		subproblem.addConformationProcessor(this);
		addChildProcessors(subproblem, eFunc);

		energyFunction = eFunc;
		sourceProblem = subproblem;
		int subproblemLocalConfs = subproblem.getTotalMConformations().intValue();
		lambdaHeaps = new ArrayList<>();
		templateHeaps = new ArrayList<>();
		while(lambdaHeaps.size() <= subproblemLocalConfs)
			lambdaHeaps.add(null);
		while(templateHeaps.size() <= subproblemLocalConfs)
			templateHeaps.add(null);
		if(leftSubproblem != null)
		{
			RightConfManager rightSide = null;
			if(rightSubproblem != null)
			{
				rightSide = new RightConfManager(rightSubproblem, sourceProblem.rightSubproblem);
			}
			childConfs = new ChildConfManager(leftSubproblem, rightSide);
		}
	}

	private void addChildProcessors (Subproblem subproblem, PartialConformationEnergyFunction eFunc) {
		if(subproblem.leftSubproblem != null)
		{
			leftSubproblem = new SubproblemConfEnumerator(subproblem.leftSubproblem, eFunc);
		}
		if(subproblem.rightSubproblem != null)
		{
			rightSubproblem  = new SubproblemConfEnumerator(subproblem.rightSubproblem, eFunc);
		}
	}

	@Override
	public void processConformation (RCTuple conformation) {
		RCTuple MAssignment = sourceProblem.extractSubproblemMAssignment(conformation);
		RCTuple lambdaAssingment = sourceProblem.extractSubproblemLambdaAssignment(conformation);
		double selfEnergy = energyFunction.computePartialEnergyGivenPriorConformation(MAssignment, lambdaAssingment);
		double leftEnergy = 0;
		if(leftSubproblem != null)
			leftEnergy = leftSubproblem.nextBestEnergy(conformation);
		double rightEnergy = 0;
		if(rightSubproblem != null)
			rightEnergy = rightSubproblem.nextBestEnergy(conformation);

		PriorityQueue<ScoredAssignment> templateHeap = getTemplateHeap(conformation);
		ScoredAssignment assignment = new ScoredAssignment(conformation, selfEnergy, leftEnergy + rightEnergy, 0);
		debugPrint("Processing "+conformation+", adding new template conf "+assignment);
		templateHeap.add(assignment);
		
	}
	


	public RCTuple nextBestConformation(RCTuple queryAssignment)
	{
		PriorityQueue<ScoredAssignment> lambdaHeap = getHeap(queryAssignment);
		debugPrint("===================== Start "+sourceProblem+" =============================================");
		debugPrint("Beginning recursion with "+queryAssignment+" at \n"+sourceProblem);
		printHeap(lambdaHeap);
		if(lambdaHeap.size() < 1)
		{
			System.err.println("Should not have polled, empty heap at \n"+sourceProblem);
		}
		ScoredAssignment previousHeapRoot = lambdaHeap.poll();
		double sanityCheckEnergy = previousHeapRoot.score();
		RCTuple curBestAssignment = previousHeapRoot.assignment;
		RCTuple outputAssignment = previousHeapRoot.assignment.copy();

		outputAssignment = outputAssignment.combineRC(queryAssignment);
		RCTuple combinedqueryAssignment = queryAssignment.combineRC(curBestAssignment);

		if(childConfs != null)
		{
			debugPrint("Processing children with "+combinedqueryAssignment);
			if(!childConfs.hasMoreConformations(combinedqueryAssignment))
			{
				System.err.println("Attempting to get new child conf but there are none.");
			}
			assert(childConfs.hasMoreConformations(combinedqueryAssignment));
			RCTuple nextBestChildConf = childConfs.getNextChildAssignment(combinedqueryAssignment);
			outputAssignment = outputAssignment.combineRC(nextBestChildConf);
			if(childConfs.hasMoreConformations(combinedqueryAssignment))
			{
				debugPrint("Re-adding "+combinedqueryAssignment+" to :\n"+sourceProblem);
				previousHeapRoot.updateLeftScore(childConfs.nextBestEnergy(combinedqueryAssignment));
				lambdaHeap.add(previousHeapRoot);
			}
			else
			{
				debugPrint("Removing "+combinedqueryAssignment+" from :\n"+sourceProblem);
				debugPrint("Remaining confs:");
				debugPrint(lambdaHeap+"-"+lambdaHeap.hashCode());
			}
		}


		debugPrint("Returning "+outputAssignment+" from \n"+sourceProblem);
		debugPrint("Heap status: ");
		printHeap(lambdaHeap);
		debugPrint("===================== End "+sourceProblem+" =============================================");
		if(lambdaHeap.size() > 0 && lambdaHeap.peek().score() < sanityCheckEnergy)
		{
			System.err.println("The next conf will have a LOWER energy: "+lambdaHeap.peek().score()+" < "+sanityCheckEnergy);
			debugPrint(lambdaHeap.peek());
			debugPrint(outputAssignment+", energy "+sanityCheckEnergy);
		}
		assert(lambdaHeap.size() < 1 || lambdaHeap.peek().score() >= sanityCheckEnergy);
		return outputAssignment;
	}
	

	
	private void checkHeap(PriorityQueue<ScoredAssignment> templateHeap)
	{
		for(ScoredAssignment conf : templateHeap)
		{
			if(conf.assignment == null)
			{
				debugPrint("Null assignment. Weird...");
			}
			if(!sourceProblem.isValidConf(conf.assignment))
			{
				debugPrint("Invalid conf in heap.");
			}
		}
	}

	private PriorityQueue<ScoredAssignment> getHeap (RCTuple queryAssignment) {
		int lambdaHeapIndex = sourceProblem.mapSubproblemConfToIndex(queryAssignment);
		String RCTupleKey = queryAssignment.toString();
		//TODO: Add a check to see if the queryAssignment belongs to a template heap.
		if(lambdaHeaps.size() <= lambdaHeapIndex || lambdaHeaps.get(lambdaHeapIndex) == null)
		{
			Map<String, PriorityQueue<ScoredAssignment>> heapMap = new HashMap<>();
			lambdaHeaps.set(lambdaHeapIndex, heapMap);
		}
		if(!lambdaHeaps.get(lambdaHeapIndex).containsKey(RCTupleKey))
			initializeHeap(queryAssignment, RCTupleKey);
		PriorityQueue<ScoredAssignment> output = lambdaHeaps.get(lambdaHeapIndex).get(RCTupleKey);
		if(output == null)
		{
			System.err.println("Heap not found.");
		}
		checkHeap(output);
		if(!output.isEmpty() && !queryAssignment.consistentWith(output.peek().assignment))
		{
			System.err.println("ERROR: Heap does not match query assignment!!");
		}
		return output;
	}
	
	private PriorityQueue<ScoredAssignment> getTemplateHeap(RCTuple MULambdaAssignment)
	{

		int templateHeapIndex = sourceProblem.mapSubproblemConfToIndex(MULambdaAssignment);
		if(templateHeaps.get(templateHeapIndex) == null)
		{
			assert(sourceProblem.isMULambdaAssignment(MULambdaAssignment));
			templateHeaps.set(templateHeapIndex, new PriorityQueue<>());
		}
		return templateHeaps.get(templateHeapIndex);
	}
	

	private void initializeHeap (RCTuple queryAssignment, String RCTupleKey) {
		int lambdaHeapIndex = sourceProblem.mapSubproblemConfToIndex(queryAssignment);
		RCTuple MTuple = sourceProblem.extractSubproblemAssignment(queryAssignment);

		Map<String, PriorityQueue<ScoredAssignment>> heapMap = lambdaHeaps.get(lambdaHeapIndex);
		PriorityQueue<ScoredAssignment> newHeap = new PriorityQueue<ScoredAssignment>();
		for(ScoredAssignment conf: getTemplateHeap(MTuple))
		{
			newHeap.add(conf.copy());
		}
		heapMap.put(RCTupleKey, newHeap);
	}
	
	public double nextBestEnergy()
	{
		return nextBestEnergy(emptyConf);
	}
	
	public double nextBestEnergy(RCTuple queryConf)
	{
		PriorityQueue<ScoredAssignment> heap = getHeap(queryConf);
		if(heap.size() < 1)
		{
			debugPrint("No confs at "+sourceProblem);
			getHeap(queryConf);
		}
		debugPrint("Getting next best energy for "+queryConf+" at "+sourceProblem);
		printHeap(heap);
		ScoredAssignment bestConf = heap.peek();
		double bestScore = bestConf.score();

		debugPrint("Next best energy is "+heap.peek().score());
		return bestScore;
	}
	
	public boolean hasMoreConformations()
	{
		return hasMoreConformations(emptyConf);
	}
	

	public boolean hasMoreConformations(RCTuple queryConf)
	{
		if(getHeap(queryConf).size() < 1)
			return false;
		if(childConfs!= null)
		{
			RCTuple topHeapConf = getHeap(queryConf).peek().assignment.copy();
			topHeapConf = topHeapConf.combineRC(queryConf);
			
			if(!childConfs.hasMoreConformations(topHeapConf) && getHeap(queryConf).size() > 0)
			{
				debugPrint("This will fail. No child confs, but still have lambda confs.");
			}
			
		}

		return getHeap(queryConf).size() > 0;
	}
	

	public RCTuple peekNextBestConformation (RCTuple queryConf) {
		PriorityQueue<ScoredAssignment> heap = getHeap(queryConf);
		if(heap.size() < 1)
		{
			debugPrint("No confs...");
			getHeap(queryConf);
		}
		ScoredAssignment bestConf = heap.peek();
		return bestConf.assignment.copy();
	}

	@Override
	public boolean recurse () {
		return false;
	}
	
	private void debugPrint(Object print)
	{
		if(debugOutput)
			System.out.println(print);
	}
	
	private void printHeap(PriorityQueue<ScoredAssignment> heap)
	{
		if(!debugOutput)
			return;
		debugPrint("Heap "+heap.hashCode()+":");
		if(heap instanceof LazyHeap)
			debugPrint("Lazy Heap default Energy: "+((LazyHeap) heap).initialRightEnergy);
		PriorityQueue<ScoredAssignment> cleanup = new PriorityQueue<>();
		int maxConfsToPrint = 10;
		int numPrinted = 0;
		while(!heap.isEmpty() && numPrinted < maxConfsToPrint)
		{
			ScoredAssignment conf = heap.poll();
			cleanup.add(conf);
			debugPrint(conf);
			numPrinted++;
		}
		while(!cleanup.isEmpty())
			heap.add(cleanup.poll());
	}
	
	
	/*** This class will encapsulate ALL of the crazy logic that 
	 * goes into maintaining the heaps required for sparse enumeration.
	 * @author Jon
	 *
	 */
	
	private class ChildConfManager
	{
		private Map<String, PriorityQueue<ScoredAssignment>> leftHeapMap = new HashMap<>();;
		private RightConfManager rightConfs;
		private SubproblemConfEnumerator leftSubproblem;
		
		public ChildConfManager(SubproblemConfEnumerator leftEnumerator, RightConfManager rightManager)
		{
			leftSubproblem = leftEnumerator;
			rightConfs = rightManager;
		}
		
		public RCTuple getNextChildAssignment(RCTuple queryAssignment)
		{
			RCTuple nextBestChildConf = null;
			if(rightConfs != null)
			{
				debugPrint("Two children found, recursing via secondary heap...");
				PriorityQueue<ScoredAssignment> leftTrackingHeap = getLeftHeapMap(queryAssignment);
				printHeap(leftTrackingHeap);

				if(leftTrackingHeap.size() < 1)
				{
					System.err.println("Tracking heap is empty??");
				}
				ScoredAssignment nextBestChildAssignment = leftTrackingHeap.poll();
				debugPrint("Polled "+nextBestChildAssignment);
				if(nextBestChildAssignment == null)
				{
					System.err.println("Null best child assignment...");
				}
				
				RCTuple bestChildConf = nextBestChildAssignment.assignment;
				debugPrint("Best left conf is "+bestChildConf+", maching to right side...");

				RCTuple nextBestRightConf = rightConfs.getNextBestRightConf(bestChildConf);
				debugPrint("RightConf returned is "+nextBestRightConf);
				updateChildAssignment(queryAssignment,nextBestChildAssignment);
				nextBestChildConf = bestChildConf.combineRC(nextBestRightConf);

				debugPrint("Left heap final status: ");
				printHeap(leftTrackingHeap);
			}
			else
			{
				debugPrint("Recursing to left child...");
				nextBestChildConf = leftSubproblem.nextBestConformation(queryAssignment);
				nextBestChildConf = nextBestChildConf.combineRC(queryAssignment);
			}
			assert(nextBestChildConf != null);
			return nextBestChildConf;

		}

		public double nextBestEnergy (RCTuple queryAssignment) {
			if(rightConfs != null)
			{
				return getLeftHeapMap(queryAssignment).peek().score();
			}
			else
			{
				return leftSubproblem.nextBestEnergy(queryAssignment);
			}
		}
		
		private PriorityQueue<ScoredAssignment> getLeftHeapMap(RCTuple queryAssignment)
		{
			if(!leftHeapMap.containsKey(queryAssignment.toString()))
			{
				RCTuple leftMAssignment = sourceProblem.leftSubproblem.extractSubproblemMAssignment(queryAssignment);
				double initialRightenergy = rightConfs.peekNextBestEnergy(queryAssignment);
				PriorityQueue<ScoredAssignment> newHeap = new LazyHeap(queryAssignment, leftSubproblem, initialRightenergy);
				if(newHeap.size() < 1)
				{
					System.err.println("Empty new heap. Should be full of template nodes.");
					leftSubproblem.getHeap(leftMAssignment);
				}
				leftHeapMap.put(queryAssignment.toString(), newHeap);
			}
			
			return leftHeapMap.get(queryAssignment.toString());
		}

		public boolean hasMoreConformations(RCTuple queryAssignment)
		{
			if(rightConfs!=null)
			{
				return getLeftHeapMap(queryAssignment).size() > 0;
			}
			else
			{
				return leftSubproblem.hasMoreConformations(queryAssignment);
			}
		}

		private void updateChildAssignment (RCTuple queryAssignment, ScoredAssignment nextBestChildAssignment) {
			RCTuple nextBestChildQueryConf = nextBestChildAssignment.assignment;
			if(rightConfs.hasMoreConformations(nextBestChildQueryConf))
			{
				debugPrint("More right confs for left assignment "+nextBestChildQueryConf+" at:\n"+sourceProblem);
				nextBestChildAssignment.updateLeftScore(rightConfs.peekNextBestEnergy(nextBestChildQueryConf));
				leftHeapMap.get(queryAssignment.toString()).add(nextBestChildAssignment);
			}
			else 
			{
				debugPrint("No more right confs for left assignment "+nextBestChildQueryConf+" at:\n"+sourceProblem);
				/*
				if(leftSubproblem.hasMoreConformations(queryAssignment))
				{
					double nextLeftEnergy = leftSubproblem.nextBestEnergy(queryAssignment);
					RCTuple nextLeftConf = leftSubproblem.nextBestConformation();
					ScoredAssignment nextLeftTrackingConf = new ScoredAssignment(nextLeftConf, nextLeftEnergy,0,0);
					leftHeapMap.get(queryAssignment.toString()).add(nextLeftTrackingConf);
				}
				*/
			}
			
		}
	
		
		
	}
	
	private class RightConfManager 
	{
		private Subproblem rightSubproblem;
		private SubproblemConfEnumerator rightSubproblemEnum;
		private Map<String, RightConf> rightConfMap = new HashMap<>();
		private Map<String, List<ScoredAssignment>> rightConfLists = new HashMap<>();

		public RightConfManager(SubproblemConfEnumerator rightSideEnum, Subproblem rightProblem)
		{
			rightSubproblemEnum = rightSideEnum;
			rightSubproblem = rightProblem;
		}
		
		private List<ScoredAssignment> getRightConfList (RCTuple queryConf)
		{
			if(!rightConfLists.containsKey(queryConf.toString()))
			{
				RCTuple templateConf = rightSubproblem.extractSubproblemMAssignment(queryConf);
				if(!rightConfLists.containsKey(templateConf.toString()))
				{
					debugPrint("Creating new conf list for "+queryConf);
					List<ScoredAssignment> newConfList = new UniqueConfList();
					double nextRightConfE = rightSubproblemEnum.nextBestEnergy(templateConf);
					debugPrint("Polling new conformation from right side with "+templateConf+"...");
					RCTuple nextRightConf = rightSubproblemEnum.nextBestConformation(templateConf);
					RCTuple rightPart = rightSubproblemEnum.sourceProblem.extractSubproblemLAssignment(nextRightConf);
					newConfList.add(new ScoredAssignment(rightPart, nextRightConfE, 0, 0));
					assert(newConfList.size() == 1);
					rightConfLists.put(templateConf.toString(), newConfList);
					debugPrint("Created new conf list for "+queryConf+": "+newConfList+";"+newConfList.hashCode());
				}
				rightConfLists.put(queryConf.toString(), rightConfLists.get(templateConf.toString()));
			}
			return rightConfLists.get(queryConf.toString());
		}
		

		private RightConf getRightConf(RCTuple queryConf)
		{
			if(!rightConfMap.containsKey(queryConf.toString()))
			{
				debugPrint("Creating new RightConf for "+queryConf);
				rightConfMap.put(queryConf.toString(), new RightConf(
						queryConf, getRightConfList(queryConf)));
			}
			return rightConfMap.get(queryConf.toString());
		}
		
		public double peekNextBestEnergy(RCTuple queryConf)
		{
			return getRightConf(queryConf).peekConfEnergy();
		}
		
		public RCTuple getNextBestRightConf(RCTuple queryConf)
		{
			RightConf rightConf = getRightConf(queryConf);
			RCTuple output = rightConf.pollCurConf();
			rightConf.updateConf(rightSubproblemEnum);
			return output;
		}
		
		public boolean hasMoreConformations(RCTuple leftConf)
		{
			return getRightConf(leftConf).hasMoreConformations();
		}
	}
	
	private class RightConf
	{
		private RCTuple queryAssignment;
		private int confListIndex = 0;
		private List<ScoredAssignment> confList;
		private SubproblemConfEnumerator rightSubproblem;
		
		public RightConf(RCTuple conf, List<ScoredAssignment> rightConfs){
			queryAssignment = conf;
			confList = rightConfs;
			checkList(confList);
		}

		public void updateConf (SubproblemConfEnumerator rightSubproblem) {
			RCTuple rightMAssignment = rightSubproblem.sourceProblem.extractSubproblemMAssignment(queryAssignment);
			if(!hasMoreConformations())
			{
				debugPrint("Used up existing "+queryAssignment+", checking to see if there are more...");
				debugPrint("Querying for more right conformations with "+rightMAssignment);
			}
			if(!hasMoreConformations() && rightSubproblem.hasMoreConformations(rightMAssignment))
			{
				debugPrint("Appending new right Conformation for "+queryAssignment);
				double nextRightConfE = rightSubproblem.nextBestEnergy(rightMAssignment);
				RCTuple nextRightConf = rightSubproblem.nextBestConformation(rightMAssignment);
				RCTuple rightPart = rightSubproblem.sourceProblem.extractSubproblemLAssignment(nextRightConf);
				confList.add(new ScoredAssignment(rightPart, nextRightConfE, 0, 0));
			}
		}

		public double peekConfEnergy () {
			return confList.get(confListIndex).score();
		}

		public boolean hasMoreConformations () {
			if(confListIndex < confList.size())
			{
				debugPrint("Assignment "+queryAssignment+" has more confs. Specifically: "+confList.get(confListIndex));
				debugPrint(confList);
			}
			return confListIndex < confList.size();
		}

		public RCTuple pollCurConf () {

			checkList(confList);

			if(confListIndex >= confList.size())
			{
				System.err.println("No confs...");
			}
			assert(confListIndex < confList.size());
			debugPrint("Using conf "+confListIndex+" from list "+confList+";"+confList.hashCode());
			RCTuple output = confList.get(confListIndex).assignment;
			confListIndex++;
			debugPrint("Polled conf "+output+", index advanced to "+confListIndex);
			return output;
		}
		
		public RCTuple peekCurConf () {
			return confList.get(confListIndex).assignment;
		}
		
		private void checkList(List<ScoredAssignment> confs)
		{
			if(confs == null)
			{
				System.err.println("No confList. Why do we exist??");
			}
				
			Set<String> confSet = new TreeSet<String>();
			for(ScoredAssignment conf : confs)
			{
				String assignmentString = conf.assignment.toString();
				if(confSet.contains(assignmentString))
				{
					System.err.println("Failure. Duplicate right conformations...");
				}
				assert(!confSet.contains(assignmentString));
				confSet.add(assignmentString);
			}
		}
	}
	
	private class LazyHeap extends PriorityQueue<ScoredAssignment>
	{
		private boolean dirty;
		private ScoredAssignment cleanAssignment = null;
		private SubproblemConfEnumerator leftSubproblem;
		private RCTuple queryAssignment;
		private double initialRightEnergy = 0;
		
		public LazyHeap(RCTuple queryConf, SubproblemConfEnumerator leftChild, double defaultEnergy)
		{
			if(defaultEnergy > 0)
			{
				debugPrint("Every one of these conformations will have a positive right energy. check this.");
			}
			queryAssignment = queryConf;
			leftSubproblem = leftChild;
			initialRightEnergy = defaultEnergy;
			addNewLeftConf();
			peek();
		}
		

		
		private void addNewLeftConf()
		{
			double nextBestEnergy = leftSubproblem.nextBestEnergy(queryAssignment);
			RCTuple nextBestLeftConf = leftSubproblem.nextBestConformation(queryAssignment);
			cleanAssignment = new ScoredAssignment(nextBestLeftConf, nextBestEnergy, initialRightEnergy, 0);

			debugPrint("Creating new left conf "+cleanAssignment);
			add(cleanAssignment);
		}
		
		public ScoredAssignment poll()
		{
			cleanHeap();
			ScoredAssignment nextBestAssignment = super.poll();
			if(nextBestAssignment == cleanAssignment)
				dirty = true;

			cleanHeap();
			return nextBestAssignment;
		}
		
		private void cleanHeap()
		{
			if((size() < 1 || dirty) && leftSubproblem.hasMoreConformations(queryAssignment))
			{
				addNewLeftConf();
				dirty = false;
			}
		}
		
		public ScoredAssignment peek()
		{
			cleanHeap();
			ScoredAssignment output = super.peek();
			return output;
		}
		
		public String toString()
		{
			return "Default Right Energy: "+initialRightEnergy+", heap: "+super.toString();
		}
		
	}
	

	public RCTuple nextBestConformation () {
		return nextBestConformation(emptyConf);
	}
	
	private class UniqueConfList extends LinkedList<ScoredAssignment>
	{
		private Set<ScoredAssignment> confSet = new TreeSet<ScoredAssignment>();
		@Override
		public boolean add(ScoredAssignment conf)
		{
			if(confSet.contains(conf))
				System.err.println("Dupe conf.");
			assert(!confSet.contains(conf));
			return super.add(conf);
		}
		
	}
}




/* Old code follows.
public void bTrackBestConfRemoveEarlyNew(RotTypeMap bestPosAARot[], int[] bestState)
{
    boolean reinsert = false;
    PriorityQueue<Conf> outHeap = getHeap(bestPosAARot, bestState, ""+this.hashCode());
    
    
    if(printHeap)
        outputInitialDebugData(bestPosAARot, outHeap);
    
    if(outHeap.size() > 1 && lambda.size() < 1)
        debugPrint("Impossibiruuuu");
    debugPrint("Begin. Polling heap...");
    Conf nextState = outHeap.poll();
    // Add lambda to the solution 
    nextState.fillRotTypeMap(bestPosAARot);
    RotTypeMap[] bestPosAARotOld = Arrays.copyOf(bestPosAARot, bestPosAARot.length);
    
    double curBest = nextState.energy;
    debugPrint("Next state is "+nextState+", energy "+nextState.energy);

    if(rightChild != null)
    {
    	debugPrint("right child is not null, performing left and right operations...");
        TreeEdge leftEdge = leftChild.getCofEdge();
        int[] leftMLambda = nextState.conformation; 
        int[] leftM = getMstateForEdgeCurState(leftMLambda, leftEdge);
        TreeEdge rightEdge = rightChild.getCofEdge();
        
        debugPrint("Getting right solutoins for "+RTMToString(bestPosAARot));
        List<RightConf> rightConfs = getRightSolutions(bestPosAARot);
        if(rightConfs.size() < 1)
        {
        	debugPrint("Generating new list, list is empty...");
            int[] rightMLambda = nextState.conformation; 
            int[] rightM = getMstateForEdgeCurState(rightMLambda, rightEdge);
            getNewRightconf(bestPosAARotOld, rightEdge, rightConfs, rightM);
        }
        debugPrint("Getting secondary heap from "+leftEdge.L+leftEdge.lambda+"...");
        LazyHeap<Conf> secondaryHeap = getSecondaryHeap(bestPosAARotOld, leftM);
        if(leftChild.getCofEdge().lambda.size() < 1 && secondaryHeap.size() > 1)
            debugPrint("IMPOSSIBIRU!?!?!!");
        if((secondaryHeap.dirty || secondaryHeap.size() < 1) && leftEdge.moreConformations(bestPosAARotOld, leftM)) 
        {
        	debugPrint("Populating heap with new conformation...");
            // Acquire new left conformation for the heap 
            RotTypeMap[] bestPosAARotLeft = Arrays.copyOf(bestPosAARotOld, bestPosAARot.length);
            
            double newLeftEnergy = leftEdge.peekEnergy(bestPosAARotLeft, leftM);
            leftEdge.bTrackBestConfRemoveEarlyNew(bestPosAARotLeft, leftM);
            Conf newLeftConf = new Conf(bestPosAARotLeft, newLeftEnergy);
            newLeftConf.updateLeftEnergy(rightConfs.get(0).energy);
            secondaryHeap.cleanNode = newLeftConf;
            secondaryHeap.add(newLeftConf);
            secondaryHeap.dirty = false;
            if(leftChild.getCofEdge().lambda.size() < 1 && secondaryHeap.size() > 1)
                debugPrint("IMPOSSIBIRU!?!?!!");
        }
        
        // Maintain cleanliness 
		PriorityQueue<Conf> copy = new PriorityQueue<Conf>();

		for(Conf c : secondaryHeap)
		{
		    copy.add(c);
		}

		while(!copy.isEmpty())
		{
		    Conf c = copy.poll();
		    debugPrint(c+"$"+c.energy+", left "+c.leftEnergy+", self "+c.selfEnergy+", right "+c.rightEnergy);
		}
        Conf leftConf = secondaryHeap.poll();
        debugPrint("Polled left conformation from secondary heap: "+leftConf+", "+leftConf.energy);

		PriorityQueue<Conf> copy2 = new PriorityQueue<Conf>();

		for(Conf c : secondaryHeap)
		{
		    copy2.add(c);
		}

		while(!copy2.isEmpty())
		{
		    Conf c = copy2.poll();
		    debugPrint(c+"$"+c.energy+", left "+c.leftEnergy+", self "+c.selfEnergy+", right "+c.rightEnergy);
		}
        leftConf.fillConf(bestPosAARot);
        
        if(leftConf == secondaryHeap.cleanNode)
            secondaryHeap.dirty = true;
        int index = getRightOffset(leftConf.toString());
        debugPrint("Right offset index is "+index);

        int[] rightMLambda = nextState.conformation; 
        int[] rightM = getMstateForEdgeCurState(rightMLambda, rightEdge);
        while(rightConfs.size() <= index + 1 && rightEdge.moreConformations(bestPosAARotOld, rightM))
        {
            debugPrint("Getting extra conformation for next run...");
            getNewRightconf(bestPosAARotOld, rightEdge, rightConfs, rightM);
        }
                    
        RightConf rightConf = rightConfs.get(index);    
        rightConf.fillRotTypeMap(bestPosAARot);        
        debugPrint("Returned right conformation is "+rightConf+", energy "+rightConf.energy);
        
        if(index + 1 < rightConfs.size() || rightEdge.moreConformations(bestPosAARotOld, rightM))
        {
            RightConf nextRightConf = rightConfs.get(index + 1);
            leftConf.updateLeftEnergy(nextRightConf.energy);
            rightSolutionOffset.put(leftConf.toString(), index + 1);
            debugPrint("Reinsert "+leftConf+" with new energy "+leftConf.energy
            		+" into secondary heap, right conf is "+nextRightConf+", energy "+nextRightConf.energy);
            secondaryHeap.add(leftConf);
            reinsert = true;
        }
        
        if((secondaryHeap.dirty || secondaryHeap.size() < 1) && leftEdge.moreConformations(bestPosAARotOld, leftM)) 
        {
            // Acquire new left conformation for the heap 
            RotTypeMap[] bestPosAARotLeft = Arrays.copyOf(bestPosAARotOld, bestPosAARot.length);
            
            double newLeftEnergy = leftEdge.peekEnergy(bestPosAARotLeft, leftM);
            leftEdge.bTrackBestConfRemoveEarlyNew(bestPosAARotLeft, leftM);
            Conf newLeftConf = new Conf(bestPosAARotLeft, newLeftEnergy);
            newLeftConf.updateLeftEnergy(rightConfs.get(0).energy);
            secondaryHeap.cleanNode = newLeftConf;
            secondaryHeap.add(newLeftConf);
            secondaryHeap.dirty = false;
            if(leftChild.getCofEdge().lambda.size() < 1 && secondaryHeap.size() > 1)
                debugPrint("IMPOSSIBIRU!?!?!!");
        }

        if(secondaryHeap.size() > 0)
        {
            nextState.updateLeftEnergy(secondaryHeap.peek().energy);
            reinsert = true;
        }
        else
            reinsert = false;
        
    }

    else if(leftChild != null)
    {
    	debugPrint("Only one child. Recursing...");
        TreeEdge leftEdge = leftChild.getCofEdge();
        int[] leftMLambda = nextState.conformation; 
        int[] leftM = getMstateForEdgeCurState(leftMLambda, leftEdge);
        if(leftM == null)
            leftM = leftMLambda;
        leftEdge.bTrackBestConfRemoveEarlyNew(bestPosAARot, leftM);        
        nextState.updateLeftEnergy(leftEdge.peekEnergy(bestPosAARot, leftM));
        debugPrint(nextState+" has new energy "+nextState.energy);
        reinsert = leftEdge.moreConformations(bestPosAARot, leftM);
    }
    if(reinsert)
        outHeap.add(nextState);
    
    checkHeap(outHeap);
    if(outHeap.size() > 1)
    {
        double curNewBest = outHeap.peek().energy;
        if(curNewBest < curBest)
        {
            System.err.println("Out of order. Terminating...");
            //System.exit(-1);
        }
    }
    if(printHeap)
        printEndDebugInfo(bestPosAARot, outHeap);

}
*/