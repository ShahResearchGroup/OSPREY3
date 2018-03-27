package edu.duke.cs.osprey.sparse;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.TreeSet;
import edu.duke.cs.osprey.confspace.RCTuple;


public class SubproblemSeqEnumerator implements ConformationProcessor {
	
	private PartialConformationEnergyFunction energyFunction;
	private ArrayList<Map<String, PriorityQueue<ScoredAssignment>>> lambdaHeaps; 
	private ArrayList<PriorityQueue<ScoredAssignment>> templateHeaps;
	private Map<String, Double> partialSums;
	private Subproblem sourceProblem;
	private SubproblemSeqEnumerator leftSubproblem;
	private SubproblemSeqEnumerator rightSubproblem;
	private ChildConfManager childConfs;
	private static RCTuple emptyConf = new RCTuple();
	public static boolean debugOutput = true;
	
	public SubproblemSeqEnumerator (Subproblem subproblem, PartialConformationEnergyFunction eFunc) {
		subproblem.addConformationProcessor(this);
		addChildProcessors(subproblem, eFunc);

		energyFunction = eFunc;
		sourceProblem = subproblem;
		int subproblemLocalConfs = subproblem.getTotalMConformations().intValue();
		lambdaHeaps = new ArrayList<>();
		templateHeaps = new ArrayList<>();
		partialSums = new HashMap<>();
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
			leftSubproblem = new SubproblemSeqEnumerator(subproblem.leftSubproblem, eFunc);
		}
		if(subproblem.rightSubproblem != null)
		{
			rightSubproblem  = new SubproblemSeqEnumerator(subproblem.rightSubproblem, eFunc);
		}
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

	@Override
	public boolean recurse () {
		// TODO Auto-generated method stub
		return false;
	}
	

	
	/***
	 * This function takes a complete sequence and returns
	 * the corresponding partial Sparse K* score.
	 * @param queryConf
	 * @return
	 */
	public double getlogPartialKStarBound(RCTuple queryConf)
	{
		/*
		 * 1. Map conf to sequence.
		 * 2. Get Sequence heap
		 * 3. For each conf in the sequence heap, get the child conf 
		 * 		partial sum.
 		 * 4. Return sum of child conf partial sums.
		 */
		KStarHeap heap = getMultiSequenceHeap(queryConf);
		double sum = childConfs.getChildKStarMultiSequenceBounds(queryConf);
		double localSum = heap.getHeapSum();
		return sum+localSum;
	}
	
	/***
	 * This function takes a complete sequence and returns
	 * the corresponding partial Sparse K* score.
	 * @param queryConf
	 * @return
	 */
	public double getlogPartialKStarScore(RCTuple queryConf)
	{
		/*
		 * 1. Map conf to sequence.
		 * 2. Get Sequence heap
		 * 3. For each conf in the sequence heap, get the child conf 
		 * 		partial sum.
 		 * 4. Return sum of child conf partial sums.
		 */
		KStarHeap heap = getSingleSequenceHeap(queryConf);
		double sum = childConfs.getChildKStarSingleSequenceSums(queryConf);
		double localSum = heap.getHeapSum();
		return sum+localSum;
	}
	

	private KStarHeap getSingleSequenceHeap(RCTuple queryConf) {
		/*
		 * 1. Map conf to sequence.
		 * 2. Get Sequence heap
 		 * 3. Return heap.
		 */
		return null;
	}
	
	/***
	 * This function takes a MUlambda assignment and returns the
	 * corresponding heap used to make the multisequence bound.
	 * @param queryConf
	 * @return
	 */
	public KStarHeap getMultiSequenceHeap(RCTuple queryConf)
	{
		return null;
	}

	@Override
	public void processConformation (RCTuple conformation) {
		/*
		 * Conformation processing pseudocode:
		 * 1. Map conformation to its M and lambda sequences
		 * 2. Add conformation energy to its respective M Conf MS bound heap (one MS bound
		 *       per M conf per lambda sequence)
		 */
		String sequence = sourceProblem.getSequenceForConf(conformation);
		System.out.println("Processing "+conformation+", sequence ["+sequence+"]");
		RCTuple MAssignment = sourceProblem.extractSubproblemMAssignment(conformation);
		RCTuple lambdaAssignment = sourceProblem.extractSubproblemLambdaAssignment(conformation);
		double energy = energyFunction.computePartialEnergyGivenPriorConformation(MAssignment,lambdaAssignment);
		double runningPartialSum = partialSums.get(sequence);
		partialSums.put(sequence, runningPartialSum+energy);
		if(leftSubproblem == null && rightSubproblem == null)
		{
			/* Process partial score per conformation */
		}

		/*
		 * Sometime after this, we need to:
		 * 1. For each M sequence, contruct the heap of lambda sequences
		 * 2. For each lambda sequence for each M conf, add the MS bounds for each lambda
		 *      sequence to compute the overall MS bound for each lambda sequence
		 *      given an M sequence.
		 * 3. Insert these sequences into the M Sequence lambda sequence heap.
		 */
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
	
	/*** This class will encapsulate ALL of the crazy logic that 
	 * goes into maintaining the heaps required for sparse enumeration.
	 * @author Jon
	 *
	 */
	
	private class ChildConfManager
	{
		private Map<String, PriorityQueue<ScoredAssignment>> leftHeapMap = new HashMap<>();;
		private RightConfManager rightConfs;
		private SubproblemSeqEnumerator leftSubproblem;
		
		public ChildConfManager(SubproblemSeqEnumerator leftEnumerator, RightConfManager rightManager)
		{
			leftSubproblem = leftEnumerator;
			rightConfs = rightManager;
		}
		
		public double getChildKStarMultiSequenceBounds(RCTuple queryConf) {
			// TODO Auto-generated method stub
			return 0;
		}

		public double getChildKStarSingleSequenceSums(RCTuple queryConf) {
			// TODO Auto-generated method stub
			return 0;
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
	
	private class KStarHeap extends PriorityQueue<ScoredAssignment>
	{
		public double getHeapSum()
		{
			return -1;
		}
	}
	
	private class RightConfManager 
	{
		private Subproblem rightSubproblem;
		private SubproblemSeqEnumerator rightSubproblemEnum;
		private Map<String, RightConf> rightConfMap = new HashMap<>();
		private Map<String, List<ScoredAssignment>> rightConfLists = new HashMap<>();

		public RightConfManager(SubproblemSeqEnumerator rightSideEnum, Subproblem rightProblem)
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
		private SubproblemSeqEnumerator rightSubproblem;
		
		public RightConf(RCTuple conf, List<ScoredAssignment> rightConfs){
			queryAssignment = conf;
			confList = rightConfs;
			checkList(confList);
		}

		public void updateConf (SubproblemSeqEnumerator rightSubproblem) {
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
	
	private class ScoredSeq implements Comparable<ScoredSeq>
	{
		
		private double cumulativeScore = 0;
		private double selfScore = 0;
		private double leftScore = 0;
		private double rightScore = 0;
		private boolean debugOutput = false;

		@Override
		public int compareTo (ScoredSeq arg0) 
		{
			return Double.compare(cumulativeScore, arg0.cumulativeScore);
		}
		
	}
	
	private class LazyHeap extends PriorityQueue<ScoredAssignment>
	{
		private boolean dirty;
		private ScoredAssignment cleanAssignment = null;
		private SubproblemSeqEnumerator leftSubproblem;
		private RCTuple queryAssignment;
		private double initialRightEnergy = 0;
		
		public LazyHeap(RCTuple queryConf, SubproblemSeqEnumerator leftChild, double defaultEnergy)
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
		private Set<ScoredAssignment> confSet = new TreeSet<>();
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
