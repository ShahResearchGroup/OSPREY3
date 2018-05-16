package edu.duke.cs.osprey.sparse.sequence;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.TreeSet;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.sparse.ConformationProcessor;
import edu.duke.cs.osprey.sparse.PartialConformationEnergyFunction;
import edu.duke.cs.osprey.sparse.Subproblem;


public class SubproblemSeqEnumerator implements ConformationProcessor {
	
	private PartialConformationEnergyFunction energyFunction;
	private ArrayList<Map<String, PriorityQueue<UpperBoundedSeq>>> lambdaHeaps;
	private ArrayList<Map<String, PriorityQueue<MultiSequenceBound>>> lambdaMSBoundHeaps;
	private Map<RCTuple, Map<Sequence, UpperBoundedSeq>> partialSparseScores = new HashMap<>();
	private Map<String, Set<RCTuple>> MSequenceConformations = new HashMap<>();
	private Map<String, Set<RCTuple>> MULambdaConformationsForMULambdaSequence = new HashMap<>();
	private ArrayList<MultiSequenceBound> templateHeaps;
	private Subproblem sourceProblem;
	private SubproblemSeqEnumerator leftSubproblem;
	private SubproblemSeqEnumerator rightSubproblem;
	private ChildConfManager childConfs;
	private static Sequence emptyConf = new Sequence();
	public static boolean debugOutput = true;
	
	public SubproblemSeqEnumerator (Subproblem subproblem, PartialConformationEnergyFunction eFunc) {
		subproblem.addConformationProcessor(this);
		addChildProcessors(subproblem, eFunc);

		energyFunction = eFunc;
		sourceProblem = subproblem;
		int subproblemLocalConfs = subproblem.getTotalMConformations().intValue();
		lambdaHeaps = new ArrayList<>();
		templateHeaps = new ArrayList<>();

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
	
	private void checkHeap(PriorityQueue<UpperBoundedSeq> output)
	{
		for(UpperBoundedSeq seq : output)
		{
			if(seq.assignment == null)
			{
				debugPrint("Null assignment. Weird...");
			}
			if(!sourceProblem.isValidSequence(seq.assignment))
			{
				debugPrint("Invalid conf in heap.");
			}
		}
	}
	
	private PriorityQueue<UpperBoundedSeq> getHeap (RCTuple queryAssignment) {
		return getHeap(sourceProblem.getSequenceForConf(queryAssignment));
	}
	
	private PriorityQueue<UpperBoundedSeq> getHeap (Sequence queryAssignment) {
		int lambdaHeapIndex = sourceProblem.mapSubproblemSeqToIndex(queryAssignment);
		String RCTupleKey = queryAssignment.toString();
		//TODO: Add a check to see if the queryAssignment belongs to a template heap.
		while(lambdaHeaps.size() <= lambdaHeapIndex)
			lambdaHeaps.add(null);
		if(lambdaHeaps.get(lambdaHeapIndex) == null)
		{
			Map<String, PriorityQueue<UpperBoundedSeq>> heapMap = new HashMap<>();
			lambdaHeaps.set(lambdaHeapIndex, heapMap);
		}
		if(!lambdaHeaps.get(lambdaHeapIndex).containsKey(RCTupleKey))
			initializeHeap(queryAssignment, RCTupleKey);
		PriorityQueue<UpperBoundedSeq> output = lambdaHeaps.get(lambdaHeapIndex).get(RCTupleKey);
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
	
	/**
	 * Return the template heap of multisequence bounds that corresponds to the 
	 * query MULambdaAssignment
	 * @param MULambdaAssignment
	 * @return
	 */
	private MultiSequenceBound getTemplateHeap(Sequence MULambdaAssignment)
	{
		int templateHeapIndex = sourceProblem.mapSubproblemSeqToIndex(MULambdaAssignment);
		while(templateHeaps.size() <= templateHeapIndex)
		{
			templateHeaps.add(null);
		}
		if(templateHeaps.get(templateHeapIndex) == null)
		{
			assert(sourceProblem.isMULambdaAssignment(MULambdaAssignment));
			templateHeaps.set(templateHeapIndex, new MultiSequenceBound());
		}
		return templateHeaps.get(templateHeapIndex);
	}
	
	private void initializeHeap (Sequence queryAssignment, String RCTupleKey) {
		int lambdaHeapIndex = sourceProblem.mapSubproblemSeqToIndex(queryAssignment);
		Sequence MTuple = sourceProblem.extractSubproblemMAssignment(queryAssignment);
		// Heap = getTemplateHeap(MSequence)
		// copy = copy(Heap)
		// HeapMap.put(RCTupleKey, copy) 
	}

	public Sequence nextBestSequence(Sequence queryAssignment)
	{
		
		/**
		 * Algorithm: Select sequence which contributes the most to a given
		 * M-sequence's MS bound.
		 * 1. Map input sequence to current node's M sequence.
		 * 2. Map M sequence to MS bound heap, which is sorted by
		 * 		sequence contribution to MS bound score.
		 * 3. Pop top lambda sequence off of MS bound heap
		 * 4. Query child nodes for next best lambda sequence
		 * 5. Update next best lambda sequence score for top lambda sequence,
		 * 6. Add/Update next best lambda sequence to the lambda sequence heaps for each M-conformation
		 * 7. Compute new contribution of lambda sequence to MS bound for M sequence
		 * 8. Reinsert lambda sequence into MS bound heap.
		 */
		
		PriorityQueue<UpperBoundedSeq> lambdaHeap = getHeap(queryAssignment);
		debugPrint("===================== Start "+sourceProblem+" =============================================");
		debugPrint("Beginning recursion with "+queryAssignment+" at \n"+sourceProblem);
		printHeap(lambdaHeap);
		if(lambdaHeap.size() < 1)
		{
			System.err.println("Should not have polled, empty heap at \n"+sourceProblem);
		}
		UpperBoundedSeq previousHeapRoot = lambdaHeap.poll();
		double sanityCheckEnergy = previousHeapRoot.score();
		Sequence curBestAssignment = previousHeapRoot.assignment;
		Sequence outputAssignment = previousHeapRoot.assignment.copy();

		outputAssignment = outputAssignment.combineSequence(queryAssignment);
		Sequence combinedqueryAssignment = queryAssignment.combineSequence(curBestAssignment);

		if(childConfs != null)
		{
			debugPrint("Processing children with "+combinedqueryAssignment);
			if(!childConfs.hasMoreConformations(combinedqueryAssignment))
			{
				System.err.println("Attempting to get new child conf but there are none.");
			}
			assert(childConfs.hasMoreConformations(combinedqueryAssignment));
			Sequence nextBestChildConf = childConfs.getNextChildAssignment(combinedqueryAssignment);
			outputAssignment = outputAssignment.combineSequence(nextBestChildConf);
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
	
	public Sequence nextBestSequence () {
		return nextBestSequence(emptyConf);
	}
	
	private void debugPrint(Object print)
	{
		if(debugOutput)
			System.out.println(print);
	}
	
	private void printHeap(PriorityQueue<UpperBoundedSeq> heap)
	{
		if(!debugOutput)
			return;
		debugPrint("Heap "+heap.hashCode()+":");
		if(heap instanceof LazyHeap)
			debugPrint("Lazy Heap default Energy: "+((LazyHeap) heap).initialRightEnergy);
		PriorityQueue<UpperBoundedSeq> cleanup = new PriorityQueue<>();
		int maxConfsToPrint = 10;
		int numPrinted = 0;
		while(!heap.isEmpty() && numPrinted < maxConfsToPrint)
		{
			UpperBoundedSeq conf = heap.poll();
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
	 * This function takes a MUlambda assignment and returns the
	 * corresponding heap used to make the multisequence bound.
	 * @param queryConf
	 * @return
	 */
	public UpperBoundedSeq getMultiSequenceBound(RCTuple queryConf)
	{
		return null;
	}
	
	public double getSparseKStarScoreForSequence(RCTuple querySeq)
	{
		Sequence sequence = sourceProblem.getSequenceForConf(querySeq);
		double sum = 0;
		for(RCTuple lambdaConf : getLambdaConfsForMULambdaSeq(sequence))
		{
			RCTuple queryConf = lambdaConf.combineRC(querySeq);
			if(leftSubproblem != null)
				sum += leftSubproblem.getPartialKStarBound(queryConf);		
			if(rightSubproblem != null)
				sum += rightSubproblem.getPartialKStarBound(queryConf);
		}
		return sum;
	}
	
	private Set<RCTuple> getMConfsForMSeq (Sequence MSeq) {
		String MSeqString = MSeq.toString();
		if(!MSequenceConformations.containsKey(MSeqString))
			MSequenceConformations.put(MSeqString, new HashSet<>());
		return MSequenceConformations.get(MSeqString);
	}
	
	private Set<RCTuple> getLambdaConfsForMULambdaSeqFromConf (RCTuple conf) {
		assert(sourceProblem.isValidConf(conf));
		conf = sourceProblem.extractSubproblemAssignment(conf);
		Sequence sequence = sourceProblem.getSequenceForConf(conf);
		if(!MSequenceConformations.containsKey(sequence))
		{
			System.err.println("This sequence wasn't queried through an existing" 
					+" sequence, or an RCTuple. Invalid use.");
			throw new NullPointerException();
		}
		return getLambdaConfsForMULambdaSeq(sequence);
	}
	
	private Set<RCTuple> getLambdaConfsForMULambdaSeq (Sequence sequence) {
		String sequenceString = sequence.toString();
		if(!MULambdaConformationsForMULambdaSequence.containsKey(sequenceString))
			MULambdaConformationsForMULambdaSequence.put(sequenceString, new HashSet<>());
		return MULambdaConformationsForMULambdaSequence.get(sequenceString);
	}
	

	/***
	 * This function takes a complete sequence and returns
	 * the corresponding partial Sparse K* score.
	 * @param queryConf
	 * @return
	 */
	public double getPartialKStarScore(RCTuple queryConf, Sequence seq)
	{
		/*
		 * 1. Extract M conformation.
		 * 2. For each lambda conformation in the provided sequence:
		 * 		a. Compute the energy of lambda | M
		 * 		b. Compute the left subtree's energy with respect to M U lambda for sequence
		 * 		c. Compute the right subtree's energy with respect to M U lambda for sequence
		 * 		d. Add values from a, b, and c into one sum, and add that to the ExponentialSum 
		 * 			for all lambda conformations.
		 * 3. Return the resulting ExponentialSum.
		 */
		RCTuple MULambdaConf = sourceProblem.extractSubproblemAssignment(queryConf);
		RCTuple MConf = sourceProblem.extractSubproblemMAssignment(queryConf);
		RCTuple queryLambdaConf = sourceProblem.extractSubproblemLambdaAssignment(queryConf);
		ExponentialSum partialBound = new ExponentialSum();
		for(RCTuple lambdaConf : getLambdaConfsForMULambdaSeqFromConf(MULambdaConf))
		{
			RCTuple queryConfUlambdaConf = queryConf.combineRC(lambdaConf);
			double leftScore = 0;
			double rightScore = 0;
			if(leftSubproblem != null)
				leftScore = leftSubproblem.getPartialKStarScore(queryConfUlambdaConf, seq);
			if(rightSubproblem != null)
				rightScore = rightSubproblem.getPartialKStarScore(queryConfUlambdaConf, seq);
			double score = energyFunction.computePartialEnergyGivenPriorConformation(MConf, lambdaConf);
			partialBound.add(leftScore+rightScore+score);
		}
		return partialBound.evaluate();
	}
	
	public double getPartialKStarBound(RCTuple queryConf)
	{
		/*
		 * 1. Extract M conformation.
		 * 2. For each lambda conformation in the provided sequence s:
		 * 		a. Compute the energy of lambda | M
		 * 		b. Compute the left optimal sequence with respect to M U lambda for s
		 * 		c. Compute the right optimal sequence with respect to M U lambda for s
		 * 		d. Add values from a, b, and c into one sum, and add that to the ExponentialSum 
		 * 			for all lambda conformations.
		 * 3. Return the resulting ExponentialSum.
		 */
		RCTuple MULambdaConf = sourceProblem.extractSubproblemAssignment(queryConf);
		RCTuple MConf = sourceProblem.extractSubproblemMAssignment(queryConf);
		RCTuple queryLambdaConf = sourceProblem.extractSubproblemLambdaAssignment(queryConf);
		ExponentialSum partialBound = new ExponentialSum();
		for(RCTuple lambdaConf : getLambdaConfsForMULambdaSeqFromConf(MULambdaConf))
		{
			RCTuple queryConfUlambdaConf = queryConf.combineRC(lambdaConf);
			double leftScore = 0;
			double rightScore = 0;
			if(leftSubproblem != null)
				leftScore = leftSubproblem.getPartialKStarBound(queryConfUlambdaConf);
			if(rightSubproblem != null)
				rightScore = rightSubproblem.getPartialKStarBound(queryConfUlambdaConf);
			double score = energyFunction.computePartialEnergyGivenPriorConformation(MConf, lambdaConf);
			partialBound.add(leftScore+rightScore+score);
		}
		return partialBound.evaluate();
	}

	@Override
	public void processConformation (RCTuple conformation) {
		/*
		 * Conformation processing pseudocode:
		 * 1. Update MULambda sequence sum
		 * 2. Map conformation to M-conformation lambda-sequence heaps
		 * 3. Update M-conformation lambda-sequence heap. 
		 */
		
		// Step 1: Separate out and score the conformation
		Sequence seq = sourceProblem.getSequenceForConf(conformation);
		System.out.println("Processing "+conformation);
		
		RCTuple MConf = sourceProblem.extractSubproblemMAssignment(conformation);
		Sequence MSeq = sourceProblem.getSequenceForConf(MConf);
		
		RCTuple lambdaConf = sourceProblem.extractSubproblemLambdaAssignment(conformation);
		Sequence lambdaSeq = sourceProblem.getSequenceForConf(lambdaConf);

		double selfScore = energyFunction.computePartialEnergyGivenPriorConformation(MConf, lambdaConf);
		
		double childScores = peekNextBestMultiSequenceBound(conformation);
		updateBoundForSequence(MConf, lambdaSeq);
		
		PriorityQueue<MultiSequenceBound> heap = getMSBoundHeap(MSeq);
		double contribution = computeContributionForSequence(MConf, lambdaSeq);
		heap.peek().addScore(MConf, lambdaSeq, contribution);

		//Bookkeeping: store MConfs.
		Set<RCTuple> MConfs = getMConfsForMSeq(MSeq);
		MConfs.add(MConf.copy());
		
		//Bookkeeping: store necessary parts of multisequence bound for MConf U lambdaSequence.
		
		//Bookeeping: If leaf, precompute partial K* score of the leaf.
		
		// Update the MULambda sequence multisequence bound
		Map<RCTuple, UpperBoundedSeq> lambdaSequenceMSBounds = getMSBoundsForMConf(MConf);
		double score = selfScore;
		if(leftSubproblem != null)
			score+=leftSubproblem.peekNextBestBoundForConf(conformation);
		if(rightSubproblem != null)
			score+=rightSubproblem.peekNextBestBoundForConf(conformation);
		
		
		if(leftSubproblem == null && rightSubproblem == null)
		{
			/* Process partial score per conformation */
		}
	}
	
	private MultiSequenceBound getTemplateMSBoundHeap(Sequence MULambdaSeq)
	{
		int templateHeapIndex = sourceProblem.mapSubproblemSeqToIndex(MULambdaSeq);
		while(templateHeaps.size() <= templateHeapIndex)
		{
			templateHeaps.add(null);
		}
		if(templateHeaps.get(templateHeapIndex) == null)
		{
			assert(sourceProblem.isMULambdaAssignment(MULambdaSeq));
			templateHeaps.set(templateHeapIndex, new MultiSequenceBound());
		}
		return templateHeaps.get(templateHeapIndex);
	}
	
	private PriorityQueue<MultiSequenceBound> getMSBoundHeap (Sequence MULambdaSeq) {
		int lambdaHeapIndex = sourceProblem.mapSubproblemSeqToIndex(MULambdaSeq);
		String RCTupleKey = MULambdaSeq.toString();
		//TODO: Add a check to see if the queryAssignment belongs to a template heap.
		while(lambdaHeaps.size() <= lambdaHeapIndex)
			lambdaHeaps.add(null);
		if(lambdaHeaps.get(lambdaHeapIndex) == null)
		{
			Map<String, PriorityQueue<MultiSequenceBound>> heapMap = new HashMap<>();
			lambdaMSBoundHeaps.set(lambdaHeapIndex, heapMap);
		}
		if(!lambdaMSBoundHeaps.get(lambdaHeapIndex).containsKey(RCTupleKey))
			initializeHeap(MULambdaSeq, RCTupleKey);
		PriorityQueue<MultiSequenceBound> output = lambdaMSBoundHeaps.get(lambdaHeapIndex).get(RCTupleKey);
		if(output == null)
		{
			System.err.println("Heap not found.");
		}
		return output;
	}

	private double computeContributionForSequence (RCTuple mConf, Sequence lambdaSeq) {
		// TODO Auto-generated method stub
		return 0;
	}

	private void updateBoundForSequence (RCTuple mConf, Sequence lambdaSeq) {
		// TODO Auto-generated method stub
		
	}

	private double peekNextBestMultiSequenceBound (RCTuple conformation) {
		double output = 0;
		if(leftSubproblem!=null)
			output+=leftSubproblem.peekNextBestBoundForConf(conformation);
		if(rightSubproblem!=null)
			output+=rightSubproblem.peekNextBestBoundForConf(conformation);
		return output;
	}

	private Map<RCTuple, UpperBoundedSeq> getMSBoundsForMConf (RCTuple mConf) {
		// TODO Auto-generated method stub
		return null;
	}

	private boolean doneProcessing (Sequence mSeq) {
		// TODO Auto-generated method stub
		return false;
	}

	private UpperBoundedSeq getMultiSequenceBoundForMConformation (RCTuple MConf, RCTuple lambdaConf) {
		Sequence lambdaConfSequence = sourceProblem.getSequenceForConf(lambdaConf);
		if(!partialSparseScores.containsKey(MConf))
			partialSparseScores.put(MConf, new HashMap<>());
		Map<Sequence, UpperBoundedSeq> partialSparseScoresForConf = partialSparseScores.get(MConf);
		if(!partialSparseScoresForConf.containsKey(MConf))
			partialSparseScoresForConf.put(lambdaConfSequence, new UpperBoundedSeq(lambdaConfSequence));
		return partialSparseScoresForConf.get(lambdaConfSequence);
	}
	

	

	public double nextBestEnergy()
	{
		return nextBestEnergy(emptyConf);
	}
	
	public double nextBestEnergy(Sequence emptyConf2)
	{
		PriorityQueue<UpperBoundedSeq> heap = getHeap(emptyConf2);
		if(heap.size() < 1)
		{
			debugPrint("No confs at "+sourceProblem);
			getHeap(emptyConf2);
		}
		debugPrint("Getting next best energy for "+emptyConf2+" at "+sourceProblem);
		printHeap(heap);
		UpperBoundedSeq bestConf = heap.peek();
		double bestScore = bestConf.score();

		debugPrint("Next best energy is "+heap.peek().score());
		return bestScore;
	}
	
	private double peekNextBestBoundForConf(RCTuple queryAssignment)
	{
		RCTuple MConf = sourceProblem.extractSubproblemMAssignment(queryAssignment);
		Sequence MSeq = sourceProblem.getSequenceForConf(MConf);
		PriorityQueue<UpperBoundedSeq> sequenceHeap = getHeap(MConf);
		return sequenceHeap.peek().score();
	}

	public boolean hasMoreConformations()
	{
		return hasMoreConformations(emptyConf);
	}
	

	public boolean hasMoreConformations(Sequence rightMAssignment)
	{
		if(getHeap(rightMAssignment).size() < 1)
			return false;
		if(childConfs!= null)
		{
			Sequence topHeapConf = getHeap(rightMAssignment).peek().assignment.copy();
			topHeapConf = topHeapConf.combineSequence(rightMAssignment);
			
			if(!childConfs.hasMoreConformations(topHeapConf) && getHeap(rightMAssignment).size() > 0)
			{
				debugPrint("This will fail. No child confs, but still have lambda confs.");
			}
			
		}

		return getHeap(rightMAssignment).size() > 0;
	}
	

	public Sequence peekNextBestConformation (Sequence queryConf) {
		PriorityQueue<UpperBoundedSeq> heap = getHeap(queryConf);
		if(heap.size() < 1)
		{
			debugPrint("No confs...");
			getHeap(queryConf);
		}
		UpperBoundedSeq bestConf = heap.peek();
		return bestConf.assignment.copy();
	}
	
	/*** This class will encapsulate ALL of the crazy logic that 
	 * goes into maintaining the heaps required for sparse enumeration.
	 * @author Jon
	 *
	 */
	
	private class ChildConfManager
	{
		private Map<String, PriorityQueue<UpperBoundedSeq>> leftHeapMap = new HashMap<>();;
		private RightConfManager rightConfs;
		private SubproblemSeqEnumerator leftSubproblem;
		
		public ChildConfManager(SubproblemSeqEnumerator leftEnumerator, RightConfManager rightManager)
		{
			leftSubproblem = leftEnumerator;
			rightConfs = rightManager;
		}
		
		public Sequence getNextChildAssignment(Sequence combinedqueryAssignment)
		{
			Sequence nextBestChildConf = null;
			if(rightConfs != null)
			{
				debugPrint("Two children found, recursing via secondary heap...");
				PriorityQueue<UpperBoundedSeq> leftTrackingHeap = getLeftHeapMap(combinedqueryAssignment);
				printHeap(leftTrackingHeap);

				if(leftTrackingHeap.size() < 1)
				{
					System.err.println("Tracking heap is empty??");
				}
				UpperBoundedSeq nextBestChildAssignment = leftTrackingHeap.poll();
				debugPrint("Polled "+nextBestChildAssignment);
				if(nextBestChildAssignment == null)
				{
					System.err.println("Null best child assignment...");
				}
				
				Sequence bestChildConf = nextBestChildAssignment.assignment;
				debugPrint("Best left conf is "+bestChildConf+", maching to right side...");

				Sequence nextBestRightConf = rightConfs.getNextBestRightConf(bestChildConf);
				debugPrint("RightConf returned is "+nextBestRightConf);
				updateChildAssignment(combinedqueryAssignment,nextBestChildAssignment);
				nextBestChildConf = bestChildConf.combineSequence(nextBestRightConf);

				debugPrint("Left heap final status: ");
				printHeap(leftTrackingHeap);
			}
			else
			{
				debugPrint("Recursing to left child...");
				nextBestChildConf = leftSubproblem.nextBestSequence(combinedqueryAssignment);
				nextBestChildConf = nextBestChildConf.combineSequence(combinedqueryAssignment);
			}
			assert(nextBestChildConf != null);
			return nextBestChildConf;

		}

		public double nextBestEnergy()
		{
			return nextBestEnergy(emptyConf);
		}
		
		public double nextBestEnergy(Sequence combinedqueryAssignment)
		{
			PriorityQueue<UpperBoundedSeq> heap = getHeap(combinedqueryAssignment);
			if(heap.size() < 1)
			{
				debugPrint("No confs at "+sourceProblem);
				getHeap(combinedqueryAssignment);
			}
			debugPrint("Getting next best energy for "+combinedqueryAssignment+" at "+sourceProblem);
			printHeap(heap);
			UpperBoundedSeq bestConf = heap.peek();
			double bestScore = bestConf.score();

			debugPrint("Next best energy is "+heap.peek().score());
			return bestScore;
		}
		
		
		private PriorityQueue<UpperBoundedSeq> getLeftHeapMap(Sequence combinedqueryAssignment)
		{
			if(!leftHeapMap.containsKey(combinedqueryAssignment.toString()))
			{
				Sequence leftMAssignment = sourceProblem.leftSubproblem.extractSubproblemMAssignment(combinedqueryAssignment);
				double initialRightenergy = rightConfs.peekNextBestEnergy(combinedqueryAssignment);
				PriorityQueue<UpperBoundedSeq> newHeap = new LazyHeap(combinedqueryAssignment, leftSubproblem, initialRightenergy);
				if(newHeap.size() < 1)
				{
					System.err.println("Empty new heap. Should be full of template nodes.");
					leftSubproblem.getHeap(leftMAssignment);
				}
				leftHeapMap.put(combinedqueryAssignment.toString(), newHeap);
			}
			
			return leftHeapMap.get(combinedqueryAssignment.toString());
		}

		public boolean hasMoreConformations(Sequence combinedqueryAssignment)
		{
			if(rightConfs!=null)
			{
				return getLeftHeapMap(combinedqueryAssignment).size() > 0;
			}
			else
			{
				return leftSubproblem.hasMoreConformations(combinedqueryAssignment);
			}
		}

		private void updateChildAssignment (Sequence combinedqueryAssignment, UpperBoundedSeq nextBestChildAssignment) {
			Sequence nextBestChildQueryConf = nextBestChildAssignment.assignment;
			if(rightConfs.hasMoreConformations(nextBestChildQueryConf))
			{
				debugPrint("More right confs for left assignment "+nextBestChildQueryConf+" at:\n"+sourceProblem);
				nextBestChildAssignment.updateLeftScore(rightConfs.peekNextBestEnergy(nextBestChildQueryConf));
				leftHeapMap.get(combinedqueryAssignment.toString()).add(nextBestChildAssignment);
			}
			else 
			{
				debugPrint("No more right confs for left assignment "+nextBestChildQueryConf+" at:\n"+sourceProblem);
				/*
				if(leftSubproblem.hasMoreConformations(queryAssignment))
				{
					double nextLeftEnergy = leftSubproblem.nextBestEnergy(queryAssignment);
					RCTuple nextLeftConf = leftSubproblem.nextBestConformation();
					UpperBoundedSeq nextLeftTrackingConf = new UpperBoundedSeq(nextLeftConf, nextLeftEnergy,0,0);
					leftHeapMap.get(queryAssignment.toString()).add(nextLeftTrackingConf);
				}
				*/
			}
			
		}
	
		
		
	}
	
	private class RightConfManager 
	{
		private Subproblem rightSubproblem;
		private SubproblemSeqEnumerator rightSubproblemEnum;
		private Map<String, RightConf> rightConfMap = new HashMap<>();
		private Map<String, List<UpperBoundedSeq>> rightConfLists = new HashMap<>();

		public RightConfManager(SubproblemSeqEnumerator rightSideEnum, Subproblem rightProblem)
		{
			rightSubproblemEnum = rightSideEnum;
			rightSubproblem = rightProblem;
		}
		
		private List<UpperBoundedSeq> getRightConfList (Sequence combinedqueryAssignment)
		{
			if(!rightConfLists.containsKey(combinedqueryAssignment.toString()))
			{
				Sequence templateConf = rightSubproblem.extractSubproblemMAssignment(combinedqueryAssignment);
				if(!rightConfLists.containsKey(templateConf.toString()))
				{
					debugPrint("Creating new conf list for "+combinedqueryAssignment);
					List<UpperBoundedSeq> newConfList = new UniqueConfList();
					double nextRightConfE = rightSubproblemEnum.nextBestEnergy(templateConf);
					debugPrint("Polling new conformation from right side with "+templateConf+"...");
					Sequence nextRightConf = rightSubproblemEnum.nextBestSequence(templateConf);
					Sequence rightPart = rightSubproblemEnum.sourceProblem.extractSubproblemLAssignment(nextRightConf);
					newConfList.add(new UpperBoundedSeq(rightPart, nextRightConfE, 0, 0));
					assert(newConfList.size() == 1);
					rightConfLists.put(templateConf.toString(), newConfList);
					debugPrint("Created new conf list for "+combinedqueryAssignment+": "+newConfList+";"+newConfList.hashCode());
				}
				rightConfLists.put(combinedqueryAssignment.toString(), rightConfLists.get(templateConf.toString()));
			}
			return rightConfLists.get(combinedqueryAssignment.toString());
		}
		

		private RightConf getRightConf(Sequence combinedqueryAssignment)
		{
			if(!rightConfMap.containsKey(combinedqueryAssignment.toString()))
			{
				debugPrint("Creating new RightConf for "+combinedqueryAssignment);
				rightConfMap.put(combinedqueryAssignment.toString(), new RightConf(
						combinedqueryAssignment, getRightConfList(combinedqueryAssignment)));
			}
			return rightConfMap.get(combinedqueryAssignment.toString());
		}
		
		public double peekNextBestEnergy(Sequence combinedqueryAssignment)
		{
			return getRightConf(combinedqueryAssignment).peekConfEnergy();
		}
		
		public Sequence getNextBestRightConf(Sequence bestChildConf)
		{
			RightConf rightConf = getRightConf(bestChildConf);
			Sequence output = rightConf.pollCurConf();
			rightConf.updateConf(rightSubproblemEnum);
			return output;
		}
		
		public boolean hasMoreConformations(Sequence nextBestChildQueryConf)
		{
			return getRightConf(nextBestChildQueryConf).hasMoreConformations();
		}
	}
	
	private class RightConf
	{
		private Sequence queryAssignment;
		private int confListIndex = 0;
		private List<UpperBoundedSeq> confList;
		private SubproblemSeqEnumerator rightSubproblem;
		
		public RightConf(Sequence combinedqueryAssignment, List<UpperBoundedSeq> rightConfs){
			queryAssignment = combinedqueryAssignment;
			confList = rightConfs;
			checkList(confList);
		}

		public void updateConf (SubproblemSeqEnumerator rightSubproblem) {
			Sequence rightMAssignment = rightSubproblem.sourceProblem.extractSubproblemMAssignment(queryAssignment);
			if(!hasMoreConformations())
			{
				debugPrint("Used up existing "+queryAssignment+", checking to see if there are more...");
				debugPrint("Querying for more right conformations with "+rightMAssignment);
			}
			if(!hasMoreConformations() && rightSubproblem.hasMoreConformations(rightMAssignment))
			{
				debugPrint("Appending new right Conformation for "+queryAssignment);
				double nextRightConfE = rightSubproblem.nextBestEnergy(rightMAssignment);
				Sequence nextRightConf = rightSubproblem.nextBestSequence(rightMAssignment);
				Sequence rightPart = rightSubproblem.sourceProblem.extractSubproblemLAssignment(nextRightConf);
				confList.add(new UpperBoundedSeq(rightPart, nextRightConfE, 0, 0));
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

		public Sequence pollCurConf () {

			checkList(confList);

			if(confListIndex >= confList.size())
			{
				System.err.println("No confs...");
			}
			assert(confListIndex < confList.size());
			debugPrint("Using conf "+confListIndex+" from list "+confList+";"+confList.hashCode());
			Sequence output = confList.get(confListIndex).assignment;
			confListIndex++;
			debugPrint("Polled conf "+output+", index advanced to "+confListIndex);
			return output;
		}
		
		public Sequence peekCurConf () {
			return confList.get(confListIndex).assignment;
		}
		
		private void checkList(List<UpperBoundedSeq> confs)
		{
			if(confs == null)
			{
				System.err.println("No confList. Why do we exist??");
			}
				
			Set<String> confSet = new TreeSet<String>();
			for(UpperBoundedSeq conf : confs)
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
	

		

	
	private class SplitConfScore
	{
		private double leftScore;
		private double rightScore;
		private double selfScore;
		private double sum;
		private boolean sumIsStale;
		
		public double getSum()
		{
			assert(sum == leftScore+rightScore+selfScore);
			return sum;
		}
		
		public void updateLeftScore(double newLeftScore)
		{
			leftScore = newLeftScore;
			sum = selfScore + rightScore + leftScore;
		}
		
		public void updateRightScore(double newRightScore)
		{
			rightScore = newRightScore;
			sum = selfScore + rightScore + leftScore;
		}
	}
	
	
	private class LazyHeap extends PriorityQueue<UpperBoundedSeq>
	{
		private boolean dirty;
		private UpperBoundedSeq cleanAssignment = null;
		private SubproblemSeqEnumerator leftSubproblem;
		private Sequence queryAssignment;
		private double initialRightEnergy = 0;
		
		public LazyHeap(Sequence combinedqueryAssignment, SubproblemSeqEnumerator leftChild, double defaultEnergy)
		{
			if(defaultEnergy > 0)
			{
				debugPrint("Every one of these conformations will have a positive right energy. check this.");
			}
			queryAssignment = combinedqueryAssignment;
			leftSubproblem = leftChild;
			initialRightEnergy = defaultEnergy;
			addNewLeftConf();
			peek();
		}
		

		
		private void addNewLeftConf()
		{
			double nextBestEnergy = leftSubproblem.nextBestEnergy(queryAssignment);
			Sequence nextBestLeftSeq = leftSubproblem.nextBestSequence(queryAssignment);
			cleanAssignment = new UpperBoundedSeq(nextBestLeftSeq);

			debugPrint("Creating new left conf "+cleanAssignment);
			add(cleanAssignment);
		}
		
		public UpperBoundedSeq poll()
		{
			cleanHeap();
			UpperBoundedSeq nextBestAssignment = super.poll();
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
		
		public UpperBoundedSeq peek()
		{
			cleanHeap();
			UpperBoundedSeq output = super.peek();
			return output;
		}
		
		public String toString()
		{
			return "Default Right Energy: "+initialRightEnergy+", heap: "+super.toString();
		}
		
	}
	
	


	
	private class UniqueConfList extends LinkedList<UpperBoundedSeq>
	{
		private Set<UpperBoundedSeq> confSet = new TreeSet<UpperBoundedSeq>();
		@Override
		public boolean add(UpperBoundedSeq conf)
		{
			if(confSet.contains(conf))
				System.err.println("Dupe conf.");
			assert(!confSet.contains(conf));
			return super.add(conf);
		}
		
	}

}
