package edu.duke.cs.osprey.sparse;

import edu.duke.cs.osprey.confspace.RCTuple;

public class ScoredAssignment implements Comparable<ScoredAssignment>
	{
		private double cumulativeScore = 0;
		private double selfScore = 0;
		private double leftScore = 0;
		private double rightScore = 0;
		private boolean debugOutput = false;
		RCTuple assignment;
		public ScoredAssignment (RCTuple conformation, double selfEnergy,
				double leftBestEnergy, double rightBestEnergy) {
			assignment = conformation;
			if(assignment == null)
			{
				System.err.println("Created a conformation with no actual conformation.");
			}
			debugPrint("Creating conf from "+conformation+" with self energy "
					+selfEnergy+", left energy "+leftBestEnergy+", right energy "+rightBestEnergy);
			selfScore = selfEnergy;
			leftScore = leftBestEnergy;
			if(leftScore == -16.36759805536081)
				debugPrint("where's this coming from?");
			rightScore = rightBestEnergy;
			cumulativeScore = leftScore+rightScore+selfScore;
		}
		
		public double score()
		{
			return cumulativeScore;
		}
		
		public void updateScore(double leftBestEnergy, double rightBestEnergy)
		{
			leftScore = leftBestEnergy;
			rightScore = rightBestEnergy;
			cumulativeScore = leftScore+rightScore+selfScore;
		}
		
		public void updateLeftScore(double leftBestEnergy)
		{
			leftScore = leftBestEnergy;
			cumulativeScore = leftScore+rightScore+selfScore;
		}
		
		public void updateRightScore(double rightBestEnergy)
		{
			rightScore = rightBestEnergy;
			cumulativeScore = leftScore+rightScore+selfScore;
		}
		@Override
		public int compareTo (ScoredAssignment arg0) {
			return Double.compare(cumulativeScore, arg0.cumulativeScore);
		}
		
		public boolean equals(ScoredAssignment other)
		{
			if(other.cumulativeScore != cumulativeScore)
				return false;
			if(other.assignment.size() != assignment.size())
				return false;
			if(!other.assignment.pos.containsAll(assignment.pos))
				return false;
			if(!other.assignment.isSameTuple(assignment))
				return false;
			return true;
		}
		
		public ScoredAssignment copy()
		{
			return new ScoredAssignment(assignment, selfScore, leftScore, rightScore);
		}
		
		public String toString()
		{
			return "{{" + assignment+":"+cumulativeScore+", self:"+selfScore+", left: "+leftScore+", right:"+rightScore+"}}";
		}
		
		private void debugPrint(Object print)
		{
			if(debugOutput)
				System.out.println(print);
		}
		
	}
