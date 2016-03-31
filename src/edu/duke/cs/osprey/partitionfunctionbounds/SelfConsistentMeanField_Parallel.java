/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds;

import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.tools.ExpFunction;
import java.math.BigDecimal;
import java.util.ArrayList;

/**
 *
 * @author hmn5
 */
public class SelfConsistentMeanField_Parallel implements InferenceCalculator {

    BigDecimal partitionFunction;

    //Nodelist defines the graphical model, since nodes keep track of their neighbors
    ArrayList<MRFNode> nodeList;
    //Emat defines the potential functions for the Probabilistic Graphical Model
    EnergyMatrix emat;
    //interactionGraph determines if two nodes are neighbors
    boolean[][] interactionGraph;

    //threshold for convergence (change in beliefs must be less than threshold)
    double threshold = 1e-8;
    final int minNumberIterations = 1000;
    final int maxNumberIterations = 10000;
    
    double lambda = 0.70;
    
    double scmfTemp = PoissonBoltzmannEnergy.constRT;
    public ExpFunction ef = new ExpFunction();

    double freeEnergy;
    boolean hasCalculatedFreeEnergy;

    public SelfConsistentMeanField_Parallel(MarkovRandomField mrf) {
        this.nodeList = mrf.nodeList;
        this.emat = mrf.emat;
        this.interactionGraph = mrf.interactionGraph;
    }

    public void run() {
        int iter = 0;
        boolean hasConverged = false;
        initializeBeliefs();
        //set the temperature to be high
        this.scmfTemp = this.scmfTemp;
        //while we haven't converged and are below maxNumberIteration;
        while (!hasConverged && (iter < this.maxNumberIterations)) {
            //update all beliefs

            hasConverged = updateAllBeliefs();
            
            if (iter == 300){
                this.freeEnergy  = calcFreeEnergy();
                this.scmfTemp = this.scmfTemp*1000;
                lambda = 1.0;
            }
            //update old Belief to new belief;
            changeBeliefs();
            checkBeliefs();
            //lower the temperature
            lowerTemperature();
            //Make sure we are at least beyond the minimum number of iterations
            hasConverged = (iter > this.minNumberIterations);
            iter++;
        }
        this.freeEnergy = Math.min(this.freeEnergy,calcFreeEnergy());
        printMaxLabel();
        resetTemperature();
    }

    //Initialize the beliefs of each label in a node to be uniform among all possible
    //labels that belong to the node
    private void initializeBeliefs() {
        for (MRFNode node : this.nodeList) {
            int numLabels = node.labelList.size();
            for (MRFLabel label : node.labelList) {
                label.currentBelief = 1.0 / (double) numLabels;
                label.oldBelief = 1.0 / (double) numLabels;
            }
        }
    }
    private void changeBeliefs(){
        for (MRFNode node : this.nodeList) {
            for (MRFLabel label : node.labelList) {
                label.oldBelief = lambda*label.currentBelief + (1-lambda)*label.oldBelief;
            }
        }        
    }

    //Returns the MeanFieldEnergy for a label
    private double getMeanFieldEnergy(MRFNode node, MRFLabel label) {
        double meanFieldE = 0.0;
        for (MRFNode neighbor : node.neighborList) {
            //get the average energy between neighbor and our label
            double averageE = getAverageEnergy(node, label, neighbor);
            meanFieldE += averageE;
        }
        return meanFieldE;
    }

    //Get average energy of label1 from node1, with respect to all labels of node2
    //this is used as a subroutine in getMeanField()
    private double getAverageEnergy(MRFNode node1, MRFLabel label1, MRFNode node2) {
        double averageE = 0.0;
        for (MRFLabel label2 : node2.labelList) {
            double E = this.emat.getPairwise(node1.nodeNum, label1.labelNum, node2.nodeNum, label2.labelNum);
            averageE += E * label2.oldBelief;
        }
        return averageE;
    }

    //updates all beliefs and returns true if 
    //max difference in beliefs is less than threshold
    private boolean updateAllBeliefs() {
        boolean allNodesConverged = true;
        for (MRFNode node : this.nodeList) {
            //We converge if and only if every node has converged
            boolean nodeConverged = updateNodeBeliefs(node);
            if(!(allNodesConverged && nodeConverged)){
                allNodesConverged = false;
            }
        }
        return allNodesConverged;
    }

    //update all of the beliefs belonging to node
    //returns true if max difference in beleifs is less than threshold
    private boolean updateNodeBeliefs(MRFNode node) {
        //create a normalizing constant to normalize beliefs
        BigDecimal partFunction = new BigDecimal(0.0);
        //keep track of difference between beliefs for convergence;
        double maxEpsilon = 0.0;
        //keep track of unnormalized Beliefs so we don't need to recompute
        ArrayList<BigDecimal> unNormalizedBeliefs = new ArrayList<>();
        ArrayList<Double> logUnnormalizedBeliefs = new ArrayList<>();
        //We keep track of best energy to avoid numerical innacuracy
        double bestE = Double.NEGATIVE_INFINITY;
        //iterate over labels to get partition function value
        for (int labelIndex=0; labelIndex< node.labelList.size(); labelIndex++) {
            MRFLabel label = node.labelList.get(labelIndex);
            double oneBodyE = this.emat.getOneBody(node.nodeNum, label.labelNum);
            double meanFieldE = getMeanFieldEnergy(node, label);
            //unnormalized updateBelief
            double logUnnormalizedBelief = -(oneBodyE + meanFieldE)/scmfTemp;
            logUnnormalizedBeliefs.add(logUnnormalizedBelief);
            if (logUnnormalizedBelief > bestE){
                bestE = logUnnormalizedBelief;
            }
        }
        //TODO: This needs to be optimized better by caching results in first for loop
        for (int labelIndex=0; labelIndex< node.labelList.size(); labelIndex++) {
            MRFLabel label = node.labelList.get(labelIndex);
            double logUnnormalizedBelief = logUnnormalizedBeliefs.get(labelIndex);
            double rescaleE = -bestE +  logUnnormalizedBelief;
            BigDecimal updateBelief = this.ef.exp(rescaleE);
            //update partition function
            partFunction = partFunction.add(updateBelief);
            //store unnormalizedBeliefs
            unNormalizedBeliefs.add(updateBelief);
        }
        //now we update the beliefs using our partFunction (normalizing constant)
        for (int labelIndex = 0; labelIndex < node.labelList.size(); labelIndex++) {
            MRFLabel label = node.labelList.get(labelIndex);
            //keep track of old belief
            double oldBelief = label.oldBelief;
            //if partFunction is 0.0, all labels had really bad energies and we will make them uniform
            if (partFunction.doubleValue() == 0.0) {
                label.currentBelief = 1.0 / node.labelList.size();
            } //otherwise we update the beliefs normally
            else {
                label.currentBelief = unNormalizedBeliefs.get(labelIndex).divide(partFunction, ExpFunction.mc).doubleValue();
            }
            //update maxEpsilon
            if (Math.abs(label.currentBelief - oldBelief) > maxEpsilon) {
                maxEpsilon = Math.abs(label.currentBelief - oldBelief);
            }
        }
        boolean hasConverged = false;
        if (maxEpsilon < threshold) {
            hasConverged = true;
        }
        return hasConverged;
    }

    //lower the temperature, never getting below the true value;
    private void lowerTemperature() {
        this.scmfTemp = this.scmfTemp * 0.98;
        if (this.scmfTemp < PoissonBoltzmannEnergy.constRT) {
            this.scmfTemp = PoissonBoltzmannEnergy.constRT;
        }
    }

    private void resetTemperature() {
        this.scmfTemp = PoissonBoltzmannEnergy.constRT;
    }

    private double getSingleNodeEnthalpy(MRFNode node) {
        double enthalpy = 0.0;
        for (MRFLabel label : node.labelList) {
            double E = this.emat.getOneBody(node.nodeNum, label.labelNum);
            enthalpy += E * label.currentBelief;
        }
        return enthalpy;
    }

    private double getPairwiseNodeEnthalpy(MRFNode node1, MRFNode node2) {
        double enthalpy = 0.0;
        for (MRFLabel label1 : node1.labelList) {
            for (MRFLabel label2 : node2.labelList) {
                double E = emat.getPairwise(node1.nodeNum, label1.labelNum, node2.nodeNum, label2.labelNum);
                enthalpy += E * label1.currentBelief * label2.currentBelief;
            }
        }
        return enthalpy;
    }

    private double getEnthalpy() {
        double enthalpy = 0.0;
        for (int i = 0; i < this.nodeList.size(); i++) {
            MRFNode node1 = nodeList.get(i);
            enthalpy += getSingleNodeEnthalpy(node1);
            for (int j = 0; j < i; j++) {
                MRFNode node2 = nodeList.get(j);
                if (this.interactionGraph[node1.nodeNum][node2.nodeNum]) {
                    enthalpy += getPairwiseNodeEnthalpy(node1, node2);
                }
            }
        }
        return enthalpy;
    }

    private double getSingleNodeEntropy(MRFNode node) {
        double entropy = 0.0;
        for (MRFLabel label : node.labelList) {
            if (label.currentBelief == 0.0) {
                entropy += 0.0;
            } else {
                entropy += (-1.0) * label.currentBelief * Math.log(label.currentBelief);
            }
        }
        return entropy;
    }

    private double getEntropy() {
        double entropy = 0.0;
        for (MRFNode node : this.nodeList) {
            entropy += getSingleNodeEntropy(node);
        }
        return entropy;
    }

    private double calcFreeEnergy() {
        double enthalpy = getEnthalpy();
        double entropy = getEntropy();

        double freeEnergy = enthalpy - this.scmfTemp * entropy;

        return freeEnergy;
    }

    @Override
    public BigDecimal calcPartitionFunction() {
        if (!hasCalculatedFreeEnergy){
            this.freeEnergy = calcFreeEnergy();
        }
        BigDecimal partitionFunction = this.ef.exp(-((freeEnergy+this.emat.getConstTerm()) / this.scmfTemp));
        return partitionFunction;
    }
    
     //Calculates the natural log of the partition function
    public double calcLBLogZ(){
        return -(calcFreeEnergy() + this.emat.getConstTerm())/this.scmfTemp;
    }
    
    //Calculates log_10 of the partition function
    public double calcLBLog10Z(){
        return (Math.log10(Math.E))*calcLBLogZ();
    }
    
    
    public boolean checkBeliefs(){
        boolean beliefsGood = true;
        for (MRFNode node : this.nodeList){
            double beliefSum = 0.0;
            for (MRFLabel label : node.labelList){
                beliefSum += label.currentBelief;
            }
            if (!equals( beliefSum, 1.0, 0.001)){
                beliefsGood = false;
            }
        }
        return beliefsGood;
    }
    
    public void printMaxLabel(){
        String maxLabels = "";
        for(MRFNode node: nodeList){
            double max = 0.0;
            int maxLabelNum = -1;
            for (MRFLabel label : node.labelList){
                if (label.currentBelief > max){
                    max = label.currentBelief;
                    maxLabelNum = label.labelNum;
                }
            }
            maxLabels = maxLabels +maxLabelNum+" ";
        }
        System.out.println(maxLabels);
    }
    
    public boolean equals(double a , double b, double epsilon){
        return a == b ? true : Math.abs(a - b ) < epsilon;
    }
}