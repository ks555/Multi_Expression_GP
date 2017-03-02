package core;
import core.GpRun;
import programElements.Operator;
import programElements.Terminal;
import utils.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import java.math.*;

public class MIndividual extends Individual {

	private static final long serialVersionUID = 7L;
	
	protected  static int numPrograms=2;

	protected static long nextId;

	protected long id;
	protected double k;
	protected double kunseen;
	protected double w;
	protected ArrayList<Individual> programs;
	protected int depth;
	protected double [] ratios;
	protected double [] ratiosK;
	protected double [] ratiosW;
	protected double[][] trainingErrorVectors =  new double[numPrograms][];
	protected double [][] unseenErrorVectors =  new double[numPrograms][];
	protected double [] unseenOutputAverages =  new double[numPrograms];
	protected double [] trainingOutputAverages =  new double[numPrograms];
	protected double [] unseenOutputSD =  new double[numPrograms];
	protected double [] trainingOutputSD =  new double[numPrograms];
	//protected int[] reorderedVectors = new int[numPrograms];
	protected double trainingTheta,unseenTheta;
	protected double reconTrainingError,reconUnseenError;
	protected double minDistance;

	protected int evaluateIndex;
	protected int maximumDepthAchieved;
	protected int depthCalculationIndex;
	protected int printIndex;
	
	protected Individual reconstructed;
	protected MIndividual mp1;
	protected MIndividual mp2;
	protected boolean sizeOverride;
	protected int computedSize;

	public MIndividual() {
		programs = new ArrayList <Individual>();
		id = getNextId();	
		if(numPrograms==2){
			w=2;
		}
	}

	protected static long getNextId() {
		return nextId++;
	}
	
	public void evaluate(Data data) {
		evaluateTrainingTheta(data);
		evaluateUnseenTheta(data);
		double[][] trainingData = data.getTrainingData();
		double[][] unseenData = data.getUnseenData();
		double[] reconstructedTrainingOutput = reconstructTrainingSemantics(trainingData);
		double[] reconstructedUnseenOutput = reconstructUnseenSemantics(unseenData);
		outputRMS();
		reconTrainingError=calculateRMSE(trainingData, reconstructedTrainingOutput);
		reconUnseenError=calculateRMSE(unseenData, reconstructedUnseenOutput);
	}	
	
	public void outputRMS(){
		double sumTrain;
		double sumUnseen;
		double sumTrain2;
		double sumUnseen2;
		double[] outputsTrain;
		double[] outputsUnseen;
		for (int i=0;i<numPrograms;i++){
			sumTrain=0;
			sumUnseen=0;
			sumTrain2=0;
			sumUnseen2=0;
			outputsTrain=getProgram(i).getTrainingDataOutputs();
			outputsUnseen=getProgram(i).getUnseenDataOutputs();
			//calc mean training output
			for (int j = 0; j < outputsTrain.length; j++) {
				
				sumTrain +=outputsTrain[j];
			}
			trainingOutputAverages[i]=sumTrain / outputsTrain.length;

			
			//calc sd of training output
			for (int j = 0; j < outputsTrain.length; j++) {
				
				sumTrain2 += Math.pow(outputsTrain[j]-trainingOutputAverages[i], 2.0);
			}
			trainingOutputSD[i]= Math.sqrt(sumTrain2 / outputsTrain.length);
			
			//calc mean unseen output
			for (int j = 0; j < outputsUnseen.length; j++) {
				
				sumUnseen += outputsUnseen[j];
			}
			unseenOutputAverages[i]=sumUnseen / outputsUnseen.length;
			
			//calc unseen standard deviation
			for (int j = 0; j < outputsUnseen.length; j++) {
				
				sumUnseen2 += Math.pow(outputsUnseen[j]-unseenOutputAverages[i], 2.0);
			}
			unseenOutputSD[i]=Math.sqrt(sumUnseen2 / outputsUnseen.length);
		}
	}
	
	public void evaluateTrainingTheta(Data data) {
		//for each program, add error vector to trainingErrorVectors
		//calc theta using trainingErrorVectors		
		for(int i=0;i<numPrograms;i++){			
			trainingErrorVectors[i]=this.getProgram(i).evaluateErrorVectorOnTrainingData(data);			
			//System.out.print(getProgram(i).getId()+" "+trainingErrorVectors[i][0]+" ");			
		}
		//System.out.println(); 
		evaluateTrainingTheta(trainingErrorVectors);		
	}
	
	public void evaluateUnseenTheta(Data data) {
		//for each program, add error vector to trainingErrorVectors
		//calc theta using trainingErrorVectors			
		for(int i=0;i<numPrograms;i++){
			//calc theta
			unseenErrorVectors[i]=this.getProgram(i).evaluateErrorVectorOnUnseenData(data);
		}	

		evaluateUnseenTheta(unseenErrorVectors);		
	}
	
	public void evaluateTrainingTheta(double[][] errorVector) {
		
		double dotProd;
		double normA;
		double normB;
		double[][] normalized=new double[numPrograms][];
		double tempTheta;
		trainingTheta=360;
		
		if (numPrograms == 2){
			dotProd = dot(errorVector[0], errorVector[1]);			
			normA = magnitude(errorVector[0]); 			
			normB = magnitude(errorVector[1]);
			trainingTheta = Math.acos(dotProd/(normA*normB));
		}
		else if (numPrograms == 3){
			
			//this is for if we are checking all combos of expressions for all planes...
			//if use this, change planeVector to a vector of same length as errorVector and make it a property of the class		

			double[][] expressionOrder=new double[numPrograms][];
			int [] reorderedVectors=new int[numPrograms];
			int index;
			int[] order=new int[numPrograms];
			for(int i=0;i<numPrograms;i++){		
				index=0;
				for (int j=0;j<numPrograms;j++){	
					expressionOrder[numPrograms-1]=errorVector[i];
					order[numPrograms-1]=i;					
				
					if(j!=i){
						expressionOrder[index]=errorVector[j];	
						order[index]=j;						
						index++;					
					}					
				}
				
				
				normalized[0]=scalar(expressionOrder[0],1/magnitude(expressionOrder[0]));
				normalized[1]=scalar(minus(expressionOrder[1],scalar(normalized[0],dot(normalized[0],expressionOrder[1]))),1/magnitude(minus(expressionOrder[1],scalar(normalized[0],dot(normalized[0],expressionOrder[1])))));
				normalized[2]=scalar(minus(minus(expressionOrder[2],scalar(normalized[0],dot(normalized[0],expressionOrder[2]))),scalar(normalized[1],dot(normalized[1],expressionOrder[2]))),1/magnitude(minus(minus(expressionOrder[2],scalar(normalized[0],dot(normalized[0],expressionOrder[2]))),scalar(normalized[1],dot(normalized[1],expressionOrder[2])))));
				tempTheta=Math.asin(dot(expressionOrder[2],normalized[2])/Math.sqrt(dot(expressionOrder[2],expressionOrder[2])));
				if(tempTheta<trainingTheta){
					trainingTheta=tempTheta;
					for (int j=0;j<numPrograms;j++){
						reorderedVectors[j]=order[j];
					}
				}
				
				//System.out.println("Theta "+ tempTheta + " with ");
						
//				for(int j=0;j<numPrograms;j++){	
//					System.out.print(getProgram(order[j]).getId()+" "+expressionOrder[j][0]+" ");			
//				}
//				System.out.println(); 
			}			
			reorderMindividual(reorderedVectors);
		}		
		//equation for angle of w to a plane containing u,v
//		    #V1 = Normalize[v];
//		    #U1 = Normalize[u - (V1.u) V1];
//		    #W1 = Normalize[w - (w.V1) V1 - (w.U1) U1];
//		    #theta3 = asin(dot(w,W1)/sqrt(dot(w,w)))
		trainingTheta =Math.toDegrees(trainingTheta);
		//trainingTheta =Math.abs(Math.toDegrees(trainingTheta));
		//if(trainingTheta>90)trainingTheta=180-Math.abs(trainingTheta);

	}
	
	public void evaluateUnseenTheta(double[][] errorVector) {
		
		double dotProd;
		double normA;
		double normB;
		double[][] normalized=new double[numPrograms][];
		
		if (numPrograms == 2){
			
			dotProd = dot(errorVector[0], errorVector[1]);			
			normA = magnitude(errorVector[0]); 			
			normB = magnitude(errorVector[1]);
			unseenTheta = Math.acos(dotProd/(normA*normB));

		}	
		else if (numPrograms == 3){
			
			normalized[0]=scalar(errorVector[0],1/magnitude(errorVector[0]));
			normalized[1]=scalar(minus(errorVector[1],scalar(normalized[0],dot(normalized[0],errorVector[1]))),1/magnitude(minus(errorVector[1],scalar(normalized[0],dot(normalized[0],errorVector[1])))));
			normalized[2]=scalar(minus(minus(errorVector[2],scalar(normalized[0],dot(normalized[0],errorVector[2]))),
					scalar(normalized[1],dot(normalized[1],errorVector[2]))),1/magnitude(minus(minus(errorVector[2],scalar(normalized[0],dot(normalized[0],errorVector[2]))),
							scalar(normalized[1],dot(normalized[1],errorVector[2])))));
			unseenTheta=Math.asin(dot(errorVector[2],normalized[2])/Math.sqrt(dot(errorVector[2],errorVector[2])));

		}		
		unseenTheta =Math.toDegrees(unseenTheta);
		//unseenTheta =Math.abs(Math.toDegrees(unseenTheta));
		//if(unseenTheta>90)unseenTheta=180-Math.abs(unseenTheta);
	}
		
    public double[] scalar(double[] vectorOne, double scalar) {
        double[] solution =new double[vectorOne.length];
        for (int i = 0; i < vectorOne.length; i++)
            solution[i] = vectorOne[i]*scalar;
        return solution;
      }
      
      public double[] minus(double[] vectorOne, double[] vectorTwo) {
      	 if (vectorOne.length != vectorTwo.length)
               throw new IllegalArgumentException("Dimensions disagree for angle calculation");
          double[] solution =new double[vectorOne.length];
          for (int i = 0; i < vectorOne.length; i++)
              solution[i] = vectorOne[i]-vectorTwo[i];
          return solution;
        }

	public void addProgramAtIndex(Individual program, int index) {
		programs.add(index, program);
	}

	public void removeProgramAtIndex(int index) {
		programs.remove(index);
	}
	public void replaceProgramAtIndex(Individual program, int index) {
		programs.set(index, program);
	}

//		//writes errors, theta etc of each individual in population
	protected void output2(int currentGeneration, File file2)throws IOException{
		   //"CurrentRun,Generation,ID,trainingTheta,UnseenTheta,reconTrainingError,reconUnseenError,LowestDistance"
		      FileWriter writer = new FileWriter(file2,true); 
		      
		      if (currentGeneration==0||mp1==null&&mp2==null){
			      // Writes the content to the file
			      writer.write("\n"+Main.CURRENTRUN+","+currentGeneration+","+getId()+","+","+
			    		  ","+getProgram(0).getTrainingError()+","+getProgram(0).getUnseenError()+","+getProgram(1).getTrainingError()+","+ getProgram(1).getUnseenError()+
			    		  ","+ getTrainingTheta()+","+getUnseenTheta()+","+getReconTrainingError()+
			    		  ","+getReconUnseenError()+","+minDistance);
			      writer.flush();
			      writer.close();	  
		      }
		      else if(mp2==null){
		    	 
			      // Writes the content to the file
			      writer.write("\n"+Main.CURRENTRUN+","+currentGeneration+","+getId()+","+mp1.getId()+","+
			    		  ","+getProgram(0).getTrainingError()+","+getProgram(0).getUnseenError()+","+getProgram(1).getTrainingError()+","+ getProgram(1).getUnseenError()+
			    		  ","+ getTrainingTheta()+","+getUnseenTheta()+","+getReconTrainingError()+
			    		  ","+getReconUnseenError()+","+minDistance);
			      writer.flush();
			      writer.close();	
		      }
		      else{
			      // Writes the content to the file
			      writer.write("\n"+Main.CURRENTRUN+","+currentGeneration+","+getId()+","+mp1.getId()+","+mp2.getId()+
			    		  ","+getProgram(0).getTrainingError()+","+getProgram(0).getUnseenError()+","+getProgram(1).getTrainingError()+","+ getProgram(1).getUnseenError()+
			    		  ","+getTrainingTheta()+","+getUnseenTheta()+","+getReconTrainingError()+
			    		  ","+getReconUnseenError()+","+minDistance);
			      writer.flush();
			      writer.close();	
		      }
		   
	}
//	// Writes the expression training error vectors and training output vectors  to  file
//	public void printVectors(int currentGeneration,File file) throws IOException{
//		
//		for(int i=0;i<numPrograms;i++){
//			FileWriter writerOut = new FileWriter(file,true);
//			writerOut.write("\n"+Main.CURRENTRUN+"_"+currentGeneration+","+getId()+
//				","+i+","+ Arrays.toString(trainingErrorVectors[i]).replace("[","").replace("]", "")+","+Arrays.toString(getProgram(i).getTrainingDataOutputs()).replace("[","").replace("]", ""));
//	      writerOut.flush();
//	      writerOut.close();
//		}   		
//	}
	// Writes the training output vectors  to  file
	public void printVectors(int currentGeneration,File file) throws IOException{
		
		for(int i=0;i<numPrograms;i++){
			FileWriter writerOut = new FileWriter(file,true);
			writerOut.write("\n"+Main.CURRENTRUN+"_"+currentGeneration+","+getId()+
				","+i+","+ Arrays.toString(getProgram(i).getTrainingDataOutputs()).replace("[","").replace("]", ""));
	      writerOut.flush();
	      writerOut.close();
		}   		
	}
    //checks the distances of all the expressions in a MIndividual from each other, returns the min distance found.
    public double calcDistances(){

        double [] ioutputs = this.getProgram(0).getTrainingDataOutputs();
        double [] joutputs = this.getProgram(1).getTrainingDataOutputs();
        double distance=calculateEuclideanDistance(ioutputs, joutputs);
        
        for(int i=2; i<numPrograms; i++){
        	  for(int j=i + 1; j<numPrograms; j++){
        		  ioutputs = this.getProgram(i).getTrainingDataOutputs();
        		  joutputs = this.getProgram(j).getTrainingDataOutputs();
  	    		double temp = calculateEuclideanDistance(ioutputs, joutputs);
  	    		if (temp<distance){
  	    			distance=temp;
  	    			}
        	  }
        	}        	
        minDistance =  distance;
        
		return distance;
	}
    //overloaded method
    //checks distances of a new expression from the expressions already within the MIndividual. Returns the minimum distance
   public double calcDistances(Individual newInd, int j){

	        double [] joutputs = newInd.getTrainingDataOutputs();
	        double [] koutputs = this.getProgram(0).getTrainingDataOutputs();
	        double distance = calculateEuclideanDistance(koutputs, joutputs);
	        
	    	for (int k=1;k<j;k++){	
	    		
	    		koutputs = this.getProgram(k).getTrainingDataOutputs();    		
	    		double temp = calculateEuclideanDistance(koutputs, joutputs);
	    		if (temp<distance){
	    			distance=temp;
	    			}	    		
	    	}  
	    minDistance = distance;    	
   	return distance;
   }
   
   //calcs ratios between outputs of the two expressions.
   public void calcRatios(double[]oneSemantics,double[]twoSemantics){

		double[] ratiosTemp=new double[oneSemantics.length];
					
		for (int i=0;i<oneSemantics.length;i++){
			if(twoSemantics[i]!=0){
				ratiosTemp[i]=oneSemantics[i]/twoSemantics[i];
			}
		}
		//remove zeros from the array
//	    int j = 0;
//	    for( int i=0;  i<ratiosTemp.length;  i++ )
//	    {
//	        if (ratiosTemp[i] != 0)
//	        	ratiosTemp[j++] = ratiosTemp[i];
//	    }
//	    double [] newArray = new double[j];
//	    System.arraycopy( ratiosTemp, 0, newArray, 0, j );
//	    Arrays.sort(newArray)
		Arrays.sort(ratiosTemp);;
		ratiosK=ratiosTemp;
		
	}
   //calc ratios between outputs of the three expressions
   public void calcRatios(double[]p1s,double[]p2s,double[]p3s, double[][] trainingData){

	   int sum=0;
		for (int i=0;i<p1s.length-1;i++){
			for(int j=i+1;j<p1s.length;j++){
				sum++;
			}
		}
		double[] ratiosTemp=new double[sum];
		int index = 0;
		for (int i=0;i<p1s.length-1;i++){
			for(int j=i+1;j<p1s.length;j++){
				double t1 = trainingData[i][trainingData[0].length - 1];
				double t2 = trainingData[j][trainingData[0].length - 1];
				if((p1s[i]*p3s[j]-p1s[j]*p3s[i]-p2s[i]*p3s[j]+p2s[j]*p3s[i]-p1s[i]*t2+p1s[j]*t1+p2s[i]*t2-p2s[j]*t1)!=0){
					ratiosTemp[index]=(p1s[i]*p2s[j]-p1s[j]*p2s[i]-p1s[i]*t2+p1s[j]*t1+p2s[i]*t2-p2s[j]*t1)/
						(p1s[i]*p3s[j]-p1s[j]*p3s[i]-p2s[i]*p3s[j]+p2s[j]*p3s[i]-p1s[i]*t2+p1s[j]*t1+p2s[i]*t2-p2s[j]*t1);
				}
				
				index++;
			}			
		}
		//remove zeros from the array
//	    int j = 0;
//	    for( int i=0;  i<ratiosTemp.length;  i++ )
//	    {
//	        if (ratiosTemp[i] != 0)
//	        	ratiosTemp[j++] = ratiosTemp[i];
//	        	//System.out.println(ratiosTemp[i]);
//	    }
//	    double [] newArray = new double[j];
//	    System.arraycopy( ratiosTemp, 0, newArray, 0, j );
	    Arrays.sort(ratiosTemp);
		ratiosK=ratiosTemp;		
		
	}
   //checks ratios between outputs of the two expressions.
   public boolean checkRatios(){
		double[] programOneSemantics = getProgram(0).getTrainingErrorVector();
		double[] programTwoSemantics = getProgram(1).getTrainingErrorVector();
		calcRatios(programOneSemantics, programTwoSemantics);
	    boolean check=false;
		for (int i=0;i<ratios.length;i++){
			
			if (Math.abs(ratios[i])>=1.2||Math.abs(ratios[i])<=0.98){
				
			}
			else{
				check=true;
			}
		}
		
		return check;
	} 
   //checks ratios between outputs of the two expressions.
   //overloaded method, checks ratio of new expression to current expressions in MIndividual
   public boolean checkRatios(Individual ind, int j){
		double[] programOneSemantics = getProgram(0).getTrainingErrorVector();
		//WHEN MORE THAN TWO EXPRESSIONS, THIS WILL BE DIFFERENT - GET OUTPUTS OF EACH EXPRESSION (UP TO INDEX J-1), PLUS OUT PUT OF IND
		double[] programTwoSemantics = ind.getTrainingErrorVector();
		calcRatios(programOneSemantics, programTwoSemantics);
	    boolean check=false;
		for (int i=0;i<ratios.length;i++){
			
			if (Math.abs(ratios[i])>=1.2||Math.abs(ratios[i])<=0.98){
				
			}
			else{
				check=true;
			}
		}
		
		return check;
	} 
	protected double calculateEuclideanDistance(double[] koutputs, double[] joutputs) {
		double sum = 0.0;
		for (int i = 0; i < koutputs.length; i++) {
			
			sum += Math.pow(koutputs[i] - joutputs[i], 2.0);
		
		}
		
		return Math.sqrt(sum);
	}

	
	public double[] reconstructTrainingSemantics(double[][] trainingData){
		//!!!change these functions to deal with N expressions	
		double[] reconstructedTrainingSemantics = new double[trainingErrorVectors[0].length];		
		
		if (numPrograms==2){
			double[] programOneSemantics = getProgram(0).getTrainingDataOutputs();
			double[] programTwoSemantics = getProgram(1).getTrainingDataOutputs();	
			double[] programOneError = getProgram(0).getTrainingErrorVector();
			double[] programTwoError = getProgram(1).getTrainingErrorVector();
			k = calculateK(trainingData,programOneError,programTwoError);
			
			for (int i = 0; i < programOneSemantics.length; i++) {
				reconstructedTrainingSemantics[i] =  1/(1-k)*programOneSemantics[i]-k/(1-k)*programTwoSemantics[i];
			}
		}
		else{
			//calc k and w using median from each index
			//create array variable to hold k, one for w
			//use function to get medians
			//reconstruct
			double[] programOneSemantics = getProgram(0).getTrainingDataOutputs();
			double[] programTwoSemantics = getProgram(1).getTrainingDataOutputs();
			double[] programThreeSemantics = getProgram(2).getTrainingDataOutputs();
			double[] programOneError = getProgram(0).getTrainingErrorVector();
			double[] programTwoError = getProgram(1).getTrainingErrorVector();
			double[] programThreeError = getProgram(2).getTrainingErrorVector();
			w = calculateW(trainingData);
			k = calculateK(trainingData,programOneError,programTwoError,programThreeError);
			for (int i = 0; i < programOneSemantics.length; i++) {
				reconstructedTrainingSemantics[i] =  1/((1-k)*(1-w))*programOneSemantics[i]-w/((1-k)*(1-w))*programTwoSemantics[i]
						-k/(1-k)*programThreeSemantics[i];
			}					
			
		}
		return reconstructedTrainingSemantics;
	}	
	
	public double[] reconstructUnseenSemantics(double[][] unseenData){
		
		//!!!change these functions to deal with N expressions
		double[] reconstructedUnseenSemantics = new double[unseenErrorVectors[0].length];
		
		if (numPrograms==2){	
			double[] programOneSemantics = getProgram(0).getUnseenDataOutputs();
			double[] programTwoSemantics = getProgram(1).getUnseenDataOutputs();
			double[] programOneError = getProgram(0).getUnseenErrorVector();;
			double[] programTwoError = getProgram(1).getUnseenErrorVector();
			//kunseen = calculateK(unseenData,programOneError,programTwoError);
			for (int i = 0; i < programOneSemantics.length; i++) {
				reconstructedUnseenSemantics[i] = 1/(1-k)*programOneSemantics[i]-k/(1-k)*programTwoSemantics[i];
			}	
		}
		else{
			double[] programOneSemantics = getProgram(0).getUnseenDataOutputs();
			double[] programTwoSemantics = getProgram(1).getUnseenDataOutputs();
			double[] programThreeSemantics = getProgram(2).getUnseenDataOutputs();
			//using k and w calculated from training data
			for (int i = 0; i < programOneSemantics.length; i++) {
				reconstructedUnseenSemantics[i] =  1/((1-k)*(1-w))*programOneSemantics[i]-w/((1-k)*(1-w))*programTwoSemantics[i]
						-k/(1-k)*programThreeSemantics[i];
			}
		}
		return reconstructedUnseenSemantics;
	}	
	
	//method for calculating K of an individual with all expressions added
	public double calculateK(double[][] data,double[] programOneError,double[] programTwoError){
		double ktemp;
		calcRatios(programOneError, programTwoError);
		//get median
		if (ratiosK.length % 2 == 0)
		   ktemp = ((double)ratiosK[ratiosK.length/2] + (double)ratiosK[ratiosK.length/2 - 1])/2;
			//k = (double)ratiosK[ratiosK.length/2];
		else
		   ktemp = (double) ratiosK[ratiosK.length/2];
				
		return ktemp;
	}
	
	//method for calculating K of an individual with all expressions added - three expressions
	public double calculateK(double[][] data,double[] p1s,double[] p2s,double[] p3s){

			calcRatios(p1s,p2s,p3s,data);
			
			if (ratiosK.length==0){
				k=1;
			}
			
			else if (ratiosK.length % 2 == 0){
				k = ((double)ratiosK[ratiosK.length/2] + (double)ratiosK[ratiosK.length/2 - 1])/2;
					//k = (double)ratios[ratios.length/2];
			}
			else{
				k = (double) ratiosK[ratiosK.length/2];
			}
		
		return k;
	}
	
	//overloaded for calculating what W would be before adding the next individual
	public double calculateW(Individual ind, double[][] trainingData){
		//!!seperate this out to match the K calc methods - a calcRatiosW and calcW separate
		double[] p1s = getProgram(0).getTrainingErrorVector();
		double[] p2s = getProgram(1).getTrainingErrorVector();
		double[] p3s = ind.getTrainingErrorVector();
		int sum=0;
		for (int i=0;i<p1s.length-1;i++){
			for(int j=i+1;j<p1s.length;j++){
				sum++;
			}
		}
		ratiosW=new double[sum];
		int index = 0;
		for (int i=0;i<p1s.length-1;i++){
			for(int j=i+1;j<p1s.length;j++){
				
					double t1 = trainingData[i][trainingData[0].length - 1];
					double t2 = trainingData[j][trainingData[0].length - 1];
					if((p2s[i]*p3s[j]-p2s[j]*p3s[i]-p2s[i]*t2+p2s[j]*t1+p3s[i]*t2-p3s[j]*t1)!=0){
					ratiosW[index]=(p1s[i]*p3s[j]-p1s[j]*p3s[i]-p1s[i]*t2+p1s[j]*t1+p3s[i]*t2-p3s[j]*t1)/
							(p2s[i]*p3s[j]-p2s[j]*p3s[i]-p2s[i]*t2+p2s[j]*t1+p3s[i]*t2-p3s[j]*t1);
					}
					index++;
				}			
		}
		Arrays.sort(ratiosW);
		if (ratiosW.length % 2 == 0)
			   w = ((double)ratiosW[ratiosW.length/2] + (double)ratiosW[ratiosW.length/2 - 1])/2;
		else
			    w = (double) ratiosW[ratiosW.length/2];
		return w;
	}
	
	public double calculateW(double[][] trainingData){
		//!!!fix to accept N expressions
		double[] p1s = getProgram(0).getTrainingErrorVector();
		double[] p2s = getProgram(1).getTrainingErrorVector();
		double[] p3s = getProgram(2).getTrainingErrorVector();
		
		int sum=0;
		for (int i=0;i<p1s.length-1;i++){
			for(int j=i+1;j<p1s.length;j++){
				sum++;
			}
		}
		ratiosW=new double[sum];
		int index = 0;
		for (int i=0;i<p1s.length-1;i++){
			for(int j=i+1;j<p1s.length;j++){
				
					double t1 = trainingData[i][trainingData[0].length - 1];
					double t2 = trainingData[j][trainingData[0].length - 1];
					if((p2s[i]*p3s[j]-p2s[j]*p3s[i]-p2s[i]*t2+p2s[j]*t1+p3s[i]*t2-p3s[j]*t1)!=0){
						ratiosW[index]=(p1s[i]*p3s[j]-p1s[j]*p3s[i]-p1s[i]*t2+p1s[j]*t1+p3s[i]*t2-p3s[j]*t1)/
							(p2s[i]*p3s[j]-p2s[j]*p3s[i]-p2s[i]*t2+p2s[j]*t1+p3s[i]*t2-p3s[j]*t1);
					}
					index++;
				}			
		}
		Arrays.sort(ratiosW);
		if (ratiosW.length % 2 == 0)
			   w = ((double)ratiosW[ratiosW.length/2] + (double)ratiosW[ratiosW.length/2 - 1])/2;
			else
			    w = (double) ratiosW[ratiosW.length/2];
		return w;
	}
	
//	public double calculateRMSE(double[][] data, double[] outputs) {
//		double errorSum = 0.0;
//		for (int i = 0; i < outputs.length; i++) {
//			double target = data[i][data[0].length - 1];
//			errorSum += Math.pow(outputs[i] - target, 2.0);
//		}
//		return Math.sqrt(errorSum / data.length);
//	}
	
    // return the inner product of vectors a and b
    public double dot(double[] vectorOne, double[] vectorTwo) {
//        if (vectorOne.length != vectorTwo.length)
//            throw new IllegalArgumentException("Dimensions disagree for angle calculation");
        double sum = 0.0;
        for (int i = 0; i < vectorOne.length; i++)
            sum = sum + (vectorOne[i] * vectorTwo[i]);
        return sum;
    }
    
    // return the Euclidean norm of this Vector
    public double magnitude(double[] vector) {
        return Math.sqrt(dot(vector,vector));
    }
    
    public void reorderMindividual(int[] order){
//    	//copy mindividual to new
    	MIndividual temp=new MIndividual();
    	
    	for(int i=0;i<numPrograms;i++){
    		temp.addProgramAtIndex(this.getProgram(i), i);
    		//System.out.println(this.getProgram(i).getId());
    	}
    	//iterate through new, assign new[i] to mindividual[order[i]]
    	for(int i=0;i<numPrograms;i++){   
    		//System.out.println("reordered "+temp.getProgram(i).getId()+" "+order[i]);
    		replaceProgramAtIndex(temp.getProgram(order[i]), i);
    		//System.out.println("reordered "+temp.getProgram(i).getId()+" to "+ order[i]); 
    	}    	 	   	
    }
    
	public void print() {
		Individual p;
		for(int i=0;i<numPrograms;i++){
			p=getProgram(i);
			if (sizeOverride == true) {
				System.out.println(" [Individual not constructed]");
			} else {
				printIndex = 0;
				System.out.print("Expression "+(i+1)+" =>");
				printInner(p);
				System.out.println();
				
			}			
		}
	}

	protected void printInner(Individual p) {
		if (p.program.get(printIndex) instanceof Terminal) {
			System.out.print(" " + p.program.get(printIndex));
		} else {
			System.out.print(" (");
			System.out.print(p.program.get(printIndex));
			Operator currentOperator = (Operator) p.program.get(printIndex);
			for (int i = 0; i < currentOperator.getArity(); i++) {
				printIndex++;
				printInner(p);
			}
			System.out.print(")");
		}
		
	}
	
	public void writeToObjectFile() {
		if (sizeOverride == true) {
			System.out.println("Individual not constructed!");
		} else {
			Utils.writeObjectToFile(this, Main.resultsFolderName+"/"+Main.DATA_FILENAME+"/"+Main.modelName+"/individuals/final_individual_"+Main.CURRENTRUN+".obj");

		}
	}

//	// ##### get's and set's from here on #####
	
	public double getK(){
		
		return k;
	}
	
	public double getKunseen(){
		
		return kunseen;
	}
	
	public double getW(){
		
		return w;
	}
	
	public double getTrainingTheta() {
		
		return trainingTheta;
	}
	
	public double getUnseenTheta() {
		
		return unseenTheta;
	}

	public long getId() {
			return id;
		}
	

	public  int getNumPrograms() {
		return numPrograms;
	}

	public Individual getProgram(int index) {
		Individual program = programs.get(index);
		return program;
	}
	
	public double getReconTrainingError(){
		return reconTrainingError;
	}
	
	public double getReconUnseenError(){
		return reconUnseenError;
	}
	
	public double getReconDepth(){
		return reconstructed.getDepth();
	}
	
	public double getReconSize(){
		return reconstructed.getSize();
	}
	
	public double[] getUnseenOutputAvg(){
		return unseenOutputAverages;
	}
	
	public double[] getTrainingOutputAvg(){
		return trainingOutputAverages;
	}
	
	public double[] getUnseenOutputSD(){
		return unseenOutputSD;
	}
	
	public double[] getTrainingOutputSD(){
		return trainingOutputSD;
	}
	
	public double getTrainingOutputSD(int expression){
		return trainingOutputSD[expression];
	}
	
	public double getTrainingErrorSum(){
		double trainingErrorSum=0;
		for(int i=0;i<trainingErrorVectors.length;i++){
			trainingErrorSum+=this.getProgram(i).getTrainingError();
		}
		return trainingErrorSum;
	}
	public Individual getReconstructedIndividual(){
		return reconstructed;
	}
	
	public void setSizeOverride(boolean sizeOverride) {
		this.sizeOverride = sizeOverride;
	}

	public void setMp1(MIndividual p1) {
		this.mp1 = p1;
	}
	public void setMp2(MIndividual p2) {
		this.mp2 = p2;
	}
	public void setReconstructedIndividual(Individual ind) {
		this.reconstructed = ind;
	}
}
