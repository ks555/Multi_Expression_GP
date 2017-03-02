package core;
//GSGP crossover and mutation used
//Also a method for switching expressions. can use instead of crossover - need to comment out/in two 
//places to implement. Also check crossover probability on GPrun.

import java.io.IOException;
import java.util.ArrayList;

import core.MIndividual;
import programElements.Addition;
import programElements.Constant;
import programElements.Multiplication;
import programElements.ProtectedDivision;
import programElements.Subtraction;

import java.io.File;
import java.io.FileWriter;

public class EsgsgpRun extends GsgpRun {

	private static final long serialVersionUID = 7L;

	protected double mutationStep;
	protected boolean boundedMutation;
	//protected boolean buildIndividuals;
	protected ArrayList<Population> populations;
	protected MPopulation mpopulation;
	protected int numPrograms;
	protected MIndividual currentMBest;
	protected double minDistance;
	protected double minKDiff;
	protected double minUnseen;
	
	protected Population reconstructedPopulation;
	//protected File file;
	protected File file2;
	//protected File OutputsFile;

	public EsgsgpRun(Data data) throws IOException {
		super(data);				
	}

	protected void initialize() throws IOException {
		System.out.println("Master");

		populations = new ArrayList <Population>();
		mpopulation = new MPopulation();
		numPrograms = MIndividual.numPrograms;
		minDistance=50;
		minKDiff=0.02;
		file2 = new File("results/population.txt");
		//OutputsFile = new File("results/outputs.txt");
		
		
		
		for (int i = 0; i <  numPrograms; i++) {
			
			super.initialize();
			
			populations.add(population);
			
		}
		for (int i = 0; i < population.getSize(); i++) {
			
			//build the Mindividuals using the individuals from the initialized populations
			MIndividual mindividual = new MIndividual();
			
			for (int j = 0; j < numPrograms; j++) {
				Individual ind = populations.get(j).getIndividual(i);
				//check if any errors of the program is > maxError . if not, add to Mindividual
				ind.evaluate(data);	
				mindividual.addProgramAtIndex(ind,j);						
			}

		//if k is too small, replace ind with a new program
			double k;
			double w;
			k = mindividual.getK();
			if (numPrograms>2){
				w=mindividual.getW();	
			}
			else{
				w=2;
			}
			
			while((Math.abs(k)<1+minKDiff&&Math.abs(k)>1-minKDiff)||(Math.abs(w)<1+minKDiff&&Math.abs(w)>1-minKDiff)||
					(Math.abs(w)<0+minKDiff&&Math.abs(w)>0-minKDiff)||(Math.abs(k)<0+minKDiff&&Math.abs(k)>0-minKDiff)){
				//if k is too small, replace ind with a new program
				for (int j = 0; j < numPrograms; j++) {
					Individual ind = grow(this.getMaximumDepth());
					//check if any errors of the program is > maxError . if not, add to Mindividual
					ind.evaluate(data);	
					mindividual.replaceProgramAtIndex(ind,j);
				}
				mindividual.evaluate(data);
				k = mindividual.getK();
				if (numPrograms>2){
					w=mindividual.getW();	
				}
				else{
					w=2;
				}
			}

			mindividual.evaluate(data);
			reconstructIndividual(mindividual);
			mpopulation.addIndividual(mindividual);			
			
		}	
		
		updateCurrentMBest();	
		//initial min unseen as the unseen of best ind at gen 0
		minUnseen=currentMBest.getReconUnseenError();
		currentMBest.writeToObjectFile();
		printMPopState(0);
		
		try {
			createOutputFile();
			output(0);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

//	public void printFirstGen(){
//		for (int i = 0; i < mpopulation.getSize(); i++) {
//			//print to each mindividual output2, (population.txt)
//			try {
//				mpopulation.getMIndividual(i).output2(0, file2);
//				output();
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//			
//		}
//	}
	
	public void evolve(int numberOfGenerations) throws IOException {
		//first tournament selection - selectLowestError for num inds in population
		//set this to population
		
		// evolve for a given number of generations
		while (currentGeneration <= numberOfGenerations) {
			MPopulation offspring = new MPopulation();

			// generate a new offspring population
			while (offspring.getSize() < population.getSize()) {
				MIndividual mp1, newIndividual;
				//mp1 = nestedSelectMParent();
				mp1 = selectMParent();
				
				// apply crossover to parents selected by tournament selection
				if (randomGenerator.nextDouble() < crossoverProbability) {
					//MIndividual mp2 = nestedSelectMParent();
					MIndividual mp2 = selectMParent();
					newIndividual = applyStandardCrossover(mp1, mp2);
					//newIndividual = applyExpressionSwitch(mp1, mp2);
					boolean flag=false;
					//check errors of each program
					for(int i=0; i<newIndividual.numPrograms;i++){
						if(newIndividual.getProgram(i).checkMaxError()){
							flag=true;
						}
					}
					//check distances of each program
					double distance=newIndividual.calcDistances();
//					double k = newIndividual.calculateK(data.getTrainingData());
//					double w = newIndividual.calculateW(data.getTrainingData());
					newIndividual.evaluate(data);
					double k =  newIndividual.getK();
					double w = newIndividual.getW();
					while((Math.abs(k)<1+minKDiff&&Math.abs(k)>1-minKDiff)||(Math.abs(w)<1+minKDiff&&Math.abs(w)>1-minKDiff)||(Math.abs(w)<0+minKDiff&&Math.abs(w)>0-minKDiff)||(Math.abs(k)<0+minKDiff&&Math.abs(k)>0-minKDiff)){
						
//						mp1 = nestedSelectMParent();
//						mp2 = nestedSelectMParent();
						mp1 = selectMParent();
						mp2 = selectMParent();
						newIndividual = applyStandardCrossover(mp1, mp2);
						//newIndividual = applyExpressionSwitch(mp1, mp2);
						flag=false;
						for(int i=0; i<newIndividual.numPrograms;i++){
							if(newIndividual.getProgram(i).checkMaxError()){
								flag=true;
							}
						}
						distance=newIndividual.calcDistances();
//						k = newIndividual.calculateK(data.getTrainingData());
//						w = newIndividual.calculateW(data.getTrainingData());
						newIndividual.evaluate(data);
						k =  newIndividual.getK();
						w = newIndividual.getW();
					}
					
			
					newIndividual.setMp1(mp1);					
					newIndividual.setMp2(mp2);
					
				}
				// apply mutation
				else {
					
					newIndividual = applyStandardMutation(mp1);
					
					//newIndividual=applyConstantMultiplicationSemantics(newIndividual);
					newIndividual.evaluate(data);
					double k =  newIndividual.getK();
					double w = newIndividual.getW();
					if(numPrograms>2){
												
						w = newIndividual.getW();
					}
					else{
						w=2;
					}
					while((Math.abs(k)<1+minKDiff&&Math.abs(k)>1-minKDiff)||(Math.abs(w)<1+minKDiff&&Math.abs(w)>1-minKDiff)||(Math.abs(w)<0+minKDiff&&Math.abs(w)>0-minKDiff)||(Math.abs(k)<0+minKDiff&&Math.abs(k)>0-minKDiff)){

						//mp1 = nestedSelectMParent();	
						mp1 = selectMParent();
						//flag=false;
						newIndividual = applyStandardMutation(mp1);
						//newIndividual=applyConstantMultiplicationSemantics(newIndividual);
						newIndividual.evaluate(data);
						k =  newIndividual.getK();
						
						if(numPrograms>2){
							w = newIndividual.getW();
							//w = newIndividual.calculateW(data.getTrainingData());
						}
					}
					newIndividual.setMp1(mp1);
				}
				
				/*
				 * add the new individual to the offspring population if its
				 * depth is not higher than the maximum (applicable only if the
				 * depth limit is enabled)
				 */
				boolean flag = false;
				for (int i=0; i<numPrograms; i++){
					if(newIndividual.getProgram(i).getDepth()>maximumDepth){
						flag=true;
					}
				}
				if (applyDepthLimit && flag) {
					newIndividual = mp1;
					
				} else {
					for (int i=0; i<numPrograms; i++){
						//evaluate each program in the individual
						newIndividual.getProgram(i).evaluate(data);
						
					}					
					
					newIndividual.evaluate(data);
					//newIndividual.printVectors(currentGeneration,OutputsFile);
				}
				reconstructIndividual(newIndividual);
				offspring.addIndividual(newIndividual);
			}
			
			mpopulation = selectSurvivors(offspring);
			//setReconPopulation();
			updateCurrentMBest();
			currentMBest.printVectors(currentGeneration, OutputsFile);
			if(currentMBest.getReconUnseenError()<minUnseen){
				minUnseen=currentMBest.getReconUnseenError();
				currentMBest.writeToObjectFile();
			}
			printMPopState();
			output();
			
			currentGeneration++;
		}
	}
	
	protected void printMPopState() {
		if (printAtEachGeneration&&numPrograms==2) {
			System.out.println("\nGeneration:\t\t" + currentGeneration);
			System.out.printf("Training Theta (deg):\t\t%.2f\nUnseen Theta (deg):\t\t%.2f\nReconstructed Training Error:"
					+ "\t\t%.2f\nReconstructed Unseen Error:\t\t%.2f\nTraining Error Program 1:\t\t%.2f\nTraining Error Program 2:\t\t%.2f\n"
					+ "Id:\t\t%d\nK:\t\t%.2f\n",
					currentMBest.getTrainingTheta(), currentMBest.getUnseenTheta(), currentMBest.getReconTrainingError(),
					currentMBest.getReconUnseenError(),currentMBest.getProgram(0).getTrainingError(),currentMBest.getProgram(1).getTrainingError(),
					currentMBest.getId(),currentMBest.getK());
			
			
		}
		else if(printAtEachGeneration&&numPrograms==3) {
			System.out.println("\nGeneration:\t\t" + currentGeneration);
			System.out.printf("Training Theta (deg):\t\t%.2f\nUnseen Theta (deg):\t\t%.2f\nReconstructed Training Error:"
					+ "\t\t%.2f\nReconstructed Unseen Error:\t\t%.2f\nTraining Error Program 1:\t\t%.2f\nTraining Error Program 2:\t\t%.2f\nTraining Error Program 3:\t\t%.2f\n"
					+ "Id:\t\t%d\nK:\t\t%.2f\nW:\t\t%.2f\n",
					currentMBest.getTrainingTheta(), currentMBest.getUnseenTheta(), currentMBest.getReconTrainingError(),
					currentMBest.getReconUnseenError(),currentMBest.getProgram(0).getTrainingError(),currentMBest.getProgram(1).getTrainingError(),
					currentMBest.getProgram(2).getTrainingError(),currentMBest.getId(),currentMBest.getK(),currentMBest.getW());
			
			
		}
	}
	//overloaded for printing Generation 0 - currentGeneration is already incremented when print is called the first time
	protected void printMPopState(int generation) {
		if (printAtEachGeneration) {
			System.out.println("\nGeneration:\t\t0");
			System.out.printf("Training Theta:\t\t%.2f\nUnseen Theta:\t\t%.2f\nReconstructed Training Error:"
					+ "\t\t%.2f\nReconstructed Unseen Error:\t\t%.2f\nTraining Error Program 1:\t\t%.2f\nTraining Error Program 2:\t\t%.2f\n",
					currentMBest.getTrainingTheta(), currentMBest.getUnseenTheta(), currentMBest.getReconTrainingError(),
					currentMBest.getReconUnseenError(),currentMBest.getProgram(0).getTrainingError(),currentMBest.getProgram(1).getTrainingError());
			
		}
	}
	
	protected MIndividual nestedSelectMParent() {		
		MPopulation tournamentPopulationOne;
		MPopulation tournamentPopulationTwo;
		MPopulation tournamentPopulationThree;
		MPopulation tournamentPopulationFour;
		tournamentPopulationFour=new MPopulation();
		int tournamentSize = (int) (0.10 * population.getSize());
		
		//from tournamentPopulation, select lowest error
		//repeat N times where N = tournament size, this is tournamentPopulation2
		
		for (int m=0;m<tournamentSize;m++){	
			tournamentPopulationFour=new MPopulation();
				for (int l=0;l<tournamentSize;l++){
					tournamentPopulationThree= new MPopulation();
					for (int k=0;k<tournamentSize;k++){
						tournamentPopulationTwo = new MPopulation();
						//fill tournamentPopulationTwo with winners of diversity outputs tournaments
						for (int j=0;j<tournamentSize;j++){
							tournamentPopulationOne = new MPopulation();				
							//get random individuals to fill first tournament pop
							for (int i = 0; i < tournamentSize; i++) {
								int index = randomGenerator.nextInt((mpopulation.getSize()));					
								tournamentPopulationOne.addIndividual(mpopulation.getMIndividual(index));	
							}				
							tournamentPopulationTwo.addIndividual(tournamentPopulationOne.getMostDiverseM());
							//System.out.println(tournamentPopulationOne.getMostDiverseM().getId());
						}	
						tournamentPopulationThree.addIndividual(tournamentPopulationTwo.getHighestK());			
					}
				tournamentPopulationFour.addIndividual(tournamentPopulationThree.getLowestErrorM());
			}
		}
		return tournamentPopulationFour.getLowestAngleM();
	}
	
	
	protected MIndividual selectMParent() {
		MPopulation tournamentPopulation = new MPopulation();
		
		int tournamentSize = (int) (0.05 * population.getSize());
		
		for (int i = 0; i < tournamentSize; i++) {
			int index = randomGenerator.nextInt((mpopulation.getSize()));
			
			tournamentPopulation.addIndividual(mpopulation.getMIndividual(index));

		}
		//System.out.println("Tournament get best");
		return tournamentPopulation.getBestM();
	}	
	
	protected MIndividual applyExpressionSwitch(MIndividual mp1, MIndividual mp2) {
		MIndividual offspring = new MIndividual();
		if (buildIndividuals == false) {
			offspring.setSizeOverride(true);
		}
		//!!!for two expressions only
			
		Individual p1 =  mp1.getProgram(randomGenerator.nextInt(numPrograms));
	
		Individual p2 = mp2.getProgram(randomGenerator.nextInt(numPrograms));				
		offspring.addProgramAtIndex(p1,0);
		offspring.addProgramAtIndex(p2,1);
		
		return offspring;
	}
	
	
	protected MIndividual applyStandardCrossover(MIndividual mp1, MIndividual mp2) {
		MIndividual offspring = new MIndividual();
		if (buildIndividuals == false) {
			offspring.setSizeOverride(true);
		}
		//for each program in numPrograms
		//choose random program from mp1 and mp2 as 'parents' of the new program.
		//add each program 'offspring' to the final offspring (an MIndividual)
		for(int i=0; i<numPrograms; i++){
			
			int rand = randomGenerator.nextInt(numPrograms);
			Individual p1 =  mp1.getProgram(rand);
			rand =  randomGenerator.nextInt(numPrograms);
			Individual p2 = mp2.getProgram(rand);
			Individual pOffspring;
				
			if (buildIndividuals) {
				pOffspring = buildCrossoverIndividual(p1, p2);
				offspring.addProgramAtIndex(pOffspring, i);
			}
			else {
				pOffspring = buildCrossoverSemantics(p1, p2);
				offspring.addProgramAtIndex(pOffspring,i);
			}	
/*			//applying standard GP crossover, no semantic crossover, as above
			pOffspring=applyStandardCrossover(p1,p1);
			offspring.addProgramAtIndex(pOffspring,i);
			pOffspring.evaluate(data);*/
		}	
		//print
//		System.out.println("offspring ID: " +offspring.getId());
//		System.out.println("P1 training error " + mp1.getTrainingError());
//		System.out.println("P2 training error " + mp2.getTrainingError());
		
		return offspring;
	}	
	
	protected MIndividual applyStandardMutation(MIndividual mp) {
		
		MIndividual offspring = new MIndividual();		
		
		for(int i=0; i<mp.getNumPrograms(); i++){
			Individual pOffspring = new Individual();
			if (buildIndividuals) {
				pOffspring=buildMutationIndividual(mp.getProgram(i));
				offspring.addProgramAtIndex(pOffspring,i);
				
			} else {				
				pOffspring=buildMutationSemantics(mp.getProgram(i));
				offspring.addProgramAtIndex(pOffspring,i);				
			}			
		}
		if(!buildIndividuals)offspring.setSizeOverride(true);
		return offspring;
	}
	
	protected MIndividual applyConstantMultiplicationSemantics(MIndividual newMindividual) {
		
		int rand = randomGenerator.nextInt(numPrograms);
		Individual p=newMindividual.getProgram(rand);
		
		//get random constant
		double constant = randomGenerator.nextInt(20)+1;		
		
		Individual multiplied = new Individual();
		double[] multipliedTrainingSemantics=new double[p.getTrainingDataOutputs().length] ;
		double[] multipliedUnseenSemantics=new double[p.getUnseenDataOutputs().length] ;

		//update training output
		for (int i = 0; i < p.getTrainingDataOutputs().length; i++) {
			multipliedTrainingSemantics[i]=p.getTrainingDataOutputs()[i]*constant;
		}
		multiplied.setTrainingDataOutputs(multipliedTrainingSemantics);
		//update unseen output
		for (int i = 0; i < p.getUnseenDataOutputs().length; i++) {
			multipliedUnseenSemantics[i]=p.getUnseenDataOutputs()[i]*constant;
		}
		multiplied.setUnseenDataOutputs(multipliedUnseenSemantics);
		
		
		multiplied.setDepth(p.getDepth()+1);
	
		multiplied.setComputedSize(p.getSize()+2);
		multiplied.setSizeOverride(true);
		newMindividual.replaceProgramAtIndex(multiplied, rand);		
		
		return newMindividual;
	}
	protected MIndividual applyConstantMultiplication(MIndividual newMindividual) {
		
		int rand = randomGenerator.nextInt(numPrograms);
		Individual p=newMindividual.getProgram(rand);
		
		//get random constant
		double constant = randomGenerator.nextInt(20)+1;		
		
		Individual multiplied = new Individual();
		multiplied.addProgramElement(new Multiplication());

		// copy original expression to multiplied version
		for (int i = 0; i < p.getSize(); i++) {
			multiplied.addProgramElement(p.getProgramElementAtIndex(i));
		}
		
		multiplied.addProgramElement(new Constant(constant));
		multiplied.calculateDepth();
		multiplied.evaluate(data);
		newMindividual.replaceProgramAtIndex(multiplied, rand);
		//System.out.println("Original one " + p.getTrainingDataOutputs()[0]+" multiplied by "+constant+" "+multiplied.getTrainingDataOutputs()[0]);
		//System.out.println("Original two " + p.getTrainingDataOutputs()[1]+" multiplied by "+constant+" "+multiplied.getTrainingDataOutputs()[1]);

		//return		
		
		return newMindividual;
	}

	
	// keep the best overall + all the remaining offsprings
	protected MPopulation selectSurvivors(MPopulation newIndividuals) {
		MPopulation survivors = new MPopulation();
		
		MIndividual bestParent = mpopulation.getBestM();
		System.out.println(bestParent.getTrainingTheta());
		
		MIndividual bestNewIndividual = newIndividuals.getBestM();
		MIndividual bestOverall;			

		// the best overall is in the current population
		//CHANGE FOR FITNESS
//		if (bestParent.getTrainingTheta() < bestNewIndividual.getTrainingTheta()) {
		if (bestParent.getReconTrainingError() < bestNewIndividual.getReconTrainingError()) {	
			bestOverall = bestParent;
		}
		// the best overall is in the offspring population
		else {
			bestOverall = bestNewIndividual;
			//System.out.println("Best  new" +bestNewIndividual.getTrainingError());
		}

		survivors.addIndividual(bestOverall);
		try {
			bestOverall.output2(currentGeneration, file2);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		for (int i = 0; i < newIndividuals.getSize(); i++) {
			if (newIndividuals.getMIndividual(i).getId() != bestOverall.getId()) {
				survivors.addIndividual(newIndividuals.getMIndividual(i));
				//print to new individual to output2, population
				try {
					newIndividuals.getMIndividual(i).output2(currentGeneration, file2);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

			}
		}

		return survivors;
	}
	
	protected void updateCurrentMBest() {
		
		currentMBest = mpopulation.getBestM();
	}

	protected MIndividual reconstructIndividual(MIndividual newMindividual) {
		//!!! need to fix for N programs
		Individual e1=newMindividual.getProgram(0);
		Individual e2=newMindividual.getProgram(1);
		Individual reconstructed=new Individual();
		if (buildIndividuals){
			reconstructed.addProgramElement(new Subtraction());
			reconstructed.addProgramElement(new Multiplication());
			reconstructed.addProgramElement(new ProtectedDivision());
			reconstructed.addProgramElement(new Constant(1.0));
			reconstructed.addProgramElement(new Subtraction());
			reconstructed.addProgramElement(new Constant(1.0));		
			reconstructed.addProgramElement(new Constant(newMindividual.getK()));		
	
			// copy E1 to reconstructed
			for (int i = 0; i < e1.getSize(); i++) {
				reconstructed.addProgramElement(e1.getProgramElementAtIndex(i));
			}
			
			reconstructed.addProgramElement(new Multiplication());
			reconstructed.addProgramElement(new ProtectedDivision());
			reconstructed.addProgramElement(new Constant(newMindividual.getK()));
			reconstructed.addProgramElement(new Subtraction());
			reconstructed.addProgramElement(new Constant(1.0));
			reconstructed.addProgramElement(new Constant(newMindividual.getKunseen()));	
			
			// copy E2 to reconstructed
			for (int i = 0; i < e2.getSize(); i++) {
				reconstructed.addProgramElement(e2.getProgramElementAtIndex(i));
			}
			
			reconstructed.calculateDepth();
			reconstructed.evaluate(data);
			
			newMindividual.setReconstructedIndividual(reconstructed);
		}
		
		return newMindividual;
	}
	
   protected void createOutputFile()throws IOException{
	      
	      // creates the file
	      file.createNewFile();
	      file2.createNewFile();
	      OutputsFile.createNewFile();
	   // creates a FileWriter Object
	      FileWriter writer = new FileWriter(file); 
	      // Writes the content to the file
	      writer.write("CurrentRun,Generation,NumGens,NumRuns,popSize,trainingTheta,UnseenTheta,"
	      		+ "reconTrainError,reconUnseenError,depth,depthLimitOn,maxDepth,crossoverProb,k,ID,kunseen,size,E1TrainError,"
	      		+ "E2TrainError,E1UnseenError,E2UnseenError,E1TrainOutAvg,E2TrainOutAvg,E1UnseenOutAvg,E2UnseenOutAvg,E1TrainOutSD,E2TrainOutSD,E1UnseenOutSD,E2UnseenOutSD"); 
	      writer.flush();
	      writer.close();
	   // creates a FileWriter Object
	      FileWriter writer2 = new FileWriter(file2); 
	      // Writes the content to the file
	      writer2.write("CurrentRun,Generation,ID,p1,p2,P1Training,P1Unseen,P2Training,P2Unseen,trainingTheta,UnseenTheta,reconTrainingError,reconUnseenError,LowestDistance"); 
	      writer2.flush();
	      writer2.close();
	   // creates a FileWriter Object
	      FileWriter writer3 = new FileWriter(OutputsFile); 
	      // Writes the content to the file
	      writer3.write("CurrentRun,ID,ExpressionNum,TrainingErrorVector0,"); 
	      for(int i=1;i<data.getTrainingData().length;i++){
	    	  writer3.write("TrainingErrorVector"+i+",");
	      }
	      writer3.write("TrainingOutputVector,");
	      for(int i=1;i<data.getTrainingData().length;i++){
	    	  writer3.write("TrainingOutputVector"+i+",");
	      }
	      writer3.flush();
	      writer3.close();
	   }	
	
   //updates results file
   protected void output()throws IOException{

	      FileWriter writer = new FileWriter(file,true); 
	      // Writes the content to the file
	      writer.write("\n"+Main.CURRENTRUN+","+currentGeneration+","+Main.NUMBER_OF_GENERATIONS+","+Main.NUMBER_OF_RUNS+
	    		  ","+populationSize+","+currentMBest.getTrainingTheta()+","+currentMBest.getUnseenTheta()+","+currentMBest.getReconTrainingError()+
	    		  ","+currentMBest.getReconUnseenError()+
	    		  ","+currentMBest.getDepth() + ","+applyDepthLimit+","+maximumDepth+
	    		  ","+crossoverProbability+","+currentMBest.getK()+","+currentMBest.getId()+
	    		  ","+currentMBest.getKunseen()+","+currentMBest.getSize()+
	    		  ","+currentMBest.getProgram(0).getTrainingError()+ ","+currentMBest.getProgram(1).getTrainingError()+
	    		  ","+currentMBest.getProgram(0).getUnseenError()+ ","+currentMBest.getProgram(1).getUnseenError()+
	    		  ","+currentMBest.getTrainingOutputAvg()[0]+","+currentMBest.getTrainingOutputAvg()[1]+
	      		  ","+currentMBest.getUnseenOutputAvg()[0]+","+currentMBest.getUnseenOutputAvg()[1]+
	      	   	","+currentMBest.getTrainingOutputSD()[0]+","+currentMBest.getTrainingOutputSD()[1]+
	      	  ","+currentMBest.getUnseenOutputSD()[0]+","+currentMBest.getUnseenOutputSD()[1]
	      				  
	      				  ); 
	      
	      
	      writer.flush();
	      writer.close();	      
	   
//	      writer.write("CurrentRun,Generation,NumGens,NumRuns,popSize,trainingTheta,UnseenTheta,"
//		      		+ "reconTrainError,reconUnseenError,depth,depthLimitOn,maxDepth,crossoverProb,k,ID,kunseen,size,E1TrainError,"
//		      		+ "E2TrainError,E1UnseenError,E2UnseenError,
	      //		+ "E1TrainOutRMS,E2TrainOutRMS,E1UnseenOutRMS,E2UnseenOutRMS,E1TrainSD,E2TrainSD,E1UnseenSD,E2UnseenSD"); 
   }
   //Overloaded method for printing results of first generation (0), because currentGeneration variable has been incremented already
   protected void output(int generation)throws IOException{

	      FileWriter writer = new FileWriter(file,true); 
	      // Writes the content to the file
	      writer.write("\n"+Main.CURRENTRUN+","+generation+","+Main.NUMBER_OF_GENERATIONS+","+Main.NUMBER_OF_RUNS+
	    		  ","+populationSize+","+currentMBest.getTrainingTheta()+","+currentMBest.getUnseenTheta()+","+currentMBest.getReconTrainingError()+
	    		  ","+currentMBest.getReconUnseenError()+
	    		  ","+currentMBest.getDepth() + ","+applyDepthLimit+","+maximumDepth+
	    		  ","+crossoverProbability+","+currentMBest.getK()+","+currentMBest.getId()+
	    		  ","+currentMBest.getKunseen()+","+currentMBest.getSize()+
	    		  ","+currentMBest.getProgram(0).getTrainingError()+ ","+currentMBest.getProgram(1).getTrainingError()+
	    		  ","+currentMBest.getProgram(0).getUnseenError()+ ","+currentMBest.getProgram(1).getUnseenError()+
	    		  ","+currentMBest.getTrainingOutputAvg()[0]+","+currentMBest.getTrainingOutputAvg()[1]+
	      		  ","+currentMBest.getUnseenOutputAvg()[0]+","+currentMBest.getUnseenOutputAvg()[1]+
	      	   	","+currentMBest.getTrainingOutputSD()[0]+","+currentMBest.getTrainingOutputSD()[1]+
	      	  ","+currentMBest.getUnseenOutputSD()[0]+","+currentMBest.getUnseenOutputSD()[1]
	      				  
	      				  ); 
	      writer.flush();
	      writer.close();	      
	   
}

	public void printBest(){
		 currentMBest.print();
	}
	public void writeBestToObjectFile(){
		 currentBest.writeToObjectFile();
	}
		// ##### get's and set's from here on #####


		public double getBestTrainingError(){
			return currentMBest.getReconTrainingError();
		}
		public double getBestUnseenError(){
			return currentMBest.getReconUnseenError();
		}
		public double getBestSize(){
			return currentMBest.getReconSize();
		}
		public double getBestDepth(){
			return currentMBest.getReconDepth();
		}
	public MIndividual getCurrentBest() {
		return currentMBest;
	}
	public double getMutationStep() {
		return mutationStep;
	}
	public File getOutputsFileName() {
		return OutputsFile;
	}
	public Population getReconPop(){
		
		return reconstructedPopulation;
	}	
	public void setReconPopulation() {
		reconstructedPopulation=new Population();
		for(int i=0;i<mpopulation.getSize();i++){
			
			reconstructedPopulation.addIndividual(mpopulation.getMIndividual(i).getReconstructedIndividual());
		}	
	}
	public void setMutationStep(double mutationStep) {
		this.mutationStep = mutationStep;
	}

	public void setBoundedMutation(boolean boundedMutation) {
		this.boundedMutation = boundedMutation;
	}

	public void setBuildIndividuals(boolean buildIndividuals) {
		this.buildIndividuals = buildIndividuals;
	}
}
