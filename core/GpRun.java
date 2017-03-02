package core;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Random;

import programElements.Addition;
import programElements.Constant;
import programElements.InputVariable;
import programElements.Multiplication;
import programElements.Operator;
import programElements.ProgramElement;
import programElements.ProtectedDivision;
import programElements.Subtraction;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class GpRun implements Serializable {

	private static final long serialVersionUID = 7L;

	// ##### parameters #####
	protected Data data;
	protected ArrayList<ProgramElement> functionSet, terminalSet, fullSet;
	protected int populationSize;
	protected boolean applyDepthLimit;
	protected int maximumDepth;
	protected double crossoverProbability;
	protected boolean printAtEachGeneration;
	protected double minUnseen;

	// ##### state #####
	protected Random randomGenerator;
	protected int currentGeneration;
	protected Population population;
	protected Individual currentBest;
	
	File file = new File(Main.resultsFolderName+"/"+Main.DATA_FILENAME+"/"+Main.modelName+"/results.txt");
	File OutputsFile = new File(Main.resultsFolderName+"/"+Main.DATA_FILENAME+"/"+"/"+Main.modelName+"/outputs.txt");
	
	public GpRun(Data data) throws IOException {
		this.data = data;
		initialize();
	}
	
	public GpRun(Data data,Population reconPop) throws IOException {
		this.data = data;
		initialize(reconPop);
	}

	protected void initialize() throws IOException {

		// adds all the functions to the function set
		functionSet = new ArrayList<ProgramElement>();
		functionSet.add(new Addition());
		functionSet.add(new Subtraction());
		functionSet.add(new Multiplication());
		functionSet.add(new ProtectedDivision());

		// adds all the constants to the terminal set
		terminalSet = new ArrayList<ProgramElement>();
		double[] constants = { -1.0, -0.75, -0.5, -0.25, 0.25, 0.5, 0.75, 1.0 };
		for (int i = 0; i < constants.length; i++) {
			terminalSet.add(new Constant(constants[i]));
		}

		// adds all the input variables to the terminal set
		for (int i = 0; i < data.getDimensionality(); i++) {
			terminalSet.add(new InputVariable(i));
		}

		// creates the set which contains all the program elements
		fullSet = new ArrayList<ProgramElement>();
		for (ProgramElement programElement : functionSet) {
			fullSet.add(programElement);
		}
		for (ProgramElement programElement : terminalSet) {
			fullSet.add(programElement);
		}

		populationSize = 100;
		applyDepthLimit = true;
		maximumDepth = 17;
		crossoverProbability = 0.0;
		printAtEachGeneration = true;

		randomGenerator = new Random();
		currentGeneration = 0;

		// initialize and evaluate population
		rampedHalfAndHalfInitialization();
		//else
		
		for (int i = 0; i < populationSize; i++) {
			population.getIndividual(i).evaluate(data);
		}
		if(Main.CURRENTRUN==1)createOutputFile();
		//Check added by KS for ESGSGP acceptable printing
		if(!(this instanceof EsgpRun)&&!(this instanceof EsgsgpRun)){	
			updateCurrentBest();
			minUnseen=currentBest.getUnseenError();
			currentBest.writeToObjectFile();
			printState();
			
			try {				
				output();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		currentGeneration++;
	}

	protected void initialize(Population reconPop) throws IOException {

		// adds all the functions to the function set
		functionSet = new ArrayList<ProgramElement>();
		functionSet.add(new Addition());
		functionSet.add(new Subtraction());
		functionSet.add(new Multiplication());
		functionSet.add(new ProtectedDivision());

		// adds all the constants to the terminal set
		terminalSet = new ArrayList<ProgramElement>();
		double[] constants = { -1.0, -0.75, -0.5, -0.25, 0.25, 0.5, 0.75, 1.0 };
		for (int i = 0; i < constants.length; i++) {
			terminalSet.add(new Constant(constants[i]));
		}

		// adds all the input variables to the terminal set
		for (int i = 0; i < data.getDimensionality(); i++) {
			terminalSet.add(new InputVariable(i));
		}

		// creates the set which contains all the program elements
		fullSet = new ArrayList<ProgramElement>();
		for (ProgramElement programElement : functionSet) {
			fullSet.add(programElement);
		}
		for (ProgramElement programElement : terminalSet) {
			fullSet.add(programElement);
		}

		populationSize = 100;
		applyDepthLimit = true;
		maximumDepth = 17;
		crossoverProbability = 0.0;
		printAtEachGeneration = true;

		randomGenerator = new Random();
		currentGeneration = 0;

		// initialize and evaluate population
		population=reconPop;
		System.out.println(population.individuals.size());
		System.out.println(population.getSize());
		//else
		
		for (int i = 0; i < population.getSize(); i++) {
			population.getIndividual(i).evaluate(data);
		}
		if(Main.CURRENTRUN==1)createOutputFile();
		//Check added by KS for ESGSGP acceptable printing
		if(!(this instanceof EsgpRun)&&!(this instanceof EsgsgpRun)){	
			updateCurrentBest();
			minUnseen=currentBest.getUnseenError();
			currentBest.writeToObjectFile();
			printState();
			
			try {				
				output();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		currentGeneration=4;
		
	}
	
	protected void rampedHalfAndHalfInitialization() {
		int maximumInitialDepth = 6;
		/*
		 * depth at the root node is 0. this implies that the number of
		 * different depths is equal to the maximumInitialDepth
		 */
		int individualsPerDepth = populationSize / maximumInitialDepth;
		int remainingIndividuals = populationSize % maximumInitialDepth;
		population = new Population();
		int fullIndividuals, growIndividuals;

		for (int depth = 1; depth <= maximumInitialDepth; depth++) {
			if (depth == maximumInitialDepth) {
				fullIndividuals = (int) Math.floor((individualsPerDepth + remainingIndividuals) / 2.0);
				growIndividuals = (int) Math.ceil((individualsPerDepth + remainingIndividuals) / 2.0);
			} else {
				fullIndividuals = (int) Math.floor(individualsPerDepth / 2.0);
				growIndividuals = (int) Math.ceil(individualsPerDepth / 2.0);
			}

			for (int i = 0; i < fullIndividuals; i++) {
				population.addIndividual(full(depth));
			}
			for (int i = 0; i < growIndividuals; i++) {
				population.addIndividual(grow(depth));
			}
		}
	}

	protected Individual full(int maximumTreeDepth) {
		Individual individual = new Individual();
		fullInner(individual, 0, maximumTreeDepth);
		individual.setDepth(maximumTreeDepth);
		return individual;
	}

	protected void fullInner(Individual individual, int currentDepth, int maximumTreeDepth) {
		if (currentDepth == maximumTreeDepth) {
			ProgramElement randomTerminal = terminalSet.get(randomGenerator.nextInt(terminalSet.size()));
			individual.addProgramElement(randomTerminal);
		} else {
			Operator randomOperator = (Operator) functionSet.get(randomGenerator.nextInt(functionSet.size()));
			individual.addProgramElement(randomOperator);
			for (int i = 0; i < randomOperator.getArity(); i++) {
				fullInner(individual, currentDepth + 1, maximumTreeDepth);
			}
		}
	}

	protected Individual grow(int maximumTreeDepth) {
		Individual individual = new Individual();
		growInner(individual, 0, maximumTreeDepth);
		individual.calculateDepth();
		return individual;
	}

	protected void growInner(Individual individual, int currentDepth, int maximumTreeDepth) {
		if (currentDepth == maximumTreeDepth) {
			ProgramElement randomTerminal = terminalSet.get(randomGenerator.nextInt(terminalSet.size()));
			individual.addProgramElement(randomTerminal);
		} else {
			// equal probability of adding a terminal or an operator
			if (randomGenerator.nextBoolean()) {
				Operator randomOperator = (Operator) functionSet.get(randomGenerator.nextInt(functionSet.size()));
				individual.addProgramElement(randomOperator);
				for (int i = 0; i < randomOperator.getArity(); i++) {
					growInner(individual, currentDepth + 1, maximumTreeDepth);
				}
			} else {
				ProgramElement randomTerminal = terminalSet.get(randomGenerator.nextInt(terminalSet.size()));
				individual.addProgramElement(randomTerminal);
			}
		}
	}

	public void evolve(int numberOfGenerations) throws IOException {

		// evolve for a given number of generations
		while (currentGeneration <= numberOfGenerations) {
			Population offspring = new Population();
			System.out.println("pop size before evolving "+ population.getSize());
			// generate a new offspring population
			while (offspring.getSize() < population.getSize()) {
				Individual p1, newIndividual;
				p1 = selectParent();
				// apply crossover
				if (randomGenerator.nextDouble() < crossoverProbability) {
					Individual p2 = selectParent();
					newIndividual = applyStandardCrossover(p1, p2);
				}
				// apply mutation
				else {
					newIndividual = applyStandardMutation(p1);
				}

				/*
				 * add the new individual to the offspring population if its
				 * depth is not higher than the maximum (applicable only if the
				 * depth limit is enabled)
				 */
				if (applyDepthLimit && newIndividual.getDepth() > maximumDepth) {
					newIndividual = p1;
				} else {
					newIndividual.evaluate(data);
				}
				offspring.addIndividual(newIndividual);
			}

			population = selectSurvivors(offspring);
			System.out.println("evolved ind size " +population.individuals.size());
			System.out.println("evolved pop size " +population.getSize());
			updateCurrentBest();
			currentBest.printVectors(currentGeneration, OutputsFile);
			//currentBest.printVectors(currentGeneration, OutputsFile);
			if(currentBest.getUnseenError()<minUnseen){
				minUnseen=currentBest.getUnseenError();
				currentBest.writeToObjectFile();
			}			
			printState();
			output();
			currentGeneration++;
		}
	}

	protected void printState() {
		if (printAtEachGeneration) {
			System.out.println("\nGeneration:\t\t" + currentGeneration);
			System.out.printf("Training error:\t\t%.2f\nUnseen error:\t\t%.2f\nSize:\t\t\t%d\nDepth:\t\t\t%d\n",
					currentBest.getTrainingError(), currentBest.getUnseenError(), currentBest.getSize(),
					currentBest.getDepth());
		}
	}

	// tournament selection
	protected Individual selectParent() {
		Population tournamentPopulation = new Population();
		int tournamentSize = (int) (0.05 * population.getSize());
		for (int i = 0; i < tournamentSize; i++) {
			int index = randomGenerator.nextInt(population.getSize());
			tournamentPopulation.addIndividual(population.getIndividual(index));
		}
		return tournamentPopulation.getBest();
	}

	protected Individual applyStandardCrossover(Individual p1, Individual p2) {

		int p1CrossoverStart = randomGenerator.nextInt(p1.getSize());
		int p1ElementsToEnd = p1.countElementsToEnd(p1CrossoverStart);
		int p2CrossoverStart = randomGenerator.nextInt(p2.getSize());
		int p2ElementsToEnd = p2.countElementsToEnd(p2CrossoverStart);

		Individual offspring = p1.selectiveDeepCopy(p1CrossoverStart, p1CrossoverStart + p1ElementsToEnd - 1);

		// add the selected tree from the second parent to the offspring
		for (int i = 0; i < p2ElementsToEnd; i++) {
			offspring.addProgramElementAtIndex(p2.getProgramElementAtIndex(p2CrossoverStart + i), p1CrossoverStart + i);
		}

		offspring.calculateDepth();
		return offspring;
	}

	protected Individual applyStandardMutation(Individual p) {

		int mutationPoint = randomGenerator.nextInt(p.getSize());
		int parentElementsToEnd = p.countElementsToEnd(mutationPoint);
		Individual offspring = p.selectiveDeepCopy(mutationPoint, mutationPoint + parentElementsToEnd - 1);
		int maximumDepth = 6;
		Individual randomTree = grow(maximumDepth);

		// add the random tree to the offspring
		for (int i = 0; i < randomTree.getSize(); i++) {
			offspring.addProgramElementAtIndex(randomTree.getProgramElementAtIndex(i), mutationPoint + i);
		}

		offspring.calculateDepth();
		return offspring;
	}

	// keep the best overall + all the remaining offsprings
	protected Population selectSurvivors(Population newIndividuals) {
		Population survivors = new Population();
		Individual bestParent = population.getBest();
		Individual bestNewIndividual = newIndividuals.getBest();
		Individual bestOverall;
		// the best overall is in the current population
		if (bestParent.getTrainingError() < bestNewIndividual.getTrainingError()) {
			bestOverall = bestParent;
		}
		// the best overall is in the offspring population
		else {
			bestOverall = bestNewIndividual;
		}

		survivors.addIndividual(bestOverall);

		for (int i = 0; i < newIndividuals.getSize(); i++) {
			if (newIndividuals.getIndividual(i).getId() != bestOverall.getId()) {
				survivors.addIndividual(newIndividuals.getIndividual(i));
			}
		}
		return survivors;
	}

	protected void updateCurrentBest() {
		currentBest = population.getBest();
	}

   protected void createOutputFile()throws IOException, IOException{
		
	      // creates the file
	      file.createNewFile();
	      OutputsFile.createNewFile();
	      // creates a FileWriter Object - results file
	      FileWriter writer = new FileWriter(file); 
	      // Writes the content to the file
	      writer.write("CurrentRun,Generation,NumGens,NumRuns,popSize,trainingTheta,UnseenTheta,"
		      		+ "reconTrainError,reconUnseenError,depth,depthLimitOn,maxDepth,crossoverProb,k,ID,kunseen,size,E1TrainError,"
		      		+ "E2TrainError,E1UnseenError,E2UnseenError,E1TrainOutAvg,E2TrainOutAvg,E1UnseenOutAvg,E2UnseenOutAvg,E1TrainOutSD,E2TrainOutSD,E1UnseenOutSD,E2UnseenOutSD"); 
	      writer.flush();
	      writer.close();
	   // creates a FileWriter Object
//	      FileWriter writer2 = new FileWriter(file2); 
//	      // Writes the content to the file
//	      writer2.write("CurrentRun,Generation,ID,p1,p2,P1Training,P1Unseen,P2Training,P2Unseen,trainingTheta,UnseenTheta,reconTrainingError,reconUnseenError,LowestDistance"); 
//	      writer2.flush();
//	      writer2.close();
	   // creates a FileWriter Object - outputs file
	      FileWriter writer3 = new FileWriter(OutputsFile); 
	      // Writes the content to the file
	      writer3.write("CurrentRun,ID,ExpressionNum,"); 
//	      for(int i=1;i<data.getTrainingData().length;i++){
//	    	  writer3.write("TrainingErrorVector"+i+",");
//	      }
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
	    		  ","+populationSize+","+""+","+""+","+currentBest.getTrainingError()+","+currentBest.getUnseenError()+
	    		  ","+currentBest.getDepth() + ","+applyDepthLimit+","+maximumDepth+
	    		  ","+crossoverProbability+","+""+","+currentBest.getId()+","+""+","+currentBest.getSize()+
	    		  ","+""+ ","+""+
	    		  ","+""+ ","+""+
	    		  ","+currentBest.getTrainingOutputAvg()[0]+","+""+
	      		  ","+currentBest.getUnseenOutputAvg()[0]+","+""+
	      	   	","+currentBest.getTrainingOutputSD()[0]+","+""+
	      	  ","+currentBest.getUnseenOutputSD()[0]+","+""); 
	      writer.flush();
	      writer.close();	      
	      
   }
   public void printBest(){
	   currentBest.print();
   }
   public void writeBestToObjectFile(){
	   currentBest.writeToObjectFile();
   }
	// ##### get's and set's from here on #####


	public double getBestTrainingError(){
		return currentBest.getTrainingError();
	}
	public double getBestUnseenError(){
		return currentBest.getUnseenError();
	}
	public double getBestSize(){
		return currentBest.getSize();
	}
	public double getBestDepth(){
		return currentBest.getDepth();
	}	
	public Individual getCurrentBest() {
		return currentBest;
	}

	public ArrayList<ProgramElement> getFunctionSet() {
		return functionSet;
	}

	public ArrayList<ProgramElement> getTerminalSet() {
		return terminalSet;
	}

	public ArrayList<ProgramElement> getFullSet() {
		return fullSet;
	}

	public boolean getApplyDepthLimit() {
		return applyDepthLimit;
	}

	public int getMaximumDepth() {
		return maximumDepth;
	}

	public double getCrossoverProbability() {
		return crossoverProbability;
	}

	public int getCurrentGeneration() {
		return currentGeneration;
	}

	public Data getData() {
		return data;
	}

	public Population getPopulation() {
		return population;
	}

	public int getPopulationSize() {
		return populationSize;
	}

	public Random getRandomGenerator() {
		return randomGenerator;
	}

	public boolean getPrintAtEachGeneration() {
		return printAtEachGeneration;
	}

	public void setFunctionSet(ArrayList<ProgramElement> functionSet) {
		this.functionSet = functionSet;
	}

	public void setTerminalSet(ArrayList<ProgramElement> terminalSet) {
		this.terminalSet = terminalSet;
	}

	public void setFullSet(ArrayList<ProgramElement> fullSet) {
		this.fullSet = fullSet;
	}

	public void setApplyDepthLimit(boolean applyDepthLimit) {
		this.applyDepthLimit = applyDepthLimit;
	}

	public void setMaximumDepth(int maximumDepth) {
		this.maximumDepth = maximumDepth;
	}

	public void setCrossoverProbability(double crossoverProbability) {
		this.crossoverProbability = crossoverProbability;
	}

	public void setPrintAtEachGeneration(boolean printAtEachGeneration) {
		this.printAtEachGeneration = printAtEachGeneration;
	}
}
