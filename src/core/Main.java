package core;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import utils.Utils;
import java.math.*;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;;


public class Main {
	
	//use for files split into runs (test1.txt, train1.txt etc.)
	public static final String[] DATA_FILENAMES = {"Dataset_bioav/","Dataset_cemento/","Dataset_Dummy/","Dataset_Energy_Buildings/","Dataset_instanbul_stock/","Dataset_Kj4/","Dataset_Parkinson_motor/","Dataset_Parkinson_total/","Dataset_PPB/","Dataset_Tox/","Dataset_Vl1/"};
//	public static final String[] DATA_FILENAMES = {"Dataset_Dummy/","Dataset_Energy_Buildings/","Dataset_instanbul_stock/","Dataset_Kj4/","Dataset_Parkinson_motor/","Dataset_Parkinson_total/","Dataset_PPB/","Dataset_Tox/","Dataset_Vl1/"};
//	public static final String[] DATA_FILENAMES = {"Dataset_Proteine_3d_structure/"};
//	public static final String[] DATA_FILENAMES = {"Dataset_Kj4/","Dataset_Energy_Buildings/","Dataset_instanbul_stock/","Dataset_PPB/","Dataset_Tox/"};
	public static String DATA_FILENAME;
	public static String modelName="gp_2_errorFitness";
//	public static String resultsFolderName="results_Jan_2017";
	public static String resultsFolderName="results";
	public static final boolean SHUFFLE_AND_SPLIT = false;
	public static final boolean INITIALIZE_POP = true;

//	//data in one file, random split for each run (bioav.txt etc.)
//	public static final String DATA_FILENAME = "bioav";
//	public static final boolean SHUFFLE_AND_SPLIT = true;
	
	public static final int NUMBER_OF_RUNS = 3;
	public static final int NUMBER_OF_GENERATIONS = 100;
	public static int CURRENTRUN;	

	
	public static void main(String[] args) throws IOException {

		
	System.out.println("Running master");


	for(int j=0;j<DATA_FILENAMES.length;j++){
		DATA_FILENAME=DATA_FILENAMES[j];
		System.out.println(DATA_FILENAMES[j]);
		// if the directory for this dataset and model does not exist, create it
		File theDir = new File(resultsFolderName+"/"+DATA_FILENAME+"/"+modelName+"/individuals");
		if (!theDir.exists()) {
		    System.out.println("creating directory: " + theDir.getName());
		    boolean result = false;
		
		    try{
		        theDir.mkdirs();
		        result = true;
		    } 
		    catch(SecurityException se){
		        //handle it
		    }        
		    if(result) {    
		        System.out.println("DIR created");  
		    }
		}	
		// run current model for given number of runs
		double[][] resultsPerRun = new double[4][NUMBER_OF_RUNS];
		for (int i = 0; i < NUMBER_OF_RUNS; i++) {
			CURRENTRUN=i+1;
			System.out.println("current run: "+CURRENTRUN);
			Data data = loadData(DATA_FILENAME);
				
			System.out.printf("\n\t\t##### Run %d #####\n", i + 1);
			System.out.println("current run: "+CURRENTRUN);
			//GpRun gp = new GpRun(data);			
			//GsgpRun gp = new GsgpRun(data);
			//EsgsgpRun gp = new EsgsgpRun(data);
			EsgpRun gp = new EsgpRun(data);			
			//CHANGE OUTPUT FOLDER NAME ABOVE!!!!!
			
			//if regular, update results per run
			if(INITIALIZE_POP){
				gp.evolve(50);
				Population evolvedPop=gp.getReconPop();
				//create new gsgp instance - with pop as parameter, need overloaded constructor
				GsgpRun newgp = new GsgpRun(data,evolvedPop,50);
				newgp.evolve(NUMBER_OF_GENERATIONS-50);
				resultsPerRun[0][i] = newgp.getBestTrainingError();
				resultsPerRun[1][i] = newgp.getBestUnseenError();
				resultsPerRun[2][i] = newgp.getBestSize();
				resultsPerRun[3][i] = newgp.getBestDepth();
				System.out.println("\nBest");
				gp.printBest();
				System.out.println();			
				System.out.println();
			}
			else{				
				gp.evolve(NUMBER_OF_GENERATIONS);
				resultsPerRun[0][i] = gp.getBestTrainingError();
				resultsPerRun[1][i] = gp.getBestUnseenError();
				resultsPerRun[2][i] = gp.getBestSize();
				resultsPerRun[3][i] = gp.getBestDepth();
				System.out.println("\nBest");
				gp.printBest();
				System.out.println();			
				System.out.println();
			// write individual to object file
			//not doing this, writing to file during evolution if unseen error not increasing
			//gp.writeBestToObjectFile();
			}
		}	
			
		 //present average results
		System.out.printf("\n\t\t##### Results #####\n\n");
		System.out.printf("Average training error:\t\t%.2f\n", Utils.getAverage(resultsPerRun[0]));
		System.out.printf("Average unseen error:\t\t%.2f\n", Utils.getAverage(resultsPerRun[1]));
		System.out.printf("Average size:\t\t\t%.2f\n", Utils.getAverage(resultsPerRun[2]));
		System.out.printf("Average depth:\t\t\t%.2f\n", Utils.getAverage(resultsPerRun[3]));
		}
	}
	
	public static Data loadData(String dataFilename) {
		double[][] trainingData, unseenData;

		if (SHUFFLE_AND_SPLIT) {
			double[][] allData = readData(dataFilename + ".txt");
			List<Integer> instances = Utils.shuffleInstances(allData.length);
			int trainingInstances = (int) Math.floor(0.7 * allData.length);
			int unseenInstances = (int) Math.ceil(0.3 * allData.length);

			trainingData = new double[trainingInstances][];
			unseenData = new double[unseenInstances][];

			for (int i = 0; i < trainingInstances; i++) {
				trainingData[i] = allData[instances.get(i)];
			}

			for (int i = 0; i < unseenInstances; i++) {
				unseenData[i] = allData[instances.get(trainingInstances + i)];
			}
		} else {
			trainingData = readData(dataFilename+"train"+Main.CURRENTRUN +".txt");
			unseenData = readData(dataFilename +"test"+Main.CURRENTRUN +".txt");
		}
		return new Data(trainingData, unseenData);
	}

	public static double[][] readData(String filename) {
		double[][] data = null;
		List<String> allLines = new ArrayList<String>();
		try {
			BufferedReader inputBuffer = new BufferedReader(new FileReader(filename));
			String line = inputBuffer.readLine();
			while (line != null) {
				allLines.add(line);
				line = inputBuffer.readLine();
			}
			inputBuffer.close();
		} catch (Exception e) {
			System.out.println(e);
		}

		StringTokenizer tokens = new StringTokenizer(allLines.get(2).trim());
		int numberOfColumns = tokens.countTokens();
		data = new double[allLines.size()-2][numberOfColumns];
		System.out.println("columns " + numberOfColumns );
		System.out.println("rows " + data.length );
		int index=0;
		for (int i = 2; i < data.length; i++) {
			tokens = new StringTokenizer(allLines.get(i).trim());
			for (int k = 0; k < numberOfColumns; k++) {
				data[index][k] = Double.parseDouble(tokens.nextToken().trim());							
			}
			index++;
		}
		return data;
	}
}
