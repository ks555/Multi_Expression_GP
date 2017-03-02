package utils;

import core.Individual;
import core.MIndividual;
import core.Main;

public class FinalEvaluation {

	public static void main(String[] args) {
		// read individual
		MIndividual finalIndividual = (MIndividual) Utils.readObjectFromFile("individuals/final_individual_1.obj");
		Individual reconIndividual =finalIndividual.getReconstructedIndividual();
		reconIndividual.print();

		// evaluate individual on final data
		// String evaluationData = "dataset.txt";
		String evaluationData = "Dataset_Energy_Buildings/test1.txt";
		double[][] data = Main.readData(evaluationData);
		double[] outputs = reconIndividual.evaluate(data);
		double finalError = reconIndividual.calculateRMSE(data,outputs);
		System.out.println("Calculated finalError: "+ finalError);
		System.out.println("Stored recon individual error: "+ reconIndividual.getUnseenError());
		System.out.println("Stored calculated recon Error: "+ finalIndividual.getReconUnseenError());
		
	}
}
