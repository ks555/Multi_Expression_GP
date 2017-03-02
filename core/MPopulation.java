package core;


import java.util.ArrayList;

public class MPopulation extends Population {

	private static final long serialVersionUID = 7L;

	protected ArrayList<MIndividual> mindividuals;

	public MPopulation() {
		mindividuals = new ArrayList<MIndividual>();		
		
	}

	public MIndividual getMostDiverseM() {
		//System.out.println(mindividuals.get(getMostDiverseIndex()).getId());

		return mindividuals.get(getMostDiverseIndex());
	}
	
	public int getMostDiverseIndex() {
		int bestIndex = 0;
		double bestSDEx1 = mindividuals.get(bestIndex).getTrainingOutputSD(0);
		double bestSDEx2 = mindividuals.get(bestIndex).getTrainingOutputSD(1);
		for (int i = 1; i < mindividuals.size(); i++) {
			if (mindividuals.get(i).getTrainingOutputSD(0) > bestSDEx1&&mindividuals.get(i).getTrainingOutputSD(1) > bestSDEx2) {
				bestSDEx1 = mindividuals.get(i).getTrainingOutputSD(0);
				bestSDEx2 = mindividuals.get(i).getTrainingOutputSD(1);
				bestIndex = i;
			}
		}
		
		return bestIndex;
	}
	
	public MIndividual getHighestK() {
		return mindividuals.get(getHighestKIndex());
	}
	
	public int getHighestKIndex() {
		int bestIndex = 0;
		double bestK = mindividuals.get(bestIndex).getK();
	
		for (int i = 1; i < mindividuals.size(); i++) {
			if (mindividuals.get(i).getK() > bestK) {
				bestK = mindividuals.get(i).getK();
				bestIndex = i;
			}
		}
		return bestIndex;
	}
	
	public MIndividual getLowestErrorM() {
		return mindividuals.get(getLowestErrorIndex());
	}

	public int getLowestErrorIndex() {
		int bestIndex = 0;
		double bestTrainingError;
		bestTrainingError=mindividuals.get(bestIndex).getTrainingErrorSum();
		for (int i = 0; i < individuals.size(); i++) {
			if (mindividuals.get(i).getTrainingErrorSum() < bestTrainingError) {
				bestTrainingError = mindividuals.get(i).getTrainingErrorSum();
				bestIndex = i;
			}
		}
		return bestIndex;
	}	
	
	public MIndividual getLowestAngleM() {
		
		return mindividuals.get(getLowestAngleIndex());
	}
	

	public int getLowestAngleIndex() {
		int bestIndex = 0;
		double bestTrainingTheta = mindividuals.get(bestIndex).getTrainingTheta();
		for (int i = 1; i < mindividuals.size(); i++) {
			if (mindividuals.get(i).getTrainingTheta() < bestTrainingTheta) {
				bestTrainingTheta = mindividuals.get(i).getTrainingTheta();
				//System.out.println("get best: id " + mindividuals.get(i).getId() + " training theta " + mindividuals.get(i).getTrainingTheta() );
				bestIndex = i;
			}
		}
		return bestIndex;
	}
	public MIndividual getBestM() {
		
		return mindividuals.get(getBestIndex());
	}
	
	public MIndividual getBest() {
		
		return mindividuals.get(getBestIndex());
	}
	
//for using reconstructed error as fitness
	public int getBestIndex() {
		int bestIndex = 0;
		double bestReconTrainingError = mindividuals.get(bestIndex).getReconTrainingError();
		for (int i = 1; i < mindividuals.size(); i++) {
			if (mindividuals.get(i).getReconTrainingError()< bestReconTrainingError) {
				bestReconTrainingError = mindividuals.get(i).getReconTrainingError();
				//System.out.println("get best: id " + mindividuals.get(i).getId() + " recon error " + mindividuals.get(i).getTrainingTheta() );

				bestIndex = i;
			}
		}
		return bestIndex;
	}

//	public int getBestIndex() {
//		int bestIndex = 0;
//		double bestTrainingTheta = mindividuals.get(bestIndex).getTrainingTheta();
//		for (int i = 1; i < mindividuals.size(); i++) {
//			if (mindividuals.get(i).getTrainingTheta() < bestTrainingTheta) {
//				bestTrainingTheta = mindividuals.get(i).getTrainingTheta();
//				//System.out.println("get best: id " + mindividuals.get(i).getId() + " training theta " + mindividuals.get(i).getTrainingTheta() );
//				bestIndex = i;
//			}
//		}
//		return bestIndex;
//	}

	public void addIndividual(MIndividual mindividual) {
		mindividuals.add(mindividual);
	}
	
//	public void removeIndividual(int index) {
//		individuals.remove(index);
//	}
//
	public int getSize() {
		return mindividuals.size();
	}

	public MIndividual getMIndividual(int index) {
		return mindividuals.get(index);
	}

}
