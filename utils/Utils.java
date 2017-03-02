package utils;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStream;
import java.io.ObjectInput;
import java.io.ObjectInputStream;
import java.io.ObjectOutput;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.StringTokenizer;

import core.Data;
import core.Main;

public class Utils {

	public static Data loadData(String dataFilename) {
//		double[][] trainingData = Utils.readData(dataFilename+"_training.txt");
//		double[][] unseenData = Utils.readData(dataFilename +"_unseen.txt");
		double[][] trainingData = Utils.readData(dataFilename+"train"+Main.CURRENTRUN +".txt");
		double[][] unseenData = Utils.readData(dataFilename +"test"+Main.CURRENTRUN +".txt");
//		double[][] trainingData = Utils.readData(dataFilename+"train3"+".txt");
//		double[][] unseenData = Utils.readData(dataFilename);
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
		System.out.println(numberOfColumns);
		data = new double[allLines.size()][numberOfColumns];
		for (int i = 2; i < data.length; i++) {
			tokens = new StringTokenizer(allLines.get(i).trim());
			
			for (int k = 0; k < numberOfColumns; k++) {				
				data[i][k] = Double.parseDouble(tokens.nextToken().trim());
				if (k==0){
					//System.out.println(data[i][k]);
				}
			}
		}
		return data;
	}
	
	public static void writeObjectToFile(Object object, String filename) {
		try {
			OutputStream outputStream = new BufferedOutputStream(new FileOutputStream(filename));
			ObjectOutput objectOutput = new ObjectOutputStream(outputStream);
			try {
				objectOutput.writeObject(object);
			} finally {
				objectOutput.close();
			}
		} catch (Exception e) {
			System.out.println(e);
			System.exit(0);
		}
	}
	
	public static Object readObjectFromFile(String filename) {
		Object object = null;
		try {
			InputStream inputStream = new BufferedInputStream(new FileInputStream(filename));
			ObjectInput objectInput = new ObjectInputStream(inputStream);
			try {
				object = objectInput.readObject();
			} finally {
				objectInput.close();
			}
		} catch (Exception e) {
			System.out.println(e);
			System.exit(0);
		}
		return object;
	}
	
	public static List<Integer> shuffleInstances(int end) {
		List<Integer> list = new ArrayList<Integer>();
		for (int i = 0; i < end; i++) {
			list.add(i);
		}
		Collections.shuffle(list);
		return list;
	}
	public static double getAverage(double[] values) {
		double sum = 0.0;
		for (int i = 0; i < values.length; i++) {
			sum += values[i];
		}
		return sum / values.length;
	}

	public static double logisticFunction(double x) {
		return 1.0 / (1.0 + Math.exp(-x));
	}
}
