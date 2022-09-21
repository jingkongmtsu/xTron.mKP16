import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 * this class is used to perform CCSD calculation for the curve scanning by qchem
 * @author fenglai
 */
class CCSD {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

		// check
		if (args.length != 4) {
			System.out.println("Invalid input parameter for running the program");
			System.out.println("input should be mono1 name, mono2 name, basis set name and step length");
			System.exit(1);
		}

		// read
		String m1 = args[0];
		String m2 = args[1];
		String b    = args[2];
		double s  = Double.valueOf(args[3]);

		// build CCSD
		CCSD ccsd = new CCSD(m1,m2,b,s);
		ccsd.doCCSD();
	}

	//
	// constructor, initilize the member data
	//
	public CCSD(String mono_1, String mono_2, String basisName, double steplength) {
		basis     = basisName;
		mono1 = mono_1;
		mono2 = mono_2;
		step      = steplength;
		Emono1 = 0.0;
		Emono2 = 0.0;
		setupBeginEndDistance();
	}

	//
	// do CCSD calculation
	//
	public void doCCSD() {

		// mono1
		runQchem(0.0,"mono1");
		String line = "mono1 energy: " + Double.toString(Emono1);
		System.out.println(line);

		// mono1
		runQchem(0.0,"mono2");
		line = "mono2 energy: " + Double.toString(Emono2);
		System.out.println(line);

		// now for cluster
		for(double dist=dist1; dist<dist2; dist = dist + step) {
			double diff = runQchem(dist,"cluster");
			String distString   = String.format("%.3f", dist);
			String line1 =  distString  + "    " + Double.toString(diff);
			System.out.println(line1);
		}
	}

	//
	// set up the dist1 and dist2 etc. according to the monomers given
	//
	private void setupBeginEndDistance() {

		// lower the case
		String m1 = mono1.toLowerCase();
		String m2 = mono2.toLowerCase();

		// now let's set up data
		if (m1.equals("ne") && m2.equals("ne")) {
			dist1 = 2.0;
			dist2 = 9.0;
		}else if (m1.equals("he") && m2.equals("he")) {
			dist1 = 1.70;
			dist2 = 9.00;
		}else if (m1.equals("ar") && m2.equals("ar")) {
			dist1 = 3.00;
			dist2 = 11.00;
		}else {
			System.out.println("the input monomers are not supported for parameter optimization");
			System.out.println("please define the beginning distance as well as ending distance for the monomer");
			System.exit(1);
		}
	}


	//
	// run qchem for CCSD energy calculation
	//
	private double runQchem(double dist, String status) {

		// form input file name
		String input = null;
		String output = null;
		if (status.equals("mono1")) {
			input = "mono1.in";
			output = "mono1.out";
		}else if (status.equals("mono2")) {
			input = "mono2.in";
			output = "mono2.out";
		}else {
			String distName    = String.format("%.2f", dist);
			String body  = mono1 + "_" + mono2 + "_" + distName;
			input = body + ".in";
			output = body + ".out";
		}

		// generate input file
		inputGenerate(dist,input,status);

		// now run qchem program
		// qchem should be in the PATH
		try{
			Runtime t   = Runtime.getRuntime();
			String args = "qchem " + input + " " + output;
			Process proc0 = t.exec(args);
			proc0.waitFor();
		}catch (IOException e){
			System.out.println("Problem running qchem program.");
			System.out.println("Be sure to include qchem into the PATH env.");
			e.printStackTrace();
		}catch (InterruptedException err){
			System.out.println("Problem running qchem, internal error for the process");
			err.printStackTrace();
		}

		// now get the energy
		double energy = getEnergy(output);
		//System.out.println("the file name is ");
		//System.out.println(input);
		//System.out.println(Double.toString(energy));

		// now let's calculate the enegy difference
		if (status.equals("cluster")) {
			double diff = Math.abs(energy-Emono1-Emono2);
			diff = diff*27.212*23.061;
			return diff;
		}else if (status.equals("mono1")) {
			Emono1 = energy;
		}else if (status.equals("mono2")) {
			Emono2 = energy;
		}

		// return the energy
		return energy;
	}

	//
	// after running qchem program, we will try to get the energy data from the file
	//
	private double getEnergy(String out) {

		// now let's read in emul output file
		BufferedReader bufferReader = null;
		try {
			bufferReader = new BufferedReader(new FileReader(out));
		} catch (FileNotFoundException e) {
			System.out.println("the qchem output file can not be set up for reading");
			e.printStackTrace();
		}

		// get the energy value
		double energy = 0.0;
		boolean gotit = false;
		while(true) {

			// read in line
			String line = null;
			try {
				line = bufferReader.readLine();
			} catch (IOException e) {
				System.out.println("error in reading the content from qchem output file");
				e.printStackTrace();
			}

			// stop here
			if(line == null) break;

			// do we see the result?
			if (line.indexOf("CCSD total energy") >= 0 || line.indexOf("CCSD Total Energy")>=0) {
				String line2= line.trim();
				String [ ] contents = line2.split("\\s+");
				if (contents.length == 5) {
					energy = Double.valueOf(contents[4]);
					gotit = true;
					break;
				}else{
					System.out.println("unable to fatch the energy in the qchem output file: ");
					System.out.println(out);
					System.exit(1);
				}
			}
		}

		//
		// do we got the energy?
		//
		if (! gotit) {
			System.out.println("fail to get the energy in the qchem output file: ");
			System.out.println(out);
			System.exit(1);
		}

		// now close file
		try {
			bufferReader.close();
		} catch (IOException e) {
			System.out.println("error in closing qchem output file");
			e.printStackTrace();
		}

		// return it
		return energy;
	}


	//
	// this function will print out the input files based on the information
	// we have
	//
	private void inputGenerate(double dist, String outputFile, String status) {

		// the first section will be the geometry
		// currently we only have atom
		// open the output file
		FileWriter fileWriter = null;
		try {
			fileWriter = new FileWriter(outputFile);
		} catch (IOException e) {
			System.out.println("the input file for running with qchem can not be set up for writting");
			e.printStackTrace();
		}
		BufferedWriter bufferWriter = new BufferedWriter(fileWriter);

		// firstly, do the molecule section
		try {
			bufferWriter.newLine();
			bufferWriter.newLine();
			bufferWriter.newLine();
		    String line  = "$molecule";
			bufferWriter.write(line);
			bufferWriter.newLine();
			bufferWriter.write("0  1");
			bufferWriter.newLine();
			if (status.equals("mono1") || status.equals("cluster")) {
				line  = mono1 + " 0.0  0.0  0.0";
				bufferWriter.write(line);
				bufferWriter.newLine();
			}
			if (status.equals("mono2")) {
				line  = mono2 + " 0.0  0.0  0.0";
				bufferWriter.write(line);
				bufferWriter.newLine();
			}else if (status.equals("cluster")) {
				String distance = Double.toString(dist);
				line = mono2 + " 0.0  0.0 " + distance;
				bufferWriter.write(line);
				bufferWriter.newLine();
			}
			line = "$end";
			bufferWriter.write(line);
			bufferWriter.newLine();
			bufferWriter.newLine();
		} catch (IOException e) {
			System.out.println("error in write the qchem input file");
			e.printStackTrace();
		}

		// now close file
		try {
			bufferWriter.close();
			fileWriter.close();
		} catch (IOException e) {
			System.out.println("error in closing input file writing");
			e.printStackTrace();
		}

		// finally it's the key word section
		keywordPrint(outputFile);
	}


	//
	// this function will print out key words setting for running qchem program with CCSD
	//
	private void keywordPrint(String outputFile) {

		// open the output file in appending mode
		FileWriter fileWriter = null;
		try {
			fileWriter = new FileWriter(outputFile, true);
		} catch (IOException e) {
			System.out.println("the output file can not be set up for writting");
			e.printStackTrace();
		}
		BufferedWriter bufferWriter = new BufferedWriter(fileWriter);

		// now begin writing things into the output file
		try {

			// xcfunc section
			bufferWriter.newLine();
			String line  = "$rem";
			bufferWriter.write(line);
			bufferWriter.newLine();
			line = "jobtype                 SP";
			bufferWriter.write(line);
			bufferWriter.newLine();
			line = "basis                     " + basis;
			bufferWriter.write(line);
			bufferWriter.newLine();
			line = "exchange              HF";
			bufferWriter.write(line);
			bufferWriter.newLine();
			line = "correlation           CCSD";
			bufferWriter.write(line);
			bufferWriter.newLine();
			line = "symmetry             FALSE";
			bufferWriter.write(line);
			bufferWriter.newLine();
			line = "sym_ignore          TRUE";
			bufferWriter.write(line);
			bufferWriter.newLine();
			line = "$end";
			bufferWriter.write(line);
			bufferWriter.newLine();
			bufferWriter.newLine();

		} catch (IOException e) {
			System.out.println("error in writing the keywords into files");
			e.printStackTrace();
		}

		// now close the file
		try {
			bufferWriter.close();
			fileWriter.close();
		} catch (IOException e) {
			System.out.println("error in closing file after writing keywords");
			e.printStackTrace();
		}
	}

	// member data
	private  String basis;                  // basis set name
	private String mono1;                // this is the first monomer for ccsd
	private String mono2;               // this is the second monomer for  ccsd
	private double  dist1;               // starting distance between the two monomers
	private double dist2;                // the final distance between the two monomers
	private double step;                 // step length in the calclation
	private double Emono1;          // mono1 CCSD energy
	private double Emono2;          // mono2 CCSD energy

}
