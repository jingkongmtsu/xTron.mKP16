import java.io.*;

/**
 * this class is used to perform curve scanning job
 * @author fenglai
 */
class ScanCurve {

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		// read in infor file
		String in = null;
		if (args.length == 1) {
			in = args[0];
		}else {
			System.out.println("Invalid input parameter for running the program");
			System.out.println("input should be only the infor file, and it must be provided");
			System.exit(1);
		}

		// create the object
		ScanCurve scan = new ScanCurve(in);
		scan.doScanning();
	}


	//
	// constructor
	// by given a input file, we will try to set up all of information
	//
	public ScanCurve(String infor) {

		// set up the default value if it's not defined in the infor file
		step        = 0.1;
		exName      = "pw86x";
		ecName      = "pbec";

		// set the following to wrong value
		// they should be read in from infor file
		dist1            = -1.0;
		dist2            = -1.0;
		xdm1             = -1.0;
		xdm2             = -1.0;

		// firstly, parse the infor file and set up the basic information
		inforParse(infor);

		// finally check the parameters
		if (xdm1 < 0 || xdm2 < 0) {
			System.out.println("invalid parameters setting, xdm paramters should be initialized");
			System.exit(1);
		}
		if (dist1 < 0 || dist2 < 0) {
			System.out.println("invalid parameters setting, distance begin/end should be initialized");
			System.exit(1);
		}
	}


	//
	// this is the driver routine to perform curve scanning
	//
	private void doScanning() {

		// perform scanning
		for(double dist=dist1; dist<=dist2; dist = dist + step) {
			double e = runEmul(dist);
			String d   = String.format("%.3f", dist);
			String line = d  + "  " + Double.toString(e);
			System.out.println(line);
		}
	}


	//
	// this function will read in and parse the input data for opt process
	//
	private void inforParse(String input) {

		// now let's read in the basis.txt content
		// now keyword section - open the file
		// we note that the name is fixed
		BufferedReader bufferReader = null;
		try {
			bufferReader = new BufferedReader(new FileReader(input));
		} catch (FileNotFoundException e) {
			System.out.println("the input file can not be set up for parsing the curve scanning");
			e.printStackTrace();
		}

		// read input, and parse
		while(true) {

			// read in line
			String line = null;
			try {
				line = bufferReader.readLine();
			} catch (IOException e) {
				System.out.println("error in reading the input file for parsing");
				e.printStackTrace();
			}

			// step out
			if(line == null) break;

			// comment line?
			if (line.indexOf("#") >= 0) continue;

			// empty line?
			if (line.isEmpty()) continue;

			// trim the line
			line = line.trim();

			// now let's split the string into a subset
			String [ ]  keyList = line.split("\\s+");
			if (keyList.length != 2) {
				System.out.println("the input line is not valid");
				System.out.println(line);
				System.out.println("please correct error inside");
				System.exit(1);
			}
			String key = keyList[0];
			String val = keyList[1];
			key.toLowerCase();

			// now let's see what it is
			if (key.equals("exchange")) {
				// exchange functional name
				exName = val;
				exName.toLowerCase();
			} else if (key.equals("correlation")) {
				// correlation functional name
				ecName = val;
				ecName.toLowerCase();
			} else if (key.equals("mono1")) {
				mono1 = val;
			} else if (key.equals("mono2")) {
				mono2 = val;
			} else if (key.equals("distance_step")) {
				step  = Double.valueOf(val);
			} else if (key.equals("xdm1")) {
				xdm1  = Double.valueOf(val);
			} else if (key.equals("xdm2")) {
				xdm2  = Double.valueOf(val);
			} else if (key.equals("distance_begin")) {
				dist1 = Double.valueOf(val);
			} else if (key.equals("distance_end")) {
				dist2 = Double.valueOf(val);
			}else {
				System.out.println(key);
				System.out.println(val);
				System.out.println("invalid line reading in, can not parse the keyword from this line");
				System.exit(1);
			}
		}

		// now close file
		try {
			bufferReader.close();
		} catch (IOException e) {
			System.out.println("error in closing input file for parsing");
			e.printStackTrace();
		}
	}


	//
	// run emul for single point energy calculation
	//
	private double runEmul(double dist) {

		// form input file name
		String distName    = String.format("%.2f", dist);
		String body  = mono1 + "_" + mono2 + "_" + distName;
		String input = body + ".in";
		String output = body + ".out";

		// generate input file
		// we have three sections, two monomers and cluster
		inputGenerate(dist,input,"mono1");
		inputGenerate(dist,input,"mono2");
		inputGenerate(dist,input,"cluster");

		// now run emul program
		// emul should be in the PATH
		try{
			Runtime t   = Runtime.getRuntime();
			String args = "emul " + input + " " + output;
			Process proc0 = t.exec(args);
			proc0.waitFor();
		}catch (IOException e){
			System.out.println("Problem running emul program.");
			System.out.println("Be sure to include emul into the PATH env.");
			e.printStackTrace();
		}catch (InterruptedException err){
			System.out.println("Problem running emul, internal error for the process");
			err.printStackTrace();
		}

		// now get the energy
		double energy = getEnergy(output);
		return energy;
		//System.out.println("the file name is ");
		//System.out.println(input);
		//System.out.println(Double.toString(energy));
	}


	//
	// after running emul program, we will try to get the energy data from the file
	//
	private double getEnergy(String out) {

		// now let's read in emul output file
		BufferedReader bufferReader = null;
		try {
			bufferReader = new BufferedReader(new FileReader(out));
		} catch (FileNotFoundException e) {
			System.out.println("the emul output file can not be set up for reading");
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
				System.out.println("error in reading the content from emul output file");
				e.printStackTrace();
			}

			// stop here
			if(line == null) break;

			// is it error we meet?
			if (line.indexOf("Fatal error occurs") >= 0) {
				System.out.println("fatal error in the emul output file: ");
				System.out.println(out);
				System.exit(1);
			}

			// do we see the result?
			if (line.indexOf("kcal") >= 0) {
				String [ ] contents = line.split("\\s+");
				if (contents.length >= 4) {
					energy = Double.valueOf(contents[4]);
					gotit = true;
					break;
				}else{
					System.out.println("unable to fatch the energy in the emul output file: ");
					System.out.println(out);
					System.exit(1);
				}
			}
		}

		//
		// do we got the energy?
		//
		if (! gotit) {
			System.out.println("fail to get the energy in the emul output file: ");
			System.out.println(out);
			System.exit(1);
		}

		// now close file
		try {
			bufferReader.close();
		} catch (IOException e) {
			System.out.println("error in closing emul output file");
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

		// do we do appending file, or start a new file?
		boolean doAppend = true;
		if (status.equals("mono1")) doAppend = false;

		// the first section will be the geometry
		// currently we only have atom
		// open the output file
		FileWriter fileWriter = null;
		try {
			if (doAppend) {
				fileWriter = new FileWriter(outputFile,true);
			}else{
				fileWriter = new FileWriter(outputFile);
			}
		} catch (IOException e) {
			System.out.println("the input file for running with emul can not be set up for writting");
			e.printStackTrace();
		}
		BufferedWriter bufferWriter = new BufferedWriter(fileWriter);

		// firstly, do the molecule section
		try {
			bufferWriter.newLine();
			bufferWriter.newLine();
			bufferWriter.newLine();
			String line  = "%molecule";
			if (status.equals("cluster")) {
				line  = "%cluster";
			}
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
			line = "%end";
			bufferWriter.write(line);
			bufferWriter.newLine();
			bufferWriter.newLine();
		} catch (IOException e) {
			System.out.println("error in write the emul input file");
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

		// next section is basis set file
		basisSetPrint(outputFile);

		// finally it's the key word section
		keywordPrint(outputFile,status);
	}

	//
	// this function will print out basis set section for running emul program
	//
	private void basisSetPrint(String outputFile) {

		// open the output file
		FileWriter fileWriter = null;
		try {
			fileWriter = new FileWriter(outputFile, true);
		} catch (IOException e) {
			System.out.println("the output file can not be set up for writting");
			e.printStackTrace();
		}
		BufferedWriter bufferWriter = new BufferedWriter(fileWriter);

		// now let's read in the basis.txt content
		// now keyword section - open the file
		// we note that the name is fixed
		BufferedReader bufferReader = null;
		try {
			bufferReader = new BufferedReader(new FileReader("basis.txt"));
		} catch (FileNotFoundException e) {
			System.out.println("the template basis set file can not be set up for reading");
			e.printStackTrace();
		}

		// append the file to the input
		while(true) {

			// read in line
			String line = null;
			try {
				line = bufferReader.readLine();
			} catch (IOException e) {
				System.out.println("error in reading the basis set files");
				e.printStackTrace();
			}

			// print it to output
			if(line == null) break;
			try {
				bufferWriter.write(line);
				bufferWriter.newLine();
			} catch (IOException e) {
				System.out.println("error in write the basis set data");
				e.printStackTrace();
			}
		}

		// now close file
		try {
			bufferWriter.close();
			bufferReader.close();
			fileWriter.close();
		} catch (IOException e) {
			System.out.println("error in closing basis set files");
			e.printStackTrace();
		}
	}

	//
	// this function will print out key words setting for running emul program
	//
	private void keywordPrint( String outputFile, String status) {

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
			String line  = "%xcfunc";
			bufferWriter.write(line);
			bufferWriter.newLine();
			line = "name " + exName + " " + ecName;
			bufferWriter.write(line);
			bufferWriter.newLine();
			line = "%end";
			bufferWriter.write(line);
			bufferWriter.newLine();
			bufferWriter.newLine();

			// xcints section
			bufferWriter.newLine();
			line = "%xcints";
			bufferWriter.write(line);
			bufferWriter.newLine();
			line = "grid_points  128  302";
			bufferWriter.write(line);
			bufferWriter.newLine();
			line = "%end";
			bufferWriter.write(line);
			bufferWriter.newLine();
			bufferWriter.newLine();

			// scf section
			bufferWriter.newLine();
			line = "%scf";
			bufferWriter.write(line);
			bufferWriter.newLine();
			line = "keep_result_in_memory = true";
			bufferWriter.write(line);
			bufferWriter.newLine();
			if (status.equals("cluster")) {
				line = "use_perturbed_vdw        = true";
				bufferWriter.write(line);
				bufferWriter.newLine();
				line = "scf_method                      = FRAGMENT_PERTURBED";
				bufferWriter.write(line);
				bufferWriter.newLine();
			}
			line = "%end";
			bufferWriter.write(line);
			bufferWriter.newLine();
			bufferWriter.newLine();

			// xdm section
			// only cluster need this
			bufferWriter.newLine();
			if (status.equals("cluster")) {
				line = "%xdm";
				bufferWriter.write(line);
				bufferWriter.newLine();
				line = "parameters  " + Double.toString(xdm1) + " " +  Double.toString(xdm2);
				bufferWriter.write(line);
				bufferWriter.newLine();
				line = "%end";
				bufferWriter.write(line);
				bufferWriter.newLine();
				bufferWriter.newLine();
			}

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


	//
	// member data for curve scanning
	//
	private String mono1;                // this is the first monomer for opt
	private String mono2;                // this is the second monomer for  opt
	private double  dist1;                // starting distance between the two monomers
	private double dist2;                 // the final distance between the two monomers
	private double step;                  // step length in the calclation
	private double xdm1;                //   the result XDM parameter 1
	private double xdm2;                //   the result XDM parameter 2

	//
	// functional name
	//
	private String exName;             // exchange functional name
	private String ecName;             // correlation functional name
}

