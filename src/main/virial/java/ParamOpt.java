import java.io.*;

/**
 * this class is used to perform parameter optimization job
 * currently we only run it with XDM parameter optimization
 * @author fenglai
 */
class ParamOpt {

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
		ParamOpt paramOpt = new ParamOpt(in);
		if (paramOpt.inScanning()) {
			paramOpt.xdmScanning();
		}else{
			paramOpt.xdmSinglePointCal();
		}

		// output the result
		System.out.println("best XDM parameter 1                    :  " + paramOpt.getXDM1());
		System.out.println("best XDM parameter 2                    :  " + paramOpt.getXDM2());
		System.out.println("error for the best XDM parameter:  " + paramOpt.getError());
	}


	//
	// constructor
	// by given a input file, we will try to set up all of information
	//
	public ParamOpt(String infor) {

		// set up the default value if it's not defined in the infor file
		optMethod = "single_point";
		step             = 0.1;
		exName      = "pw86x";
		ecName      = "pbec";
		xdmStep    = 0.01;
		standardRatio = 0.02;
		cleaningFiles  = true;

		// set the following to wrong value
		// they should be read in from infor file
		dist1            = -1.0;
		dist2            = -1.0;
		xdm1Begin = -1.0;
		xdm1End    = -1.0;
		xdm2Begin = -1.0;
		xdm2End    = -1.0;

		// initilize the result
		finalXDM1  = 0.0;
		finalXDM2  = 0.0;
		rmsError     = 10000.0;

		// firstly, parse the infor file and set up the basic information
		inforParse(infor);

		// now set up the begin and end distance in terms of the monomer
		setupBeginEndDistance();

		// set up the energy and R for bottom of curve
		setupBottomData();

		// finally check
		if (inScanning()) {

			// distance should be also properly set
			if (dist1<0 || dist2<0) {
				System.out.println("the distance range value is not properly set");
				System.exit(1);
			}
		}

		// check XDM range
		// xdm range should be properly set
		if (xdm1Begin < 0 || xdm1End < 0 || xdm2Begin < 0 || xdm2End < 0) {
			System.out.println("the xdm range value is not properly set");
			System.exit(1);
		}
	}

	//
	// do we use the method of scanning?
	//
	public boolean inScanning() {
		if (optMethod.equals("scanning")) return true;
		return false;
	}

	//
	// do we use the method of scanning?
	//
	public boolean singlePointCal() {
		if (optMethod.equals("single_point")) return true;
		return false;
	}

	//
	// get the result error
	//
	public String getError() {
		return Double.toString(rmsError);
	}

	//
	// get the final xdm paramter 1
	//
	public String getXDM1() {
		return Double.toString(finalXDM1);
	}

	//
	// get the final xdm paramter 2
	//
	public String getXDM2() {
		return Double.toString(finalXDM2);
	}

	//
	// this is the driver routine to perform optimization based on xddm paramters
	// by comparing with bottom energy
	//
	public void xdmSinglePointCal() {
		for(double xdm1=xdm1Begin; xdm1<=xdm1End; xdm1 = xdm1 + xdmStep) {
			for(double xdm2=xdm2Begin; xdm2<=xdm2End; xdm2 = xdm2 + xdmStep) {

				// perform single point calculation
				double dist = bottomR;
				double err = runEmul(xdm1,xdm2,dist);

				// print out the value
				String xdm1String   = String.format("%.3f", xdm1);
				String xdm2String   = String.format("%.3f", xdm2);
				String line = "xdm1: " + xdm1String + " xdm2: " + xdm2String  + " abs error: " + Double.toString(err);
				System.out.println(line);

				// now let's refresh the result
				if (err<rmsError) {
					rmsError = err;
					finalXDM1 = xdm1;
					finalXDM2 = xdm2;
				}
			}
		}
	}


	//
	// this is the driver routine to perform optimization based on xdm parameters
	//
	public void xdmScanning() {
		for(double xdm1=xdm1Begin; xdm1<=xdm1End; xdm1 = xdm1 + xdmStep) {
			for(double xdm2=xdm2Begin; xdm2<=xdm2End; xdm2 = xdm2 + xdmStep) {

				// perform scanning
				double rms = 0.0;
				int npoints  = 0;
				for(double dist=dist1; dist<=dist2; dist = dist + step) {
					double err = runEmul(xdm1,xdm2,dist);
					if (err<0) {

						// in this case, it means that we can stop the scanning now
						if (dist>bottomR) {
							break;
						}else{
							continue;
						}
					}
					rms = rms + err*err;
					npoints++;
				}

				// now let's get the real value
				rms = rms/npoints;
				rms = Math.sqrt(rms);

				// print out the value
				String xdm1String   = String.format("%.3f", xdm1);
				String xdm2String   = String.format("%.3f", xdm2);
				String line = "xdm1: " + xdm1String + " xdm2: " + xdm2String  + " rms error: " + Double.toString(rms);
				System.out.println(line);

				// now let's refresh the result
				if (rms<rmsError) {
					rmsError = rms;
					finalXDM1 = xdm1;
					finalXDM2 = xdm2;
				}
			}
		}
	}


	//
	// set up the bottom R and bottom E in terms of the monomers given
	// energry is in kcal/mol
	// distance is in angstrom
	//
	// the data here is either taken from experiments, or from high level
	// QC calculation (CCSD with large basis sets)
	//
	private void setupBottomData() {
		String m1 = mono1.toLowerCase();
		String m2 = mono2.toLowerCase();
		if (m1.equals("ne") && m2.equals("ne")) {
			bottomE = -0.0840862875;
			bottomR = 3.09;
		}else if (m1.equals("he") && m2.equals("he")) {
			bottomE = 0.02183733436;
			bottomR = 2.97;
		}else if (m1.equals("ar") && m2.equals("ar")) {
			bottomE = -0.285;
			bottomR = 3.76;
		}else {
			System.out.println("the input monomers are not supported for parameter optimization");
			System.out.println("please define the bottom energy as well as bottom distance for the monomer");
			System.exit(1);
		}
	}

	//
	// set up the dist1 and dist2 etc. according to the monomers given
	//
	private void setupBeginEndDistance() {

		// lower the case
		String m1 = mono1.toLowerCase();
		String m2 = mono2.toLowerCase();

		// check that whether we use curve scanning method
		// if not, we do not need to do it
		if (singlePointCal()) return;

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
			System.out.println("the input file can not be set up for parsing the optimization process");
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
			} else if (key.equals("opt")) {
				optMethod = val;
				optMethod.toLowerCase();

				// additional check for opt method
				if (! optMethod.equals("single_point") && ! optMethod.equals("scanning")) {
					System.out.println("the input opt method is invalid");
					System.out.println("the possible choices are: single_point or scanning");
					System.exit(1);
				}
			} else if (key.equals("distance_step")) {
				step  = Double.valueOf(val);
			} else if (key.equals("xdm_p1_begin")) {
				xdm1Begin = Double.valueOf(val);
			} else if (key.equals("xdm_p2_begin")) {
				xdm2Begin = Double.valueOf(val);
			} else if (key.equals("xdm_p1_end")) {
				xdm1End    = Double.valueOf(val);
			} else if (key.equals("xdm_p2_end")) {
				xdm2End   = Double.valueOf(val);
			} else if (key.equals("xdm_step")) {
				xdmStep   = Double.valueOf(val);
			} else if (key.equals("standard_ratio")) {
				standardRatio   = Double.valueOf(val);
			} else if (key.equals("clean_file")) {
				String v = val.toLowerCase();
				if (v.equals("true")) {
					cleaningFiles = true;
				}else if (v.equals("false")) {
					cleaningFiles = false;
				}else {
					System.out.println("invalid paramter set for clean_file, only true or false allowed");
					System.exit(1);
				}
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
	private double runEmul(double xdm1, double xdm2, double dist) {

		// form input file name
		String distName    = String.format("%.2f", dist);
		String xdm1Name = String.format("%.2f", xdm1);
		String xdm2Name = String.format("%.2f", xdm2);
		String body  = mono1 + "_" + mono2 + "_" + distName + "_" + xdm1Name + "_" + xdm2Name;
		String input = body + ".in";
		String output = body + ".out";

		// generate input file
		// we have three sections, two monomers and cluster
		inputGenerate(dist,xdm1,xdm2,input,"mono1");
		inputGenerate(dist,xdm1,xdm2,input,"mono2");
		inputGenerate(dist,xdm1,xdm2,input,"cluster");

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
		//System.out.println("the file name is ");
		//System.out.println(input);
		//System.out.println(Double.toString(energy));

		// if it's on scanning, we need to figure out whether we account in
		// this energy value
		double error = 0.0;
		if (inScanning()) {
			if (dist>=bottomR) {

				// for the curve scanning in front of the bottom point,
				// the important thing is to check whether it's large than 0
				// if so, we return some negative value to indicate that
				// this value is not useful
				if (energy >0) error =  -1.0;

			}else{

				// firstly, if energy >0, we do not need it
				if (energy >0) error =   -1.0;

				// for the curve scanning after the bottom point,
				// we just want to see whether we should stop here
				double ratio = Math.abs(energy)/bottomE;
				if (ratio<standardRatio) error =  -1.0;
			}
		}

		// compare the energy with the standard one
		double standEnerrgy = getStandardEnergy(dist);

		// now let's calculate the error for this value
		error = Math.abs(standEnerrgy-energy);

		// do we clean the files?
		if (cleaningFiles) {
			File inf    = new File(input);
			File outf = new File(output);
			inf.delete();
			outf.delete();
		}

		// return the abs error
		return error;
	}


	//
	// get the standard energy result either from bottom data, or from the CCSD curve
	//
	private double getStandardEnergy(double dist) {

		// check the opt method
		// for single point calculation, we only do the bottom E
		double energy = bottomE;
		if (singlePointCal()) {
			return energy;
		}

		// now let's read in ccsd file
		String ccsdFile = "ccsd.txt";
		BufferedReader bufferReader = null;
		try {
			bufferReader = new BufferedReader(new FileReader(ccsdFile));
		} catch (FileNotFoundException e) {
			System.out.println("the CCSD data file can not be set up for reading");
			e.printStackTrace();
		}

		// get the energy value
		boolean gotit = false;
		while(true) {

			// read in line
			String line = null;
			try {
				line = bufferReader.readLine();
			} catch (IOException e) {
				System.out.println("error in reading the data from CCSD file");
				e.printStackTrace();
			}

			// stop here
			if(line == null) break;

			// comment line?
			if (line.indexOf("#") >= 0) continue;

			// empty line?
			if (line.isEmpty()) continue;

			// do we see the result?
			String [ ] contents = line.split("\\s+");
			if (contents.length == 2) {
				double d = Double.valueOf(contents[0]);
				double e = Double.valueOf(contents[1]);
				if (Math.abs(d-dist) <= 0.00001) {
					energy = e;
					gotit = true;
					break;
				}
			}else{
				System.out.println("wrong file format for CCSD data ");
				System.exit(1);
			}
		}

		//
		// do we got the energy?
		//
		if (! gotit) {
			System.out.println("fail to get the energy in CCSD file: ");
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
	private void inputGenerate(double dist, double xdm1, double xdm2, String outputFile, String status) {

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
		keywordPrint(xdm1,xdm2,outputFile,status);
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
	private void keywordPrint(double xdm1, double xdm2, String outputFile, String status) {

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
	// member data
	// for parameter optimization, we only need two monomers
	// these member data are constant during the calculation
	//
	private String mono1;                // this is the first monomer for opt
	private String mono2;               // this is the second monomer for  opt
	private String optMethod;       // which opt method we use, scanning the curve or single point calculation?
	private double  dist1;               // starting distance between the two monomers
	private double dist2;                // the final distance between the two monomers
	private double step;                 // step length in the calclation
	private double bottomR;         // the distance between two monomers for bottom energy
	private double bottomE;         // bottom energy
	private double standardRatio;  // this value is the limit for abs(E)/bottomE. if it's smaller than this value,
	                                                          // it indicates that E is near 0 and the E is useless for Virial calculation

	//
	// functional name
	//
	private String exName;             // exchange functional name
	private String ecName;             // correlation functional name

	//
	// optimization data
	// in terms of XDM parameters
	//
	private double xdm1Begin;    //  the beginning value for range of xdm1
	private double xdm1End;       //  the end value for range of xdm1
	private double xdm2Begin;    //  the beginning value for range of xdm2
	private double xdm2End;       //  the end value for range of xdm2
	private double xdmStep;        //  step size for the xdm range

	//
	// after running emul program, do we clean the input and output files?
	//
	private boolean cleaningFiles;

	//
	// error accounting
	// also result output
	//
	private double finalXDM1;    //   the result XDM parameter 1
	private double finalXDM2;    //   the result XDM parameter 2
	private double rmsError;       //    the rms error for the result
}
