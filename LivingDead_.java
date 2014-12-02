/**
LivingDead - Classifies red and green fluorescent Synechocystis cells 
for the analysis of vitality.
A full description of the technique can be found in "K Schulze, D A Lopez, 
U M Tillich, and M Frohme, A Simple Live / Dead Analysis for unicellular
Cyanobacteria using a new Autofluorescence Assay, automated Microscopy, 
and ImageJ, BMC Biotechnology"

Copyright (C) 2011  Katja Schulze (kschulze@th-wildau.de)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

import ij.*;
import ij.measure.ResultsTable;
import ij.process.*;
import ij.gui.*;

import java.awt.image.ColorModel;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Vector;

import javax.swing.JFileChooser;

import ij.plugin.*;
import ij.plugin.frame.*;


public class LivingDead_ implements PlugIn {

	public void run(String arg) {
		
		//Ask for parameters:
        GenericDialog gd = new GenericDialog("Set Calibration Parameters");
        gd.addNumericField("Image height in mm:", 0.273, 4);
        gd.addNumericField("Image width in mm:", 0.364, 4);
        gd.addNumericField("Height of used counting chamber in mm:", 0.02, 4);
        gd.addNumericField("Dilution of the sample:", 2, 4);
        gd.addCheckbox("Show images?", true);
        gd.addCheckbox("Save data?", true);
        gd.showDialog();
        if (gd.wasCanceled()) return;
        
        double height_mm = (double) gd.getNextNumber();
        double width_mm = (double) gd.getNextNumber();
        double height_ch = (double) gd.getNextNumber();
        double dilution = (double) gd.getNextNumber();
        boolean process = gd.getNextBoolean();
        boolean save = gd.getNextBoolean();
        final double volumen = (height_mm * width_mm * height_ch)/1000;
		
        //Choose directory 
        JFileChooser chooser = new JFileChooser(); 
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.setDialogTitle("Select target directory");
        chooser.setAcceptAllFileFilterUsed(false);
        File directory = null;
        if (chooser.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) {
        	directory = chooser.getSelectedFile();
          } else {
        	IJ.showMessage("No Selection."); return;
        }
        
		//Filter for bmp and tif and jpg
		FilenameFilter filter  = new FilenameFilter()
		{
				public boolean accept(File dir, String name) {
					return name.toLowerCase().endsWith(".bmp") || name.toLowerCase().endsWith(".tif") || name.toLowerCase().endsWith(".tiff") ||
						   name.toLowerCase().endsWith(".jpg") || name.toLowerCase().endsWith(".jpeg") ;
					
		}};

		
		if (directory.listFiles().length == 0 ){IJ.showMessage("Directory contains no files.");return;}
		
		File[] imageList = directory.listFiles(filter);
		if (imageList.length == 0 ){IJ.showMessage("Directory contains no suited images (bmp, tif or jpg)."); return;}
		
		java.util.Arrays.sort(imageList);
		int countRed = 0; 
		int countGreen = 0;
		Vector<Double> yValues = new Vector<Double>();
		RoiManager rm = null ;
		int imagecount = 0; 
		
		// Analyze images in folder
		for (int k = 0; k<imageList.length; k++){
			ImagePlus imp = new ImagePlus(imageList[k].getAbsolutePath());
			if (process == true) imp.show(); 
			
			// analyzes only RGB images
			if (imp.getType() != ImagePlus.COLOR_RGB) continue;
						
			rm = RoiManager.getInstance();
			if (rm==null){
				rm = new RoiManager();
				if (process == false)rm.hide();
			}
			rm.runCommand("reset");
			
			ColorProcessor cp = (ColorProcessor)imp.getProcessor();

			// exclude images which contain no fluorescent signal
			if(imp.getProcessor().getStatistics().max < 15.0){
				if (process == true){
					imp.changes = false;
					imp.close();
				}
				continue; 
			}else{ imagecount++;}
			
			// threshold RGB chanels
			byte[] R = new byte[imp.getWidth()*imp.getHeight()];
			byte[] G = new byte[imp.getWidth()*imp.getHeight()];
			byte[] B = new byte[imp.getWidth()*imp.getHeight()];
			cp.getRGB(R, G, B);
			
			ByteProcessor rp = new ByteProcessor(imp.getWidth(), imp.getHeight(), R, ColorModel.getRGBdefault());
			ByteProcessor gp = new ByteProcessor(imp.getWidth(), imp.getHeight(), G, ColorModel.getRGBdefault());
			
			AutoThresholder adjust = new AutoThresholder();
			int globalThreshold = adjust.getThreshold("MaxEntropy", rp.getHistogram());
			rp.threshold(globalThreshold);
			rp.invert();
			
			globalThreshold = adjust.getThreshold("MaxEntropy", gp.getHistogram());
			if (globalThreshold < gp.getStatistics().mode){
				gp.setColor(256);
				gp.fill();
			}
			else{
			gp.threshold(globalThreshold);
			gp.invert();
			}
			
			R = (byte[]) rp.getPixels();
			G = (byte[]) gp.getPixels();
			
			// Combine thresholded chanels
			for (int a = 0; a <R.length; a++){
				if (G[a] == 0) R[a] = 0; 
			}
			
			//Register particles
			ImagePlus particle = NewImage.createByteImage("particle", imp.getWidth(), imp.getHeight(),1, NewImage.FILL_BLACK);
			particle.getProcessor().setPixels(R);
			IJ.run(particle, "Analyze Particles...", "size=5-Infinity circularity=0.00-1.00 show=Nothing add");

			// Classify red and green particles
			cp.getRGB(R, G, B);
			ByteProcessor gint = new ByteProcessor(imp.getWidth(), imp.getHeight(),G, ColorModel.getRGBdefault());
			double mean = 0.0; 
			Roi[] rois =  rm.getRoisAsArray();
			
			for (int i = 0 ; i < rois.length; i++){
				gint.setRoi(rois[i]);
				mean = gint.getStatistics().mean;
				if (mean < 50){
					countRed++;
				}
				else{
					countGreen++;
				}
				yValues.add(mean);
			}
			
			if (process == true){
				rm.runCommand("Show All");
				IJ.wait(500);
				imp.changes = false;
				imp.close();
			}
		}
		
		// Create histogram
		double[] histX = new double[257];
		double[] histY = new double[257];
		for(int i = 0 ; i<257; i++){
			histY[i]=0;
		}
		
		for(int i = 1 ; i<257; i++){
			for (int j = 0; j<yValues.size(); j++){
				double value = yValues.get(j);
				if ((int) value == i)histY[i]++; 
			}
			histX[i]=i;
		}
		
		int total = countRed+countGreen; 
	
		//Normalize histogram
		for(int i = 1 ; i<257; i++){
				histY[i]=histY[i]/total; 
		}
		
		// Plot and save histogram
		Plot plot = new Plot(directory.getAbsolutePath(), "Green-Mean of Roi (Cyano) ", "Frequency of Green-Mean", histX, histY ); 
		plot.show();
		ImagePlus imgplot = new ImagePlus("Plot",plot.getImagePlus().getImage());
		if (save == true) IJ.saveAs(imgplot, "Tiff", directory.getAbsolutePath()+"/plot.tif");
		
		// Save plot data
		ResultsTable plotdata = new ResultsTable();
		for (int i = 0; i<histX.length; i++){
			plotdata.incrementCounter();
			plotdata.addValue("x", histX[i]);
			plotdata.addValue("y", histY[i]);
		}
		if (save == true){
			try {
				plotdata.saveAs(directory+"/plotdata.csv");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		plotdata.reset();
		
		double totalvol = volumen * imagecount;
		double concentration = (total / totalvol)*dilution; 
		double percent_r = (double) countRed/total*100; 
		double percent_g = (double) countGreen/total*100 ; 
		
		
		// Show and save results
		if (WindowManager.getFrame("Log") == null){
			IJ.log("directory; total cell count; red (live) count; green (dead) count; red (live) in %; green (dead) in %; concentration (cells/ml)");
		}
		IJ.log(directory+"; "+total+"; "+countRed+"; "+countGreen+"; "+percent_r+"; "+percent_g+"; "+concentration);
		
		if (save == true){
			ResultsTable countdata = new ResultsTable();
			countdata.incrementCounter();
			countdata.addValue("total cell count", total);
			countdata.addValue("red (live) count", countRed);
			countdata.addValue("green (dead) count", countGreen);
			countdata.addValue("red (live) in %", percent_r);
			countdata.addValue("green (dead) in %", percent_g);
			countdata.addValue("concentration (cells/ml)", concentration);
			try {
				countdata.saveAs(directory+"/countdata.csv");
			} catch (IOException e) {
				e.printStackTrace();
			}
			countdata.reset();
		}
		
		rm.close();
	}
}
