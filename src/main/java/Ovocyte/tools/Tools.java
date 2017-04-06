/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Ovocyte.tools;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.Duplicator;
import ij.plugin.filter.Analyzer;
import ij.process.ImageProcessor;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.formats.FormatException;
import loci.plugins.util.ImageProcessorReader;

/**
 *
 * @author phm
 */
public class Tools {
    
    static public final Calibration cal = new Calibration();
    // spatial calibration (micron)
    static final public double pixelWidth = 0.1135;
    static final public String pixel_unit = "microns";
    //Frame interval (sec)
    static final public double frame_time = 0.5;
    static final public String time_unit = "sec";
    static final private int theta = 10;
    static double nucXcent, nucYcent;
    // Object min size 50 microns 
    static final public int minSize = 100;
    
    /**
     * Create out folder
     * 
     * @return 
     */
    static public String createDir(String imageDir) {
        String OutDir = imageDir + File.separator + "Out" + File.separator;
        File outDir = new File(OutDir);
        if (!outDir.isDirectory())
            outDir.mkdirs();
        return(OutDir);
    }
    
    /**
     * 2D distance
     * @param Xc firstpoint
     * @param Yc firstpoint
     * @param Xp last point 
     * @param Yp last point
     * @return dist 2D distance
     * 
     */
    static public double calculateDistance(double Xc, double Yc, double Xp, double Yp) {
        double dist =  Math.sqrt( Math.pow( (Xc - Xp), 2) + Math.pow( (Yc - Yp), 2));
        return(dist);
    }
    
    /**
     * Open
     * @param r
     * @param imageName
     * @return
     * @throws FormatException
     * @throws IOException 
     */
    static public ImagePlus OpenImage(ImageProcessorReader r, String imageName) throws FormatException, IOException {
        // Read image
       r.setSeries(0);
       ImageStack stack = new ImageStack(r.getSizeX(), r.getSizeY());
       for (int n = 0; n < r.getSizeT(); n++) {
            ImageProcessor ip = null;
            try {
                ip = r.openProcessors(n)[0];
            } catch (FormatException ex) {
                Logger.getLogger(Tools.class.getName()).log(Level.SEVERE, null, ex);
            }
            IJ.showProgress(n, r.getImageCount());
            stack.addSlice("" + (n + 1), ip);
        }
        ImagePlus imgStack = new ImagePlus(imageName, stack);
        imgStack.setDimensions(1, 1, r.getSizeT());
        imgStack.setCalibration(cal);
        return imgStack;
    }
        
    static public void crop_nucleus(ImagePlus img) {
        ImagePlus imp = new Duplicator().run(img);
        // move to middle stack
        img.setT(imp.getNFrames()/2);
        // gaussian filter
        IJ.run(imp, "Gaussian Blur...", "sigma=10 slice");
        // Threshold nucleus
        IJ.setAutoThreshold(imp, "Default dark");
        IJ.run(imp, "Create Selection", "");
        IJ.run(imp, "To Bounding Box", "");
        IJ.run(imp, "Enlarge...", "enlarge=5");
        Roi rect = imp.getRoi();
        img.setRoi(rect);
        IJ.run(img, "Crop", "");
        img.deleteRoi();
        imp.flush();
        imp.close();
    }
    
        /**
     * Compute centroid distance to border for teta angle
     * @param perimeter, Xc, Yc
     * @return double []
     */
     
    static private ArrayList<Double> centerToBorder(ImagePlus img, double Xc, double Yc) {
        ArrayList distance = new ArrayList();
        for (int i = 0; i < 360; i+=theta) {
            double rad = i*(Math.PI/180);            
            for (int p = 0; p < img.getWidth(); p++) {
                int ptX = (int)(cal.getRawX(Xc) + Math.cos(rad)*p);
                int ptY = (int)(cal.getRawY(Yc) - Math.sin(rad)*p);
                if (img.getProcessor().getPixelValue(ptX, ptY) == 0) {
                    distance.add(calculateDistance((float)Xc, (float)Yc, (float)cal.getX(ptX), (float)cal.getY(ptY)));
                    p =  img.getWidth(); 
                }
            }
        }
        return(distance);
    }
    
        /** find centroid of the Nucleus  
    /*
    */
    static private void getCentroid(ImagePlus imgBin, int t) {
        ResultsTable table = new ResultsTable();
        Analyzer measure = new Analyzer(imgBin,Measurements.CENTROID,table);
        imgBin.setT(t);
        // measure parameters
        measure.measure();
	nucXcent = table.getValue("X",table.size()-1);
	nucYcent = table.getValue("Y",table.size()-1);
        table.reset();
    }
    
     static public void computeDistance (ImagePlus img, String imageName, String OutDir) throws IOException {
        ResultsTable distResults = new ResultsTable();
        distResults.setPrecision(3);
        for (int t = 1; t <= img.getNFrames(); t++) {
           nucXcent = 0;
           nucYcent = 0;
           int angle = 0;
           getCentroid(img, t);
           ArrayList<Double> centerBorderDist = centerToBorder(img, nucXcent, nucYcent);
           String header = "t"+t;
           for (int i = 0; i < centerBorderDist.size(); i++) { 
               distResults.setValue("tetha", i, angle);
               angle += theta;
               distResults.setValue(header, i,centerBorderDist.get(i));
            }     
        }
        // write table
        distResults.saveAs(OutDir+imageName+"_distances.xls");
        distResults.reset();
    }
     
     
}
