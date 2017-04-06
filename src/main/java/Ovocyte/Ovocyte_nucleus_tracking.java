package Ovocyte;

import Ovocyte.tools.MTrack2;
import Ovocyte.tools.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.plugin.Duplicator;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.util.ImageProcessorReader;
import ome.units.quantity.Length;
import ome.units.quantity.Time;


/**
 *
 * @author phm
 */
public class Ovocyte_nucleus_tracking implements PlugIn {
    
    private String OutDir;
    private String imageName; 
    

/** Track centroid of nucleus and nucleolus  
    /* 
    */
    private void trackObject2(ImagePlus img, boolean nuc) throws IOException {
        // Mtrack2 tracks black objects
        ImagePlus imgNucleus = new Duplicator().run(img);
        ImagePlus imgNucleolus = new ImagePlus();
        ImageCalculator ic = new ImageCalculator();
        IJ.run(imgNucleus, "Fill Holes", "stack");
        if(nuc) {
            imgNucleolus = ic.run("Subtract create stack", imgNucleus, img);
            FileSaver imgNucleolusSave = new FileSaver(imgNucleolus);
            imgNucleolusSave.saveAsTiff(OutDir+imageName+"_NucleolusMask2.tif");
            MTrack2.track(imgNucleolus,Tools.minSize, Integer.MAX_VALUE,20,OutDir,imageName+"Nucleolus.xls"); 
            imgNucleolus.close();
        }
        else {
            FileSaver imgNucleusSave = new FileSaver(imgNucleus);
            imgNucleusSave.saveAsTiff(OutDir+imageName+"_NucleusMask2.tif");
            MTrack2.track(imgNucleus,Tools.minSize, Integer.MAX_VALUE,20,OutDir,imageName+"Nucleus.xls"); 
        }
            
    }
    
    @Override
   public void run(String arg) {
           
        try {
            String imageDir = IJ.getDirectory("Choose Directory Containing TIF Files...");
            if (imageDir == null) return;
            OutDir = Tools.createDir(imageDir);
            File inDir = new File(imageDir);
            String [] imageFile = inDir.list();
            if (imageFile == null) return;
        // Open tif files
        // create OME-XML metadata store of the latest schema version
            ServiceFactory factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            for (int f = 0; f < imageFile.length; f++) {
                if (imageFile[f].endsWith(".tif")) {
                    imageName = imageFile[f].substring(0, imageFile[f].indexOf(".tif"));
                    String tifFile = inDir + File.separator + imageFile[f];
                    reader.setId(tifFile);
                    reader.getSeriesCount();
                    int series =  reader.getSeriesCount() - 1;
                    reader.setSeries(series);
                    // Read time set interval
                    Time st = meta.getPixelsTimeIncrement(series);
                    if (st != null) {
                        Tools.cal.frameInterval = st.value().doubleValue();
                        Tools.cal.setTimeUnit(st.unit().toString());
                    }
                    else {
                        Tools.cal.frameInterval = Tools.frame_time;
                        Tools.cal.setTimeUnit(Tools.time_unit);
                    }
                    
                // Read spatial calibration
                    Length px = meta.getPixelsPhysicalSizeX(series);
                    if (px != null) {
                        Tools.cal.pixelHeight = Tools.cal.pixelWidth = px.value().doubleValue();
                        Tools.cal.setUnit(Tools.pixel_unit);
                    }
                    else {
                        Tools.cal.pixelHeight = Tools.cal.pixelWidth = Tools.pixelWidth;
                        Tools.cal.setUnit(Tools.pixel_unit);
                    }
                    
                // write headers
                    ImagePlus img = Tools.OpenImage(reader, imageName);

                // crop nucleus
                    Tools.crop_nucleus(img);
                    
                // Nucleolus
                    ImagePlus imgNucleolus = new Duplicator().run(img);
                    imgNucleolus.setDimensions(1, 1, reader.getSizeT());
                    imgNucleolus.setCalibration(Tools.cal);
                    IJ.run(imgNucleolus,"Gaussian Blur...", "sigma=4 stack");
                    
                    imgNucleolus.setSlice(imgNucleolus.getNSlices()/2);
                    imgNucleolus.show();
                    IJ.run("Threshold...", ""); 
                    IJ.setAutoThreshold(imgNucleolus, "Default dark");
                    new WaitForUserDialog("Choose threshold, click Apply and Ok").show();
                    imgNucleolus.hide();
                    // stack registration
                    IJ.showStatus("Stack registration ...");
                    IJ.run(imgNucleolus, "StackReg", "transformation=[Rigid Body]");
                    Prefs.blackBackground = false;
                    IJ.run(imgNucleolus, "Convert to Mask", "method=Default background=Dark");
                    img.close();
                    // Track nucleolus 
                    
                    trackObject2(imgNucleolus, true);
                    // Track nucleus 
                    
                    trackObject2(imgNucleolus, false);
            }
        }
        IJ.showStatus("Process done");
        } catch (DependencyException | ServiceException | FormatException | IOException ex) {
            Logger.getLogger(Ovocyte_nucleus_tracking.class.getName()).log(Level.SEVERE, null, ex);
        }
    } 
}

