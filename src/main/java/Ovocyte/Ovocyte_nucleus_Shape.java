package Ovocyte;


import Ovocyte.tools.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.io.FileSaver;
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
public class Ovocyte_nucleus_Shape implements PlugIn {

   
   public void run(String arg) {
           
        try {
            String imageDir = IJ.getDirectory("Choose Directory Containing TIF Files...");
            String OutDir = Tools.createDir(imageDir);
            if (imageDir == null) return;
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
                    String imageName = imageFile[f].substring(0, imageFile[f].indexOf(".tif"));
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
                        Tools.cal.setUnit(px.unit().toString());
                    }
                    else {
                        Tools.cal.pixelHeight = Tools.cal.pixelWidth = Tools.pixelWidth;
                        Tools.cal.setUnit(Tools.pixel_unit);
                    }
                    
                    // write headers
                    ImagePlus img = Tools.OpenImage(reader, imageName);

                    // crop nucleus
                    Tools.crop_nucleus(img);
                    // run pure denoise filter
                    IJ.run(img, "PureDenoise ...", "parameters='1 1' estimation='Auto Global' ");
                    WindowManager.getWindow("Log").setVisible(false);
                    img.close();
                    img.flush();
                    
                    ImagePlus imgFilter = WindowManager.getCurrentImage();
                    imgFilter.setDimensions(1, 1, reader.getSizeT());
                    imgFilter.setCalibration(Tools.cal);
                    imgFilter.hide();
                    IJ.run(imgFilter,"8-bit","");
                    IJ.run(imgFilter,"Gaussian Blur...", "sigma=2 stack");
                    // Threshold nucleus
                    Prefs.blackBackground = false;
                    IJ.run(imgFilter, "Convert to Mask", "method=Default background=Dark calculate");
                    //Close arround centroid
                    IJ.run(imgFilter, "Options...", "interations=2 count=1 do=Close stack");
                    // fill hole
                    IJ.run(imgFilter, "Fill Holes", "stack");
                    // stack registration
                    IJ.run(imgFilter, "StackReg", "transformation=[Rigid Body]");
                    FileSaver imgSave = new FileSaver(imgFilter);
                    imgSave.saveAsTiff(OutDir+imageName+"_Mask.tif");
                    // for each time compute distance center to border print all values
                    Tools.computeDistance(imgFilter, imageName, OutDir);
                    imgFilter.close();
                    imgFilter.flush();
                }
            }
            IJ.showStatus("Process done");
        } catch (DependencyException | ServiceException | FormatException | IOException ex) {
            Logger.getLogger(Ovocyte_nucleus_Shape.class.getName()).log(Level.SEVERE, null, ex);
        }
       } 
    }
