/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Ovocyte.tools;

import java.util.*;
import java.io.*;
import ij.*;
import ij.plugin.filter.ParticleAnalyzer;
import ij.measure.*;
 

/**
	Uses ImageJ's particle analyzer to track the movement of
	multiple objects through a stack. 
	Based on the Object Tracker plugin filter by Wayne Rasband

	Based on Multitracker, but should be quite a bit more intelligent
	Nico Stuurman, Vale Lab, UCSF/HHMI, May,June 2003
*/
public class MTrack2  {


	int		nParticles;
	float[][]	ssx;
	float[][]	ssy;
	static int 	minTrackLength = 20;
	static boolean 	bSaveResultsFile = true;
        static float   	maxVelocity = 50;
	static int 	maxColumns=75;


	static public class particle {
		float	x;
		float	y;
		int	z;
		int	trackNr;
		boolean inTrack=false;
                boolean flag=false;

		public void copy(particle source) {
			this.x=source.x;
			this.y=source.y;
			this.z=source.z;
			this.inTrack=source.inTrack;
			this.flag=source.flag;
		}

		public float distance (particle p) {
			return (float) Math.sqrt(sqr(this.x-p.x) + sqr(this.y-p.y));
		}
	}
        
	static public void track(ImagePlus imp, int minSize, int maxSize, float maxVelocity, String directory, String filename) {
		int nFrames = imp.getStackSize();
		if (nFrames<2) {
			IJ.showMessage("Tracker", "Stack required");
			return;
		}

		ImageStack stack = imp.getStack();
		int options = 0; // set all PA options false
		int measurements = ParticleAnalyzer.CENTROID;

		// Initialize results table
		ResultsTable rt = new ResultsTable();
		rt.reset();

		// create storage for particle positions
		List[] theParticles = new ArrayList[nFrames];
		int trackCount=0;

		// record particle positions for each frame in an ArrayList
		for (int iFrame=1; iFrame<=nFrames; iFrame++) {
			theParticles[iFrame-1]=new ArrayList();
			rt.reset();
			ParticleAnalyzer pa = new ParticleAnalyzer(options, measurements, rt, minSize, maxSize);
			pa.analyze(imp, stack.getProcessor(iFrame));
			float[] sxRes = rt.getColumn(ResultsTable.X_CENTROID);				
			float[] syRes = rt.getColumn(ResultsTable.Y_CENTROID);
			if (sxRes==null)
				continue;

			for (int iPart=0; iPart<sxRes.length; iPart++) {
				particle aParticle = new particle();
				aParticle.x=sxRes[iPart];
				aParticle.y=syRes[iPart];
				aParticle.z=iFrame-1;
				theParticles[iFrame-1].add(aParticle);
			}
			IJ.showProgress((double)iFrame/nFrames);
		}

		// now assemble tracks out of the particle lists
		// Also record to which track a particle belongs in ArrayLists
		List theTracks = new ArrayList();
		for (int i=0; i<=(nFrames-1); i++) {
			IJ.showProgress((double)i/nFrames);
			for (ListIterator j=theParticles[i].listIterator();j.hasNext();) {
				particle aParticle=(particle) j.next();
				if (!aParticle.inTrack) {
					// This must be the beginning of a new track
					List aTrack = new ArrayList();
					trackCount++;
					aParticle.inTrack=true;
					aParticle.trackNr=trackCount;
					aTrack.add(aParticle);
					// search in next frames for more particles to be added to track
					boolean searchOn=true;
					particle oldParticle=new particle();
					particle tmpParticle=new particle();
					oldParticle.copy(aParticle);
					for (int iF=i+1; iF<=(nFrames-1);iF++) {
						boolean foundOne=false;
						particle newParticle=new particle();
						for (ListIterator jF=theParticles[iF].listIterator();jF.hasNext() && searchOn;) { 
							particle testParticle =(particle) jF.next();
							float distance = testParticle.distance(oldParticle);
							// record a particle when it is within the search radius, and when it had not yet been claimed by another track
							if ( (distance < maxVelocity) && !testParticle.inTrack) {
								// if we had not found a particle before, it is easy
								if (!foundOne) {
									tmpParticle=testParticle;
									testParticle.inTrack=true;
									testParticle.trackNr=trackCount;
									newParticle.copy(testParticle);
									foundOne=true;
								}
								else {
									// if we had one before, we'll take this one if it is closer.  In any case, flag these particles
									testParticle.flag=true;
									if (distance < newParticle.distance(oldParticle)) {
										testParticle.inTrack=true;
										testParticle.trackNr=trackCount;
										newParticle.copy(testParticle);
										tmpParticle.inTrack=false;
										tmpParticle.trackNr=0;
										tmpParticle=testParticle;
									}
									else {
										newParticle.flag=true;
									}	
								}
							}
							else if (distance < maxVelocity) {
							// this particle is already in another track but could have been part of this one
							// We have a number of choices here:
							// 1. Sort out to which track this particle really belongs (but how?)
							// 2. Stop this track
							// 3. Stop this track, and also delete the remainder of the other one
							// 4. Stop this track and flag this particle:
								testParticle.flag=true;
							}
						}
						if (foundOne)
							aTrack.add(newParticle);
						else
							searchOn=false;
						oldParticle.copy(newParticle);
					}
					theTracks.add(aTrack);
				}
			}
		}

		// Create the column headings based on the number of tracks
		// with length greater than minTrackLength
		// since the number of tracks can be larger than can be accomodated by Excell, we deliver the tracks in chunks of maxColumns
		// As a side-effect, this makes the code quite complicated
		String strHeadings = "Frame";
		int frameCount=1;
		for (ListIterator iT=theTracks.listIterator(); iT.hasNext();) {
			List bTrack=(ArrayList) iT.next();
			if (bTrack.size() >= minTrackLength) {
				if (frameCount <= maxColumns)
                                    strHeadings += "\tX" + frameCount + "\tY" + frameCount;
				frameCount++;
			}
		}

		String contents="";
		boolean writefile=false;
		if (filename != null) {
			File outputfile=new File (directory,filename);
			if (!outputfile.canWrite()) {
				try {
					outputfile.createNewFile();
				}
				catch (IOException e) {
					IJ.showMessage ("Error", "Could not create "+directory+filename);
				}
			}
			if (outputfile.canWrite())
				writefile=true;
			else
				IJ.showMessage ("Error", "Could not write to " + directory + filename);
		}

		int repeat=(int) ( (frameCount/maxColumns) );
		float reTest = (float) frameCount/ (float) maxColumns;
		if (reTest > repeat)
			repeat++;
		
		// write to file
		if (writefile) {
			try {
				File outputfile=new File (directory,filename);
				BufferedWriter dos= new BufferedWriter (new FileWriter (outputfile));
				dos.write(strHeadings+"\n",0,strHeadings.length()+1);
				for (int j=1; j<=repeat;j++) {
					int to=j*maxColumns;
					if (to > frameCount-1)
						to=frameCount-1;
					//String stLine="Tracks " + ((j-1)*maxColumns+1) +" to " +to;
					//dos.write(stLine + "\n",0,stLine.length()+1);
					for (int i=0; i<=(nFrames-1); i++) {
						String strLine = "" + (i+1);
						int trackNr=0;
						int listTrackNr=0;
						for (ListIterator iT=theTracks.listIterator(); iT.hasNext();) {
							trackNr++;
							List bTrack=(ArrayList) iT.next();
							boolean particleFound=false;
							if (bTrack.size() >= minTrackLength) {
								listTrackNr++;
								if ( (listTrackNr>((j-1)*maxColumns)) && (listTrackNr<=(j*maxColumns))) {
									for (ListIterator k=theParticles[i].listIterator();k.hasNext() && !particleFound;) {
										particle aParticle=(particle) k.next();
										if (aParticle.trackNr==trackNr) {
											particleFound=true;
											strLine+="\t" + aParticle.x + "\t" + aParticle.y;
										}
									}
									if (!particleFound)
										strLine+="\t \t ";
								}
							}
						}
						dos.write(strLine + "\n",0,strLine.length()+1);
					}
				}
				dos.close();
			}
			catch (IOException e) {
				if (filename != null)
					IJ.error ("An error occurred writing the file. \n \n " + e);
			}
		}
	}

	// Utility functions
	static double sqr(double n) {return n*n;}
	
	int doOffset (int center, int maxSize, int displacement) {
		if ((center - displacement) < 2*displacement) {
			return (center + 4*displacement);
		}
		else {
			return (center - displacement);
		}
	}

 	double s2d(String s) {
		Double d;
		try {d = new Double(s);}
		catch (NumberFormatException e) {d = null;}
		if (d!=null)
			return(d.doubleValue());
		else
			return(0.0);
	}

}


