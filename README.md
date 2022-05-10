# BeadScan
Analyzes tiled 2p-digiFLIM images of GSBs (Koveal et al, Nature Communications 2022)

Current compatibility is 
   ThorimageLS v4.2.2020.3061
   MATLAB R2019b

Input files are .dFLIM files from a digiFLIM-enabled version of ThorimageLS, collected as X-Y tiled sets. The same tiles are imaged after each change of solution. The experiment.xml file is used automatically to determine image field size and tile numbers, and the fileset is read in automatically. 

BeadSurveyorT reads the files from specified folder location, displays them, and analyzes them into circular ROIs, based on specified radius limits (in pixels) and thresholds.  Sets should be read one at a time; autonumbering is \Set\, \Set_001\, \Set_002\, etc. If run during acquisition, the program will monitor file creation and read images when they are ready.  To use a different microscope file format, the user will need to import the data into the format of dFLIM_Image objects (as shown in the dFLIM_Classes folder).

First enter the target folder for the set to be read in, enter the conditions (e.g. concentration, pH) under which that set was acquired, then press the Collect button.  If parameter adjustment is needed, they can be optimized using Test Segmentation button; then reset the (auto-incremented) folder name and the set number back to one, and Collect again.  BeadSurveyor saves cumulative beaddata files, and the last one can be read by BeadGazer2.

BeadGazer2 is for analyzing the circular ROI data collected by BeadSurveyorT. The inputs are the last beaddata file (for the last set, since saving is cumulative), and optionally the roiimages.mat file saved (for the first set).  If the images are included, the individual data points are colored green for images that meet a neural-network-trained set of criteria.  First use the Read File button, then the Replot button.  Default X value is brightness and Y value is maximum delta-tau.  These can be changed by entering formulae into the Alternate X and Alternate Y fields and checking the adjacent box; sample formulae are shown (s# is the lifetime value for that ROI in the specified set number; n# is the photon count).

Once the scatter plot is created, the mouse can be used to select data points by "brushing" (choosing a box).  The corresponding ROI numbers are collected in a list to the right and plotted, using the options in the righthand part of the window.

