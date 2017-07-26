#Data log and important points for capelin prediction project


## General Notes
E00 files converted to RData files are organized as a list with a dataframe (ice attributes) and polygons that spatially represent these attributes.  
* There are two subsets that must be performed:
  + First subset removes rows of the dataframe that do not meet criteria for concentration (Cx) and stage of development (Sx).  This represents “Among polygons”
  + Second subset removes columns based on the stage of development and recalculates areas of the polygon.  Completed.  Calculations look correct but should confirm.
  
Not clear if these subsets removes polygons or not for mapping?  Check this as well.

## e00 files
* Notes on the e00 files (resolved):
  +	On the egg diagrams, there are specific protocols for stages of ice development.  Dots to the left of a number indicate that number, and all others to the left, have a dot, i.e., 74. Really means 7.4.  However, it appears that in the e00 file, this has been accurately  recorded, i.e., a 7.4. really is a 7.4.
  +	1 file downloaded in error (one day after the regular file) – wrote code to eliminate this file.
	

* Notes on the e00 files (unresolved):
  +	Note that ~9+ is a code in E_CS for “strips and patches of ice” (see egg metadata).  Is this an issue?  If yes, will need to deal with a ~9+
  +	What is the difference between “egg” and “remote egg”
  +	Egg attribute code is given for “Remote egg” and “egg” but not for “Fast Ice”.
  +	About 1% of all e00 files do not convert to *.Rdata.  
    - Not clear why this is the case and both Paul Regular and I have tried to figure it out.  The code creates a avc_data folder with appropriate subfolders but only creates a few files.  
    - Dates are as follows: 19730528, 19780122, 19780423, 19830102, 19840313, 19900226, 19931213, 20110718, 20111114, 20111121, 20120723, 20120730, 20121126, 20131223, 20150216, 20151109, 20160208, 20161205, 20170213
    - Try to convert in QGIS on recommendations from Pete.

## Notes on Area calculations 
* Unresolved:
  +	Is area being calculated properly?  I have done some tests and it appears that it is.  However, a systematic examination is required.
  + Get new data
    - ice (done)
    - capelin
    - calanus
    - phytoplankton
  + create ice hypotheses
  + check out ggmap
## To Do:
+	Create a markdown file for ouputs of ice-capelin-update

