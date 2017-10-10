#Data log and important points for capelin prediction project


## General Notes
E00 files converted to RData files are organized as a list with a dataframe (ice attributes) and polygons that spatially represent these attributes.  
* There are two subsets that must be performed:
  + First subset removes rows of the dataframe that do not meet criteria for concentration (Cx) and stage of development (Sx).  This represents “Among polygons”
  + Second subset removes columns based on the stage of development and recalculates areas of the polygon.  Completed.  
  
  Metadata for e00 files: http://www.ec.gc.ca/glaces-ice/default.asp?lang=En&n=503E8E74-1
  
  This has attribute information AND projection information
  
  Projection Information

Projection: Lambert Conformal Conic
Units: Metres
Spheroid: Clarke 1866 
Datum: NAD27 
1st Standard Parallel: 49 0 0.000
2nd Standard Parallel: 77 0 0.000
Central Meridian: -100 0 0.000
Latitude of Projections Origin: 40 0 0.000
False Easting: 0.000
False Northing: 0.00


## e00 files
* Notes on the e00 files (resolved):
  +	On the egg diagrams, there are specific protocols for stages of ice development.  Dots to the left of a number indicate that number, and all others to the left, have a dot, i.e., 74. Really means 7.4.  However, it appears that in the e00 file, this has been accurately  recorded, i.e., a 7.4. really is a 7.4.
  +	1 file downloaded in error (one day after the regular file) – wrote code to eliminate this file.
  +	Egg attribute code is given for “Remote egg” and “egg” and for “Fast Ice” but values do not appear in the partial concentrations and stage of developments.  
    - The code extracts the Egg_Attr data
	
* Notes on the e00 files (unresolved):
  +	Note that ~9+ is a code in E_CS for “strips and patches of ice” (see egg metadata).  Is this an issue?  If yes, will need to deal with a ~9+
    - No action
  +	What is the difference between “egg” and “remote egg”. 
    - Have searched and cannot find.  May need to contact EC
  +	About 1% of all e00 files do not convert to *.Rdata.  
    - Not clear why this is the case and both Paul Regular and I have tried to figure it out.  The code creates a avc_data folder with appropriate subfolders but only creates a few files.  
    - Dates are as follows: 19730528, 19780122, 19780423, 19830102, 19840313, 19900226, 19931213, 20110718, 20111114, 20111121, 20120723, 20120730, 20121126, 20131223, 20150216, 20151109, 20160208, 20161205, 20170213
    - Solution (thanks to Pete):
      - Load files into QGIS.  Layer-> add Layer -> Add Vector Layer -> select file - select all
      - Save teh PAL (polygon layer) as a shapefile.  Select all 4 layers in teh layer pannel -> right click "Save As"Layer -> name file and hit ok (use the defaults)
      - Import to R.  Now, it is a SPDF rather than a SLineDF.  Dataframe is there and can be querried.
      - Problems: doesn't plot well although all data is there. Have explored but not sure why
  + queries for ice volume
    - no action
  + errors:
   there appear to be two main types of errors---------------
 # conversion errors - see spreadsheet - these don't really matter
 1) 4 errors for one file
Error in e00toavc(e00file, file.path(avcdir, "bin")) : 
  ERROR 4: Unable to create coverage directory: avc_data/19730528/bin.

Error in get.arcdata(avcdir, "bin") : 
  ERROR 3: Attempt to read past EOF in avc_data/19730528\bin\arc.adf.

Error in get.paldata(avcdir, "bin") : 
  ERROR 4: Failed to open file avc_data/19730528\bin\pal.adf

Error in get.tabledata(file.path(avcdir, "info"), "bin.PAT") : 
  ERROR 3: Attempt to read past EOF in avc_data/19730528/info\arc.dir.
 these vary from file to file and sometimes, the error is different, e.g.
 Error in get.tabledata(file.path(avcdir, "info"), "bin.PAT") : 
 ERROR 3: Attempt to read past EOF in avc_data/19780423/info\arc.dir.


 2) This error is less important - we have gone through the work and downloaded teh files but it would be ncie to solve for the future
Error in file(con, "r") : cannot open the connection

In addition: Warning message:
 In file(con, "r") :
  cannot open file 'e00_data/19730604.e00': Too many open files

calculation errors  
  3) [1] "1982-03-18 NST"
Error in RGEOSBinTopoFunc(spgeom1, spgeom2, byid, id, drop_lower_td, unaryUnion_if_byid_false,  : 
  TopologyException: Input geom 0 is invalid: Self-intersection at or near point -57.942815797715888 54.922564573680816 at -57.942815797715888 54.922564573680816
  
  [1] "1984-02-26 NST"
Error in RGEOSBinTopoFunc(spgeom1, spgeom2, byid, id, drop_lower_td, unaryUnion_if_byid_false,  : 
  TopologyException: Input geom 0 is invalid: Self-intersection at or near point -50.184149377353116 48.265060927600707 at -50.184149377353116 48.265060927600707
  
  [1] "1988-03-06 NST"
Error in RGEOSBinTopoFunc(spgeom1, spgeom2, byid, id, drop_lower_td, unaryUnion_if_byid_false,  : 
  TopologyException: Input geom 0 is invalid: Self-intersection at or near point -53.98627133956596 49.656478912159606 at -53.98627133956596 49.656478912159606
  
  [1] "1992-02-20 NST"
Error in RGEOSBinTopoFunc(spgeom1, spgeom2, byid, id, drop_lower_td, unaryUnion_if_byid_false,  : 
  TopologyException: Input geom 0 is invalid: Self-intersection at or near point -63.997875049438377 46.957036130392055 at -63.997875049438377 46.957036130392055
In addition: Warning message:
In iceVolume(x@data) : NAs introduced by coercion

[1] "2012-01-30 NST"
Error in RGEOSBinTopoFunc(spgeom1, spgeom2, byid, id, drop_lower_td, unaryUnion_if_byid_false,  : 
  TopologyException: Input geom 0 is invalid: Self-intersection at or near point -58.466159741596755 57.632328446064101 at -58.466159741596755 57.632328446064101

  [1] "within1"
[1] "2012-04-02 NDT" (note that error and warning only occur in this month)
Error in RGEOSBinTopoFunc(spgeom1, spgeom2, byid, id, drop_lower_td, unaryUnion_if_byid_false,  : 
  TopologyException: Input geom 0 is invalid: Self-intersection at or near point -62.102833099999998 57.543979810000003 at -62.102833099999998 57.543979810000003

In addition: Warning messages:
1: In gBuffer(egg, byid = TRUE, width = 0) :
  Spatial object is not projected; GEOS expects planar coordinates
2: In gBuffer(egg, byid = TRUE, width = 0) :
  Polygons object missing comment attribute ignoring hole(s). See function createSPComment.
3: In gBuffer(egg, byid = TRUE, width = 0) :
  Polygons object missing comment attribute ignoring hole(s). See function createSPComment.
4: In gBuffer(egg, byid = TRUE, width = 0) :
  Polygons object missing comment attribute ignoring hole(s). See function createSPComment.
5: In gBuffer(egg, byid = TRUE, width = 0) :
  Polygons object missing comment attribute ignoring hole(s). See function createSPComment.
6: In gBuffer(egg, byid = TRUE, width = 0) :
  Polygons object missing comment attribute ignoring hole(s). See function createSPComment.
  
  4) probably 1986-12-30 (one earlier than this too)

This appears at the end of all calculations
  [1] "within1"
Warning messages:
1: In gBuffer(egg, byid = TRUE, width = 0) :
  Spatial object is not projected; GEOS expects planar coordinates
2: In gBuffer(egg, byid = TRUE, width = 0) :
  Polygons object missing comment attribute ignoring hole(s). See function createSPComment.
3: In gBuffer(egg, byid = TRUE, width = 0) :
  Polygons object missing comment attribute ignoring hole(s). See function createSPComment.
4: In gBuffer(egg, byid = TRUE, width = 0) :
  Polygons object missing comment attribute ignoring hole(s). See function createSPComment.
5: In gBuffer(egg, byid = TRUE, width = 0) :
  Polygons object missing comment attribute ignoring hole(s). See function createSPComment.

However, these have no visible effect on the calculations aside from skewing the medians somewhat

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
  + note that NULL is not working in the calcAreaVolume function.  Not sure why but have a hack - just created a list with all possible values.  This works fine for now.
  Not clear if these subsets removes polygons or not for mapping?  Check this as well.
## To Do:
+	Create a markdown file for ouputs of ice-capelin-update

## Notes on minlat
The e00 data comes in LCC format which makes perfect sense for Canadian waters (see metadata).  But its hard to visualize.  Further, Paul made the filters in WGS84.  While its like adding an extra set, the calcAreaVol() function converts data to WGS84 to match the filters and then converts it back.

To keep minlat in degrees and WGS84, calculate iceTiming() before Area and Volume when the projection is WGS84.  
* Problems:
  + Area is 0 but volume is not and minlat has a value.  
    - this is because ice area is not 0 but ~ 0.  Therefore, there is a polygon and values for volume/minlat
  + if there is 1 feature in the SPDF, it usually results in a null feature in teh sub.egg.  
    - The default minlat for these has been set to 55.
  + HOwever, there are polygons with 0 features.  These result in a sub.egg with no features.  This has been corrected by changing the function to make the default value 55 lat for is.null and "try-error" issues

## Predictive Model
  + what additional covariates
  + what hypotheses for subsetting ice
    - all ice v. Gary's proposed subset v. ?????
  
    

