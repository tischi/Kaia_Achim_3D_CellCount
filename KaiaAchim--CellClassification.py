from ij.io import OpenDialog
from ij.io import Opener
from ij.gui import GenericDialog
from ij.plugin import ZProjector, RGBStackMerge, SubstackMaker, Concatenator
from ij import IJ, ImagePlus, ImageStack
from ij.plugin import Duplicator
from ij.process import StackStatistics
from ij.plugin import ImageCalculator
from ij.measure import ResultsTable
import os
import os.path
import re
from jarray import array
from ij.process import ImageConverter
import math
from ij.macro import MacroRunner

from loci.plugins import BF
from loci.common import Region
from loci.plugins.in import ImporterOptions

from automic.table import TableModel			# this class stores the data for the table
from automic.table import ManualControlFrame 	#this class visualises TableModel via GUI
from java.io import File

# import my analysis function collection
import os, sys, inspect
this_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if this_folder not in sys.path:
  print this_folder
  sys.path.insert(0, this_folder)
import ct_analysis_functions as af
reload(af)

all_results = []

# Structure for batch analysis:
# 
# - main 
#  - parameters = get_analysis_parameters()
#  - folder = get_folder()
#  - table = init_results_table()
#  - data_info = get_data_info(folder) 
#  - batch_analyze(folder, parameters, table, dada) 
#    - iDataSet = 0
#    - for iDataSet in 'data_info' 
#       - imp, filename = return_imp(iDataSet, data_info) # might be different for different projects
#       - analyze(imp, folder, filename, parameters, table)
#         - do analysis on imp
#         - write results to table(iDataSet)
#         - write segmentation overlay images to filepath_output = f(folder, filename)
# - show interactive results table
# - save results table

# Improved structure:

# Structure for batch analysis:
#  
# - main 
#  - parameters = get_analysis_parameters()
#  - folder = get_folder()
#  - table = init_results_table()
#  - get_data_info(folder, table) 
#  - batch_analyze(parameters, table) 
#    - for row in table 
#       - analyze(row, table, parameters)
#         - imp = load_imp(table, row)
#         - write results to table(row)
#         - write segmentation overlay images (use from table(row))


def analyze(iDataSet, tbModel, p, output_folder):

  #
  # LOAD FILES
  #

  # get path to some file in the directory for the "Image Sequence..."
  file_path = tbModel.getFileAPth(iDataSet, "DAPI", "IMG")
  file_path = os.path.join(os.path.split(file_path)[0], os.listdir(os.path.split(file_path)[0])[1])
  
  # Trans
  file_id = tbModel.getFileName(iDataSet, "TRANS", "IMG")
  parameters = "open="+file_path+" file=(^"+file_id+"--W.*Trans.tif) sort"
  IJ.run("Image Sequence...", parameters);
  imp = IJ.getImage()
  output_file = str(file_id)+"--Trans--Segmentation.tif"
  IJ.saveAs(imp, "TIFF", os.path.join(output_folder, output_file))
  tbModel.setFileAPth(output_folder, output_file, iDataSet, "TRANS","IMG")
  imp.close()
  
  # Dapi
  file_id = tbModel.getFileName(iDataSet, "DAPI", "IMG")
  parameters = "open="+file_path+" file=(^"+file_id+"--W.*Dapi.tif) sort"
  IJ.run("Image Sequence...", parameters);
  imp = IJ.getImage()
  output_file = str(file_id)+"--Dapi--Segmentation.tif"
  IJ.saveAs(imp, "TIFF", os.path.join(output_folder, output_file))
  tbModel.setFileAPth(output_folder, output_file, iDataSet, "DAPI","IMG")

  #
  # INIT
  #
  IJ.run("Options...", "iterations=1 count=1"); 


  #
  # FLOATING POINT CONVERSION
  #

  IJ.run(imp, "32-bit", "");
  
  #
  # SCALING
  #

  # IJ.run(imp, "Scale...", "x=0.5 y=0.5 z=1.0 width=768 height=768 depth=9 interpolation=Bilinear average process"); imp.setRoi(392,386,750,762);
  IJ.run(imp, "Scale...", "x=0.25 y=0.25 z=1.0 width=384 height=384 depth=9 interpolation=Bilinear average process"); imp.setRoi(578,576,380,382);
  IJ.run(imp, "Crop", "");
  
  #
  # BACKGROUND SUBTRACTION
  #
  # IJ.run(imp, "Subtract...", "value=32768 stack");
  
  #
  # ENHANCEMENT
  #
  imp = af.filter_diff_gauss_3d(imp, dx1=1, dy1=1, dz1=1, dx2=5, dy2=5, dz2=3)
  
  #
  # SEGMENTATION
  # 

  IJ.setMinAndMax(imp, 0, 4095); IJ.run(imp, "16-bit", "");
  IJ.run(imp, "3D Spot Segmentation", "seeds_threshold="+str(p["threshold"])+" local_background=0 radius_0=2 radius_1=20 radius_2=22 weigth=0.5 radius_max=10 sd_value=1 local_threshold=[Local Mean] seg_spot=Block volume_min=1 volume_max=1000000 seeds=Automatic spots="+imp.getTitle()+"radius_for_seeds=2 output=[Label Image]");
  # get number of objects in 3D
  imp_label = IJ.getImage()
  num_objects = StackStatistics(imp_label).max
  imp_bw = af.threshold(imp_label, 1)

  
  '''
  IJ.run(imp, "3D Fast Filters","filter=MaximumLocal radius_x_pix=3.0 radius_y_pix=3.0 radius_z_pix=3.0 Nb_cpus=4");
  imp_bw = IJ.getImage()
  imp_bw = af.threshold(imp_bw, 20)
  '''

  '''
  IJ.run(imp, "Extended Min & Max 3D", "operation=[Extended Maxima] dynamic=60 connectivity=6");
  imp_bw = IJ.getImage()
  imp_bw = af.threshold(imp_bw, 1)
  '''

  
  '''
  imp_bw = af.find_maxima_stack(imp, 50)
  imp_bw = af.threshold(imp_bw, 1)
  '''
  
  '''
  imp_bw = imp.duplicate()
  # maybe i can find the maximum in the image for the high threshold?
  IJ.run(imp_bw, "3D Hysteresis Thresholding", "high=60 low=30");
  imp_bw = af.threshold(imp_bw, 1)
  '''

  '''
  imp_bw = af.threshold(imp, p["threshold"]); imp_bw.show()
  '''
  
  '''
  IJ.run(imp, "3D Iterative Thresholding", "min_vol_pix=10 max_vol_pix=10000 min_threshold=0 criteria_method=ELONGATION threshold_method=STEP value_method=1 starts");
  imp_bw = IJ.getImage()
  imp_bw = af.extract_channel_frame(imp_bw, 1, 1, title="")
  imp_bw = af.threshold(imp_bw, 1); imp_bw.show()
  '''
    
  '''
  # auto global threshold
  stats = StackStatistics(imp); #print stats.mean, stats.stdDev, stats.mean + 3 *stats.stdDev
  threshold = stats.mean + 8 * stats.stdDev
  imp_bw = af.threshold(imp, threshold); #af.show_max(imp_bw)
  
  # auto local threshold
  #imp_8bit = af.make_8bit_using_min_max(imp)
  #imp_8bit.show()
  #imp_bw = imp_8bit.duplicate()
  #IJ.run(imp_bw, "Auto Local Threshold", "method=Niblack radius=80 parameter_1=2 parameter_2=0 white stack");
  #fff
  '''
  
  #
  # MEASURE
  #
  
  IJ.run("Set Measurements...", "area centroid mean stack skewness kurtosis shape redirect=None decimal=2");
  IJ.run(imp_bw, "Analyze Particles...", "size="+str(p["minimal_particle_area"])+"-1000 show=[Bare Outlines] clear stack");
  rt = ResultsTable.getResultsTable()
  imp_outlines = IJ.getImage(); IJ.run(imp_outlines, "Invert", "stack");
  
  '''
  from ij.plugin.frame import RoiManager;
  RoiManager rm=RoiManager.getInstance();
  rm.runCommand("Save", "absolutePath\fileName.zip"); 
  '''
  
  #imp_outlines.getProcessor().invertLut()
  #imp_outlines.show()
  #ddd
  #IJ.run(imp_outlines, "Invert LUT", "");
  
  
  
  #
  # EVALUATE MEASUREMENTS
  #

  #print rt.getColumnHeadings()
  
  issues = ""
  tbModel.setNumVal(num_objects, iDataSet, "Segmented_Particles")
  
  if num_objects == 0:
    issues = issues + " no_particles "  
  elif num_objects > 1:
    issues = issues + " multiple_particles " 
  else:    
    # not round
    rc  = rt.getColumn(rt.getColumnIndex("Round"))
    tbModel.setNumVal(min(rc), iDataSet, "Minimal_Roundness")
    if min(rc) < p["minimal_roundness"]:
      issues = issues + " elongated "      
    # shifting
    xy  = zip(rt.getColumn(rt.getColumnIndex("X")), rt.getColumn(rt.getColumnIndex("Y")))
    shifts = [math.sqrt((v[0]-xy[0][0])**2 + (v[1]-xy[0][1])**2) for v in xy]
    tbModel.setNumVal(max(shifts), iDataSet, "Max_Shift")
    if max(shifts) > p["maximal_shift"]:
      issues = issues + " shifting "
 
  tbModel.setValue(issues, iDataSet, "Issues")
  tbModel.setBoolVal(len(issues)==0, iDataSet, "OK")
  
  # make and save outline image
  imp_merged = af.merge_images(imp, imp_outlines, force_same_bit_depth=True)
  # todo: determine LUT
  # todo: save ROIs instead
  imp_merged = af.change_composite_display(imp_merged, [["white", 1, StackStatistics(imp).max*1.5],["red",1,255]])
  imp_merged.show()
  output_file = output_file + "--fiji-out.tif"
  IJ.saveAs(imp_merged, "TIFF", os.path.join(output_folder, output_file))
  tbModel.setFileAPth(output_folder, output_file, iDataSet, "SEG","IMG")
  

#
# ANALYZE INPUT FILES
#
def determine_input_files(folder_name, tbModel):
  
  pattern = re.compile('(.*)--W(.*)--P(.*)--Z(.*)--T(.*)--.*')  
  ids = []
  for root, directories, filenames in os.walk(folder_name):
	for filename in filenames:
	   if filename == "Thumbs.db":
	     continue
	   match = re.search(pattern, filename)
	   if match.group(1) == None:
	     print "filename does not match:",filename
	     ddd
	   ids.append(match.group(1))
  unique_ids = map(str,sorted(map(int,set(ids))))
  #print unique_ids
 
  file_path = os.path.join(folder_name, filename)
      
  for i, file_id in enumerate(unique_ids):

    tbModel.addRow()

    # Dapi
    tbModel.setFileAPth(folder_name, file_id, i, "DAPI","IMG")
    
    # Trans   
    tbModel.setFileAPth(folder_name, file_id, i, "TRANS","IMG")
    
  return(tbModel)


  
def get_parameters(p, num_data_sets):
  gd = GenericDialog("Correct 3D Drift Options")

  gd.addMessage("Found "+str(num_data_sets)+" data sets")
  gd.addStringField("analyse", "all");

  gd.addMessage("Image analysis parameters:")
  for k in p.keys():
    gd.addNumericField(k, p[k], 2);
  gd.showDialog()
  if gd.wasCanceled():
    return

  to_be_analyzed = gd.getNextString()
  for k in p.keys():
    p[k] = gd.getNextNumber()
    
  return to_be_analyzed, p

    
if __name__ == '__main__':
  
  #
  # GET INPUT FOLDER
  #
  od = OpenDialog("Click on one of the image files in the folder to be analysed", None	)
  input_folder = od.getDirectory()
  #input_folder = "C:/Users/almf/Desktop/kaia/data/"
  if input_folder is None:
    fff

  #
  # MAKE OUTPUT FOLDER
  #
  output_folder = input_folder[:-1]+"--fiji"
  if not os.path.isdir(output_folder):
    os.mkdir(output_folder)

 
  #
  # INIT OUTPUT TABLE
  #
  tbModel = TableModel(output_folder)
  tbModel.addFileColumns('TRANS','IMG')
  tbModel.addFileColumns('DAPI','IMG')
  tbModel.addFileColumns('SEG','IMG')
  tbModel.addValColumn("Minimal_Roundness", "NUM")
  tbModel.addValColumn("Segmented_Particles", "NUM")
  tbModel.addValColumn("Max_Shift", "NUM")
  tbModel.addColumn("Issues")
  tbModel.addValColumn("OK", "BOOL")
  frame=ManualControlFrame(tbModel)
  frame.setVisible(True)

  #
  # DETERMINE INPUT FILES
  #
  tbModel = determine_input_files(input_folder, tbModel)

  #
  # GET PARAMETERS
  #
  p = dict()
  scale = 0.25
  p["minimal_particle_area"] = 0 * scale *scale
  p["threshold"] = 30
  p["minimal_roundness"] = 0.6
  p["maximal_shift"] = 20 * scale
    
  to_be_analyzed, p = get_parameters(p, tbModel.getRowCount())
  
  #
  # ANALYZE
  #
  
  if to_be_analyzed is not "all":
    af.close_all_image_windows()
    analyze(int(to_be_analyzed)-1, tbModel, p, output_folder)
  else:
    for i in range(tbModel.getRowCount()):
      af.close_all_image_windows()
      analyze(i, tbModel, p, output_folder)

  af.close_all_image_windows()