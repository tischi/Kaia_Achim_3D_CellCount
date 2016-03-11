from ij.io import OpenDialog
from ij.io import Opener
from ij.plugin import ZProjector, RGBStackMerge
from ij import IJ, ImagePlus, ImageStack, WindowManager
from ij.plugin import Duplicator, Filters3D, GaussianBlur3D, HyperStackConverter, FolderOpener, ImageCalculator
from ij.process import StackStatistics, ImageProcessor, ImageConverter, LUT
from ij.plugin.filter import MaximumFinder
from ij.measure import ResultsTable, Calibration
#from mcib_plugins import Watershed_Split3D
#from mcib_plugins.analysis import simpleMeasure
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from loci.formats import ImageReader
from sc.fiji.analyzeSkeleton import AnalyzeSkeleton_
#from sc.fiji.skeletonize3D import Skeletonize3D_
from features import TubenessProcessor
#from skeleton_analysis import AnalyzeSkeleton_, Skeletonize3D_
from java.awt import Color

#from ij3d import ImageJ3DViewer

import time
import os


def StepMaker():
  def StepReader(func):
    def MakeStep(*args, **kwargs):
      start_time = time.time()
      output = func(*args, **kwargs)
      elapsed_time = time.time() - start_time
      print func.__name__,"-- time:",round(elapsed_time,3)
      return output
    return MakeStep
  return StepReader



@StepMaker()
def get_lif_series_length(fpath):
  reader = ImageReader()
  reader.setId(fpath)
  return reader.getSeriesCount();

@StepMaker()
def open_from_lif_series(fpath, i):
  options = ImporterOptions()
  options.setId(fpath)
  options.clearSeries()
  options.setSeriesOn(i, True)
  imps = BF.openImagePlus(options)
  return imps[0]

@StepMaker()
def find_maxima_stack(imp, tolerance):
  stack = imp.getImageStack()
  stack_out = ImageStack(imp.width, imp.height)
  for i in range(1, imp.getNSlices() + 1):
    ip = stack.getProcessor(i)
    # findMaxima(ImageProcessor ip, double tolerance, double threshold, int outputType, boolean excludeOnEdges, boolean isEDM)
    ip_max = MaximumFinder().findMaxima(ip, tolerance, ImageProcessor.NO_THRESHOLD, MaximumFinder.IN_TOLERANCE, False, False)
    stack_out.addSlice(str(i), ip_max)
    #segip.invert()
  return ImagePlus("Maxima", stack_out)

@StepMaker()
def close_non_image_window(window_name):
  # forcefully closes one specific non image window
  frames = WindowManager.getNonImageWindows()
  for i in range(len(frames)):
    #IJ.log(frames[i].getTitle())
    if (frames[i].getTitle() == window_name):
      frames[i].dispose()

@StepMaker()
def show_3D(imp):
  imp.show()
  viewer = ImageJ3DViewer()
  viewer.setCoordinateSystem("false")
  viewer.add(imp.getTitle(), "None", imp.getTitle(), "0", "true", "true", "true", "1", "0")

  
@StepMaker()
def change_composite_display(imp, luts):
  '''
  luts = [["green",0,1],["red",0,1]]
  '''
  # make Luts
  grayLut = LUT.createLutFromColor(Color.white)
  #grayLut.min = luts[0][1]; grayLut.max = luts[0][2]
  #print grayLut
  
  greenLut = LUT.createLutFromColor(Color.green)
  #greenLut.min = luts[1][1]; greenLut.max = luts[1][2]
  #print greenLut
  # change display mode (otherwise the Lut changes do not work)
  imp.setDisplayMode(IJ.COLOR)
  # set Luts
  print "imp.isComposite()",imp.isComposite()
  imp.setC(1); imp.setChannelLut(grayLut)
  imp.setDisplayRange(luts[0][1], luts[0][2])
  imp.setC(2); imp.setChannelLut(greenLut)
  #IJ.run(imp, "Fire", "");
  #IJ.setMinAndMax(imp, 100, 255);
  imp.setDisplayRange(luts[1][1], luts[1][2])
  #imp.getProcessor().invertLut()
  # check that it worked
  print imp.getLuts()
  imp.setDisplayMode(IJ.COMPOSITE)
  return imp
  

@StepMaker()
def clear_results_table_if_exists():
  rt = ResultsTable.getResultsTable()
  if rt!=None:
    rt.reset()
  
@StepMaker()
def close_all_image_windows():
  # forcefully closes all open images windows
  ids = WindowManager.getIDList();
  if (ids==None):
    return
  for i in ids:
     imp = WindowManager.getImage(i)
     if (imp!=None):
       win = imp.getWindow()
       if (win!=None):
         imp.changes = False # avoids the "save changes" dialog
         win.close()
         

@StepMaker()
def keepSlices(imp, first, last):
  stack = imp.getImageStack()
  # remove front part
  for i in range(1,first+1):
    stack.deleteSlice(1)
  # remove back part
  n = stack.getSize()
  for i in range(last-first,n+1):
    stack.deleteSlice(last-first)
  imp.setStack(stack)
  return imp


@StepMaker()
def geometrical_measure_3D(imp_label):
  # remove results table if exists:
  close_non_image_window("Results")
  # measure
  #imp_label.show()
  rt = ResultsTable.getResultsTable()
  IJ.run(imp_label, "3D Geometrical Measure", "");
  rt = ResultsTable.getResultsTable()
  return rt


@StepMaker()
def intensity_measure_3D(imp_intens, imp_label):
  imp_intens.show()
  imp_label.show()
  IJ.run(imp_intens, "3D Intensity Measure", "objects=[%s] signal=[%s]" % (imp_label.getTitle(), imp_intens.getTitle()) );
  imp_intens.hide()
  imp_label.hide()
  rt = ResultsTable.getResultsTable()
  return rt

@StepMaker()
def enhance_tubes(imp, tube_radius):
  tp = TubenessProcessor(tube_radius, True)
  imp_out = tp.generateImage(imp)
  return imp_out

# todo
@StepMaker()
def intensity_measure_3D_direct(imp_intens, imp_label):
  # todo: handle the arraylist output
  sm = simpleMeasure(imp_label)
  arraylist = sm.getMeasuresStats(imp_intens)
  for i in arraylist:
    print i
  
@StepMaker()
def save_hyperstack_as_image_sequence(imp, path):
  # make dir
  if not os.path.isdir(path): 
    os.mkdir(path)
  #imp.show()
  ssize = imp.getStackSize()
  titleext = imp.getTitle()
  #print titleext
  title = os.path.splitext(titleext)[0]
  dimA = imp.getDimensions()
  print "Title",title,"     -----   Dimensions",dimA
  for c in range(dimA[2]):   
    for z in range(dimA[3]):
        for t in range(dimA[4]):
            # imagej is one-based and we save zero-based
            numberedtitle = \
            title + "--C" + IJ.pad(c, 2) + \
            "--Z" + IJ.pad(z, 4) + \
            "--T" + IJ.pad(t, 4) + ".tif"
            stackindex = imp.getStackIndex(c+1, z+1, t+1) # imagej is one-based
            aframe = ImagePlus(numberedtitle, imp.getStack().getProcessor(stackindex))
            IJ.log("saving: " + os.path.join(path, numberedtitle))
            IJ.saveAs(aframe, "TIFF", os.path.join(path, numberedtitle))
            
@StepMaker()
def load_image_sequence_from_path_to_hyperstack(path, nc, nz, nt):
  imp = FolderOpener().openFolder(path)
  imp2 = HyperStackConverter.toHyperStack(imp, nc, nz, nt, "default", "Color");
  return imp2

 
@StepMaker()
def project_z(imp, method):
  zp = ZProjector(imp)
  if method=="max":
    zp.setMethod(ZProjector.MAX_METHOD)
  if method=="sum":
    zp.setMethod(ZProjector.SUM_METHOD)
    
  zp.setStopSlice(imp.getNSlices())
  zp.doHyperStackProjection(True)
  zpimp = zp.getProjection()
  
  #IJ.run(imp, "Z Project...", "projection=[Max Intensity] all");
  #impout = IJ.getImage()
  
  return zpimp

@StepMaker()
def show_max(imp):
  imp_max = project_z(imp, "max")
  imp_max.show()
  

@StepMaker()
def make_16bit_using_minmax(_imp):
  # ????? this does not make sense, or ???
  imp = _imp.duplicate()
  stats = StackStatistics(imp)
  IJ.setMinAndMax(imp, stats.min, stats.max)
  IJ.run(imp, "16-bit", "")
  IJ.setMinAndMax(imp, 0, 65535)
  return imp

@StepMaker()
def append_stack(imp, imp_append):
  stack_append = imp_append.getImageStack()
  stack = imp.getImageStack()
  for i in range(1, stack_append.getSize() + 1):
    stack.addSlice("", stack_append.getProcessor(i))
  return ImagePlus(imp.getTitle(), stack)
  
@StepMaker()
def skeletonize(imp):
  imp_skel = imp.duplicate()
  skel = Skeletonize3D_()
  skel.setup("", imp_skel)
  skel.run(None)
  return imp_skel

  
@StepMaker()
def analyze_skeleton(imp_skel):
  skel = AnalyzeSkeleton_()
  skel.setup("",imp_skel)
  skelResult = skel.run()
  #print "slab voxels",skeleton.getListOfSlabVoxels()
  #print "starting slab voxels",skeleton.getListOfStartingSlabVoxels() 
  #print "slabs",skeleton.getSlabs() 
  #print "endpoints",skeleton.getEndPoints() 
  graph = skelResult.getGraph()
  print "num trees:",skelResult.getNumOfTrees()
  print "len graph:",len(graph)
  #print graph
  #print graph[0]
  #print graph[0].getEdges()
  '''
  https://github.com/fiji/AnalyzeSkeleton/blob/master/src/main/java/sc/fiji/analyzeSkeleton/AnalyzeSkeleton_.java
  http://javadoc.imagej.net/Fiji/skeleton_analysis/Graph.html
  http://javadoc.imagej.net/Fiji/skeleton_analysis/Edge.html
  http://javadoc.imagej.net/Fiji/skeleton_analysis/Vertex.html
  '''
  for g in graph:
    print "NEW GRAPH"
    for e in g.getEdges():
      print e.getLength()
      print e.getSlabs()
      voxels = e.getSlabs()
      #print voxels
      print "slab voxel",voxels[0]
      print "slab voxel",voxels[0].x
      
      
  # add to new image
  # stack.setVoxel(int x, int y, int z, double value) 
  return skeleton

@StepMaker()
def make_image_from_graphs(graphs):
  skel = AnalyzeSkeleton_()
  skel.setup("",imp_skel)
  skelResult = skel.run()
  #print "slab voxels",skeleton.getListOfSlabVoxels()
  #print "starting slab voxels",skeleton.getListOfStartingSlabVoxels() 
  #print "slabs",skeleton.getSlabs() 
  #print "endpoints",skeleton.getEndPoints() 
  graph = skelResult.getGraph()
  print "num trees:",skelResult.getNumOfTrees()
  print "len graph:",len(graph)
  #print graph
  #print graph[0]
  #print graph[0].getEdges()
  '''
  https://github.com/fiji/AnalyzeSkeleton/blob/master/src/main/java/sc/fiji/analyzeSkeleton/AnalyzeSkeleton_.java
  http://javadoc.imagej.net/Fiji/skeleton_analysis/Graph.html
  http://javadoc.imagej.net/Fiji/skeleton_analysis/Edge.html
  http://javadoc.imagej.net/Fiji/skeleton_analysis/Vertex.html
  '''
  for g in graph:
    print "NEW GRAPH"
    for e in g.getEdges():
      print e.getLength()
      print e.getSlabs()
      #stack.setVoxel(int x, int y, int z, double value) 
 
  # add to new image
  return skeleton

  


@StepMaker()
def get_skeleton_slabs(imp):
  skel = AnalyzeSkeleton_()
  skel.setup("",imp)
  skel.run()
  # imp_labelled_skel = ImagePlus("", skel.getLabeledSkeletons()) # label image of skeleton
  imp_tagged_skel = ImagePlus("", skel.getResultImage(False)) # label image of skeleton
  #imp_tagged_skel.show()
  imp_slabs = gate(imp_tagged_skel,126,128)
  return imp_slabs



@StepMaker()
def measureSkeletonTotalLength(imp):
  IJ.run(imp,"Analyze Skeleton (2D/3D)", "prune=none")
  totalLength = 0
  nBranches = 0
  rt = ResultsTable.getResultsTable()
  avgLengths = rt.getColumn(rt.getColumnIndex("Average Branch Length"))
  for i in range(len(avgLengths)):
    totalLength = totalLength + rt.getValue("# Branches", i) * rt.getValue("Average Branch Length", i)
    nBranches = nBranches + rt.getValue("# Branches", i)	
  print root,",",imp.getTitle(),",",totalLength,",",nBranches

@StepMaker()
def openCroppedImageUsingBF(filepath,x,y,xw,yw,zs,ze, crop):
  # read in and display ImagePlus(es) with arguments
  options = ImporterOptions()
  options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE)
  if crop:
    options.setCrop(True)
    options.setCropRegion(0, Region(x,y,xw,yw))
    options.setTBegin(0, zs)
    options.setTEnd(0, ze)
  options.setId(filepath)
  imps = BF.openImagePlus(options)
  return imps[0]

@StepMaker()
def scale(imp,x,y,z):
  #imp = _imp.duplicate()
  #impout.setTitle("scaled")
  IJ.run(imp, "Scale...", "x=%s y=%s z=%s interpolation=Bilinear average process create" % (x, y, z));
  imp_out = IJ.getImage(); imp_out.setTitle("Scaled"); imp_out.hide()
  return imp_out


@StepMaker()
def make_outline_image(imp):
  
  #imp.show() # one needs to display it, otherwise below function cannot find it
  IJ.run(imp, "3D Fast Filters", "filter=Minimum radius_x_pix=1.0 radius_y_pix=1.0 radius_z_pix=1.0");
  #imp.hide()
  
  imp_minimum = WindowManager.getImage("3D_Minimum")
  imp_outlines = imp.duplicate()
  ic = ImageCalculator()
  imp_outlines = ic.run("Subtract stack", imp_outlines, imp_minimum)
  # clean up
  imp_minimum.close()
  imp_outlines.setTitle("Outlines")  
  return imp_outlines

@StepMaker()
#def imp_math(imp2,imp1,op,bit_depth_result):
#  # private static String[] operators = {"Add","Subtract","Multiply","Divide", "AND", "OR", "XOR", "Min", "Max", "Average", "Difference", "Copy", "Transparent-zero"};
#  
#  imp_out = imp2.duplicate()
#  
#  if bit_depth_result == 32:
#    IJ.run(imp_out, "32-bit", "")
#  elif bit_depth_result == 16:
#    IJ.run(imp_out, "16-bit", "")
#  elif bit_depth_result == 8:
#    IJ.run(imp_out, "8-bit", "")
#  calc = ImageCalculator()
#  if op=="add":
#    imp_out = ic.run(op+" stack", imp_out, imp1)
#  elif op=="subtract":
#    imp_out = ic.run(op+" stack", imp_out, imp1)
#  return imp_out

  
  
@StepMaker()
def merge_images(_imp1, _imp2, force_same_bit_depth):
  imp1 = _imp1.duplicate()
  imp2 = _imp2.duplicate()
  imp1.show(); imp2.show();
  if(force_same_bit_depth):
    bd = max(imp1.getBitDepth(),imp2.getBitDepth())
    if bd==16:
      convert_to_16bit(imp1)
      convert_to_16bit(imp2)
    if bd==32:
      convert_to_32bit(imp1)
      convert_to_32bit(imp2)
  IJ.run(imp1, "Merge Channels...", "c1=[%s] c2=[%s] create" % (imp1.getTitle(), imp2.getTitle()))
  
  imp_out = WindowManager.getImage("Composite") 
  if not imp_out:
    imp_out = WindowManager.getImage("Merged")
  imp_out.hide()
  imp1.hide(); imp2.hide()
  return imp_out


@StepMaker()
def objects_counter_3d(imp_bw, imp_intens, measurements, threshold, vol_min, vol_max):
  # set options
  options = "%s dots_size=5 font_size=10 show_numbers white_numbers redirect_to=%s" % (measurements, imp_intens.getTitle())
  IJ.run("3D OC Options", options);
  # make sure the output image does not exist already as this leasds to weird behaviors....
  outlines_title = "Surface map of %s redirect to %s" % (imp_bw.getTitle(),imp_intens.getTitle())
  if WindowManager.getImage(outlines_title) != None: 
    WindowManager.getImage(outlines_title).close()
  # analyze
  imp_bw.show(); imp_intens.show();
  IJ.run(imp_bw, "3D Objects Counter", "threshold=%s min.=%s max.=%s exclude_objects_on_edges surfaces statistics" % (threshold, vol_min, vol_max))
  imp_out = IJ.getImage(); imp_out.setTitle("outlines"); imp_out.hide()
  return imp_out

@StepMaker()
def set_xyz_calibration(imp, dx, dy, dz, unit):
  cal = Calibration(imp)
  cal.setUnit(unit)
  cal.pixelWidth = dx
  cal.pixelHeight = dy
  cal.pixelDepth = dz
  imp.setCalibration(cal)
  return imp


@StepMaker()
def convert_to_8bit(imp, vmin, vmax):
  IJ.setMinAndMax(imp, vmin, vmax)
  ImageConverter(imp).convertToGray8() 
  return imp

@StepMaker()
def convert_to_16bit(imp):
  ImageConverter(imp).convertToGray16() 
  return imp

@StepMaker()  
def make_8bit_using_min_max(_imp):
  imp = _imp.duplicate()
  stats = StackStatistics(imp)
  IJ.setMinAndMax(imp, stats.min, stats.max)
  ImageConverter(imp).convertToGray8() 
  return imp

@StepMaker()
def convert_to_32bit(imp):
  IJ.run(imp, "32-bit", "")
  return imp

@StepMaker()
def filter_gauss_3d(imp,dx,dy,dz):
  imp_out = imp.duplicate() 
  GaussianBlur3D().blur(imp_out, dx, dy, dz)
  return imp_out


@StepMaker()
def filter_diff_gauss_3d(imp,dx1,dy1,dz1,dx2,dy2,dz2):
  
  imp_1 = imp.duplicate() 
  GaussianBlur3D().blur(imp_1, dx1, dy1, dz1)
  
  imp_2 = imp.duplicate() 
  GaussianBlur3D().blur(imp_2, dx2, dy2, dz2)
  
  ic = ImageCalculator()
  imp_diff = ic.run("Subtract stack create", imp_1, imp_2)  
  imp_diff.show()
  
  return imp_diff


@StepMaker()
def filter_stack_max_3d(imp, dx ,dy, dz):
  stack_filtered = Filters3D.filter(imp.getStack(), Filters3D.MAX, dx, dy, dz)
  return ImagePlus("max filtered", stack_filtered)
  
@StepMaker()
def filter_open_gray_3d(imp,dx,dy,dz):
  imp.show()
  IJ.run("3D Fast Filters","filter=OpenGray radius_x_pix=%s radius_y_pix=%s radius_z_pix=%s" % (dx,dy,dz))
  imp_out = WindowManager.getImage("3D_OpenGray")
  imp_out.hide()
  return imp_out

"""
@StepMaker()
def filter_open_gray_3d_direkt(imp,dx,dy,dz):
  imp.show()
  IJ.run("3D Fast Filters","filter=OpenGray radius_x_pix=%s radius_y_pix=%s radius_z_pix=%s" % (dx,dy,dz))
  imp_out = WindowManager.getImage("3D_OpenGray")
  imp_out.hide()
  return imp_out

# public static ImageStack filterImageStack(ImageStack stackorig, int filter, float vx, float vy, float vz, int nbcpus, boolean showstatus) 

               if ((depth == 8) || (depth == 16)) {
                    res = FastFilters3D.filterIntImageStack(stack, filter, voisx, voisy, voisz, nbcpus, true);
                } else if (imp.getBitDepth() == 32) {
                    //if ((filter != FastFilters3D.SOBEL)) {
                    res = FastFilters3D.filterFloatImageStack(stack, filter, voisx, voisy, voisz, nbcpus, true);
//                    } else {
//                        if (debug) {
//                            IJ.log("Not implemented for 32-bits images");
//                        }
                    // }
                    // }
                } else {
                    IJ.log("Does not wotk with stack with bitDepth " + depth);
                }
                if (res != null) {
                    ImagePlus plus = new ImagePlus("3D_" + filters[filter], res);
                    plus.setCalibration(calibration);
                    plus.show();
""" 
  
@StepMaker()
def extract_channel_frame(imp, nChannel, nFrame, title=""):
  """ Extract a stack for a specific color channel and timeframe """
  # todo: preserve spatial calibration
  stack = imp.getImageStack()
  ch = ImageStack(imp.width, imp.height)
  for i in range(1, imp.getNSlices() + 1):
    index = imp.getStackIndex(nChannel, i, nFrame)
    ch.addSlice(str(i), stack.getProcessor(index))
  if title:
    imp_out = ImagePlus(title, ch)
  else:
    imp_out = ImagePlus("Channel " + str(nChannel), ch)
  return imp_out

  
@StepMaker()
def watershed_3d(imp, radius):
  imp.show() # one needs to display it, otherwise below function cannot find it
  IJ.run("3D Watershed Split", "binary=%s seeds=Automatic radius=%s" % (imp.getTitle(), radius))
  imp.hide()
  # get the output image
  imp_out = IJ.getImage()
  imp_out.hide()
  # remove the EDT image, which is also displayed
  WindowManager.getImage("EDT").close()
  return imp_out


@StepMaker()  
def fill_holes_3d(imp):
  imp_out = imp.duplicate() 
  IJ.run(imp_out, "3D Fill Holes", "");
  imp_out.setTitle("fill_holes_3d")
  return imp_out

@StepMaker()  
def fill_holes_2d_stack(imp):
  imp_out = imp.duplicate() 
  IJ.run(imp_out, "Fill Holes", "stack")
  imp_out.setTitle("fill_holes_2d_stack")
  return imp_out

      
@StepMaker()  
def attenuation_correction(imp, intensity_zmin, intensity_zmax):
  # check: http://rsb.info.nih.gov/ij/plugins/stack-contrast/index.htm
  stack = imp.getImageStack()
  stack_out = ImageStack(imp.width, imp.height)
  a = intensity_zmin
  b = intensity_zmax
  nz = imp.getNSlices()
  for iz in range(1, nz + 1):
    ip_out = stack.getProcessor(iz).duplicate()
    correction = (a+b)/2 * (1/(a-(a-b)*(iz-1)/(nz-1)))
    #print iz, correction
    ip_out.multiply(correction) 
    stack_out.addSlice(str(iz), ip_out)
    #segip.invert()
  return ImagePlus("Attenuation_corrected", stack_out)

@StepMaker()  
def multiply(imp, value):
  # check: http://rsb.info.nih.gov/ij/plugins/stack-contrast/index.htm
  stack = imp.getImageStack()
  stack_out = ImageStack(imp.width, imp.height)
  nz = stack.getSize()
  for iz in range(1, nz + 1):
    ip_out = stack.getProcessor(iz).duplicate()
    #print iz, correction
    ip_out.multiply(value) 
    stack_out.addSlice(str(iz), ip_out)
  #segip.invert()
  return ImagePlus("Multiplied", stack_out)

@StepMaker()
def measureSumIntensity3D(imp):
  stats = StackStatistics(imp)
  return stats.mean * stats.pixelCount

@StepMaker()
def auto_threshold(imp, method):
  impout = imp.duplicate() 
  IJ.run(impout, "Auto Threshold", "method=" + method + " white stack use_stack_histogram");
  impout.setTitle("Auto_threshold")
  impout.hide()
  return impout
  
@StepMaker()
def threshold(_img, threshold):
  imp = Duplicator().run(_img)
  #imp.show(); time.sleep(0.2)
  #IJ.setThreshold(imp, mpar['lthr'], mpar['uthr'])
  IJ.setThreshold(imp, threshold, 1000000000)
  IJ.run(imp, "Convert to Mask", "stack")
  imp.setTitle("threshold");
  #IJ.run(imp, "Divide...", "value=255 stack");
  #IJ.setMinAndMax(imp, 0, 1);
  return imp

@StepMaker()
def gate(_img, lower, higher):
  imp = Duplicator().run(_img)
  #imp.show(); time.sleep(0.2)
  #IJ.setThreshold(imp, mpar['lthr'], mpar['uthr'])
  IJ.setThreshold(imp, lower, higher)
  IJ.run(imp, "Convert to Mask", "stack")
  imp.setTitle("gated");
  #IJ.run(imp, "Divide...", "value=255 stack");
  #IJ.setMinAndMax(imp, 0, 1);
  return imp

