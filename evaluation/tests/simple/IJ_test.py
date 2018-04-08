from ij import IJ, ImagePlus, ImageStack
from script.imglib.analysis import DoGPeaks  
from script.imglib.color import Red  
from script.imglib.algorithm import Scale2D  
from script.imglib.math import Compute  
from script.imglib import ImgLib  
from ij3d import Image3DUniverse  
# from javax.vecmath import Color3f, Point3f  
  
# Load a stack of images: a fly brain, in RGB  
imp = IJ.openImage("http://imagej.nih.gov/ij/images/flybrain.zip")  
stack = imp.getImageStack()
  
print "number of slices:", imp.getNSlices()  

# Show the peaks as spheres in 3D, along with orthoslices:  
univ = Image3DUniverse(512, 512)  
# univ.addIcospheres(ps, Color3f(1, 0, 0), 2, cell_diameter/2, "Cells").setLocked(True)  
# univ.addOrthoslice(imp).setLocked(True)  
univ.show()

# A list of green slices  
greens = []  
  
# Iterate each slice in the stack  
for i in xrange(1, imp.getNSlices()+1):  
  # Get the ColorProcessor slice at index i  
  cp = stack.getProcessor(i)  
  # Get its green channel as a FloatProcessor  
  fp = cp.toFloat(1, None)  
  # ... and store it in a list  
  greens.append(fp)  
  
# Create a new stack with only the green channel  
stack2 = ImageStack(imp.width, imp.height)  
for fp in greens:  
  stack2.addSlice(None, fp)  
  
# Create a new image with the stack of green channel slices  
imp2 = ImagePlus("Green channel", stack2)  
# Set a green look-up table:  
IJ.run(imp2, "Green", "")  
imp2.show()
