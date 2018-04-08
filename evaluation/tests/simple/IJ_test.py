from ij import IJ, ImagePlus, ImageStack  
  
# Load a stack of images: a fly brain, in RGB  
imp = IJ.openImage("http://imagej.nih.gov/ij/images/flybrain.zip")  
stack = imp.getImageStack()  
  
print "number of slices:", imp.getNSlices()  
  
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