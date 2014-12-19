thrufoc_rssimage.py

    You will need to edit hardwired directory paths at line 61 and 64 of the .py 
depending on where you want to put the distortion model data imgdist.txt. 

The sample .txt shows you what I get when I run it on the 20141124 6645 filter 
thrufocus data from 1124. On my machine, this takes about 2 minutes. Basically, 
the program looks at the middle image in the run to make an updated distortion 
map using a SExtraction of the grid spots, then SExtracts those spots from all 
the images, forming a cube of SExtractor data versus row and column. "spots" is 
the number of spots found. There are actually 330 of them, so the spurious ones 
are thrown out, if necessary. Then it cubic-spline interpolates in the fwhm 
cube to calculate the rss mean of fwhm over a plane, and varies focus, tip, and 
tilt to minimize it. This is for two cases, one using grid points in a cross on 
the x and y axis, for longslit work, and the other using the whole grid to get 
the imaging optimum. (Presumably you want longslit, since that mode is most 
sensitive to focus). You can see that it advises a tilt adjustment of about 4 
turns, in the negative sense, so the right side of the detector should come 
out. One would like to get the suggested adjustments down to less than a turn 
or so. 
