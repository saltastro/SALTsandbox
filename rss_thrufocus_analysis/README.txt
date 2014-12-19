thrufoc_rssimage.py

    You will need to edit hardwired directory path "datadir" at line 22 of the
thrufoc_rssimage.py. Then put qred.sex and spectrograph/imgdist.txt there. 

The sample .txt shows you what I get when I run it on the 20141212 6645 filter 
thrufocus data. On my machine, this takes about 2 minutes. Basically, the 
program looks at the middle image in the run to make an updated distortion map 
using a SExtraction of the grid spots, then SExtracts those spots from all the 
images, forming a cube of SExtractor data versus row and column. "spots" is the 
number of spots found. There are actually 330 of them, so the spurious ones are 
thrown out, if necessary. Then it cubic-spline interpolates in the fwhm cube to 
calculate the rss mean of fwhm over a plane, and varies focus, tip, and tilt to 
minimize it. This is for three cases, one using grid points in a cross on the x 
and y axis, for longslit work, one using the whole grid, and one just using the 
spots at the rim of the mask. (Right now, I prefer "maskrim", since there is 
evidence that the N99 mask is not flat, so the spots at the center give the 
wrong focus). You can see that it advises a tilt adjustment of about 2 turns, 
in the negative sense, so the right side of the detector should come out. One 
would like to get the suggested adjustments down to less than a turn or so. 
