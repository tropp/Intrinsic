#!/usr/bin/env python
"""
maptest.py: Compute phase maps from streams of images.
Algorithm:
1. compute (or get) period of image sequence and number of repeats
2. fold data
3. take mean across all sequences (reduction to one sequence or cycle)
4. compute phase map (peak phase) as a function of (x,y) within the map
(we also compute amplitude, as it is useful to know...)
The result is both plotted (matplotlib) and written to disk.

To do:
1. sum phase maps (allow two inputs) to get absolute delay
2. speed phase calculation
3. data reduction in image to nxn blocks in (x,y)

5/12/2010 Paul B. Manis
UNC Chapel Hill
pmanis@med.unc.edu

"""

import numpy
import pickle
import matplotlib
import pylab
from pyqtgraph import metaarray
D = []
d = []

class testAnalysis():
    def __init__(self):
        global d
        self.imageData = []
        self.times = []
        #fileName = '/Volumes/Promise Pegasus/ManisLab_Data3/IntrinsicImaging/video_004.ma'
        fileName = '/Volumes/TRoppData/2016.03.14_000/Sound_Stimulation_001/000/Camera/frames.ma'
        try:
            im = metaarray.MetaArray(file=fileName)
        except:
            print "Error loading file: %s\n" % fileName
            return
        self.times = im.axisValues('Time')
        self.nFrames = numpy.shape(im)[0]
        self.imageData = numpy.array(im).astype(numpy.float32, copy=False)
        d = numpy.random.random(self.imageData.shape).astype('float32')
        for i in range(0,8):
	        dx = i*55
	        d[:, dx:dx+50, 50:200] += numpy.sin(numpy.linspace(0, 2.0*numpy.pi*d.shape[0]/39.0, d.shape[0])-i*numpy.pi/8)[:, numpy.newaxis, numpy.newaxis]
        self.imageData = d
        print "Image Loaded"

    def Analysis_FourierMap(self):
        global D
        # print "times: ", self.times # self.times has the actual frame times in it. 
        # first squeeze the image to 3d if it is 4d
        sh = numpy.shape(self.imageData);
        if len(sh) == 4:
           self.imageData = numpy.squeeze(self.imageData)
           sh = numpy.shape(self.imageData)
        print 'image shape: ', sh
        self.imagePeriod = 8.0 # image period in seconds.
        w = 2.0 * numpy.pi * self.imagePeriod
        # identify an interpolation for the image for one cycle of time
        dt = numpy.mean(numpy.diff(self.times)) # get the mean dt
        maxt = numpy.amax(self.times) # find last image time
        n_period = int(numpy.floor(maxt/self.imagePeriod)) # how many full periods in the image set?
        n_cycle = int(numpy.floor(self.imagePeriod/dt)); # estimate image points in a stimulus cycle
        ndt = self.imagePeriod/n_cycle
        i_times = numpy.arange(0, n_period*n_cycle*ndt, ndt) # interpolation times
        n_times = numpy.arange(0, n_cycle*ndt, ndt) # just one cycle
        print "dt: %f maxt: %f # images %d" % (dt, maxt, len(self.times))
        print "# full per: %d  pts/cycle: %d  ndt: %f #i_times: %d" % (n_period, n_cycle, ndt, len(i_times))
        B = numpy.zeros([sh[1], sh[2], n_period, n_cycle])
        #for i in range(0, sh[1]):
#            for j in range(0, sh[2]):
#                B[i,j,:] = numpy.interp(i_times, self.times, self.imageData[:,i,j])
        B = self.imageData[range(0, n_period*n_cycle),:,:]
        print 'new image shape: ', numpy.shape(self.imageData)
        print "B shape: ", numpy.shape(B)
        C = numpy.reshape(B, (n_period, n_cycle, sh[1], sh[2]))
        print 'C: ', numpy.shape(C)
        D = numpy.mean(C, axis=0)
        print "D: ", numpy.shape(D)
        sh = numpy.shape(D)
        A = numpy.zeros((sh[0], 2), float)
        print "A: ", numpy.shape(A)
        A[:,0] = numpy.sin(w*n_times)
        A[:,1] = numpy.cos(w*n_times)
        #
 
        sparse = 4
        resshape1 = sh[1]/sparse+1
        resshape2 = sh[2]/sparse+1
        rs = (resshape1, resshape2)
        print 'resshape:', rs
        self.phaseImage1 = numpy.zeros(rs)
        self.amplitudeImage1 = numpy.zeros(rs)
#        ik = 0
#        for i in range(0, sh[1], sparse):
#            jk = 0
#            for j in range(0, sh[2], sparse):
#                (p, residulas, rank, s) = numpy.linalg.lstsq(A, D[:,i,j])
#                self.amplitudeImage1[ik,jk] = numpy.hypot(p[0],p[1])
#                self.phaseImage1[ik, jk] = numpy.arctan2(p[1],p[0])
#                jk = jk + 1
#            ik = ik + 1

        DF = numpy.fft.fft(D, axis = 0)
        da = numpy.abs(DF[1,:,:])
        self.amplitudeImage = da
        dp = numpy.angle(DF[1,:,:])
        self.phaseImage = dp


        f = open('img_phase.dat', 'w')
        pickle.dump(self.phaseImage, f)
        f.close()
        f = open('img_amplitude.dat', 'w')
        pickle.dump(self.amplitudeImage, f)
        f.close()

        pylab.figure()
        pylab.subplot(2,2,1)
        pylab.title('Amplitude Map')
        pylab.imshow(self.amplitudeImage)
        pylab.colorbar()
#        pylab.subplot(2,2,3)
#        pylab.imshow(self.amplitudeImage1)
#        pylab.colorbar()
        pylab.subplot(2,2,2)
        pylab.title('Phase Map')
        pylab.imshow(self.phaseImage)
        pylab.colorbar()
#        pylab.subplot(2,2,4)
#        pylab.imshow(self.phaseImage1)
#        pylab.colorbar()

        pylab.show()



if __name__ == "__main__":
    ta=testAnalysis()
    ta.Analysis_FourierMap()
   