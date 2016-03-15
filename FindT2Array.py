import dicom
import numpy as np
import glob
from operator import itemgetter
import matplotlib.pyplot as plt
import sys


files = glob.glob("dicom/*.dcm")

#################### dicom loaders ####################

dcmobj_stack = {}
i=0
for ele in files:
	dcmobj_stack['{num:03d}'.format(num=i)] = dicom.read_file(ele); i+=1
	
	
te_vala_stack ={} # TE : value array 
for dcm,val in dcmobj_stack.iteritems():
	if val.EchoTime in te_vala_stack:
		te_vala_stack[val.EchoTime].append(val.pixel_array)
	else:
		te_vala_stack[val.EchoTime] = [val.pixel_array]

#################### average values in each (x,y) of TEs ####################
i = 0.
te_avgval_stack = {}
for te, val_lists in te_vala_stack.iteritems():
	i += 1.
	if len(val_lists) > 1:
		avg_val = np.zeros((256,256))
		#val_list[0] = [y[x] y[x] y[x] y[x]]  val_list[0][0] = [x x x] val_list[0][0][0] = x
		#va_lists[whole img][y][x]
		for y in range(0,256): # range(0,256) should be replaced with roi
			for x in range(0,256):
				avg = np.average( [ val_lists[img][y][x] for img in range(0,len(val_lists)) ] )
				avg_val[x][y] = avg
		te_avgval_stack[te] = avg_val
	print '\rAveraging pixel values : {0} %'.format(i/len(te_vala_stack)*100),

print 'done'

##################### screen the value used to calculate T2 #####################

sorted_tes = []
for key, val_array in te_avgval_stack.iteritems():
	sorted_tes.append(key)
sorted_tes = sorted(sorted_tes)

cut_te_avgval_stack = {}
for ele in sorted_tes:
	cut_te_avgval_stack[ele] = np.zeros((256,256))

for y in range(0,256): # range(0,256) should be replaced with roi
	print '\rScreening pixel values : {0} %'.format(y/256. * 100.), 
	for x in range(0,256):
		limit = np.absolute( (te_avgval_stack[sorted_tes[0]][y][x] - te_avgval_stack[sorted_tes[1]][y][x])/2 )
		if limit == 0: continue
		count = 0
		for i in range(0,len(sorted_tes)):
			if count >= 6: break
			try:
				width = np.absolute( te_avgval_stack[sorted_tes[i]][y][x] - te_avgval_stack[sorted_tes[i+1]][y][x] )
				if  width < limit:
					cut_te_avgval_stack[sorted_tes[i]][y][x] = te_avgval_stack[sorted_tes[i]][y][x]
					count += 1
				else:
					cut_te_avgval_stack[sorted_tes[i]][y][x] = te_avgval_stack[sorted_tes[i]][y][x]
					count = 0
			except IndexError:
				break
print 'TE avg # :', len(te_avgval_stack), 'screened TE avg # :', len(cut_te_avgval_stack)

# return cut_te_avgval_stack as {te: 'screened TE value'}

##################### log every values #####################

for te in sorted_tes:
	print '\rTaking logarithm : {0} %'.format(te/len(sorted_tes) * 100.),
	for y in range(0,256):
		for x in range(0,256):
			if cut_te_avgval_stack[te][y][x] > 0:
				cut_te_avgval_stack[te][y][x] = np.log(cut_te_avgval_stack[te][y][x])
			else:
				cut_te_avgval_stack[te][y][x] = 0
##################### calculate T2 #####################

t2_array = np.zeros((256,256))

teslist = np.array(sorted_tes)
x1 = np.ones_like(teslist)
x2 = -teslist

for y in range(0,256):
	print '\rFinding T2 : {0} %'.format(y/256. * 100.),
	for x in range(0,256):
		loglist = []
		for te in teslist:
			try:
				loglist.append(cut_te_avgval_stack[te][y][x])
			except IndexError:
				continue
			
		yt2 = loglist
		xt2 = np.concatenate((x1[:, np.newaxis], x2[:,np.newaxis]), 1)
		beta, _,_,_ = np.linalg.lstsq(xt2,yt2) 

		s0_ = np.exp(beta[0])
		if beta[1] <= 0:
			t2_ = 1000
		else:
			t2_ = 1/beta[1]
		t2_array[y][x] = format(t2_, '.4f')

ds = dicom.read_file('result.dcm')

for y in range(0,256):
	for x in range(0,256):
		if t2_array[y][x] > 1000:
			ds.pixel_array[y][x] = 1000
		else:
			ds.pixel_array[y][x] = t2_array[y][x]
		
ds.PixelData = ds.pixel_array.tostring()
ds.save_as("resultout.dcm")
