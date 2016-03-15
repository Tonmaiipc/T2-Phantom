import dicom
import numpy as np
import glob
from operator import itemgetter
import matplotlib.pyplot as plt


files = glob.glob("dicom/*.dcm")

## dicom loaders ##
dcm_stack = {}
i=0
for ele in files:
	dcm_stack['img{num:02d}'.format(num=i)] = dicom.read_file(ele); i+=1
	
	

######################################################	
########### collect log stack and find T2 ###########

def findT2(dcm_stack, x, y):
		
	te_val_stack, tesorigin, intsorigin = extract(dcm_stack, x, y)
		
	te_average_ints_dic = average(te_val_stack) # {[te, ints],[te, ints]}
	
	te_ints_log_stack = findlog(te_average_ints_dic)
	
	if te_average_ints_dic == {} : return 0, 0, 0, 0, 0, 0, 0
	
	teslist = [];	intslist = [];	loglist = []
	for tes,ints,log in te_ints_log_stack:
		teslist.append(tes)
		intslist.append(ints)
		loglist.append(log)
	
	teslist = np.array(teslist)
	intslist = np.array(intslist)
	loglist = np.array(loglist)
	
	x1 = np.ones_like(teslist) 
	x2 = -teslist
	y = loglist
	x = np.concatenate((x1[:, np.newaxis], x2[:,np.newaxis]), 1)
	beta, _,_,_ = np.linalg.lstsq(x,y) # give B1 & B2
	
	s0_ = np.exp(beta[0])
	t2_ = 1/beta[1]
	return s0_, t2_, teslist, intslist, loglist, tesorigin, intsorigin
	
	
######################################################	
################## data extract ##################

def extract(dcm_stack, x, y):
	te_val_stack = []
	
	for key, dicom_obj in dcm_stack.iteritems():
		intensity_stack = dicom_obj.pixel_array # screen&specitic x,y  give 
		
		te_val_stack.append((float(dicom_obj.EchoTime), intensity_stack[y][x]))
		
	tesorigin = []
	intsorigin = []
	for te, ints in te_val_stack:
		tesorigin.append(te)
		intsorigin.append(ints)
	return te_val_stack, tesorigin, intsorigin
		
	
	
######################################################	
######### extract pixel value and find log #########

def findlog(te_average_ints_dic): #file = dicom object

	te_ints_log_stack = []
	for tes, ints in te_average_ints_dic.iteritems():
		if ints == 0:
			log = 0
		else:
			log = np.log(ints)
		te_ints_log_stack.append([tes,ints,log])
	
	return te_ints_log_stack # [[tes,ints,log],[tes,ints,log]]
	
	
	
######################################################
################## value averager ##################


def average(te_val_stack):
	
	dic = {} # {te : ints, te : ints}
	#average
	for i in range(0, len(te_val_stack)):
		tes, ints = te_val_stack[i]
		if tes not in dic:
			dic[tes] = [ints]
		else:
			dic[tes].append(ints)
	for key,val in dic.iteritems():
		dic[key] = np.average(val)
	
	# cut off intensity 
	lis = []
	for te,ints in dic.iteritems():
		lis.append([te, ints])
	
	lis = sorted(lis,key=itemgetter(0)) # sort the lis of [te,ints] 
	
	cut_te_ints_lis = []
	count = 0
	border = np.absolute(lis[0][1] - lis[1][1])
	i = 0
	while count < 6:
		try:
			if  np.absolute(lis[i][1]-lis[i+1][1])/2 < border:
				cut_te_ints_lis.append([lis[i][0], lis[i][1]]) # [[te,ints],[te,ints]]
				count += 1
				i += 1
			else:
				count = 0
				i += 1
		except IndexError:
			break
			
	te_average_ints_dic = dict()
	
	for te, ints in cut_te_ints_lis:
		te_average_ints_dic[te] = ints
	
	return te_average_ints_dic
	

T2co = []
i = 0
for y in range(30,214):
	print 'calculating #:', (float(y)-30)/184*100., "%"
	for x in range(37,216):
		#print 'calculating at coordinate (x,y):', x,y
		_,t2,_,_,_,_,_ = findT2(dcm_stack, x, y)
		if t2 >= 0:
			T2co.append([(x,y), t2])
		else:
			T2co.append([(x,y), 0])
print 'done'
	

ds = dicom.read_file('result.dcm')

for (x,y), t2 in T2co:
	ds.pixel_array[y][x] = t2

ds.PixelData = ds.pixel_array.tostring()
ds.save_as("resultout.dcm")
