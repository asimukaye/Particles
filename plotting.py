import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import csv

def plot_vortices():
	for k in [1.0, 2.0]:
		file = open('time='+str(k)+'.csv', 'rb')
		reader = csv.reader(file, delimiter=',', quotechar=' ')
		blob_pos = []
		for row in reader:
			print "===%s"%row 
			blob_pos.append(float(row[0])+1.0j*float(row[1]))
			blob_str.append(float(row[2]))
		blob_str = np.array(blob_str)
		blob_pos = np.array(blob_pos)
		from IPython.core.debugger import Tracer; Tracer()()
		print k
		positive_blobs = []
		print k
		negative_blobs = []
		for i, gamma in enumerate(blob_str):
        	if gamma >= 0.0:
        		positive_blobs.append(blob_pos[i])
        	else:
        		negative_blobs.append(blob_pos[i])
    	positive_blobs = np.array(positive_blobs)
    	negative_blobs = np.array(negative_blobs)
    	max_x = abs(max(blob_pos.real))
    	max_y = abs(max(blob_pos.imag))
    	min_x = abs(min(blob_pos.real))
    	min_y = abs(min(blob_pos.imag))
    	length = max(max_x, max_y, min_x, min_y)
    	plt.figure(figsize=(17.0,9.0))
    	plt.plot(positive_blobs.real, positive_blobs.imag, '*', label='Positive strength blobs')
    	plt.plot(negative_blobs.real, negative_blobs.imag, 'r*', label='Negative strength blobs')
    	theta = np.linspace(0.0, 2.0*np.pi, 41)
    	plt.plot(np.cos(theta), np.sin(theta), label='Cylinder of radius 1.0')
    	plt.legend()
    	plt.axis([-2.0*length, 2.0*length, -length, length])
    	plt.title("Position of Vortices at time = " + str(k) + 's', fontsize=28, fontweight='bold')
    	outfile.savefig()
    	plt.close()    

if __name__ == '__main__':
	plt.ioff()
	plt.clf()
	outfile = PdfPages('report.pdf')
	plot_vortices()
	outfile.close()