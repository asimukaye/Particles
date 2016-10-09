#Assignment 1 Question 2
#Vikas Kurapati
#130010058
from matplotlib import pyplot as plt
import numpy as np
def stream(x,y):
	strm = y + (0.5/np.pi)*(np.arctan2(y,(x+1)) - np.arctan2(y,(x-1)))
	return strm
delt = 0.001
delx = 0.00001
dely = 0.00001
n = 5000
#Euler Scheme
x = np.zeros((11,n))
y = np.zeros((11,n))
x[:,0] = -2
y[:,0] = np.linspace(-2,2,11)
for i in range(1,n):
	x[:,i] = x[:,i-1] + (0.5*(1/dely)*delt*(stream(x[:,i-1],(y[:,i-1]+dely)) - stream(x[:,i-1],(y[:,i-1]-dely))))
	y[:,i] = y[:,i-1] - (0.5*(1/delx)*delt*(stream((x[:,i-1]+delx),y[:,i-1]) - stream((x[:,i-1]-delx),y[:,i-1])))
for j in range(11):
	plt.plot(x[j,:],y[j,:])
plt.title('Trajectory of the particles in Eulers scheme')
plt.show()
#Runge-Kutta Scheme
x0 = np.zeros((11,n))
y0 = np.zeros((11,n))
x0[:,0] = -2
y0[:,0] = np.linspace(-2,2,11)
for i in range(1,n): 
	kx1 = delt*(0.5*(1/dely)*(stream(x0[:,i-1],(y0[:,i-1]+dely)) - stream(x0[:,i-1],(y0[:,i-1]-dely))))
	kx2 = delt*(0.5*(1/dely)*(stream((x0[:,i-1]+kx1),(y0[:,i-1]+dely)) - stream((x0[:,i-1]+kx1),(y0[:,i-1]-dely))))
	ky1 = -delt*(0.5*(1/delx)*(stream((x0[:,i-1]+delx),y0[:,i-1]) - stream((x0[:,i-1]-delx),y0[:,i-1])))
	ky2 = -delt*(0.5*(1/delx)*(stream((x0[:,i-1]+delx),(y0[:,i-1]+ky1)) - stream((x0[:,i-1]-delx),(y0[:,i-1]+ky1))))
	x0[:,i] = x0[:,i-1] + 0.5*(kx1 + kx2)
	y0[:,i] = y0[:,i-1] + 0.5*(ky1 + ky2)
for j in range(11):
	plt.plot(x0[j,:],y0[j,:])
plt.title('Trajectory of the particles in Runge-Kutta scheme')
plt.show()