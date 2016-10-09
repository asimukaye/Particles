#Assignment1 Question3
#Vikas Kurapati
#130010058
from matplotlib import pyplot as plt
import numpy as np
delt = 0.01
n = 1000
#Euler Scheme
x1 = np.zeros((n))
y1 = np.zeros((n))
x2 = np.zeros((n))
y2 = np.zeros((n))
x1[0] = -0.5
y1[0] = 0.
x2[0] = 0.5
y2[0] = 0.
for i in range(1,n):
	x1[i] = x1[i-1] + (delt*(y1[i-1] - y2[i-1])/(((x1[i-1] - x2[i-1])*(x1[i-1] - x2[i-1])) + ((y1[i-1] - y2[i-1])*(y1[i-1] - y2[i-1]))))
	y1[i] = y1[i-1] - (delt*(x1[i-1] - x2[i-1])/(((x1[i-1] - x2[i-1])*(x1[i-1] - x2[i-1])) + ((y1[i-1] - y2[i-1])*(y1[i-1] - y2[i-1]))))
	x2[i] = x2[i-1] + (delt*(y2[i-1] - y1[i])/(((x2[i-1] - x1[i])*(x2[i-1] - x1[i])) + ((y2[i-1] - y1[i])*(y2[i-1] - y1[i]))))
	y2[i] = y2[i-1] - (delt*(x2[i-1] - x1[i])/(((x2[i-1] - x1[i])*(x2[i-1] - x1[i])) + ((y2[i-1] - y1[i])*(y2[i-1] - y1[i]))))
plt.plot(x1,y1)
plt.plot(x2,y2)
plt.title('Trajectory of the vortices in Eulers scheme with delta t = 0.01')
plt.show()
#Runge-Kutta Method
xr1 = np.zeros((n))
yr1 = np.zeros((n))
xr2 = np.zeros((n))
yr2 = np.zeros((n))
xr1[0] = -0.5
yr1[0] = 0
xr2[0] = 0.5
yr2[0] = 0
for j in range(1,n):
	kx11 = (delt*(yr1[j-1] - yr2[j-1])/(((xr1[j-1] - xr2[j-1])*(xr1[j-1] - xr2[j-1])) + ((yr1[j-1] - yr2[j-1])*(yr1[j-1] - yr2[j-1]))))
	kx12 = (delt*(yr1[j-1] - yr2[j-1])/(((xr1[j-1] +kx11 - xr2[j-1])*(xr1[j-1] +kx11 - xr2[j-1])) + ((yr1[j-1] - yr2[j-1])*(yr1[j-1] - yr2[j-1]))))
	xr1[j] = xr1[j-1] + 0.5*(kx11 + kx12)
	ky11 = (delt*(xr1[j] - xr2[j-1])/(((xr1[j] - xr2[j-1])*(xr1[j] - xr2[j-1])) + ((yr1[j-1] - yr2[j-1])*(yr1[j-1] - yr2[j-1]))))
	ky12 = (delt*(xr1[j] - xr2[j-1])/(((xr1[j] - xr2[j-1])*(xr1[j] - xr2[j-1])) + ((yr1[j-1] + ky11- yr2[j-1])*(yr1[j-1] +ky11 - yr2[j-1]))))
	yr1[j] = yr1[j-1] - 0.5*(ky11 + ky12)
	kx21 = (delt*(yr2[j-1] - yr1[j])/(((xr2[j-1] - xr1[j])*(xr2[j-1] - xr1[j])) + ((yr2[j-1] - yr1[j])*(yr2[j-1] - yr1[j]))))
	kx22 = (delt*(yr2[j-1] - yr1[j])/(((xr2[j-1] + kx21- xr1[j])*(xr2[j-1] + kx21 - xr1[j])) + ((yr2[j-1] - yr1[j])*(yr2[j-1] - yr1[j]))))
	xr2[j] = xr2[j-1] + 0.5*(kx21 + kx22)
	ky21 = (delt*(xr2[j] - xr1[j])/(((xr2[j] - xr1[j])*(xr2[j] - xr1[j])) + ((yr2[j-1] - yr1[j])*(yr2[j-1] - yr1[j]))))
	ky22 = (delt*(xr2[j] - xr1[j])/(((xr2[j] - xr1[j])*(xr2[j] - xr1[j])) + ((yr2[j-1] +ky21 - yr1[j])*(yr2[j-1] +ky21 - yr1[j]))))
	yr2[j] = yr2[j-1] - 0.5*(ky21 + ky22)
plt.plot(xr1,yr1)
plt.plot(xr2,yr2)
plt.title('Trajectory of the vortices in Runge-Kutta scheme with delta t = 0.01')
plt.show()