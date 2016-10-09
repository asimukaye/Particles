#Assignment 1 Question 1
#Vikas Kurapati
#130010058
import numpy as np
from matplotlib import pyplot as plt
x=np.linspace(-2.,2.,2000)
y=np.linspace(-2.,2.,2000)
X,Y = np.meshgrid(x,y)
z=X+1j*Y
compl = z + (np.log(z+1) - np.log(z-1))/(2*np.pi)
pot = compl.real
stream = compl.imag
plt.figure(1)
POT = plt.contour(X,Y,pot)
plt.title('Potential lines of the flow')
plt.show()
plt.figure(2)
STREAM = plt.contour(X,Y,stream)
plt.title('Streamlines of the flow')
plt.show()