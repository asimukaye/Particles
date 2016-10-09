import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def random_vortices(num, mu, gamma, t):
    mean = 0.0
    sigma = np.sqrt(2.0*mu*t)
    z = np.random.normal(mean, sigma, num) + 1.0j*np.random.normal(mean, sigma, num)
    return z

def plot_vortices():
    z = random_vortices(100, 0.1, 1.0, 1.0)
    plt.figure(figsize=(15.0,9.0))
    plt.plot(z.real, z.imag,'o')
    plt.title("Distribution of 50 Vortices by Random Normal Distribution", fontsize=22, fontweight='bold')
    outfile.savefig()    
    plt.close()

def meshing(grid_num, num_of_vort, mu, gamma, t):
    z = random_vortices(num_of_vort, mu, gamma, t)
    x, y = np.mgrid[-1.5 : 1.5 : (grid_num+1)*1j, -1.5 : 1.5 : (grid_num+1)*1j]  
    return np.ravel(x + 1.0j*y), z

class Vortex_Boundary:
    def __init__(self, grid, vort_pos, h):
        x = -1.5 + h*int((vort_pos.real + 1.5)/h)
        y = -1.5 + h*int((vort_pos.imag + 1.5)/h)
        self.z1 = x + 1.0j*y
        self.z2 = x+h + 1.0j*y
        self.z3 = x + 1.0j*(y+h)
        self.z4 = x+h + 1.0j*(y+h)
        
def vorticity_distribution(grid, vort_pos, gamma, h):
    vort_str = gamma/(len(vort_pos))
    d = {}
    for z in vort_pos:
        Vortex = Vortex_Boundary(grid, z, h)
        z1 = np.round(Vortex.z1, 6)
        z2 = np.round(Vortex.z2, 6)
        z3 = np.round(Vortex.z3, 6)
        z4 = np.round(Vortex.z4, 6)
        x = z.real - z1.real
        y = z.imag - z1.imag
        d[z1] = d[z1] + vort_str*(h-x)*(h-y)   if z1 in d.keys() else vort_str*(h-x)*(h-y)
        d[z2] = d[z2] + vort_str*x*(h-y)       if z2 in d.keys() else vort_str*x*(h-y)
        d[z3] = d[z3] + vort_str*(h-x)*y       if z3 in d.keys() else vort_str*(h-x)*y
        d[z4] = d[z4] + vort_str*x*y           if z4 in d.keys() else vort_str*x*y
    return np.array(d.values()), np.array(d.keys())

def actual_vorticity_distribution(grid, mu, t, gamma, h):
    exact_vorticity = np.exp(-(abs(grid)**2)/(4.0*mu*t))/(4.0*np.pi*mu*t)
    return gamma*exact_vorticity*h*h

def vorticity_error(grid_num=100, num_of_vort=100, mu=0.1, gamma=1.0, t=1.0):
    grid, vort_pos = meshing(grid_num, num_of_vort, mu, gamma, t)
    h = 3.0/grid_num
    numerical, new_grid = vorticity_distribution(grid, vort_pos, gamma, h)
    actual = actual_vorticity_distribution(new_grid, mu, t, gamma, h)
    error = abs(numerical - actual)
    avg_error = sum(error)/len(error) 
    return avg_error

def plot1():
    error = []
    num = [100, 150, 200, 500, 1000]
    for i in num:
        error.append(vorticity_error(num_of_vort=i))
    plt.figure(figsize=(15.0,9.0))
    plt.plot(num, error)
    plt.title("Variation of Error with number of vortices", fontsize=22, fontweight='bold')
    plt.xlabel("Number of Vortices", fontsize=18, fontweight='bold')
    plt.ylabel("Average Error", fontsize=18, fontweight='bold')
    outfile.savefig()
    plt.close()

def plot2():
    num = np.array([10, 50, 100, 200, 500])
    h = 3.0/num
    error = []
    for i in num:
        error.append(vorticity_error(grid_num=i))
    plt.figure(figsize=(15.0,9.0))
    plt.plot(h, error)
    plt.title("Variation of Error with grid size", fontsize=22, fontweight='bold')    
    plt.ylabel("Average Error", fontsize=18, fontweight='bold')
    plt.xlabel("Grid size", fontsize=18, fontweight='bold')
    outfile.savefig()
    plt.close()

if __name__ == '__main__':
    plt.ioff()
    plt.clf()
    outfile = PdfPages('report.pdf')
    plot_vortices()
    plot1()
    plot2()
    outfile.close()