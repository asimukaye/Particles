import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class mesh:
    def __init__(self, h = 0.05, x_max= 2.0, y_max= 2.0):
        self.h = h
        self.x_max = x_max
        self.y_max = y_max
        x = np.arange(-x_max, x_max+h, h)
        y = np.arange(-y_max, y_max+h, h)
        self.x_len = len(x)
        self.y_len = len(y)
        self.x, self.y = np.meshgrid(x, y)

def vorticity_actual(z, nu, t):
    a = 1/(4*nu*t)
    gamma = a* np.exp(-abs(z)**2/a)/np.pi
    return gamma

def generate_vortices(n, sigma=1.0,  mean=0.0):
    x = np.random.normal(mean, sigma, n)
    y = np.random.normal(mean, sigma, n)
    return x + 1j*y

def rvm(n, mesh, gamma_init=1.0, nu=0.1, t=1.0):
    sigma = np.sqrt(2*nu*t)
    z_vor = generate_vortices(n, sigma)
    gamma = np.zeros((mesh.x_len, mesh.y_len))
    for i in range(n):
        p = int(abs(z_vor[i].real - mesh.x[0,0])/mesh.h)
        q = int(abs(z_vor[i].imag - mesh.y[0,0])/mesh.h)
        if (p<n and q<n):
            del_x = z_vor[i].imag - mesh.x[p,q] 
            del_y = z_vor[i].imag - mesh.y[p,q]
            gamma[p, q]     += gamma_init/n*(mesh.h - del_x)*(mesh.h - del_y)
            gamma[p+1, q]   += gamma_init/n*(mesh.h + del_x)*(mesh.h - del_y)
            gamma[p, q+1]   += gamma_init/n*(mesh.h - del_x)*(mesh.h + del_y)
            gamma[p+1, q+1] += gamma_init/n*(mesh.h + del_x)*(mesh.h + del_y)
        elif(p>n and q<n):
            gamma[n, q] += gamma_init * mesh.h**2 /n
        elif(p<n and q>n):
            gamma[p, n] += gamma_init * mesh.h**2 /n
        else:
            gamma[n, n] += gamma_init * mesh.h**2 /n
    return gamma


def exact_sol(mesh, gamma_init=1.0, nu=0.1, t=1.0):
    z = mesh.x +1j*mesh.y
    gamma = gamma_init*vorticity_actual(z, nu, t)*mesh.h**2
    return gamma

def errorandplots(runs=5):
    m1 = mesh()
    gamma_exact = exact_sol(m1)
    n_array = [20, 50, 100, 200, 300, 500]
    plt.figure(1)
    CS = plt.contour(m1.x, m1.y, gamma_exact)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.title('Exact Gamma Distribution Contour Plot')
    pp.savefig()
    for i in range(runs):
        Err = []
        Gamma_rvm = []
        for j, n in enumerate(n_array):
            gamma_rvm = rvm(n, m1)
            Gamma_rvm.append(gamma_rvm)
            error = sum(sum(abs(gamma_rvm - gamma_exact)))/len(gamma_exact)**2
            Err.append(error)
        plt.figure(i+2)
        plt.scatter(m1.x, m1.y, s=2000*Gamma_rvm[5])  #2000 is just to scale the scatter size
        plt.title('RVM Gamma Distribution for Run '+ str(i+1)+' with 500 vortices')
        pp.savefig()
        plt.figure(runs+i+2)
        plt.plot(n_array, Err)
        plt.title('Error in Run '+ str(i+1))
        plt.xlabel('Number of vortices')
        plt.ylabel('Error')
        pp.savefig()
        #plt.close()

if __name__ == '__main__':
    global pp
    pp = PdfPages('a6_130010012.pdf')
    errorandplots()
    #plt.show()
    pp.close()