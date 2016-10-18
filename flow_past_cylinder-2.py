import numpy as np
import pymetabiosis.auto
import csv

def write_data(data, filename, type='complex'):
	file  = open(filename+'.csv', 'wb')
	writer = csv.writer(file, delimiter=',', quotechar=' ', quoting=csv.QUOTE_NONNUMERIC)
	for row in data:
		if type=='complex':
			row = [row.real, row.imag]
		elif type=='2d':
			row = row
		elif type=='coupled_data':
			row = [row[0].real, row[0].imag, row[1]]
		writer.writerow(row)

class Panel:
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.mid = 0.5*(start + end) 
        self.length = abs(start - end)
        self.angle = np.angle(end - start)
        ang = np.angle(end - start) + np.pi/2.0
        self.normal = np.exp(1.0j*ang)

def discretise_cylinder(n, r):
    ang_pos = np.linspace(0.0, 2.0*np.pi, n+1)[:-1]
    panel_start_pos = r*np.exp(1.0j*ang_pos)
    panel_end_pos = r*np.exp(1.0j*(ang_pos+2.0*np.pi/n))
    panels = []
    for i in range(n):
        panels.append(Panel(panel_start_pos[i], panel_end_pos[i]))
    return panels

def blob_vel(z, blob_pos, gamma, delta):
    r = abs(z-blob_pos)
    kr = r/delta if r < delta else 1.0
    kz = 1.0j*(z-blob_pos)/(2.0*np.pi*r*r) if r>0.0 else 0.0
    return kr*kz*gamma

def blob_blob_vel(blob_pos, gamma, delta):
    velocities = np.zeros_like(blob_pos)*(1.0+0.0j)
    for i,z in enumerate(blob_pos):
        for j,z1 in enumerate(blob_pos):
            velocities[i] = velocities[i] + blob_vel(z, z1, gamma[j], delta)
    return np.asarray(velocities)    

def panel_velocity(z, gamma1, gamma2, length):
    vel_gamma1 = -1.0j*gamma1*(1.0 + (z/length - 1.0)*np.log(1.0 - length/z))/(2.0*np.pi)
    vel_gamma2 = 1.0j*gamma2*(1.0 + (z/length)*np.log(1.0 - length/z))/(2.0*np.pi)
    return (vel_gamma1+vel_gamma2).conjugate()

def panel_point_vel(panels, gamma, pos):
    velocities = np.zeros_like(pos)*(1.0+0.0j)
    n = len(panels)
    for i, z in enumerate(pos):
        for j, panel in enumerate(panels):
            rotate = np.exp(1.0j*panel.angle)
            z_new = (z - panel.start)/rotate
            v = panel_velocity(z_new, gamma[j], gamma[(j+1)%n], panel.length)
            velocities[i] = velocities[i] + v*rotate
    return np.asarray(velocities)

def diffuse(num, mu, gamma, t):
    mean = 0.0
    sigma = np.sqrt(2.0*mu*t)
    z = np.random.normal(mean, sigma, num) + 1.0j*np.random.normal(mean, sigma, num)
    return z

def reflect_blobs(panel, blobs):
    reflected = []
    for z in blobs:
        slope = np.tan(panel.angle)
        x = panel.start.real
        y = panel.start.imag
        condition = (y - slope*x)*(z.real*slope - z.imag + y - slope*x) > 0
        if condition:
            rotate = np.exp(1.0j*panel.angle)
            z_new = (z - panel.mid)/rotate
            z_reflected = z_new.conjugate()*rotate + panel.mid
        else:
            z_reflected = z
        reflected.append(z_reflected)
    return np.array(reflected)

def correct_course(old_pos, new_pos, r):
    corrected = []
    for z1, z2 in zip(old_pos, new_pos):
        if abs(z2) < r:
            slope = np.tan(np.angle(z2-z1))
            a = slope*slope + 1.0
            b = 2.0*z1.imag*slope - 2.0*slope*slope*z1.real
            c = z1.imag*z1.imag + (slope*z1.real)**2 - r*r - 2.0*z1.real*z1.imag*slope
            x1 = (-b + np.sqrt(b*b - 4.0*a*c))/(2.0*a)
            y11 = np.sqrt(r*r - x1*x1) if r*r - x1*x1 > 0.0 else 0.0
            y12 = -np.sqrt(r*r - x1*x1) if r*r - x1*x1 > 0.0 else 0.0
            x2 = (-b - np.sqrt(b*b - 4.0*a*c))/(2.0*a)
            y21 =  np.sqrt(r*r - x2*x2) if r*r - x2*x2 > 0.0 else 0.0
            y22 = -np.sqrt(r*r - x2*x2) if r*r - x2*x2 > 0.0 else 0.0
            zc = []
            if abs(z1.imag - y11 + slope*(x1 - z1.real)) <= 1.0e-10:
                zc.append(x1+1.0j*y11)
            else:
                zc.append(x1+1.0j*y12)
            if abs(z1.imag - y21 + slope*(x2 - z1.real)) <= 1.0e-10:
                zc.append(x2+1.0j*y21)
            else:
                zc.append(x2+1.0j*y22) 
            zc = np.array(zc)
            i = np.argmin(abs(z1-zc))
            z_c = zc[i]
            if abs(z1-z_c) <= 1.0e-10:
                z_corrected = z1*abs(1.0 - z2/z1)
            else:
                alpha = np.pi - 2.0*(np.angle(z1-z_c) - np.angle(z_c))
                z_corrected = z_c + (z2-z_c)*np.exp(1.0j*alpha)
            corrected.append(z_corrected)
        else:
            corrected.append(z2)
    return np.array(corrected)

def matrix_a(n, r):
    panels = discretise_cylinder(n, r)
    a = np.zeros([n+1,n])
    for i in range(n):
        for j in range(n):
            rotate1 = np.exp(1.0j*panels[j%n].angle)               
            mid_new1 = (panels[i].mid - panels[j%n].start)/rotate1               
            v_ij_1 = -1.0j*(1.0 + (mid_new1/panels[j%n].length - 1.0)*np.log(1.0 - panels[j%n].length/mid_new1))/(2.0*np.pi)
            rotate2 = np.exp(1.0j*panels[(j-1)%n].angle)
            mid_new2 = (panels[i].mid - panels[(j-1)%n].start)/rotate2               
            v_ij_2 = 1.0j*(1.0 + (mid_new2/panels[(j-1)%n].length)*np.log(1.0 - panels[(j-1)%n].length/mid_new2))/(2.0*np.pi)
            v_ij = v_ij_1.conjugate()*rotate1 + v_ij_2.conjugate()*rotate2
            a[i][j] = (v_ij.conjugate()*panels[i].normal).real
    a[n][:] = 1.0
    return a

def matrix_b(n, r, v_fs, blob_pos, blob_str, delta):
    b = np.zeros([n+1, 1])
    panels = discretise_cylinder(n, r)
    for i in range(n):
        v_blob = 0.0 + 0.0j
        for j, z in enumerate(blob_pos):
            v_blob = v_blob + blob_vel(panels[i].mid, z, blob_str[j], delta)
        b[i] = ((-v_fs - v_blob)*((panels[i].normal).conjugate())).real
    return b

"""
def plot_vortices(num_panels, blob_pos, blob_str, r, t):
    positive_blobs = []
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
    theta = np.linspace(0.0, 2.0*np.pi, num_panels+1)
    plt.plot(r*np.cos(theta), r*np.sin(theta), label='Cylinder of radius ' + str(r))
    plt.legend()
    plt.axis([-2.0*length, 2.0*length, -length, length])
    plt.title("Position of Vortices at time = " + str(t) + 's', fontsize=28, fontweight='bold')
    outfile.savefig()    
    plt.close()
"""

def velocity_field(n, r, blob_pos, blob_str, panels, panel_str, delta, region):
    x, y = np.mgrid[region[0]:region[1]:n*1j, region[2]:region[3]:n*1j]
    z = x + 1.0j*y
    velocity = []
    for i in z:
        vel = []
        for j in i:
            v = 0.0+0.0j
            for k, gamma in enumerate(blob_str):
                v = v + blob_vel(j, blob_pos[k], gamma, delta)
            vel.append(abs(v))
        panel_vel = abs(panel_point_vel(panels, panel_str, i))
        velocity.append(np.array(vel) + panel_vel)
    velocity = np.array(velocity)
    for i, row in enumerate(z):
        for j, z1 in enumerate(row):
            if abs(z1) < r:
                velocity[i][j] = 0.0
    return velocity
    """
    plt.figure(figsize=(17.0, 9.0))
    plt.contour(x, y, velocity, n/2)
    plt.title("Velocity Field in Region x=%s to x=%s and y=%s y=%s" % region, fontsize=28, fontweight='bold')
    outfile.savefig()
    plt.close()
    """
        
def calculate_x_momentum(blob_pos, blob_str):
    x_momentum = sum(blob_str*blob_pos.imag)
    return x_momentum

def get_cd(t_array, x_momentum, v_fs):
    cd = []
    dt = t_array[1]-t_array[0]
    for i in range(len(x_momentum)-1):
        drag = (x_momentum[i+1]-x_momentum[i])/dt
        cd.append(drag*2.0/(abs(v_fs)**2))
    t_array = t_array[1:len(cd)+1]
    return np.array(t_array) + 1.0j*np.array(cd)
    """
    plt.figure(figsize=(17.0,9.0))
    plt.plot(t_array, cd)
    plt.title("$C_D$ vs Time", fontweight='bold', fontsize=28)
    plt.xlabel("Time in s", fontweight='bold', fontsize=22)
    plt.ylabel("$C_D$", fontweight='bold', fontsize=22)
    outfile.savefig()
    plt.close()
    """

def flow_past_cylinder(dt=0.1, tot_t=2.0001, n=50, r=1.0, v_fs=1.0+0.0j, gamma_max=0.1, delta=0.1, blob_pos=[], blob_str=[]):
    re = 1000.0
    mu = (v_fs.real)*2.0*r/re
    t=0.0
    t_array = [t]
    a = matrix_a(n, r)
    panels = discretise_cylinder(n, r)
    x_momentum = [0.0]
    while t <= tot_t:
        blob_prev_pos = blob_pos
        b = matrix_b(n, r, v_fs, blob_pos, blob_str, delta)
        panel_str = np.linalg.lstsq(a, b)[0]
        panel_str = np.ravel(panel_str)
        velocity = blob_blob_vel(blob_pos, blob_str, delta) + panel_point_vel(panels, panel_str, blob_pos)
        mid_pos = blob_pos + dt*velocity/2.0
        mid_pos = correct_course(blob_prev_pos, mid_pos, r)
        b = matrix_b(n, r, v_fs, mid_pos, blob_str, delta)
        panel_str = np.linalg.lstsq(a, b)[0]
        panel_str = np.ravel(panel_str)
        mid_velocity = blob_blob_vel(mid_pos, blob_str, delta) + panel_point_vel(panels, panel_str, mid_pos)
        blob_pos = blob_pos + dt*mid_velocity
        blob_pos = correct_course(blob_prev_pos, blob_pos, r)
        for i, panel in enumerate(panels):
            gamma = -2.0*panel.length*abs(v_fs)*np.sin(panel.angle + np.pi/2.0)
            if abs(gamma) >= gamma_max:
                num_diffused_blobs = int(abs(gamma)/gamma_max) + 1
                z = diffuse(num_diffused_blobs, mu, gamma, dt)
                z = z + panel.mid*(1 + panel.length/abs(panel.mid))
                z = reflect_blobs(panel, z)
                blob_pos = np.concatenate([blob_pos, z])
                gamma_diffused = np.zeros_like(z)
                gamma_diffused[:] = gamma/num_diffused_blobs
                blob_str = np.concatenate([blob_str, gamma_diffused])
        x_momentum.append(calculate_x_momentum(blob_pos, blob_str).real)
        if (round(t,2) in [1.0, 2.0, 3.0, 4.0, 5.0]): 
            write_data(zip(blob_pos, blob_str.real), 'time='+str(t), 'coupled_data')
        print t
        t = t + dt
        t_array.append(t)
    fullregion = velocity_field(n, r, blob_pos, blob_str, panels, panel_str, delta, (-2.0, 2.0, -2.0, 2.0))
    halfregion = velocity_field(n, r, blob_pos, blob_str, panels, panel_str, delta, (0.0, 2.0, -2.0, 2.0))
    write_data(fullregion, 'fullregion', '2d')
    write_data(halfregion, 'halfregion', '2d')
    x_momentum=np.array(x_momentum)
    x_momentum = (x_momentum[:-2] + x_momentum[1:-1] + x_momentum[2:])/3.0
    cd_data = get_cd(t_array, x_momentum, v_fs)
    write_data(cd_data, 'cd', 'complex')

if __name__ == '__main__':
    flow_past_cylinder()