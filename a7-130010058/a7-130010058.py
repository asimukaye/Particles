#Assignment 7
#Vikas Kurapati
#130010058

import numpy as np
#import pymetabiosis.auto
from matplotlib import pyplot as plt
import copy
from numpy import linalg as LA

def vortex_velocity(z, vor, gamma,delta):
	if abs(z - vor) < delta:
		if abs(z-vor) < 1.0e-10:
			compl_vel = 0.0+0.0j
		else:
			compl_vel = (-1.0j*gamma/(2.0*np.pi*(z - vor)))*(abs(z-vor)/delta)
	else:
		compl_vel = -1.0j*gamma/(2.0*np.pi*(z - vor))
	return compl_vel.conjugate()

def source_velocity(z, source, strength):
	compl_vel = strength/(2.0*np.pi*(z - source))
	return compl_vel.conjugate()

def panel_vel(z,gamma1,gamma2,l):
	if abs(z) > 1.0e-10:
		vel_gamma1 = (-1.0j*gamma1*(1. + (z/l - 1.)*np.log(1. - l/z))/(2.*np.pi)).conjugate()
		vel_gamma2 = (1.0j*gamma2*(1. + (z/l)*np.log(1. - l/z))/(2.0*np.pi)).conjugate()
		return vel_gamma1+vel_gamma2
	else:
		return 0.0+0.0j

def vortex_velocities_vortices(pos, gamma):
	vel = np.zeros_like(pos) + 1.0j*np.zeros_like(pos)
	for i, z_i in enumerate(pos):
		for j, z_j in enumerate(pos):
			if i!= j:
				vel[i] = vel[i] + vortex_velocity(z_i,z_j,gamma[j],0.1)
	return vel

def vortex_velocities_panels(pos,panels_matrix,gamma_matrix):
	vel = np.zeros_like(pos) + 1.0j*np.zeros_like(pos)
	for i,z_i in enumerate(pos):
		for j,panel_j in enumerate(panels_matrix):
			pos_r = (pos[i] - panel_j.p1)*np.exp(-1.0j*panel_j.angle)
			v = panel_vel(pos_r,gamma_matrix[j%len(panels_matrix)],gamma_matrix[(j+1)%len(panels_matrix)],panel_j.length)
			vel[i] = vel[i] + v[0]*np.exp(1.0j*panel_j.angle)			
	return vel

def point_velocities(Z,pos_vor,gamma_vor,panels_matrix,gamma_matrix):
	vel = np.zeros_like(Z) + 1.0j*np.zeros_like(Z)
	for ii,i in enumerate(Z):
		for jj,j in enumerate(i):
			for k,z_k in enumerate(pos_vor):
				vel[ii][jj] = vel[ii][jj] + vortex_velocity(j,z_k,gamma_vor[k],0.1)
			vel[ii] = vel[ii] + vortex_velocities_panels(i,panels_matrix,gamma_matrix)
	return vel

class Panel:
	"""Panel with initial,final,control points, length and angle subtends"""
	def __init__(self,p1,p2):
		self.p1 = p1
		self.p2 = p2
		self.length = abs(p1 - p2)
		self.control_point = .5*(p1+p2)
		self.angle = np.angle(p2 - p1)
		self.normal = np.exp(1.0j*np.angle(self.control_point))

def panelize(r,n):
	theta = np.linspace(0.0,2.0*np.pi,n,endpoint = False)
	complex_pos = r*np.exp(1.0j*theta)
	panels_matrix = np.empty(n,dtype = object)
	for i in range(n):
		panels_matrix[i] = Panel(complex_pos[i%n],complex_pos[(i+1)%n])
	return panels_matrix

def computeA(r=1.,n=20):
	panels_matrix = panelize(r,n)
	A = np.zeros([n+1,n],dtype = float)
	for i in range(n):
		for j in range(n):
			cp1_new = (panels_matrix[i%n].control_point - panels_matrix[j%n].p1)*np.exp(-1.0j*panels_matrix[j%n].angle)
			v_1 = -1.0j*(1. + (cp1_new/panels_matrix[j%n].length - 1.)*np.log(1. - panels_matrix[j%n].length/cp1_new))/(2.*np.pi)
			cp2_new = (panels_matrix[i%n].control_point - panels_matrix[(j-1)%n].p1)*np.exp(-1.0j*panels_matrix[(j-1)%n].angle)
			v_2 = 1.0j*(1. + (cp2_new/panels_matrix[(j-1)%n].length)*np.log(1. - panels_matrix[(j-1)%n].length/cp2_new))/(2.0*np.pi)
			v = v_1.conjugate()*np.exp(1.0j*panels_matrix[j%n].angle) + v_2.conjugate()*np.exp(1.0j*panels_matrix[(j-1)%n].angle)
			A[i][j] = (v.conjugate()*panels_matrix[i].normal).real
	A[n][:] = 1.
	return A

def computeB(r,n,pos_vor,gamma_vor,u_inf = 1. + 0.0j,v_b = 0. + 0.0j,pos_source = [0+0.0j],str_source =[0]):
	panels_matrix = panelize(r,n)
	B = np.zeros([n+1,1],dtype = float)
	for i in range(n):
		v_bn = (v_b.conjugate()*panels_matrix[i].normal).real
		v_fs = (u_inf.conjugate()*panels_matrix[i].normal).real
		v_s = 0.0 + 0.0j
		v_w = 0.0 + 0.0j
		for j, z_j in enumerate(pos_source):
			v_s += source_velocity(panels_matrix[i].control_point, z_j, str_source[j])
		v_s = (v_s.conjugate()*panels_matrix[i].normal).real
		for j,z_j in enumerate(pos_vor):
			v_w += vortex_velocity(panels_matrix[i].control_point, pos_vor[j], gamma_vor[j],0.1)
		v_w = (v_w.conjugate()*panels_matrix[i].normal).real
		B[i] = v_bn - v_fs - v_s - v_w
	return B
		
def diffusion(n,nu,tf):		#n number of vortices are already placed at the center.
	"""This returns the final position of n vortices placed at origin each with a strength of gamma/n"""
	x = np.random.normal(0.0,np.sqrt(2.0*nu*tf),n) + 1.0j*np.random.normal(0.0,np.sqrt(2.0*nu*tf),n)
	return x

def reflect_vor1(pos,r):
	pos_new = []
	for z in pos:
		if abs(z) < r:
			z = (2.0*r - abs(z))*np.exp(1.0j*np.angle(z))
			pos_new.append(z)
		else:
			pos_new.append(z)
	return np.array(pos_new)

def reflect_vor2(diff,panel):
	x = []
	for z_i in diff:
		k1 = panel.p1.real*(panel.p2.imag - panel.p1.imag) - panel.p1.imag*(panel.p2.real - panel.p1.real)
		k2 = (panel.p1.real - z_i.real)*(panel.p2.imag - panel.p1.imag) - (panel.p1.imag - z_i.imag)*(panel.p2.real - panel.p1.real)
		if  k1*k2 > 0.0:
			z_i_new = (z_i - panel.control_point)*np.exp(-1.0j*panel.angle)
			z_i = (z_i_new.conjugate())*np.exp(1.0j*panel.angle) + panel.control_point
			x.append(z_i)
		else:
			x.append(z_i)
	return np.array(x)

def cylinder_flow1(dt,tf,r,n,Re,u_inf):
	nu = (2.0*abs(u_inf)*r/Re)
	pos_vor = []
	gamma_vor = []
	vor_mom = np.array([0.0])
	t = 0.0
	T =[]
	while t <= tf:
		B = computeB(r,n,pos_vor,gamma_vor,u_inf = 1. + 0.0j,v_b = 0. + 0.0j,pos_source = [0.0+0.0j],str_source =[0.0])
		gamma_diff = LA.lstsq(A,B)[0]
		v = u_inf + vortex_velocities_vortices(pos_vor, gamma_vor) + vortex_velocities_panels(pos_vor,panels_matrix,gamma_diff)
		midpos = pos_vor + v*dt/2.
		midpos = reflect_vor1(midpos,r)
		B = computeB(r,n,midpos,gamma_vor,u_inf = 1. + 0.0j,v_b = 0. + 0.0j,pos_source = [0+0.0j],str_source =[0.0])
		gamma_matrix = LA.lstsq(A,B)[0]
		midv = u_inf + vortex_velocities_vortices(midpos, gamma_vor) + vortex_velocities_panels(midpos,panels_matrix,gamma_matrix)
		pos_vor = pos_vor + midv*dt
		pos_vor = reflect_vor1(pos_vor,r)
		for i,panel in enumerate(panels_matrix):
			strength = (gamma_diff[i] + gamma_diff[(i+1)%len(gamma_diff)])*panel.length/2.0
			if abs(strength) > 0.1:
				n_d = int(abs(strength)/0.1) + 1
				x = diffusion(n_d,nu,dt) + panel.control_point
				x = reflect_vor2(x,panel)	
				pos_vor = np.concatenate([pos_vor,x])
				small_gamma = strength/n_d
				gamma_vor = np.concatenate([gamma_vor,np.full(len(x),small_gamma)])

		for k in [1,3,5]:
			if abs(t - k) < 1e-9:
				positive = []
				negative = []
				for i,gam in enumerate(gamma_vor):
					if gam > 0:
						positive.append(pos_vor[i])
					else:
						negative.append(pos_vor[i])

				positive = np.array(positive)
				negative = np.array(negative)
				plt.figure(figsize = (17.0,9.0))
				plt.plot(positive.real,positive.imag,"o")
				plt.plot(negative.real,negative.imag,"o")
				panel_x = []
				panel_y = []
				for panel in panels_matrix:
					panel_x.append(panel.p1.real)
					panel_y.append(panel.p1.imag)
				panel_x.append(panels_matrix[0].p1.real)
				panel_y.append(panels_matrix[0].p1.imag)
				plt.plot(panel_x,panel_y,label = 'Panels')
				plt.axis([-8.0,8.0,-4.0,4.0])
				#plt.show()
				plt.savefig('method1vor'+str(k)+'.png')
				plt.close()	

		mom_ij = 0.0
		for i,z_i in enumerate(pos_vor):
			mom_ij = mom_ij + gamma_vor[i]*z_i.imag
		vor_mom = np.concatenate([vor_mom,[mom_ij]])
		T.append(t)
		t = t + dt
		print t

	X,Y = np.mgrid[-2.:2.:51j,-2.:2.:51j]
	Z = X + 1.0j*Y
	V = abs(u_inf + point_velocities(Z,pos_vor,gamma_vor,panels_matrix,gamma_matrix))
	for ii,i in enumerate(Z):
		for jj,j in enumerate(i):
			if abs(j) < r:
				V[ii][jj] = 0.0
	plt.contour(X,Y,V,50)
	plt.savefig('method1contour.png')
	plt.close()
	
	X,Y = np.mgrid[-2.:2.:51j,-2.:2.:51j]
	Z = X + 1.0j*Y
	V = abs(u_inf + point_velocities(Z,pos_vor,gamma_vor,panels_matrix,gamma_matrix))
	for ii,i in enumerate(Z):
		for jj,j in enumerate(i):
			if abs(j) < r or j.real<0:
				V[ii][jj] = 0.0
	plt.contour(X,Y,V,50)
	plt.savefig('method1contour1.png')
	plt.close()

	DRAG = (vor_mom[:-1] - vor_mom[1:])/dt
	cd = DRAG/(r*abs(u_inf)*abs(u_inf))
	plt.plot(T,cd)
	plt.savefig('method1cd.png')
	plt.close()

def cylinder_flow2(dt,tf,r,n,Re,u_inf):
	nu = (2.0*abs(u_inf)*r/Re)
	pos_vor = []
	gamma_vor = []
	vor_mom = np.array([0.0])
	t = 0.0
	T =[]
	while t <= tf:
		B = computeB(r,n,pos_vor,gamma_vor,u_inf = 1. + 0.0j,v_b = 0. + 0.0j,pos_source = [0.0+0.0j],str_source =[0.0])
		gamma_diff = LA.lstsq(A,B)[0]
		v = u_inf + vortex_velocities_vortices(pos_vor, gamma_vor) + vortex_velocities_panels(pos_vor,panels_matrix,gamma_diff)
		midpos = pos_vor + v*dt/4.
		midpos = reflect_vor1(midpos,r)
		B = computeB(r,n,midpos,gamma_vor,u_inf = 1. + 0.0j,v_b = 0. + 0.0j,pos_source = [0+0.0j],str_source =[0.0])
		gamma_matrix = LA.lstsq(A,B)[0]
		midv = u_inf + vortex_velocities_vortices(midpos, gamma_vor) + vortex_velocities_panels(midpos,panels_matrix,gamma_matrix)
		pos_vor = pos_vor + midv*dt/2.
		pos_vor = reflect_vor1(pos_vor,r)
		for i,panel in enumerate(panels_matrix):
			strength = (gamma_diff[i] + gamma_diff[(i+1)%len(gamma_diff)])*panel.length/2.0
			if abs(strength) > 0.1:
				n_d = int(abs(strength)/0.1) + 1
				x = diffusion(n_d,nu,dt) + panel.control_point
				x = reflect_vor2(x,panel)	
				pos_vor = np.concatenate([pos_vor,x])
				small_gamma = strength/n_d
				gamma_vor = np.concatenate([gamma_vor,np.full(len(x),small_gamma)])

		B = computeB(r,n,pos_vor,gamma_vor,u_inf = 1. + 0.0j,v_b = 0. + 0.0j,pos_source = [0.0+0.0j],str_source =[0.0])
		gamma_diff = LA.lstsq(A,B)[0]
		v = u_inf + vortex_velocities_vortices(pos_vor, gamma_vor) + vortex_velocities_panels(pos_vor,panels_matrix,gamma_diff)
		midpos = pos_vor + v*dt/4.
		midpos = reflect_vor1(midpos,r)
		B = computeB(r,n,midpos,gamma_vor,u_inf = 1. + 0.0j,v_b = 0. + 0.0j,pos_source = [0+0.0j],str_source =[0.0])
		gamma_matrix = LA.lstsq(A,B)[0]
		midv = u_inf + vortex_velocities_vortices(midpos, gamma_vor) + vortex_velocities_panels(midpos,panels_matrix,gamma_matrix)
		pos_vor = pos_vor + midv*dt/2.
		pos_vor = reflect_vor1(pos_vor,r)

		for k in [1,3,5]:
			if abs(t - k) < 1e-9:
				positive = []
				negative = []
				for i,gam in enumerate(gamma_vor):
					if gam > 0:
						positive.append(pos_vor[i])
					else:
						negative.append(pos_vor[i])

				positive = np.array(positive)
				negative = np.array(negative)
				plt.figure(figsize = (17.0,9.0))
				plt.plot(positive.real,positive.imag,"o")
				plt.plot(negative.real,negative.imag,"o")
				panel_x = []
				panel_y = []
				for panel in panels_matrix:
					panel_x.append(panel.p1.real)
					panel_y.append(panel.p1.imag)
				panel_x.append(panels_matrix[0].p1.real)
				panel_y.append(panels_matrix[0].p1.imag)
				plt.plot(panel_x,panel_y,label = 'Panels')
				plt.axis([-8.0,8.0,-4.0,4.0])
				#plt.show()
				plt.savefig('method2vor'+str(k)+'.png')
				plt.close()	

		mom_ij = 0.0
		for i,z_i in enumerate(pos_vor):
			mom_ij = mom_ij + gamma_vor[i]*z_i.imag
		vor_mom = np.concatenate([vor_mom,[mom_ij]])
		T.append(t)
		t = t + dt
		print t

	X,Y = np.mgrid[-2.:2.:51j,-2.:2.:51j]
	Z = X + 1.0j*Y
	V = abs(u_inf + point_velocities(Z,pos_vor,gamma_vor,panels_matrix,gamma_matrix))
	for ii,i in enumerate(Z):
		for jj,j in enumerate(i):
			if abs(j) < r:
				V[ii][jj] = 0.0
	plt.contour(X,Y,V,50)
	plt.savefig('method2contour.png')
	plt.close()
	
	X,Y = np.mgrid[-2.:2.:51j,-2.:2.:51j]
	Z = X + 1.0j*Y
	V = abs(u_inf + point_velocities(Z,pos_vor,gamma_vor,panels_matrix,gamma_matrix))
	for ii,i in enumerate(Z):
		for jj,j in enumerate(i):
			if abs(j) < r or j.real<0:
				V[ii][jj] = 0.0
	plt.contour(X,Y,V,50)
	plt.savefig('method2contour1.png')
	plt.close()

	DRAG = (vor_mom[:-1] - vor_mom[1:])/dt
	cd = DRAG/(r*abs(u_inf)*abs(u_inf))
	plt.plot(T,cd)
	plt.savefig('method2cd.png')
	plt.close()

if __name__ == '__main__':
	dt = 0.1	
	tf = 5.0
	r = 1.0
	n = 50
	Re = 1000.0
	panels_matrix = panelize(r,n)
	A = computeA(r,n) 
	u_inf = 1.0 + 0.0j
	cylinder_flow1(dt,tf,r,n,Re,u_inf)
	cylinder_flow2(dt,tf,r,n,Re,u_inf)