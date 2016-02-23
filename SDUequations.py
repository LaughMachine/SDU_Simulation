import numpy as np 
from scipy.stats import norm
import matplotlib.pyplot as plt

# Proposition 3: Fluid model, optimal allocation of nurses
def fm_b(r_I, r_S, lb, mu_C, mu_SC, p, w_C, w_SC):
	c_ratio = float(w_C)/float(w_SC)
	ratio = (r_I*p*mu_C + r_S*mu_SC)/float(r_I*mu_C)

	if c_ratio > ratio:
		b_I = min(r_I, float(lb)/mu_C )
	elif c_ratio <= ratio:
		b_I = float(r_I*r_S*mu_SC)/(r_I*mu_C*p + r_S*mu_SC)
	else:
		print 'Calculation Error'

	b_S = r_S*(1-(b_I/float(r_I)))

	#returns allocations
	return b_I, b_S

# Hazard Function used in subsequent functions
def hazard(x):

	return norm.pdf(x)/(1-norm.cdf(x))


# ---------------- THEOREM 2 (Erlang-A in Steady-State) ----------------
def E_Q_2(B, mu_C, theta):

	r = mu_C/float(theta)

	num = -1*(B*mu_C**.5)/theta + hazard(B*r**.5)*(1/theta)**.5
	denom = 1 + hazard(B*r**.5)/(hazard(-B)*r**.5)

	return num/denom


def E_I_2(B, mu_C, theta):

	r = mu_C/float(theta)

	denom = 1 + hazard(B*r**.5)/(hazard(-B)*r**.5)
	
	num = (B + hazard(-B)) / mu_C**.5

	return num*(1-denom**(-1))

#cost function
def C_2(B, mu_C, mu_SC, theta, w_SC, w_Q_C, r_S, r_I, p):

	p1 = w_Q_C*E_Q_2(B, mu_C, theta)

	p2 = B*p*mu_C**.5 + r_S*B*mu_SC/(r_I*mu_C**.5)

	p3 = (mu_SC + mu_C*p)*E_I_2(B, mu_C, theta)

	return p1+w_SC*(p2-p3)

#Bed allocation function
def B_2(B, lb, mu_C, mu_SC, r_I, r_S, p, N):

	x = lb/float(mu_C) + B*(lb/float(mu_C))**.5

	y = N*r_I*r_S*mu_SC/(r_I*mu_C*p + mu_SC*r_S)

	a = (N*r_I - lb/float(mu_C) - B*(lb/float(mu_C))**.5)*r_S/float(r_I)

	b = r_I*r_S*mu_C*p*N/(r_I*mu_C*p + r_S*mu_SC)
        
	return max(x, y), min(a, b)

#Optimal Beta function
def min_Beta_2(min1, max1, mu_C, mu_SC, theta, w_SC, w_Q_C, r_S, r_I, p):

	min_C = float('inf')
	min_B = 0

	for j in [float(i)/100 for i in range(100*min1, 100*max1+1)]:
		if C_2(j, mu_C, mu_SC, theta, w_SC, w_Q_C, r_S, r_I, p) <= min_C:
			min_C = C_2(j, mu_C, mu_SC, theta, w_SC, w_Q_C, r_S, r_I, p)
			min_B = j
	       
        
	return min_B

#All-in-one function that calculates the optimal beta and uses it to calculate bed allocation
def B_2_C(lb, mu_C, mu_SC, theta, w_SC, w_Q_C, r_I, r_S, p, N):

	max1 = int(np.ceil(100*(N*r_I/(lb/mu_C)**.5 - (lb/mu_C)**.5))/100)
	min1 = int(np.floor(-100*(lb/mu_C)**.5)/100)

	B = min_Beta_2(min1, max1, mu_C, mu_SC, theta, w_SC, w_Q_C, r_S, r_I, p)

	x = lb/float(mu_C) + B*(lb/float(mu_C))**.5

	y = N*r_I*r_S*mu_SC/(r_I*mu_C*p + mu_SC*r_S)

	a = (N*r_I - lb/float(mu_C) - B*(lb/float(mu_C))**.5)*r_S/float(r_I)

	b = r_I*r_S*mu_C*p*N/(r_I*mu_C*p + r_S*mu_SC)
        
	return max(x, y), min(a, b)

# ---------------- THEOREM 3 (ID Diffusion performance) ----------------
def E_Q_3(m, mu_C, theta, k):

	sig = (2*mu_C)**.5
	expmu = np.exp(m**2/(mu_C*sig**2))
	exptheta = np.exp(m**2/(theta*sig**2))

	NCDF1 = norm.cdf( (m/sig) * (2/theta)**.5 )
	NCDF2 = norm.cdf( (k + m/theta)*(2*theta)**.5/sig )
	NCDF3 = norm.cdf( (m/sig) * (2/mu_C)**.5 )

	num1 = np.exp( -1*(k**2 + 2*m*k/theta)*theta/sig**2 )
	num2 = (2/sig) * (np.pi/theta)**.5 * m * exptheta * (NCDF1 - NCDF2)

	num = 1 - num1 + num2

	denom1 = (1/(mu_C)**.5) * expmu * NCDF3
	denom2 = (1/(theta)**.5) * exptheta * (NCDF2 - NCDF1)

	denom = (2/sig)*np.pi**.5*(denom1 + denom2)

	return (1/(theta*mu_C**.5)) * num/denom

def E_I_3(m, mu_C, theta, k):

	sig = (2*mu_C)**.5
	expmu = np.exp(m**2/(mu_C*sig**2))
	exptheta = np.exp(m**2/(theta*sig**2))

	NCDF1 = norm.cdf( (m/sig) * (2/theta)**.5 )
	NCDF2 = norm.cdf( (k + m/theta)*(2*theta)**.5/sig )
	NCDF3 = norm.cdf( (m/sig) * (2/mu_C)**.5 )       
        
	num1 = (2/sig) * (np.pi/mu_C)**.5 * m * expmu * NCDF3

	num = (1/mu_C) * (1 + num1)

	denom1 = (1/(mu_C)**.5) * expmu * NCDF3
	denom2 = (1/(theta)**.5) * exptheta * (NCDF2 - NCDF1)

	denom = (2/sig)*np.pi**.5*(denom1 + denom2)

	return (1/mu_C**.5) * num / denom


def L_3(m, mu_C, theta, k):

	sig = (2*mu_C)**.5
	expmu = np.exp(m**2/(mu_C*sig**2))
	exptheta = np.exp(m**2/(theta*sig**2))

	NCDF1 = norm.cdf( (m/sig) * (2/theta)**.5 )
	NCDF2 = norm.cdf( (k + m/theta)*(2*theta)**.5/sig )
	NCDF3 = norm.cdf( (m/sig) * (2/mu_C)**.5 )

	num = np.exp( -1 *(k**2 + 2*m*k/theta)*theta/sig**2 )

	denom1 = (1/(mu_C)**.5) * expmu * NCDF3
	denom2 = (1/(theta)**.5) * exptheta * (NCDF2 - NCDF1)

	denom = (2/sig)*np.pi**.5*(denom1 + denom2)

	return (1/mu_C**.5) * num / denom


#cost function
def C_3(m, mu_C, mu_SC, theta, w_SC, w_Q_C, w_B_C, r_S, r_I, p, k):

	p1 = w_B_C*L_3(m, mu_C, theta, k)

	p2 = w_Q_C*E_Q_3(m, mu_C, theta, k)

	p3 = m*p / mu_C**.5 + m*r_S*mu_SC / (r_I*mu_C**1.5) - (mu_SC + mu_C*p)*E_I_3(m, mu_C, theta, k)

	return p1 + p2 + w_SC*p3

#Bed allocation function
def B_3(m, lb, mu_C, mu_SC, r_I, r_S, p, N):

	I1 = lb/mu_C + m/mu_C * (lb/mu_C)**.5
	I2 = r_I*r_S*mu_SC*N / (r_I*mu_C*p + mu_SC*r_S)


	S1 = (r_S/r_I) * (N*r_I - lb/mu_C - (m/mu_C) * (lb/mu_C)**.5 )
	S2 = r_I*r_S*mu_C*p*N / (r_I*mu_C*p + r_S*mu_SC)


	return max(I1, I2), min(S1, S2)

# Optimal m and k function
def min_m_k_3(minm, maxm, mink, maxk, mu_C, mu_SC, theta, w_SC, w_Q_C, w_B_C, r_S, r_I, p):

	min_C = float('inf')
	min_m = 0
	min_k = []

	for m in [float(i1)/100 for i1 in range(minm*100, maxm*100+1)]:
		for k in [float(i0)/10 for i0 in range(mink*10, maxk*10+1)]:

			if C_3(m, mu_C, mu_SC, theta, w_SC, w_Q_C, w_B_C, r_S, r_I, p, k) <= min_C:
				if C_3(m, mu_C, mu_SC, theta, w_SC, w_Q_C, w_B_C, r_S, r_I, p, k) == min_C:
					min_k.append(k)
				else:
					min_k = []
					min_k.append(k)
				min_C = C_3(m, mu_C, mu_SC, theta, w_SC, w_Q_C, w_B_C, r_S, r_I, p, k)
				min_m = m
        
	mk = min(min_k)    
	
	return min_m, mk

#All-in-one function that calculates the optimal m, k and uses it to calculate bed allocation
def B_3_C(lb, mu_C, mu_SC, theta, w_SC, w_Q_C, w_B_C, r_I, r_S, p, N):

	mink = 0
	maxk = 4
	minm = int(np.floor(-100*(lb/mu_C)**.5)/100)
	maxm = int(np.ceil(100*(N*r_I*mu_C/(lb/mu_C)**.5 - (lb*mu_C)**.5))/100)

	m, k = min_m_k_3(minm, maxm, mink, maxk, mu_C, mu_SC, theta, w_SC, w_Q_C, w_B_C, r_S, r_I, p)
	I1 = lb/mu_C + m/mu_C * (lb/mu_C)**.5
	I2 = r_I*r_S*mu_SC*N / (r_I*mu_C*p + mu_SC*r_S)


	S1 = (r_S/r_I) * (N*r_I - lb/mu_C - (m/mu_C) * (lb/mu_C)**.5 )
	S2 = r_I*r_S*mu_C*p*N / (r_I*mu_C*p + r_S*mu_SC)


	return max(I1, I2), min(S1, S2), k

# ---------------- Proposition 4 (balking threshold in the CD regime) ----------------
def g_4(lb, mu_SC, mu_C, r_I, r_S, p, N):

	num = N*r_I*r_S*mu_SC*mu_C 

	denom = lb*(r_I*mu_C*p + r_S*mu_SC)

	return num/denom


def d_4(B, lb, mu_C, mu_SC, r_I, r_S, p, N):

	num = -1*(N/lb)**.5 * B*r_I*mu_C*mu_SC * (r_I*r_S*p/(r_I*mu_C*p + r_S*mu_SC))**.5

	denom = r_I*mu_C*p + r_S*mu_SC

	return num/denom

#cost function
def C_4(B, mu_SC, mu_C, r_I, r_S, w_C, w_SC, p):

	c1 = mu_SC * (r_I*r_S*mu_C*p/(r_I*mu_C*p + r_S*mu_SC))**.5

	c2 = w_C*B*r_I*mu_C / (r_I*mu_C*p + r_S*mu_SC) + w_SC*hazard(-B) 

	return c1*c2 

#Bed allocation function
def B_4(B, lb, mu_C, mu_SC, r_I, r_S, p, N):

	RI = lb/mu_C

	I1 = g_4(lb, mu_SC, mu_C, r_I, r_S, p, N)*RI + d_4(B, lb, mu_C, mu_SC, r_I, r_S, p, N)*RI**.5
	I2 = r_I*N
	I3 = lb/mu_C

	B_I = min(max(I1,0), I2, I3)

	RS = B_I*mu_C*p/mu_SC

	S1 = RS + B*RS**.5
	S2 = (r_S/r_I) * (N*r_I - lb/mu_C)

	B_S = max(S1, S2)
	
	return B_I, B_S

#Optimal beta function
def min_Beta_4(min1, max1, mu_C, mu_SC, w_SC, w_C, r_S, r_I, p):

	min_C = float('inf')
	min_B = 0

	for j in [float(i)/100 for i in range(min1*100, max1*100)]:

	   if C_4(j, mu_SC, mu_C, r_I, r_S, w_C, w_SC, p) <= min_C:
	       ICU, SDU = B_4(j, lb, mu_C, mu_SC, r_I, r_S, p, N)
	       if ICU >= 0 and SDU >= 0:
	           min_C = C_4(j, mu_SC, mu_C, r_I, r_S, w_C, w_SC, p)
	           min_B = j

	return min_B

#All-in-one function that calculates the optimal beta and uses it to calculate bed allocation
def B_4_C(lb, mu_C, mu_SC, w_SC, w_C, r_I, r_S, p, N):

	max1 = 28
	min1 = -3
	B = min_Beta_4(min1, max1, mu_C, mu_SC, w_SC, w_C, r_S, r_I, p)

	RI = lb/mu_C

	I1 = g_4(lb, mu_SC, mu_C, r_I, r_S, p, N)*RI + d_4(B, lb, mu_C, mu_SC, r_I, r_S, p, N)*RI**.5
	I2 = r_I*N
	I3 = lb/mu_C

	B_I = min(max(I1,0), I2, I3)

	RS = B_I*mu_C*p/mu_SC

	S1 = RS + B*RS**.5
	S2 = (r_S/r_I) * (N*r_I - lb/mu_C)

	B_S = max(S1, S2)
	
	return B_I, B_S

# ---------------- Simulation Adjustment ----------------
#Adjusts number of beds in the SDU for simulation use
def simAdj_B_S(B_S, lb_SDU, mu_SC):

	return B_S + lb_SDU/float(mu_SC)


#___________________ Test ____________________________

mu_C = 1/4.8
mu_SC = 1/2.3
theta = 1
w_SC = 1
w_Q_C = 3
w_B_C = 15
r_S = 4
r_I = 2
p = .8
N = 20
w_C = min(w_Q_C, w_B_C)
lb = N*r_I*mu_C
T = (r_I*p*mu_C + r_S*mu_SC)/float(r_I*mu_C)

print fm_b(r_I, r_S, lb, mu_C, mu_SC, p, w_C, w_SC)
print T


# mu_C = 1/4.8
# mu_SC = 1/2.3
# theta = 1
# w_SC = 1
# w_Q_C = 15
# w_B_C = 5.5
# r_S = 4
# r_I = 2
# p = .8
# N = 20
# w_C = min(w_Q_C, w_B_C)
# lb = N*r_I*mu_C
# T = (r_I*p*mu_C + r_S*mu_SC)/float(r_I*mu_C)

# print T

# max1 = int(np.ceil(100*(N*r_I/(lb/mu_C)**.5 - (lb/mu_C)**.5))/100)
# min1 = int(np.floor(-100*(lb/mu_C)**.5)/100)

# mink = 0
# maxk = 4
# minm = int(np.floor(-100*(lb/mu_C)**.5)/100)
# maxm = int(np.ceil(100*(N*r_I*mu_C/(lb/mu_C)**.5 - (lb*mu_C)**.5))/100)
# min4 = -3
# max4 = 28


# min4 = -g_4(lb, mu_SC, mu_C, r_I, r_S, p, N)*(lb/mu_C)**.5 / d_4(1, lb, mu_C, mu_SC, r_I, r_S, p, N)


# B2 = min_Beta_2(min1, max1, mu_C, mu_SC, theta, w_SC, w_Q_C, r_S, r_I, p)
# m, k = min_m_k_3(minm, maxm, mink, maxk, mu_C, mu_SC, theta, w_SC, w_Q_C, w_B_C, r_S, r_I, p)
# B4 = min_Beta_4(min4, max4, mu_C, mu_SC, w_SC, w_C, r_S, r_I, p)

# fBI, fBS = fm_b(r_I, r_S, lb, mu_C, mu_SC, p, w_C, w_SC)
# BI2, BS2 = B_2(B2, lb, mu_C, mu_SC, r_I, r_S, p, N)
# BI3, BS3 = B_3(m, lb, mu_C, mu_SC, r_I, r_S, p, N)
# BI4, BS4 = B_4(B4, lb, mu_C, mu_SC, r_I, r_S, p, N)



# B2 = min_Beta_2(min1, max1, mu_C, mu_SC, theta, w_SC, 6.9, r_S, r_I, p)
# BI2, BS2 = B_2(B2, lb, mu_C, mu_SC, r_I, r_S, p, N)
# B4_Q = min_Beta_4(min4, max4, mu_C, mu_SC, w_SC, 6.9, r_S, r_I, p)
# BI4_Q, BS4_Q = B_4(B4_Q, lb, mu_C, mu_SC, r_I, r_S, p, N)
# B_T_Q = (BI2 + BI4_Q)/2

# m, k = min_m_k_3(minm, maxm, mink, maxk, mu_C, mu_SC, theta, w_SC, 15, 6.9, r_S, r_I, p)
# BI3, BS3 = B_3(m, lb, mu_C, mu_SC, r_I, r_S, p, N)
# B4_B = min_Beta_4(min4, max4, mu_C, mu_SC, w_SC, 6.9, r_S, r_I, p)
# BI4_B, BS4_B = B_4(B4_B, lb, mu_C, mu_SC, r_I, r_S, p, N)
# B_T_B = (BI2 + BI4_B)/2

# print B_T

# print 'Fluid solution:'
# print 'ICU Beds: ' + str(fBI)
# print 'SDU Beds: ' + str(fBS)
# print '---------------------'
# print 'ID regime, Queue Dominated:'
# print 'ICU Beds: ' + str(BI2)
# print 'SDU Beds: ' + str(BS2)
# print '---------------------'
# print 'ID regime, Balking Dominated:'
# print 'ICU Beds: ' + str(BI3)
# print 'SDU Beds: ' + str(BS3)
# print '---------------------'
# print 'CD regime, Balking Dominated:'
# print 'ICU Beds: ' + str(BI4)
# print 'SDU Beds: ' + str(BS4)

# QDI = []
# QDS = []
# BDI = []
# BDS = []
# KB = []
# xaxis = []
# xaxis1 = []

# print B_4_C(lb, mu_C, mu_SC, w_SC, w_C, r_I, r_S, p, N)[0]

# Balking Case
# for x in range(1,51):
# 	w_B_C = x*.2
# 	w_Q_C = 15
# 	w_C = min(w_Q_C, w_B_C)
# 	T = (r_I*p*mu_C + r_S*mu_SC)/float(r_I*mu_C)
# 	if w_C < T:
# 		B4 = min_Beta_4(min4, max4, mu_C, mu_SC, w_SC, w_C, r_S, r_I, p)
# 		BI4, BS4 = B_4(B4, lb, mu_C, mu_SC, r_I, r_S, p, N)
# 		BI4 = np.round(min(BI4, B_T_B))
# 		BS4 = simAdj_B_S((N-BI4)*3, lb, mu_SC)
		
# 		BDI.append(np.round(BI4))
# 		BDS.append(np.round(BS4))
		
# 		KB.append(0)
# 		xaxis.append(w_B_C)
# 	else:
# 		B2 = min_Beta_2(min1, max1, mu_C, mu_SC, theta, w_SC, w_C, r_S, r_I, p)
# 		BI2, BS2 = B_2(B2, lb, mu_C, mu_SC, r_I, r_S, p, N)
# 		m, k = min_m_k_3(minm, maxm, mink, maxk, mu_C, mu_SC, theta, w_SC, w_Q_C, w_B_C, r_S, r_I, p)
# 		BI3, BS3 = B_3(m, lb, mu_C, mu_SC, r_I, r_S, p, N)
# 		BI3 = np.round(max(BI3, B_T_B))
# 		BS3 = simAdj_B_S((N-BI3)*3, lb, mu_SC)

# 		BDI.append(np.round(BI3))
# 		BDS.append(np.round(BS3))
		
# 		KB.append(k*20**.5)
# 		xaxis.append(w_B_C)
# print BDI
# print BDS
# print KB
# print xaxis

# # Queue Case
# for x in range(1,51):
# 	w_Q_C = x*.2
# 	w_B_C = 15
# 	w_C = min(w_Q_C, w_B_C)
# 	T = (r_I*p*mu_C + r_S*mu_SC)/float(r_I*mu_C)
# 	if w_C < T:
# 		B4 = min_Beta_4(min4, max4, mu_C, mu_SC, w_SC, w_C, r_S, r_I, p)
# 		BI4, BS4 = B_4(B4, lb, mu_C, mu_SC, r_I, r_S, p, N)
# 		BI4 = np.round(min(BI4, B_T_Q))
# 		BS4 = simAdj_B_S((N-BI4)*3, lb, mu_SC)
		
# 		QDI.append(np.round(BI4))
# 		QDS.append(np.round(BS4))

# 		xaxis1.append(w_Q_C)
# 	else:
# 		B2 = min_Beta_2(min1, max1, mu_C, mu_SC, theta, w_SC, w_C, r_S, r_I, p)
# 		BI2, BS2 = B_2(B2, lb, mu_C, mu_SC, r_I, r_S, p, N)
# 		BI2 = np.round(max(BI2, B_T_Q))
		
# 		BS2 = simAdj_B_S((N-BI2)*3, lb, mu_SC)
		
# 		QDI.append(np.round(BI2))
# 		QDS.append(np.round(BS2))
		
# 		xaxis1.append(w_Q_C)

# print QDI
# print QDS
# print xaxis1


# fig = plt.figure(figsize=(16.0, 10.0))
# ax = fig.add_subplot(111)
# plt.cla()
# plt.plot(xaxis, QDI, color = 'r', marker = 'o')
# plt.plot(xaxis, QDS, color = 'b', marker = 'o')
# plt.show()