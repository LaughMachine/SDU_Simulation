#!/usr/bin/env /apps/python2.7/bin/python

# Any questions email mh3166@columbia.edu

import numpy as np
import csv
import time
import os 
import sys
import heapq

fil = open(os.getcwd() + "/SimResults_B_" + str(sys.argv[16]) + ".csv","wb")

def simulation(T, Tw, B_I, B_S, K_len, p1, p2, p3, p4, p5, mu_c, sig_c, mu_sc, sig_sc, mult, lb_C, lb_S, theta, theta_ret):
	# ============================== Variables ==============================

	# ---------- Simulation variables ----------
	TotalT = T		#Total time to run simulation
	t_warmup = Tw	#Time to let simulation run before data collection

	# ---------- Bed/Queue sizes ----------
	Beds = 30
	Ibeds = B_I 			#Beds in ICU
	Sbeds = B_S				#Beds in SDU
	K = K_len				#Length of critical patients queue
	 
	# ---------- Probabilities ----------
	p_c_sc = p1			#Probability critical patient goes to semicritical
	p_c_c_d = p2		#Probability critical patient returns as critical
	p_c_sc_d = p3		#Probability critical patient returns as semi-critical
	p_sc_c = p4			#Probability semi-critical patient returns as critical
	p_sc_sc = p5		#Probability semi-critical patient returns as semi-critical

	# ---------- Lognormal distribution ----------
	cmu = mu_c			#Target mean for critical patient service time
	csig = sig_c		#Target standard deviation for critical patient service time
	scmu = mu_sc 		#Target mean for semi-critical patient service time
	scsig = sig_sc		#Target standard deviation for semi-critical patient service time

	crit_mean = np.log(cmu**2/(csig**2+cmu**2)**.5)			#Mean service time of critical patients
	crit_std = np.log(1+csig**2/cmu**2)**.5					#Standard deviation of service time of critical patients

	semicrit_mean = np.log(scmu**2/(scsig**2+scmu**2)**.5)		#Mean service time of semi-critical patients
	semicrit_std = np.log(1+scsig**2/scmu**2)**.5				#Standard deviation service time of semi-critical patients

	SDUcrit_mean = 3.5					#Mean service time of critical patients in SDU
	SDUcrit_std = SDUcrit_mean*1.5		#Standard deviation service time of critical patients in SDU

	multiplier = mult

	# ---------- Exponential distribution ----------
	arr_rate_C = 1/float(lb_C) 		#Arrival rate of critical patients (mean)
	arr_rate_S = 1/float(lb_S) 		#Arrival rate of semi-critical patients (mean)

	aband_rate = theta 			#Abandonment rate
	ret_delay = theta_ret 		#Length of delay before return 

	patience = theta			#Abandonment rate of critical patient in SDU


	# ============================== Queues ==============================
	ICUqueue = [] 			#Holds arrival times of critical patients that are in the queue; Necessary for calculating wait time
	ICUqueueAban = []		#Holds the abandonment times of each patient in line
	ICUqueue_count = 0

	ICUreturnq = [] 		#Holds arrival times of when patients return in a critical state; calculates waiting time
	ICUreturnqAban = []		#Holds the abandonment times of each patient in return line
	ICUreturnq_count = 0

	SDUqueue = [] 			#Holds arrival time they are in the 'outside'

	RetDelayQ = [] 			#If patients return to semi-critical or critical condition after treatment they enter this queue

	# ============================== Beds ==============================
	ICUbeds_C = [] 						#Holds tuples representing critical patients
	ICUbeds_S = []						#Holds tuples representing semi-critical patients
	ICUbeds_S_queue = []				#Holds the order of arrivals of semi-critical patients into ICU
	ICUcount = 0 						#Counts number of ICU beds occupied

	SDUbeds_C = [] 						#Holds tuples representing critical patients
	SDUbeds_S = []						#Holds tuples representing semi-critical patients
	SDUcount = 0 						#Counts number of SDU beds occupied

	# ============================== Simulation ==============================

	# ------------------------------ Initialization ------------------------------
	# ________ Counters ________
	t = 0		        # Keeps track of current time in simulation
	Tmax = TotalT
	ArrivalCount = 0        # Keeps track of the number of arrivals
	BalkCount = 0			# Keeps track of the number of balks
	BumpCount = 0 			# Keeps track of the number of bumps
	AbandCount = 0 			# Keeps track of the number of abandonments
	waittime = 0 			# Keeps track of total wait time
	CritArrival = 0 		# Keeps track of the number of critical patients
	SemiCritArrival = 0     # Keeps track of the number of semi-critical patients
	OffPlaceCount = 0 		# Keeps track of offplaced patients
	OffPlaceCountC = 0 		# Keeps track of offplaced critical patients 
	OffPlaceCountS = 0 		# Keeps track of offplaced semi-critical patients
	ReturnCritCnt = 0 		# Keeps track of returns to critical
	ReturnSCritCnt = 0 		# Keeps track of returns to semi-critical
	OffPlaceCTime = 0 		# Keeps track of offplace time in critical
	OffPlaceSCTime = 0      # Keeps track of offplace time in semi-critical
	weightedLine = 0 		# Keeps track of time weighted length of line to be divided by time
	weightedTotalLine = 0 	# Keeps track of weighted average of line
	line_t = 0 				# Keeps track of the last time since line change


	served = 0
	case = ''
	servedS = 0

	# New Event variables
	newArrivalC = float('inf')
	newArrivalS = float('inf')

	Iempty = 1 				#Keeps track of if ICUbed is empty
	Icount = 0 				#Keeps track of number of beds occupied
	Sempty = 1 				#Keeps track of if ICUbed is empty
	Scount = 0 				#Keeps track of number of beds occupied
	I_index = range(0,B_I)	#Keeps track of bed index that are free for ICU
	S_index = range(0,B_S)	#Keeps track of bed index that are free for SDU


	time_case = {}
	loop_count = 0
	# ------------------------------ Loop ------------------------------
	while(t < TotalT): #Loop runs until total run time is reached, USE Tmax to end with empty hospital otherwise USE TotalT
		loop_count += 1
		run_time = time.clock()

	# ------------------------------ Events List ------------------------------
		# Generates new variables and updates event list
		if newArrivalC == float('inf') and t < TotalT:
			newArrivalC = np.random.exponential(arr_rate_C) + t
		if newArrivalS == float('inf') and t < TotalT and arr_rate_S != 0:
			newArrivalS = np.random.exponential(arr_rate_S) + t

		#Finds next abandonment in ICU queue 
		if ICUqueueAban:
			min_ICUqAban = ICUqueueAban[0][0]
		else:
			min_ICUqAban = float('inf')

		#Finds next abandonment in returning ICU queue 
		if ICUreturnqAban:
			min_ICUretqAban = ICUreturnqAban[0][0]
		else:
			min_ICUretqAban = float('inf')	

		#Finds next semi-critical patient who finishes healing in outside
		if SDUqueue:
			min_SDUq = SDUqueue[0][0]
		else:
			min_SDUq = float('inf')

		#Finds next return patient to re-enter system
		if RetDelayQ:
			min_RetQ = RetDelayQ[0][0]
		else:
			min_RetQ = float('inf')

		#Find when next ICU bed is left by critical patient
		if ICUbeds_C:
			min_ICUbed_C = ICUbeds_C[0][0]
		else:
			min_ICUbed_C = float('inf')

		#Find when next ICU bed is left by semi-critical patient
		if ICUbeds_S:
			min_ICUbed_S = ICUbeds_S[0][0]
		else:
			min_ICUbed_S = float('inf')


		#Find when next SDU bed is left by semi-critical patient
		if SDUbeds_S:
			min_SDUbed_S = SDUbeds_S[0][0]
		else:
			min_SDUbed_S = float('inf')

		#Find when next SDU bed is left by critical patient
		if SDUbeds_C:
			min_SDUbed_C = SDUbeds_C[0][0]
		else:
			min_SDUbed_C = float('inf')


	#              _________ Next Event Stage _________
		
		# Determine which event must occur next
		newEvent = min(newArrivalC, newArrivalS, min_ICUqAban, min_ICUretqAban, min_SDUq, min_RetQ, min_ICUbed_C, min_ICUbed_S, min_SDUbed_S, min_SDUbed_C)

		# Get ICU queue length
		ICUlinelength = len(ICUqueue) + len(ICUreturnq) + len(SDUbeds_C)
		prevICUlinelength = len(ICUqueue) + len(ICUreturnq) + len(SDUbeds_C)


	#	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Critical patient arrival event ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		if newEvent == newArrivalC:		#Arrival is a critical patient
			case = 'Critical Arrival'

			ArrivalCount += 1
			CritArrival += 1

			t = newEvent

	#	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ICU Bed is open ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			if ICUcount < Ibeds: #ICU Bed is open

				ICUtuple = (np.random.lognormal(crit_mean, crit_std) + t, 'C', t) #(Time service is completed, patient type, arrival time) tuple
				ICUcount += 1
				heapq.heappush(ICUbeds_C, ICUtuple)


	#	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ICU Bed has semi-critical patients to bump ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			elif ICUbeds_S and SDUcount >= Sbeds: # Check for bump opportunity

				#Bump semi-critical patient that entered the ICU earliest
				min_arr = float('inf')
				for s_pat in ICUbeds_S:
					if s_pat[1] < min_arr:
						min_arr = s_pat[1]
						tup = s_pat
				
				#Bumping should always send semi-critical patient to outside since if the SDU had an open spot it would have been filled by the semi-crit patient
				ICUbeds_S.remove(tup)

				#Adds semi-critical patient to outside
				heapq.heappush(SDUqueue, (tup[0],t,1))	#(Time service is completed, arrival time, bumped status)

				#Calculates offplace time of the bumped patient
				OffPlaceSCTime += t - tup[1]

				#Places ICU patient in ICU bed
				ICUtuple = (np.random.lognormal(crit_mean, crit_std) + t, 'C', t) #(Time service is completed, patient type, arrival time) tuple
				heapq.heappush(ICUbeds_C, ICUtuple) #put critical patient into critical section of ICU

	#	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Only SDU Bed is open ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			elif ICUcount >= Ibeds and SDUcount < Sbeds and ICUlinelength < K: 
				
				serv_time = np.random.lognormal(crit_mean, crit_std)*multiplier + t
				aban_time = np.random.exponential(patience) + t

				#Determines whether patient will abandon first or finish treatment first in the SDU
				if serv_time > aban_time:
					ICUtuple = (aban_time, serv_time, 'A','C', t) #(Time critical patient abandons SDU, Time service is completed, Aban/Service Dominated, patient type, arrival time) tuple
				else:
					ICUtuple = (serv_time, aban_time, 'S','C', t) #(Time service is completed, Time critical patient abandons SDU, Aban/Service Dominated, patient type, arrival time) tuple

				#ICU patient is added to SDU
				heapq.heappush(SDUbeds_C, ICUtuple)
				SDUcount += 1

				#Count offplacement
				OffPlaceCount += 1
				OffPlaceCountC += 1

	#	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ All Beds Full and ICU Queue is UNDER capacity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			elif ICUcount >= Ibeds and SDUcount >= Sbeds and ICUlinelength < K:  
				# print 'ICU queue'
				queue_tuple = (np.random.exponential(aband_rate) + t, t)
				heapq.heappush(ICUqueueAban, queue_tuple)
				ICUqueue.append(queue_tuple)
				ICUqueue_count += 1

	#	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ All Beds Full and ICU Queue is OVER capacity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
			elif ICUcount == Ibeds and SDUcount == Sbeds and ICUlinelength >= K: 
				BalkCount += 1
				# print 'Balk'
	
	#	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ All Beds Full and ICU Queue is OVER capacity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
			# This is techinically a repeat of the previous case
			elif ICUcount == Ibeds and ICUlinelength >= K: 
				BalkCount += 1
				# print 'Balk'

			else: 
				print 'ICU Bed Assignment Error'
				

			newArrivalC = float('inf')
			
	#	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Semi-critical patient arrival event +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		elif newEvent == newArrivalS:		#Arrival is a semi-critical patient
			
			case = 'Semi-critical arrival'

			ArrivalCount += 1
			SemiCritArrival += 1
			t = newEvent

	#	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SDU Bed is open ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			if SDUcount < Sbeds: 		
				#Add semi critical patient to Semi-critical part of SDU
				SDUtuple = (np.random.lognormal(semicrit_mean, semicrit_std) + t, t, 0)
				heapq.heappush(SDUbeds_S, SDUtuple)
				SDUcount += 1

	#	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Empty ICU queue and Only ICU bed open ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			elif SDUcount >= Sbeds and ICUcount < Ibeds and ICUlinelength == 0: 
				#Add semi-critical patient to Semi-critical part of SDU
				SDUtuple = (np.random.lognormal(semicrit_mean, semicrit_std) + t, t, 0)
				heapq.heappush(ICUbeds_S, SDUtuple)
				ICUcount += 1

				OffPlaceCount += 1
				OffPlaceCountS += 1

	#	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Only ICU bed is open and ICU queue is NOT empty ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			elif SDUcount >= Sbeds and ICUcount < Ibeds and ICUlinelength != 0: 
				#Send Semi-Critical patient to Outside
				SDUtuple = (np.random.lognormal(semicrit_mean, semicrit_std) + t, t, 0)
				heapq.heappush(SDUqueue, SDUtuple)


				print 'THIS SHOULD NOT HAPPEN'

	#	~~~~~~~~~~~~~~~ All beds full ~~~~~~~~~~~~~~~
			elif SDUcount >= Sbeds and ICUcount >= Ibeds: 
				#Send Semi-Critical patient to Outside
				SDUtuple = (np.random.lognormal(semicrit_mean, semicrit_std) + t, t, 0)
				heapq.heappush(SDUqueue, SDUtuple)

			else:
				print 'SDU Bed Assignment Error'
				print Ibeds
				print ICUcount
				print ICUlinelength
				print SDUcount
				print Sbeds

			newArrivalS = float('inf')		

	#	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Abandonment in ICU queue event +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		elif newEvent == min_ICUqAban:		#Abandonment occurs in ICUqueue
			
			case = 'Abandonment in ICU queue'

			abandoner = heapq.heappop(ICUqueueAban)
			ICUqueue.remove(abandoner)
			ICUqueue_count -= 1

			t = newEvent

			waittime += (t - abandoner[1])

			AbandCount += 1


	#	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Abandonment in returning to ICU queue event ++++++++++++++++++++++++++++++++++++++++++++
	#	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		elif newEvent == min_ICUretqAban:

			case = 'Abandonment in ICU return queue'

			abandoner = heapq.heappop(ICUreturnqAban)
			ICUreturnq.remove(abandoner)
			ICUreturnq_count -= 1

			t = newEvent

			waittime += (t - abandoner[1])

			AbandCount += 1

	#	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Patient in SDU 'outside' heals event +++++++++++++++++++++++++++++++++++++++++++++++++++
	#	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		elif newEvent == min_SDUq:

			case = 'Patient in SDU outside heals'
			t = newEvent
			SC_patient = heapq.heappop(SDUqueue) #(Service time, arrival, bump)
			SDUq_arr = SC_patient[1]
			SDUq_dep = SC_patient[0]
			bump = SC_patient[2]
			waittime += (t - SDUq_arr)

			p = np.random.random()

			#Determine if patient retuns after healing and places in delay queue
			if p < p_sc_sc:
				heapq.heappush(RetDelayQ, (t + np.random.exponential(ret_delay), 'S'))
				ReturnSCritCnt += 1
				served += 1
				servedS += 1
			elif p_sc_sc <= p < (p_sc_sc + p_sc_c):
				heapq.heappush(RetDelayQ, (t + np.random.exponential(ret_delay), 'C'))
				ReturnCritCnt += 1
				served += 1
				servedS += 1
			else: 
				# print 'Leaves System'
				if bump == 1:
					BumpCount += 1
				served += 1
				servedS += 1


	#	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Patient returns to system event +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		elif newEvent == min_RetQ:
			
			case = 'Patient returns to system'

			t = newEvent

			ArrivalCount += 1

			ret_patient = heapq.heappop(RetDelayQ)

			# If patient returns as a critical patient then they go through the same process 
			# as a arriving critical patient except there they cannot balk and are put into 
			# the returning to ICU queue
			if ret_patient[1] == 'C':
				CritArrival += 1
				if ICUcount < Ibeds: #ICU Bed is open
					ICUtuple = (np.random.lognormal(crit_mean, crit_std) + t, 'C', t) #(Time service is completed, patient type, arrival time) tuple
					ICUcount += 1
					heapq.heappush(ICUbeds_C, ICUtuple)

				elif ICUbeds_S: # Check for bump opportunity

					#Bump semi-critical patient that entered the ICU earliest
					min_arr = float('inf')
					indx = 0
					for ind, s_pat in enumerate(ICUbeds_S):
						if s_pat[1] < min_arr:
							min_arr = s_pat[1]
							tup = s_pat 
					
					#Bumping should always send semi-critical patient to outside since if the SDU had an open spot it would have been filled by the semi-crit patient
					ICUbeds_S.remove(tup)

					#Adds semi-critical patient to outside
					heapq.heappush(SDUqueue, (tup[0],t,1))	#(Time service is completed, arrival time, bumped status)

					#Calculates offplace time of the bumped patient
					OffPlaceSCTime += t - tup[1]

					#Places ICU patient in ICU bed
					ICUtuple = (np.random.lognormal(crit_mean, crit_std) + t, 'C', t) #(Time service is completed, patient type, arrival time) tuple
					heapq.heappush(ICUbeds_C, ICUtuple) #put critical patient into critical section of ICU


				elif ICUcount >= Ibeds and SDUcount < Sbeds: #Only SDU Bed is open
					serv_time = np.random.lognormal(crit_mean, crit_std)*multiplier + t
					aban_time = np.random.exponential(patience) + t

					#Determines whether patient will abandon first or finish treatment first in the SDU
					if serv_time > aban_time:
						#(Time critical patient abandons SDU, Time service is completed, Aban/Service Dominated, patient type, arrival time) tuple
						ICUtuple = (aban_time, serv_time, 'A','C', t) 
					else:
						#(Time service is completed, Time critical patient abandons SDU, Aban/Service Dominated, patient type, arrival time) tuple
						ICUtuple = (serv_time, aban_time, 'S','C', t) 

					#ICU patient is added to SDU
					heapq.heappush(SDUbeds_C, ICUtuple)
					SDUcount += 1

					#Count offplacement
					OffPlaceCount += 1
					OffPlaceCountC += 1


				elif ICUcount >= Ibeds and SDUcount >= Sbeds: # All Beds Full and ICU Queue is UNDER capacity
					# print 'ICU queue'
					abnd_t = np.random.exponential(aband_rate) + t
					ICUreturnq.append((abnd_t, t))
					heapq.heappush(ICUreturnqAban, (abnd_t, t))
					ICUreturnq_count += 1
				
				else:
					print 'Return ICU Bed Assignment Error'

			# If patient returns as a semi-critical patient then they go through the same process 
			# as a arriving semi-critical patient 
			elif ret_patient[1] == 'S':
				SemiCritArrival += 1
				if SDUcount < Sbeds: 		#SDU Bed is open
					SDUtuple = (np.random.lognormal(semicrit_mean, semicrit_std) + t, t, 0)
					heapq.heappush(SDUbeds_S, SDUtuple)
					SDUcount += 1

				elif SDUcount >= Sbeds and ICUcount < Ibeds and ICUlinelength == 0: #Only ICU bed is open and ICU queue is empty
					SDUtuple = (np.random.lognormal(semicrit_mean, semicrit_std) + t, t, 0)
					heapq.heappush(ICUbeds_S, SDUtuple)
					ICUcount += 1

					OffPlaceCount += 1
					OffPlaceCountS += 1

				elif SDUcount >= Sbeds and ICUcount < Ibeds and ICUlinelength != 0: #Only ICU bed is open and ICU queue is not empty
					SDUtuple = (np.random.lognormal(semicrit_mean, semicrit_std) + t, t, 0)
					heapq.heappush(SDUqueue, SDUtuple)

				elif SDUcount >= Sbeds and ICUcount >= Ibeds: #All Beds Full
					SDUtuple = (np.random.lognormal(semicrit_mean, semicrit_std) + t, t, 0)
					heapq.heappush(SDUqueue, SDUtuple)

				else:
					print 'SDU Bed Assignment Error 1'


			else:
				print 'Return Queue Error'


	#	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#	++++++++++++++++++++++++++++++++++++++++++++++++++ ICU bed is cleared event by critical patient ++++++++++++++++++++++++++++++++++++++++++++++++++
	#	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		elif newEvent == min_ICUbed_C or newEvent == min_ICUbed_S:
			
			case = 'ICU bed cleared'
			
			t = newEvent

			#	Determines where the critical patient goes after leaving ICU bed
			if newEvent == min_ICUbed_C:
				
				crit_tuple = heapq.heappop(ICUbeds_C)
				ICUcount -= 1
				p = np.random.random()
				if p < p_c_sc:
	
				#	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SDU Bed is open ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					if SDUcount < Sbeds:

						heapq.heappush(SDUbeds_S, (np.random.lognormal(semicrit_mean, semicrit_std) + t, t, 0))
						SDUcount += 1
						
						# print 'test'
				#	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Empty ICU queue and Only ICU bed open ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					elif SDUcount >= Sbeds and ICUcount < Ibeds and ICUlinelength == 0: 
						
						heapq.heappush(ICUbeds_S, (np.random.lognormal(semicrit_mean, semicrit_std) + t, t, 0))
						ICUcount += 1

						OffPlaceCount += 1
						OffPlaceCountS += 1
						
						# print 'test'
				#	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Only ICU bed is open and ICU queue is NOT empty ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					elif SDUcount >= Sbeds and ICUcount < Ibeds and ICUlinelength != 0: 
						
						heapq.heappush(SDUqueue, (np.random.lognormal(semicrit_mean, semicrit_std) + t, t, 1) )
						
				#	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ All beds full ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					elif SDUcount >= Sbeds and ICUcount >= Ibeds: 
						
						heapq.heappush(SDUqueue, (np.random.lognormal(semicrit_mean, semicrit_std) + t, t, 1) )
						

					else:

						print 'SDU Bed Assignment Error 2'

				
				
				else:
					p10 = np.random.random()
					if p10 < p_c_sc_d:
						heapq.heappush(RetDelayQ, (t + np.random.exponential(ret_delay), 'S'))
						served += 1
						ReturnSCritCnt += 1
					elif p_c_sc_d <= p10 < p_c_c_d + p_c_sc_d:
						heapq.heappush(RetDelayQ, (t + np.random.exponential(ret_delay), 'C'))
						ReturnCritCnt += 1
						served += 1
					else:
						# print 'Leaves System'
						served += 1

			#	Case where semicritical patient leaves
			elif newEvent == min_ICUbed_S:

				semicrit_tuple = heapq.heappop(ICUbeds_S)
				ICUcount -= 1
				OffPlaceSCTime += (t - semicrit_tuple[1])

				p10 = np.random.random()
				if p10 < p_sc_sc:
					heapq.heappush(RetDelayQ, (t + np.random.exponential(ret_delay), 'S'))
					served += 1
					ReturnSCritCnt += 1
				elif p_sc_sc <= p10 < p_sc_c + p_sc_sc:
					heapq.heappush(RetDelayQ, (t + np.random.exponential(ret_delay), 'C'))
					ReturnCritCnt += 1
					served += 1
				else:
					# print 'Leaves System'
					pass


	#	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Checks if ICU patients are waiting for ICU ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			if ICUlinelength > 0: 

				# Checks if there are offplaced patients and selects patient with the shortest service time left
				if SDUbeds_C: 
					min_arr = float('inf') 	
					SDU_ind = 0
					for pt in SDUbeds_C:
						if pt[2] == 'S':
							if pt[0] < min_arr:
								min_arr = pt[0]
								Ituple = pt
						elif pt[2] == 'A':
							if pt[1] < min_arr:
								min_arr = pt[1]
								Ituple = pt

					SDUbeds_C.remove(Ituple)
					nICUtuple = ((min_arr-t)/multiplier + t, 'C', t)
					heapq.heappush(ICUbeds_C, nICUtuple)
					ICUcount += 1

					OffPlaceCTime += (t - Ituple[4])

					# Reassign the originally offplaced bed
					# Checks return ICU queue to fill offplaced bed in SDU
					if ICUreturnq:
						
						ICU_tuple = ICUreturnq.pop(0)
						ICUreturnqAban.remove(ICU_tuple)

						waittime += (t - ICU_tuple[1])

						serv_time = np.random.lognormal(crit_mean, crit_std)*multiplier + t
						aban_time = np.random.exponential(patience) + t

						#Determines whether patient will abandon first or finish treatment first in the SDU
						if serv_time > aban_time:
							#(Time critical patient abandons SDU, Time service is completed, Aban/Service Dominated, patient type, arrival time) tuple
							ICUtuple = (aban_time, serv_time, 'A','C', t) 
						else:
							#(Time service is completed, Time critical patient abandons SDU, Aban/Service Dominated, patient type, arrival time) tuple
							ICUtuple = (serv_time, aban_time, 'S','C', t) 

						heapq.heappush(SDUbeds_C, ICUtuple)

						OffPlaceCount += 1
						OffPlaceCountC += 1
						# print '3'

					# Checks ICU queue to fill offplaced bed
					elif not ICUreturnq and ICUqueue:	

						ICU_tuple = ICUqueue.pop(0)
						ICUqueueAban.remove(ICU_tuple)

						waittime += (t - ICU_tuple[1])

						serv_time = np.random.lognormal(crit_mean, crit_std)*multiplier + t
						aban_time = np.random.exponential(patience) + t

						#Determines whether patient will abandon first or finish treatment first in the SDU
						if serv_time > aban_time:
							#(Time critical patient abandons SDU, Time service is completed, Aban/Service Dominated, patient type, arrival time) tuple
							ICUtuple = (aban_time, serv_time, 'A','C', t) 
						else:
							#(Time service is completed, Time critical patient abandons SDU, Aban/Service Dominated, patient type, arrival time) tuple
							ICUtuple = (serv_time, aban_time, 'S','C', t) 

						heapq.heappush(SDUbeds_C, ICUtuple)

						OffPlaceCount += 1
						OffPlaceCountC += 1
						# print '4'

					# ChecksSDU queue to fill offplaced bed
					elif not ICUreturnq and not ICUqueue and SDUqueue:
						
						ICU_tuple = heapq.nlargest(1, SDUqueue)[0]
						SDUqueue.remove(ICU_tuple)
						SDU_tuple = (ICU_tuple[0], t, 0)

						waittime += (t - ICU_tuple[1])

						heapq.heappush(SDUbeds_S, SDU_tuple)
					
					# If there is no demand, the bed is empty
					else:
						SDUcount -= 1


	#	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Checks ICU return queue to fill opening ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				elif not SDUbeds_C and ICUreturnq: # If there are no offplaced patients, check if there is a return queue
					
					Qtuple = ICUreturnq.pop(0)
					ICUreturnqAban.remove(Qtuple)

					waittime += (t - Qtuple[1])

					heapq.heappush(ICUbeds_C, (np.random.lognormal(semicrit_mean, semicrit_std) + t, 'C', t))
					ICUcount += 1


	#	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Checks ICU queue to fill opening ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				elif not SDUbeds_C and not ICUreturnq and ICUqueue: # If there is no return queue, then check if there is a regular queue
					
					Qtuple = ICUqueue.pop(0)
					ICUqueueAban.remove(Qtuple)

					waittime += (t - Qtuple[1])

					heapq.heappush(ICUbeds_C, (np.random.lognormal(semicrit_mean, semicrit_std) + t, 'C', t))
					ICUcount += 1



	#	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Checks for SDU patients to put in empty ICU bed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			elif ICUlinelength == 0 and SDUqueue and SDUcount >= Sbeds and ICUcount < Ibeds: # Checks if SDU patients are waiting for ICU if there are no ICU patients waiting
				
				ICU_tuple = heapq.nlargest(1, SDUqueue)[0]
				SDUqueue.remove(ICU_tuple)
				
				heapq.heappush(ICUbeds_S, (ICU_tuple[0], t, 0))
				ICUcount += 1

				OffPlaceCount += 1
				OffPlaceCountS += 1

				waittime += (t - ICU_tuple[1])


			else:
				pass




	#	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++ SDU bed is cleared event +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		elif newEvent == min_SDUbed_C or newEvent == min_SDUbed_S:

			case = 'SDU bed clears'

			t = newEvent
			

			if newEvent == min_SDUbed_C:
				S_tuple = heapq.heappop(SDUbeds_C)
				SDUcount -= 1
				OffPlaceCTime += (t - S_tuple[4])

				#Checks if bed was abandoned or not
				if S_tuple[2] == 'A':	
					AbandCount += 1

			elif newEvent == min_SDUbed_S:
				heapq.heappop(SDUbeds_S)
				SDUcount -= 1



			if ICUbeds_S and not SDUqueue and ICUlinelength == 0: # Check if there are off-placed semi-critical patients
				
				# Check for patient with the least amount of service time remaining
				min_arr = float('inf')
				for tup in ICUbeds_S:
					if tup[1] < min_arr:
						min_arr = tup[1]
						I_tuple = tup

				ICUbeds_S.remove(I_tuple)
				ICUcount -= 1

				SDU_tuple = (I_tuple[0] ,t, 0)

				heapq.heappush(SDUbeds_S, SDU_tuple)
				SDUcount += 1

				OffPlaceSCTime += (t - I_tuple[1])


			elif ICUreturnq: #Check if returning ICU queue has demand

				Q_tuple = ICUreturnq.pop(0)
				ICUreturnqAban.remove(Q_tuple)

				waittime += (t - Q_tuple[1])

				serv_time = np.random.lognormal(crit_mean, crit_std)*multiplier + t
				aban_time = np.random.exponential(patience) + t

				#Determines whether patient will abandon first or finish treatment first in the SDU
				if serv_time > aban_time:
					#(Time critical patient abandons SDU, Time service is completed, Aban/Service Dominated, patient type, arrival time) tuple
					ICUtuple = (aban_time, serv_time, 'A','C', t) 
				else:
					#(Time service is completed, Time critical patient abandons SDU, Aban/Service Dominated, patient type, arrival time) tuple
					ICUtuple = (serv_time, aban_time, 'S','C', t) 

				#ICU patient is added to SDU
				heapq.heappush(SDUbeds_C, ICUtuple)
				SDUcount += 1

				OffPlaceCount += 1
				OffPlaceCountC += 1
				# print '5'
	 
			elif ICUqueue: #Check if ICU queue has demand

				Q_tuple = ICUqueue.pop(0)
				ICUqueueAban.remove(Q_tuple)

				waittime += (t - Q_tuple[1])

				serv_time = np.random.lognormal(crit_mean, crit_std)*multiplier + t
				aban_time = np.random.exponential(patience) + t

				#Determines whether patient will abandon first or finish treatment first in the SDU
				if serv_time > aban_time:
					#(Time critical patient abandons SDU, Time service is completed, Aban/Service Dominated, patient type, arrival time) tuple
					ICUtuple = (aban_time, serv_time, 'A','C', t) 
				else:
					#(Time service is completed, Time critical patient abandons SDU, Aban/Service Dominated, patient type, arrival time) tuple
					ICUtuple = (serv_time, aban_time, 'S','C', t) 

				#ICU patient is added to SDU
				heapq.heappush(SDUbeds_C, ICUtuple)
				SDUcount += 1

				OffPlaceCount += 1
				OffPlaceCountC += 1
				# print 6

			elif SDUqueue: #Check if SDUqueue is empty
			
				Q_tuple = heapq.nlargest(1, SDUqueue)[0]
				SDUqueue.remove(Q_tuple)
				
				heapq.heappush(SDUbeds_S, (Q_tuple[0], t, 0))
				SDUcount += 1

				waittime += (t - Q_tuple[1])

				# print 7
			else:
				pass

			# else:
			# 	print 'No Demand for SDU'
			if newEvent == min_SDUbed_S or (newEvent == min_SDUbed_C and S_tuple[2] == 'S'):
				
				p = np.random.random()

				if p < p_sc_sc:
					heapq.heappush(RetDelayQ, (t + np.random.exponential(ret_delay), 'S'))
					ReturnSCritCnt += 1
					served += 1
				elif p_sc_sc <= p < (p_sc_sc + p_sc_c):
					heapq.heappush(RetDelayQ, (t + np.random.exponential(ret_delay), 'C'))
					ReturnCritCnt += 1
					served += 1
				else: 
					# print 'Leaves System'
					served += 1




	#	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	# 	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++ CLEAN UP +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		#Reset counters unless time is above warm-up time, comments are repeated from when variables were initialized 
		if t < t_warmup:
			ArrivalCount = 0	# Keeps track of the number of arrivals
			BalkCount = 0		# Keeps track of the number of balks
			BumpCount = 0 		# Keeps track of the number of bumps
			AbandCount = 0 		# Keeps track of the number of abandonments
			waittime = 0 		# Keeps track of total wait time
			CritArrival = 0 	# Keeps track of the number of critical patients
			SemiCritArrival = 0 # Keeps track of the number of semi-critical patients
			OffPlaceCount = 0 	# Keeps track of offplaced patients
			OffPlaceCountC = 0 	# Keeps track of offplaced critical patients 
			OffPlaceCountS = 0 	# Keeps track of offplaced semi-critical patients
			ReturnCritCnt = 0 	# Keeps track of returns to critical
			ReturnSCritCnt = 0 	# Keeps track of returns to semi-critical
			OffPlaceCTime = 0 	# Keeps track of offplace time in critical
			OffPlaceSCTime = 0  # Keeps track of offplace time in semi-critical
			weightedTotalLine = 0 # Keeps track of line length states multiplied by the duration of their states
			last_line_time = Tw 
			line_t = Tw  
			servedS = 0


		# Checks for errors, guarantees variables keeping track of beds are correct
		if (len(ICUbeds_S) + len(ICUbeds_C)) != ICUcount:
			print 'error'
		if (len(SDUbeds_S) + len(SDUbeds_C)) != SDUcount:
			print 'error'

		# Update Line Length and the variable 
		ICUlinelength = len(ICUqueue) + len(ICUreturnq) + len(SDUbeds_C)
		weightedTotalLine += (t-line_t)* prevICUlinelength
		line_t = t

		time_run = time.clock() - run_time
		# ~~~~~ Code below was used to test how fast the simulation would run for each case ~~~~
		# if case not in time_case:
		# 	time_case[case] = [time_run]
		# else:
		# 	time_case[case].append(time_run)

		# ---------------------------------------- LOOP END ----------------------------------------
	
	# Statistics are calculated based of simulation time minus warm-up time
	# If-else statements generally indicate cases where the parameters in the denominator are zero
	AbandonmentRate = AbandCount/float(T-t_warmup)
	BalkRate = BalkCount/float(T-t_warmup)
	BumpRate = BumpCount/float(T-t_warmup)
	EWaittime = waittime/float(ArrivalCount)
	P_return_Critical = ReturnCritCnt/float(CritArrival)
	if SemiCritArrival != 0:
		P_return_SemiCrit = ReturnSCritCnt/float(SemiCritArrival)
	else:
		P_return_SemiCrit = 0
	P_offplacementC = float(OffPlaceCountC)/ArrivalCount
	P_offplacementSC = float(OffPlaceCountS)/ArrivalCount
	if OffPlaceCountC > 0:
		E_time_offplaced_C = float(OffPlaceCTime)/OffPlaceCountC
	else:
		E_time_offplaced_C = 0
	if OffPlaceCountS > 0:
		E_time_offplaced_S = float(OffPlaceSCTime)/OffPlaceCountS
	else:
		E_time_offplaced_S = 0


#  Returns statistics as an array
	return AbandonmentRate, BalkRate, BumpRate, EWaittime, P_return_Critical, P_return_SemiCrit, P_offplacementC, P_offplacementSC, E_time_offplaced_C, E_time_offplaced_S, weightedTotalLine/float(t - Tw), servedS, ArrivalCount, CritArrival, SemiCritArrival

# Function to write table results to a file, used in the main function to write the return values to a csv file
def writeLog(fil, table):
	c1 = csv.writer(fil)

	for val in table:
		# print val
		c1.writerow(val)

def main():
	table = [] # Table containing rows that will be written to csv file.
	header = ['ICU beds', 'SDU beds' , 'ArrivalCount', 'Served Outside', 'Arrival Count', 'Critical Arrival', 'Semi-Critical Arrival', 'Balk Threshold', 'Abandonment Rate', 'Balk Rate', 'Bump Rate', 'Expected Wait Time', 'P(Return to Critical)', 'P(Return to Semi-Critical)', 'P(Offplacement of Critical patient)', 'P(Offplacement of Semi-critical patient)', 'E_time_offplaced_C', 'E_time_offplaced_S', 'Average Queue Length']
	table.append(header) 
	krange = range(0,sys.argv[17])
	krange.append(float('inf'))
	x = int(sys.argv[16])
	for k in krange:
		start_time = time.clock()
		row = []
		AR, BaR, BuR, EWT, PrC, PrS, PoC, PoS, EtoC, EtoS, aql, SS, ArrC, CArr, SCArr = simulation(sys.argv[1], sys.argv[2],  x, (20-x)*3 + 11, k, sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[8]*1.5, sys.argv[9], sys.argv[9]*sys.argv[10], sys.argv[11], sys.argv[12], sys.argv[13], sys.argv[14], sys.argv[15])
		row.append(x)
		row.append((20-x)*3 + 11)
		row.append(ArrC)
		row.append(SS)
		row.append(ArrC)
		row.append(CArr)
		row.append(SCArr)
		row.append(k)
		row.append(AR)
		row.append(BaR)
		row.append(BuR)
		row.append(EWT)
		row.append(PrC)
		row.append(PrS)
		row.append(PoC)
		row.append(PoS)
		row.append(EtoC)
		row.append(EtoS)
		row.append(aql)
		table.append(row)
		print('--- %s seconds ---' % (time.clock() - start_time))
	writeLog(fil, table)

# main()


