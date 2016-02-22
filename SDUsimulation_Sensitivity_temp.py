import ICUsimfuncO1 as ISF1
import SDUequations as SDUQ 
import csv
import time

fil_mu_sc = open("/Users/michaelhuang/Projects/Python/SDU/SimResultsO1_SDU_Balk_mu_sc_Q.csv","wb")
fil_mu_c = open("/Users/michaelhuang/Projects/Python/SDU/SimResultsO1_SDU_Balk_mu_c_Q.csv","wb")
fil_p = open("/Users/michaelhuang/Projects/Python/SDU/SimResultsO1_SDU_Balk_p.csv_Q","wb")
fil_theta = open("/Users/michaelhuang/Projects/Python/SDU/SimResultsO1_SDU_Balk_theta_Q.csv","wb")


prange = [.5, .6, .7, .8, .9, 1, 1.1, 1.2, 1.3, 1.4, 1.5]
ICU = [15,20,17]
krange = [0, float('inf')]


musc = [1, 1.04, 1.08, 1.12, 1.16, 1.20, 1.24, 1.28, 1.32, 1.36, 1.4]
muc = [2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7]
pr = [.59, .6, .61, .62, .63, .64, .65, .66, .67, .68, .69, .7, .71]
tr = [.5, .6, .7, .8, .9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2]

mu_C = 1/2.5
mu_SC = 1/1.2
theta = 1
w_SC = 1
w_Q_C = 2
w_B_C = 15
r_S = 3
r_I = 1
p = .65
N = 20
w_C = min(w_Q_C, w_B_C)
lb = N*r_I*mu_C


table = [['ICU beds', 'SDU beds' , 'Balk Threshold', 'Abandonment Rate', 'Balk Rate', 'Bump Rate', 'Expected Wait Time', 'P(Return to Critical)', 'P(Return to Semi-Critical)', 'P(Offplacement of Critical patient)', 'P(Offplacement of Semi-critical patient)', 'E_time_offplaced_C', 'E_time_offplaced_S', 'Average Queue Length', 'P-cost', 'Balk Regime Cost', 'Queue Regime Cost']]
for x in musc:
	halfcostB = 0
	nosducostB = 0
	halfcostQ = 0
	nosducostQ = 0
	for b in ICU:
		for k in krange:

			start_time = time.time()

			if b == 17 and krange == 0:
				lb = N*r_I*mu_C
				b0 = np.round(SDUQ.B_4_C(lb, mu_C, x, w_SC, w_C, r_I, r_S, p, N))
			else: 
				b0 = b

			row = []
			AR, BaR, BuR, EWT, PrC, PrS, PoC, PoS, EtoC, EtoS, aql = ISF1.simulation(901000, 1000, b0, (20-b0)*3 + 11, float('inf'), .65, .07, .07, .07, .07, 2.5, 2.5*1.5, x, x*1.5, 1.5, 9, 9, 1, 1)
			row.append(b0)
			row.append((20-b0)*3 + 11)
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
			row.append(x)
			if k == 0:
				row.append(BuR + aql*15 + BaR*2)
				if b == 17 or b == 18:
					row.append((BuR + aql*15 + BaR*2 - halfcostB)/halfcostB)
					row.append((BuR + aql*15 + BaR*2 - nosducostB)/nosducostB)
			elif k == float('inf'):
				row.append(BuR + aql*2 + BaR*15)
				if b == 17 or b == 18:
					row.append((BuR + aql*2 + BaR*15 - halfcostQ)/halfcostQ)
					row.append((BuR + aql*2 + BaR*15 - nosducostQ)/nosducostQ)

			table.append(row)

			if k == 0 and b == 15:
				halfcostB = BuR + aql*15 + BaR*2
			elif k == float('inf') and b == 15:
				halfcostQ = BuR + aql*2 + BaR*15
			elif k == 0 and b == 20:
				nosducostB = BuR + aql*15 + BaR*2
			elif k == float('inf') and b == 20:
				nosducostQ = BuR + aql*2 + BaR*15



			print x
			print b
			print k
			print b0
			print ("--- %s seconds ---" % (time.time() - start_time))

ISF1.writeLog(fil_mu_sc, table)

table = [['ICU beds', 'SDU beds' , 'Balk Threshold', 'Abandonment Rate', 'Balk Rate', 'Bump Rate', 'Expected Wait Time', 'P(Return to Critical)', 'P(Return to Semi-Critical)', 'P(Offplacement of Critical patient)', 'P(Offplacement of Semi-critical patient)', 'E_time_offplaced_C', 'E_time_offplaced_S', 'Average Queue Length', 'P-cost', 'Balk Regime Cost', 'Queue Regime Cost']]
for x in muc:
	halfcostB = 0
	nosducostB = 0
	halfcostQ = 0
	nosducostQ = 0
	for b in ICU:
		for k in krange:

			if b == 17 and krange == 0:
				lb = N*r_I*mu_C
				b0 = np.round(SDUQ.B_4_C(lb, x, mu_SC, w_SC, w_C, r_I, r_S, p, N))
			else: 
				b0 = b

			start_time = time.time()
			row = []
			AR, BaR, BuR, EWT, PrC, PrS, PoC, PoS, EtoC, EtoS, aql = ISF1.simulation(901000, 1000, b0, (20-b0)*3 + 11, float('inf'), .65, .07, .07, .07, .07, x, x*1.5, 1.2, 1.2*1.5, 1.5, 9, 9, 1, 1)
			row.append(b0)
			row.append((20-b0)*3 + 11)
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
			row.append(x)
			if k == 0:
				row.append(BuR + aql*15 + BaR*2)
				if b == 17 or b == 18:
					row.append((BuR + aql*15 + BaR*2 - halfcostB)/halfcostB)
					row.append((BuR + aql*15 + BaR*2 - nosducostB)/nosducostB)
			elif k == float('inf'):
				row.append(BuR + aql*2 + BaR*15)
				if b == 17 or b == 18:
					row.append((BuR + aql*2 + BaR*15 - halfcostQ)/halfcostQ)
					row.append((BuR + aql*2 + BaR*15 - nosducostQ)/nosducostQ)

			table.append(row)

			if k == 0 and b == 15:
				halfcostB = BuR + aql*15 + BaR*2
			elif k == float('inf') and b == 15:
				halfcostQ = BuR + aql*2 + BaR*15
			elif k == 0 and b == 20:
				nosducostB = BuR + aql*15 + BaR*2
			elif k == float('inf') and b == 20:
				nosducostQ = BuR + aql*2 + BaR*15

			print x
			print k
			print b0
			print ("--- %s seconds ---" % (time.time() - start_time))

ISF1.writeLog(fil_mu_c, table)

table = [['ICU beds', 'SDU beds' , 'Balk Threshold', 'Abandonment Rate', 'Balk Rate', 'Bump Rate', 'Expected Wait Time', 'P(Return to Critical)', 'P(Return to Semi-Critical)', 'P(Offplacement of Critical patient)', 'P(Offplacement of Semi-critical patient)', 'E_time_offplaced_C', 'E_time_offplaced_S', 'Average Queue Length', 'P-cost', 'Balk Regime Cost', 'Queue Regime Cost']]
for x in pr:
	halfcostB = 0
	nosducostB = 0
	halfcostQ = 0
	nosducostQ = 0
	for b in ICU:
		for k in krange:

			if b == 17 and krange == 0:
				lb = N*r_I*mu_C
				b0 = np.round(SDUQ.B_4_C(lb, mu_C, mu_SC, w_SC, w_C, r_I, r_S, x, N))
			else: 
				b0 = b

			start_time = time.time()
			row = []
			AR, BaR, BuR, EWT, PrC, PrS, PoC, PoS, EtoC, EtoS, aql = ISF1.simulation(901000, 1000, b0, (20-b0)*3 + 11, float('inf'), x, .07, .07, .07, .07, 2.5, 2.5*1.5, 1.2, 1.2*1.5, 1.5, 9, 9, 1, 1)
			row.append(b0)
			row.append((20-b0)*3 + 11)
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
			row.append(x)
			if k == 0:
				row.append(BuR + aql*15 + BaR*2)
				if b == 17 or b == 18:
					row.append((BuR + aql*15 + BaR*2 - halfcostB)/halfcostB)
					row.append((BuR + aql*15 + BaR*2 - nosducostB)/nosducostB)
			elif k == float('inf'):
				row.append(BuR + aql*2 + BaR*15)
				if b == 17 or b == 18:
					row.append((BuR + aql*2 + BaR*15 - halfcostQ)/halfcostQ)
					row.append((BuR + aql*2 + BaR*15 - nosducostQ)/nosducostQ)

			table.append(row)

			if k == 0 and b == 15:
				halfcostB = BuR + aql*15 + BaR*2
			elif k == float('inf') and b == 15:
				halfcostQ = BuR + aql*2 + BaR*15
			elif k == 0 and b == 20:
				nosducostB = BuR + aql*15 + BaR*2
			elif k == float('inf') and b == 20:
				nosducostQ = BuR + aql*2 + BaR*15

			print x
			print k
			print b0
			print ("--- %s seconds ---" % (time.time() - start_time))

ISF1.writeLog(fil_p, table)

table = [['ICU beds', 'SDU beds' , 'Balk Threshold', 'Abandonment Rate', 'Balk Rate', 'Bump Rate', 'Expected Wait Time', 'P(Return to Critical)', 'P(Return to Semi-Critical)', 'P(Offplacement of Critical patient)', 'P(Offplacement of Semi-critical patient)', 'E_time_offplaced_C', 'E_time_offplaced_S', 'Average Queue Length', 'P-cost', 'Cost', 'half vs diff', 'nosdu vs diff']]
for x in tr:
	halfcostB = 0
	nosducostB = 0
	halfcostQ = 0
	nosducostQ = 0
	for b in ICU:
		for k in krange:

			if b == 17 and krange == 0:
				lb = N*r_I*mu_C
				b0 = np.round(SDUQ.B_4_C(lb, mu_C, mu_SC, w_SC, w_C, r_I, r_S, p, N))
			else: 
				b0 = b

			start_time = time.time()
			row = []
			AR, BaR, BuR, EWT, PrC, PrS, PoC, PoS, EtoC, EtoS, aql = ISF1.simulation(901000, 1000, b0, (20-b0)*3 + 11, float('inf'), .65, .07, .07, .07, .07, 2.5, 2.5*1.5, 1.2, 1.2*1.5, 1.5, 9, 9, x, x)
			row.append(b0)
			row.append((20-b0)*3 + 11)
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
			row.append(x)
			if k == 0:
				row.append(BuR + aql*15 + BaR*2)
				if b == 17 or b == 18:
					row.append((BuR + aql*15 + BaR*2 - halfcostB)/halfcostB)
					row.append((BuR + aql*15 + BaR*2 - nosducostB)/nosducostB)
			elif k == float('inf'):
				row.append(BuR + aql*2 + BaR*15)
				if b == 17 or b == 18:
					row.append((BuR + aql*2 + BaR*15 - halfcostQ)/halfcostQ)
					row.append((BuR + aql*2 + BaR*15 - nosducostQ)/nosducostQ)

			table.append(row)

			if k == 0 and b == 15:
				halfcostB = BuR + aql*15 + BaR*2
			elif k == float('inf') and b == 15:
				halfcostQ = BuR + aql*2 + BaR*15
			elif k == 0 and b == 20:
				nosducostB = BuR + aql*15 + BaR*2
			elif k == float('inf') and b == 20:
				nosducostQ = BuR + aql*2 + BaR*15

			print x
			print k
			print b0
			print ("--- %s seconds ---" % (time.time() - start_time))

ISF1.writeLog(fil_theta, table)
