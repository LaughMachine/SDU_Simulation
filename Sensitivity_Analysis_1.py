import ICUsimfunc_Clean_B_E as ISF1 
import csv
import time

mode = int(sys.argv[2])

fil = open(os.getcwd() + "/SimResults_SA_1_" + str(sys.argv[2]) + ".csv","wb")

prange = [.5, .6, .7, .8, .9, 1, 1.1, 1.2, 1.3, 1.4, 1.5]
prange1 = [.5, .6, .7, .8, .9, 1]
prange2 = [1.1, 1.2, 1.3, 1.4, 1.5]


musc = [1, 1.04, 1.08, 1.12, 1.16, 1.20, 1.24, 1.28, 1.32, 1.36, 1.4]
musc1 = [1, 1.04, 1.08, 1.12, 1.16, 1.20]
musc2 = [1.24, 1.28, 1.32, 1.36, 1.4]

muc = [2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7]
muc1 = [2.2, 2.25, 2.3, 2.35, 2.4, 2.45]
muc2 = [2.5, 2.55, 2.6, 2.65, 2.7]

pr = [.59, .6, .61, .62, .63, .64, .65, .66, .67, .68, .69, .7, .71]
pr1 = [.59, .6, .61, .62, .63, .64]
pr2 = [.65, .66, .67, .68, .69, .7, .71]

tr = [.5, .6, .7, .8, .9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2]
tr1 = [.5, .6, .7, .8, .9, 1, 1.1, 1.2]
tr2 = [1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2]

if mode == 1:
	rangeset = prange1
elif mode == 2:
	rangeset = prange2
elif mode == 3:
	rangeset = prange1
elif mode == 4:
	rangeset = prange2
elif mode == 5:
	rangeset = prange1
elif mode == 6:
	rangeset = prange2
elif mode == 7: 
	rangeset = prange1
elif mode == 8:
	rangeset = prange2
elif mode == 9:
	rangeset = musc1
elif mode == 10:
	rangeset = musc2
elif mode == 11:
	rangeset = muc1
elif mode == 12:
	rangeset = muc2
elif mode == 13: 
	rangeset = pr1
elif mode == 14:
	rangeset = pr2
elif mode == 15:
	rangeset = tr1
elif mode == 16:
	rangeset = tr2

k = int(sys.argv[1])
p = int(sys.argv[2])
ICU = [15, 20, 17]
krange = [0, float('inf')]
k = krange[k]
table = [['ICU beds', 'SDU beds' ,'ArrivalCount', 'Served Outside', 'Arrival Count', 'Critical Arrival', 'Semi-Critical Arrival', 'Balk Threshold', 'Abandonment Rate', 'Balk Rate', 'Bump Rate', 'Expected Wait Time', 'P(Return to Critical)', 'P(Return to Semi-Critical)', 'P(Offplacement of Critical patient)', 'P(Offplacement of Semi-critical patient)', 'E_time_offplaced_C', 'E_time_offplaced_S', 'Average Queue Length', 'P-cost', 'Cost', 'Half', 'No SDU']]

halfcostB = 0
nosducostB = 0
halfcostQ = 0
nosducostQ = 0
for b in ICU:

	# start_time = time.time()
	AR, BaR, BuR, EWT, PrC, PrS, PoC, PoS, EtoC, EtoS, aql, SS, ArrC, CArr, SCArr = ISF1.simulation(20001000, 1000, b, (20-b)*3 + 11, k, .65, .07, .07, .07, .07, 2.5, 2.5*1.5, 1.2, 1.2*1.5, 1.5*p, 9, 9, 1, 1)
		
	row.append(b)
	row.append((20-b)*3 + 11)
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
	row.append(p)
	if k == 0:
		row.append((int(SS)/float(20001000) + aql*15 + BaR*2))
		if b == 17:
			row.append((int(SS)/float(20001000) + aql*15 + BaR*2 - halfcostB)/halfcostB)
			row.append((int(SS)/float(20001000) + aql*15 + BaR*2 - nosducostB)/nosducostB)
	elif k == float('inf'):
		row.append((int(SS)/float(20001000) + aql*2 + BaR*15))
		if b == 17:
			row.append((int(SS)/float(20001000) + aql*2 + BaR*15 - halfcostQ)/halfcostQ)
			row.append((int(SS)/float(20001000) + aql*2 + BaR*15 - nosducostQ)/nosducostQ)

	table.append(row)

	if k == 0 and b == 15:
		halfcostB = (int(SS)/float(20001000) + aql*15 + BaR*2
	elif k == float('inf') and b == 15:
		halfcostQ = (int(SS)/float(20001000) + aql*2 + BaR*15
	elif k == 0 and b == 20:
		nosducostB = (int(SS)/float(20001000) + aql*15 + BaR*2
	elif k == float('inf') and b == 20:
		nosducostQ = (int(SS)/float(20001000) + aql*2 + BaR*15

ISF1.writeLog(fil, table)



			# print 0
			# print b
			# print ("--- %s seconds ---" % (time.time() - start_time))
