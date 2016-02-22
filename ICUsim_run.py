import os

directory = os.getcwd()
for x in range(14,21):
	os.system(directory + "/ICUsimfunc_Clean_B_vary.py " + str(x))