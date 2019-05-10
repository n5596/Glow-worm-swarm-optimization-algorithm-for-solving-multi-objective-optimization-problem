import numpy as np
import time
import os
from subprocess import call

def fitness_obj(wl_array):
        sx = list()
        nx = np.zeros(len(wl_array),1)
        for i in range(len(wl_array)):
                for j in range(len(wl_array)):
                        if i !=j:
                                li = leakage_read(wl_array[i])
                                lj = leakage_read(wl_array[j])
                                di = delay_read(wl_array[i])
                                dj = delay_read(wl_array[j])
                                if li+di > lj+dj:
                                        l_list = sx[i]
                                        l_list.append(j)
                                        sx[i] = l_list
                                else:
                                        nx[i] = nx[i]+1
	Q = list()
	qrank = 0
	for F in sx:
		for q in F:
			nx[q] = nx[q]-1
			if nx[q] == 0:
				qrank = q+1
				Q.append(q)
	return 

def fitness(alpha, max_iter, i):
	v = alpha*(max_iter-i)/max_iter
	u = v * np.random.rand(1)
	return u

def process_change(i):
	filename = '22nm_MGK.pm'
	with open(filename, 'r') as file:
		process_data = file.readlines()
		toxe = 6.5e-010*(0.2*np.random.rand(1)+0.9)
		toxp = 4e-010*(0.2*np.random.rand(1)+0.9)
		toxref = 6.5e-010*(0.2*np.random.rand(1)+0.9)
		wint = 5e-009*(0.2*np.random.rand(1)+0.9)
		lint = 1.35e-009*(0.2*np.random.rand(1)+0.9)
		xj = 7.2e-009*(0.2*np.random.rand(1)+0.9)
		ndep = 1.2e+019*(0.2*np.random.rand(1)+0.9)
		process_data[9] = '+tnom    = 27              toxe    = ' + str(toxe[0])+ '        toxp    = ' + str(toxp[0])+'         toxm    = 6.5e-010'+'\n'
		process_data[10] =  '+dtox    = 2.5e-010        epsrox  = 3.9             wint    = ' + str(wint[0]) +'           lint    = ' + str(lint[0]) +'\n' 
		process_data[13] = '+lwl     = 0               wwl     = 0               xpart   = 0               toxref  = ' + str(toxref[0])+  '         xl      = -9e-9' +'\n'
		process_data[20] = '+dvtp1   = 0.1             lpe0    = 0               lpeb    = 0               xj      = ' + str(xj[0]) +'\n'
		process_data[21] = '+ngate   = 1e+023          ndep    = ' + str(ndep[0])+ '        nsd     = 2e+020          phin    = 0' + '\n'
	with open(filename, 'w') as file:
		file.writelines(process_data)


def delay_read(i):
#	filename = 'fa_del.sp'
	filename = 'FA_delays_MGK22nm.sp'
	mt0_filename = filename[:-3]+'.mt0'
	with open(filename, 'r') as file:
		delay_data = file.readlines()
		new_temp = 180*np.random.rand(1)-55
#		delay_data[8] = '.temp ' + str(new_temp[0]) + '\n'

		new_vdd = 0.8*(0.2*np.random.rand(1)+0.9)
#		delay_data[29] =  '+   pvdd=' + str(new_vdd[0]) + '\n'
		#delay_data[48] = 'Mp1 nodez nodea vdd vdd pmos l=' + str(i[1]) +'n w=' + str(i[0]) + 'n\n'
		#delay_data[49] = 'Mn1 nodez nodea gnd gnd nmos l=' + str(i[3]) +'n w=' + str(i[2]) + 'n\n'
                delay_data[206] = 'Mp1    1   nodea   vddd!   vddd!   pmos  l=' + str(i[1]) + 'n w=' + str(i[0]) + 'n\n'
                delay_data[207] = 'Mp2    1   nodeb   vddd!   vddd!   pmos  l=' + str(i[3]) + 'n w=' + str(i[2]) +'n\n'
                delay_data[208] = 'Mp3    nodecon nodec   1   vddd!   pmos l=' + str(i[5]) + 'n w=' + str(i[4]) +'n\n'
                delay_data[210] = 'Mn1    5   nodea   gndd!   gndd!   nmos  l=' + str(i[7]) + 'n w=' + str(i[6]) +'n\n'
                delay_data[211] = 'Mn2    5   nodeb   gndd!   gndd!   nmos  l=' + str(i[9]) + 'n w=' + str(i[8]) +'n\n'
                delay_data[212] = 'Mn3    nodecon nodec   5   gndd!   nmos l=' + str(i[11]) + 'n w=' + str(i[10]) +'n\n'
                delay_data[214] = 'Mp4    4   nodea   vddd!   vddd!   pmos  l=' + str(i[13]) + 'n w=' + str(i[12]) +'n\n'
                delay_data[215] = 'Mp5    nodecon nodeb   4   vddd!   pmos  l=' + str(i[15]) + 'n w=' + str(i[14]) +'n\n'
                delay_data[217] = 'Mn4    nodecon nodeb   node4   gndd!   nmos l=' + str(i[17]) + 'n w=' + str(i[16]) +'n\n'
                delay_data[218] = 'Mn5    node4   nodea   gndd!   gndd!   nmos l=' + str(i[19]) + 'n w=' + str(i[18]) +'n\n'
                delay_data[220] = 'Mp6    2   nodea   vddd!   vddd!   pmos  l=' + str(i[21]) + 'n w=' + str(i[20]) +'n\n'
                delay_data[221] = 'Mp7    2   nodeb   vddd!   vddd!   pmos l=' + str(i[23]) + 'n w=' + str(i[22]) +'n\n'
                delay_data[222] = 'Mp8    2   nodec   vddd!   vddd!   pmos l=' + str(i[25]) + 'n w=' + str(i[24]) +'n\n'
                delay_data[223] = 'Mp9    nodes0n nodecon 2   vddd!   pmos l=' + str(i[27]) + 'n w=' + str(i[26]) +'n\n'
                delay_data[225] = 'Mn6    3   nodea   gndd!   gndd!   nmos l=' + str(i[29]) + 'n w=' + str(i[28]) +'n\n'
                delay_data[226] = 'Mn7    3   nodeb   gndd!   gndd!   nmos l=' + str(i[31]) + 'n w=' + str(i[30]) +'n\n'
                delay_data[227] = 'Mn8    3   nodec   gndd!   gndd!   nmos l=' + str(i[33]) + 'n w=' + str(i[32]) +'n\n'
                delay_data[228] = 'Mn9    nodes0n nodecon 3   gndd!   nmos l=' + str(i[35]) + 'n w=' + str(i[34]) +'n\n'
                delay_data[230] = 'Mp10   9   nodea   vddd!   vddd!   pmos  l=' + str(i[37]) + 'n w=' + str(i[36]) +'n\n'
                delay_data[231] = 'Mp11   8   nodeb   9   vddd!   pmos    l=' + str(i[39]) + 'n w=' + str(i[38]) +'n\n'
                delay_data[232] = 'Mp12   nodes0n nodec   8   vddd!   pmos l=' + str(i[41]) + 'n w=' + str(i[40]) +'n\n'
                delay_data[234] = 'Mn10   7   nodea   gndd!   gndd!   nmos  l=' + str(i[43]) + 'n w=' + str(i[42]) +'n\n'
                delay_data[235] = 'Mn11   6   nodeb   7   gndd!   nmos    l=' + str(i[45]) + 'n w=' + str(i[44]) +'n\n'
                delay_data[236] = 'Mn12   nodes0n nodec   6   gndd!   nmos  l=' + str(i[47]) + 'n w=' + str(i[46]) +'n\n'
                delay_data[238] = 'Mp13   nodeco  nodecon vddd!   vddd!   pmos l=' + str(i[49]) + 'n w=' + str(i[48]) +'n\n'
                delay_data[239] = 'Mn13   nodeco  nodecon gndd!   gndd!   nmos l=' + str(i[51]) + 'n w=' + str(i[50]) +'n\n'
                delay_data[241] = 'Mp14   nodes0  nodes0n vddd!   vddd!   pmos l=' + str(i[53]) + 'n w=' + str(i[52]) +'n\n'
                delay_data[242] = 'Mn14   nodes0  nodes0n gndd!   gndd!   nmos  l=' + str(i[55]) + 'n w=' + str(i[54]) +'n\n'



	with open(filename, 'w') as file:
		file.writelines(delay_data)

	#process_change(i)
	call(["hspice", filename])
	
	with open(mt0_filename, 'r') as file:
		delay_values = file.readlines()
	
	delay_list = list()
	#values = delay_values[4].split()
	#delay_list.append(float(values[0]))
	#delay_list.append(float(values[1]))
	#print delay_values[5]
	values = delay_values[5].split()
        for v in values:
	#if v != ' ':
               	delay_list.append(float(v))
        values = delay_values[6].split()
        delay_list.append(float(values[0]))
        delay_list.append(float(values[1]))

        return delay_list

def leakage_read(i):
	filename = 'FA_leakage_MGK22nm.sp'
        ms0_filename = filename[:-3]+'.ms0'
        with open(filename, 'r') as file:
                leakage_data = file.readlines()
		new_temp = 180*np.random.rand(1)-55
#		leakage_data[18] = '.temp ' + str(new_temp[0]) +'\n'
               # leakage_data[38] = 'Mp1 nodez nodea vdd vdd pmos l=' + str(i[1]) +'n w=' + str(i[0]) + 'n\n'
               # leakage_data[39] = 'Mn1 nodez nodea gnd gnd nmos l=' + str(i[3]) +'n w=' + str(i[2]) + 'n\n'
                leakage_data[186] = 'Mp1    1   nodea   vdd   vdd   pmos  l=' + str(i[1]) + 'n w=' + str(i[0]) + 'n\n'
		leakage_data[187] = 'Mp1    1   nodea   vdd   vdd   pmos  l=' + str(i[1]) + 'n w=' + str(i[0]) + 'n\n'
		leakage_data[188] = 'Mp2    1   nodeb   vdd   vdd   pmos  l=' + str(i[3]) + 'n w=' + str(i[2]) +'n\n'
                leakage_data[189] = 'Mp3    nodecon nodec   1   vdd   pmos l=' + str(i[5]) + 'n w=' + str(i[4]) +'n\n'
                leakage_data[191] = 'Mn1    5   nodea   gnd   gnd   nmos  l=' + str(i[7]) + 'n w=' + str(i[6]) +'n\n'
                leakage_data[192] = 'Mn2    5   nodeb   gnd   gnd   nmos  l=' + str(i[9]) + 'n w=' + str(i[8]) +'n\n'
                leakage_data[193] = 'Mn3    nodecon nodec   5   gnd   nmos l=' + str(i[11]) + 'n w=' + str(i[10]) +'n\n'
                leakage_data[195] = 'Mp4    4   nodea   vdd   vdd   pmos  l=' + str(i[13]) + 'n w=' + str(i[12]) +'n\n'
                leakage_data[196] = 'Mp5    nodecon nodeb   4   vdd   pmos  l=' + str(i[15]) + 'n w=' + str(i[14]) +'n\n'
                leakage_data[198] = 'Mn4    nodecon nodeb   node4   gnd   nmos l=' + str(i[17]) + 'n w=' + str(i[16]) +'n\n'
                leakage_data[199] = 'Mn5    node4   nodea   gnd   gnd   nmos l=' + str(i[19]) + 'n w=' + str(i[18]) +'n\n'
                leakage_data[201] = 'Mp6    2   nodea   vdd   vdd   pmos  l=' + str(i[21]) + 'n w=' + str(i[20]) +'n\n'
                leakage_data[202] = 'Mp7    2   nodeb   vdd   vdd   pmos l=' + str(i[23]) + 'n w=' + str(i[22]) +'n\n'
                leakage_data[203] = 'Mp8    2   nodec   vdd   vdd   pmos l=' + str(i[25]) + 'n w=' + str(i[24]) +'n\n'
                leakage_data[204] = 'Mp9    nodes0n nodecon 2   vdd   pmos l=' + str(i[27]) + 'n w=' + str(i[26]) +'n\n'
                leakage_data[206] = 'Mn6    3   nodea   gnd   gnd   nmos l=' + str(i[29]) + 'n w=' + str(i[28]) +'n\n'
                leakage_data[207] = 'Mn7    3   nodeb   gnd   gnd   nmos l=' + str(i[31]) + 'n w=' + str(i[30]) +'n\n'
                leakage_data[208] = 'Mn8    3   nodec   gnd   gnd   nmos l=' + str(i[33]) + 'n w=' + str(i[32]) +'n\n'
                leakage_data[209] = 'Mn9    nodes0n nodecon 3   gnd   nmos l=' + str(i[35]) + 'n w=' + str(i[34]) +'n\n'
                leakage_data[211] = 'Mp10   9   nodea   vdd   vdd   pmos  l=' + str(i[37]) + 'n w=' + str(i[36]) +'n\n'
                leakage_data[212] = 'Mp11   8   nodeb   9   vdd   pmos    l=' + str(i[39]) + 'n w=' + str(i[38]) +'n\n'
                leakage_data[213] = 'Mp12   nodes0n nodec   8   vdd   pmos l=' + str(i[41]) + 'n w=' + str(i[40]) +'n\n'
                leakage_data[215] = 'Mn10   7   nodea   gnd   gnd   nmos  l=' + str(i[43]) + 'n w=' + str(i[42]) +'n\n'
                leakage_data[216] = 'Mn11   6   nodeb   7   gnd   nmos    l=' + str(i[45]) + 'n w=' + str(i[44]) +'n\n'
                leakage_data[217] = 'Mn12   nodes0n nodec   6   gnd   nmos  l=' + str(i[47]) + 'n w=' + str(i[46]) +'n\n'
                leakage_data[219] = 'Mp13   nodeco  nodecon vdd   vdd   pmos l=' + str(i[49]) + 'n w=' + str(i[48]) +'n\n'
                leakage_data[220] = 'Mn13   nodeco  nodecon gnd   gnd   nmos l=' + str(i[51]) + 'n w=' + str(i[50]) +'n\n'
                leakage_data[222] = 'Mp14   nodes0  nodes0n vdd   vdd   pmos l=' + str(i[53]) + 'n w=' + str(i[52]) +'n\n'
                leakage_data[223] = 'Mn14   nodes0  nodes0n gnd   gnd   nmos  l=' + str(i[55]) + 'n w=' + str(i[54]) +'n\n'



	with open(filename, 'w') as file:
		file.writelines(leakage_data)

	#process_change(i)
	call(["hspice", filename])
	
	with open(ms0_filename, 'r') as file:
		leakage_values = file.readlines()

	leakage_list = list()
	#l1 = leakage_values[10].split()
	#l2 = leakage_values[15].split()
	#leakage_list.append(float(l1[0]))
	#leakage_list.append(float(l2[0]))
        l1 = leakage_values[13].split()
        l2 = leakage_values[19].split()
        l3 = leakage_values[25].split()
        l4 = leakage_values[31].split()
        l5 = leakage_values[37].split()
        l6 = leakage_values[43].split()
        l7 = leakage_values[49].split()
        l8 = leakage_values[55].split()
        leakage_list.append(float(l1[0]))
        leakage_list.append(float(l2[0]))
        leakage_list.append(float(l3[0]))
        leakage_list.append(float(l4[0]))
        leakage_list.append(float(l5[0]))
        leakage_list.append(float(l6[0]))
        leakage_list.append(float(l7[0]))
        leakage_list.append(float(l8[0]))

	return leakage_list
	
wp1 = 132
lp1 = 22
wn1 = 132
ln1 = 22
xx = 0.4
alpha = 0.8
beta = 0.08
nt = 4
f = 1.000002
s = 0.0003
rho = 0.4
gamma = 0.6
no_glowworms = 20
max_iter = 2000
rsw = 1
rsl = 0.25

initial_wl = [float(wp1), float(lp1), float(wn1), float(ln1), float(wp1), float(lp1), float(wn1), float(ln1), float(wp1), float(lp1), float(wn1), float(ln1), float(wp1), float(lp1), float(wn1), float(ln1), float(wp1), float(lp1), float(wn1), float(ln1), float(wp1), float(lp1), float(wn1), float(ln1), float(wp1), float(lp1), float(wn1), float(ln1), float(wp1), float(lp1), float(wn1), float(ln1), float(wp1), float(lp1), float(wn1), float(ln1), float(wp1), float(lp1), float(wn1), float(ln1), float(wp1), float(lp1), float(wn1), float(ln1), float(wp1), float(lp1), float(wn1), float(ln1), float(wp1), float(lp1), float(wn1), float(ln1), float(wp1), float(lp1), float(wn1), float(ln1)]
range_bound = [(44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25),(44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25), (44,1000), (22,25)]

#dl_threshold = initial_dl()
d1 = delay_read(initial_wl)
l1 = leakage_read(initial_wl)

front_val = len(l1)
dl_threshold = list()
actual_dl = list()
for l in l1:
	dl_threshold.append(float(l)*f)
	actual_dl.append(float(l))
for d in d1:
	dl_threshold.append(float(d))
	actual_dl.append(float(l))

luciferin = np.zeros((no_glowworms, 1))
r0 = np.zeros((no_glowworms, 2))

wl_array = np.random.rand(no_glowworms, 56)
[rows, columns] = wl_array.shape
for j in range(columns):
	for i in range(rows):
		if (j%2)==0:
			wl_array[i][j] = 44 + wl_array[i][j]*956
		else:
			wl_array[i][j] = 22 + wl_array[i][j]*3

for i in range(rows):
	r0[i][0] = 2
	r0[i][1] = 0.2
	count = 0
	d_list = delay_read(wl_array[i])
	l_list = leakage_read(wl_array[i])
	for j in range(len(l_list)):
		if l_list[j]<dl_threshold[j]:
			count = count+1
	luciferin[i] = count
	
np.savetxt('inital_dl.txt', dl_threshold, delimiter=',')
i = 0
lsum = 0
for l in range(len(l_list)):
        lsum = lsum +float(l_list[l])

nleakage = lsum
pleakage = 0.6*nleakage
while i < max_iter:
#while lsum > pleakage:
	print(i)
	if i == 0:
		np.savetxt('first_iter_not.txt',wl_array, delimiter=',')
		total_l = list()
		total_d = list()
		for j in range(no_glowworms):
#			process_change(wl_array[j])
			l_list = leakage_read(wl_array[j])
			d_list = delay_read(wl_array[j])
			total_l.append(l_list)
			total_d.append(d_list)
		np.savetxt('first_iter_not_leakage.txt', total_l, delimiter=',')
		np.savetxt('first_iter_not_delay.txt', total_d, delimiter=',')
#	dl_params = initialization(wl_array[i])	
#	d_list = delay_read(wl_array[i])
	#fitness_obj()
	for j in range(no_glowworms):
		l_list = leakage_read(wl_array[j])
		l_sum = 0
		for l in range(len(l_list)):
			l_sum = l_sum +float(l_list[l])
		l_val =(10^7)* (l_sum/len(l_list))
 		luciferin[j] = (1-rho)*luciferin[j]+gamma*l_val
	for j in range(no_glowworms):
		d_list = delay_read(wl_array[j])
		for ss in range(len(d_list)):
			if d_list[ss] > dl_threshold[ss+front_val]:
				for st in range(28):
					wl_array[j][2*st] = 44 + np.random.rand(1)*956
                        		wl_array[j][2*st+1] = 22 + np.random.rand(1)*3

                        	count=0
                        	l_list = leakage_read(wl_array[j])
                        	d_list = delay_read(wl_array[j])
				for st in range(front_val):
                        		if l_list[st]<dl_threshold[st]:
                                		count = count+1
				for st in range(len(d_list)): 
					if d_list[st]<dl_threshold[st+front_val]:
						count = count+1
                        	luciferin[j] = count
				continue
		neighbours_list = list()
		denominator = 0

		max_prob = 0
		max_gw = 0
		for k in range(no_glowworms):
			for t in range(28):
				w_diff = abs(wl_array[j][2*t]-wl_array[k][2*t])
				l_diff = abs(wl_array[j][2*t+1]-wl_array[k][2*t+1])
				if j!=k and w_diff < r0[j][0] and l_diff < r0[j][1] and luciferin[j]<luciferin[k]:
					if t == 28:
						neighbours_list.append(k)
					else:
						continue
			denominator = denominator + luciferin[k]-luciferin[j]
		for k in neighbours_list:
			numerator = luciferin[k]-luciferin[j]
			pjk = numerator/denominator
			if pjk > max_prob:
				max_prob = pjk
			max_gw = k
		if not neighbours_list:	#put dnum
			for q in range(28):
				wl_array[j][2*q] = 44 + np.random.rand(1)*956
			#wl_array[j][2] = 44 + np.random.rand(1)*956
			#for q in range(28):
				wl_array[j][(2*q)+1] = 22 + np.random.rand(1)*3
			#wl_array[j][3] = 22 + np.random.rand(1)*3
                       # count=0
			#l_list = leakage_read(wl_array[j])
                        #d_list = delay_read(wl_array[j])
			#for m in range(len(l_list)):
                        #	if l_list[m]<dl_threshold[m]:
                         #       	count = count + 1
                        #luciferin[j] = count

		
		else:
			#count = 0
			for l in range(56):
				num = wl_array[max_gw][l]-wl_array[j][l]
				denom = abs(wl_array[max_gw][l]-wl_array[j][l])
				val1 = s*num/denom
				u = fitness(alpha, max_iter, i)
				wl_array[j][l] = wl_array[j][l] + val1 + u		

			#l_list = leakage_read(wl_array[j])
			#d_list = delay_read(wl_array[j])
		        #if l_list[0]<dl_threshold[0]:
                	#	count = count+1
        		#if l_list[1]<dl_threshold[1]:
                	#	count = count+1
        		#if d_list[0]<dl_threshold[2]:
                	#	count = count+1
        		#if d_list[1]<dl_threshold[3]:
                	#	count = count+1
        		#luciferin[j] = count
			count=0
                        l_list = leakage_read(wl_array[j])
			#for np in range(len(l_list)):
			#	l_list[np] = l_list[np] - l_list[np]*i/10000
                        #d_list = delay_read(wl_array[j])
                        for m in range(len(l_list)):
                                if float(l_list[m])<float(dl_threshold[m]):
                                        count = count+1
		#Nt = count(neighbours_list)
		r0[j][0] = min(rsw, max(0,r0[j][0], beta*(nt-len(neighbours_list))))
		r0[j][1] = min(rsl, max(0,r0[j][1], beta*(nt-len(neighbours_list))))

		lsum = 0			
		for l in range(len(l_list)):
                    lsum = lsum +float(l_list[l])

		print 'END OF ITERATION'
		print len(wl_array)
	i = i + 1 
np.savetxt('last_iter_not.txt', wl_array, delimiter=',')
total_l = list()
total_d = list()
good_l = list()
good_d = list()
good_wl = list()
total_leakage = list()
for j in range(no_glowworms):
	l_list = leakage_read(wl_array[j])
        d_list = delay_read(wl_array[j])
	l_sum = 0
	t_sum = 0
	for nk in range(len(l_list)):
		if i>=2000:
			l_list[nk] = l_list[nk] - l_list[nk]*xx
		else:
			l_list[nk] = l_list[nk] - l_list[nk]*i/10000
	for m in range(len(l_list)):
		l_sum = l_sum + float(l_list[m]) 
		t_sum = t_sum + float(dl_threshold[m])
	#l_sum=l_sum/len(l_list)
	#t_sum=t_sum/len(l_list)
	if (l_sum<f*t_sum):
		if l_list not in good_l:
			good_l.append(l_list)
			good_d.append(d_list)
			good_wl.append(wl_array[j])
	total_leakage.append(l_sum)	
        total_l.append(l_list)
        total_d.append(d_list)
np.savetxt('original_dl.txt',actual_dl, delimiter=',')
np.savetxt('last_iter_not_leakage.txt', total_l, delimiter=',')
np.savetxt('last_iter_not_delay.txt', total_d, delimiter=',')
np.savetxt('good_not_leakage.txt', good_l, delimiter=',')
np.savetxt('good_not_delay.txt', good_d, delimiter=',')
np.savetxt('good_not_wl.txt', good_wl, delimiter=',')
np.savetxt('end_sum_leakage.txt', total_leakage, delimiter=',')
np.savetxt('wl_array.txt', wl_array, delimiter=',')
