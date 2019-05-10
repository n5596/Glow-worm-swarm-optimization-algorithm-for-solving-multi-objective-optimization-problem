import numpy as np
import time
import os
from subprocess import call


def delay_read(i):
	filename = 'not_del.sp'
	mt0_filename = filename[:-3]+'.mt0'
	with open(filename, 'r') as file:
		delay_data = file.readlines()
		delay_data[48] = 'Mp1 nodez nodea vdd vdd pmos l=' + str(i[1]) +'n w=' + str(i[0]) + 'n\n'
		delay_data[49] = 'Mn1 nodez nodea gnd gnd nmos l=' + str(i[3]) +'n w=' + str(i[2]) + 'n\n'

	with open(filename, 'w') as file:
		file.writelines(delay_data)

	call(["hspice", filename])
	
	with open("not_del.mt0", 'r') as file:
		delay_values = file.readlines()
	
	delay_list = list()
	values = delay_values[4].split()
	delay_list.append(float(values[0]))
	delay_list.append(float(values[1]))
	return delay_list

def leakage_read(i):
	filename = 'not_leak.sp'
        ms0_filename = filename[:-3]+'.ms0'
        with open(filename, 'r') as file:
                leakage_data = file.readlines()
                leakage_data[38] = 'Mp1 nodez nodea vdd vdd pmos l=' + str(i[1]) +'n w=' + str(i[0]) + 'n\n'
                leakage_data[39] = 'Mn1 nodez nodea gnd gnd nmos l=' + str(i[3]) +'n w=' + str(i[2]) + 'n\n'

	with open(filename, 'w') as file:
		file.writelines(leakage_data)

	call(["hspice", filename])
	
	with open(ms0_filename, 'r') as file:
		leakage_values = file.readlines()

	leakage_list = list()
	l1 = leakage_values[10].split()
	l2 = leakage_values[15].split()
	leakage_list.append(float(l1[0]))
	leakage_list.append(float(l2[0]))
	return leakage_list

	
wp1 = 132
lp1 = 22
wn1 = 132
ln1 = 22

beta = 0.08
nt = 4
s = 0.0003
rho = 0.4
gamma = 0.6
no_glowworms = 20
max_iter = 200
rsw = 20
rsl = 5

initial_wl = [float(wp1), float(lp1), float(wn1), float(ln1)]
range_bound = [(44,1000), (22,25), (44,1000), (22,25)]

#dl_threshold = initial_dl()
d1 = delay_read(initial_wl)
l1 = leakage_read(initial_wl)

dl_threshold = [l1[0],l1[1],d1[0],d1[1]]
np.savetxt('initial_dl.txt', dl_threshold, delimiter=',')
luciferin = np.zeros((no_glowworms, 1))
r0 = np.zeros((no_glowworms, 2))

wl_array = np.random.rand(no_glowworms, 4)
[rows, columns] = wl_array.shape
for j in range(columns):
	for i in range(rows):
		if j==0 or j==2:
			wl_array[i][j] = 44 + wl_array[i][j]*956
		else:
			wl_array[i][j] = 22 + wl_array[i][j]*3

for i in range(rows):
	r0[i][0] = 4
	r0[i][1] = 0.2
	count = 0
	d_list = delay_read(wl_array[i])
	l_list = leakage_read(wl_array[i])
	if l_list[0]<dl_threshold[0]:
		count = count+1
	if l_list[1]<dl_threshold[1]:
		count = count+1
	if d_list[0]<dl_threshold[2]:
		count = count+1
	if d_list[1]<dl_threshold[3]:
		count = count+1
	luciferin[i] = count
	

i = 0
while i<max_iter:
	if i == 0:
		np.savetxt('testing_first_iter_not.txt',wl_array, delimiter=',')
		total_l = list()
		total_d = list()
		for j in range(no_glowworms):
			l_list = leakage_read(wl_array[j])
			d_list = delay_read(wl_array[j])
			total_l.append(l_list)
			total_d.append(d_list)
		np.savetxt('testing_first_iter_not_leakage.txt', total_l, delimiter=',')
		np.savetxt('testing_first_iter_not_delay.txt', total_d, delimiter=',')
#	dl_params = initialization(wl_array[i])	
#	d_list = delay_read(wl_array[i])
	for j in range(no_glowworms):
		l_list = leakage_read(wl_array[j])
		l_sum = (10^(-7))/(float(l_list[0]) + float(l_list[1]))
		d_list = delay_read(wl_array[j])
		#d_sum = float(d_list[0]) + float(d_list[1])
		#if (d_list[0]<=dl_threshold[2] and d_list[1]<=dl_threshold[3])
 		luciferin[j] = (1-rho)*luciferin[j]+gamma*l_sum
	
	for j in range(no_glowworms):
		d_list = delay_read(wl_array[j])
		if d_list[0]>dl_threshold[2] or d_list[1]>dl_threshold[3]:
			wl_array[j][0] = 44 + np.random.rand(1)*956
                        wl_array[j][2] = 44 + np.random.rand(1)*956
			wl_array[j][1] = 22 + np.random.rand(1)*3
			wl_array[j][3] = 22 + np.random.rand(1)*3

                        count=0
                        l_list = leakage_read(wl_array[j])
                        d_list = delay_read(wl_array[j])
                        if l_list[0]<dl_threshold[0]:
                                count = count+1
                        if l_list[1]<dl_threshold[1]:
                                count = count+1
                        if d_list[0]<dl_threshold[2]:
                                count = count+1
                        if d_list[1]<dl_threshold[3]:
                                count = count+1
                        luciferin[j] = count
			continue
		neighbours_list = list()
		denominator = 0

		max_prob = 0
		max_gw = -1
		for k in range(no_glowworms):
			w1_diff = abs(wl_array[j][0]-wl_array[k][0])
			l1_diff = abs(wl_array[j][1]-wl_array[k][1])
			w2_diff = abs(wl_array[j][2]-wl_array[k][2])
			l2_diff = abs(wl_array[j][3]-wl_array[k][3])
			if j!=k and w1_diff < r0[j][0] and l1_diff < r0[j][1] and w2_diff < r0[j][0] and l2_diff < r0[j][1] and luciferin[j]<luciferin[k]:
				neighbours_list.append(k)
				denominator = denominator + luciferin[k]-luciferin[j]
		for k in neighbours_list:
			numerator = luciferin[k]-luciferin[j]
			pjk = numerator/denominator
			if pjk > max_prob:
				max_prob = pjk
				max_gw = k
		if not neighbours_list:	#put dnum
			wl_array[j][0] = 44 + np.random.rand(1)*956
			wl_array[j][2] = 44 + np.random.rand(1)*956
			wl_array[j][1] = 22 + np.random.rand(1)*3
			wl_array[j][3] = 22 + np.random.rand(1)*3
                        count=0
			l_list = leakage_read(wl_array[j])
                        d_list = delay_read(wl_array[j])
                        if l_list[0]<dl_threshold[0]:
                        	count = count+1
                        if l_list[1]<dl_threshold[1]:
                                count = count+1
                        if d_list[0]<dl_threshold[2]:
                                count = count+1
                        if d_list[1]<dl_threshold[3]:
                                count = count+1
                        luciferin[j] = count

		
		else:
			#count = 0
			for l in range(4):
				num = wl_array[max_gw][l]-wl_array[j][l]
				denom = abs(wl_array[max_gw][l]-wl_array[j][l])
				val1 = s*num/denom
				wl_array[j][l] = wl_array[j][l] + val1
			l_list = leakage_read(wl_array[j])
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
		#Nt = count(neighbours_list)
		r0[j][0] = min(rsw, max(0,r0[j][0], beta*(nt-len(neighbours_list))))
		r0[j][1] = min(rsl, max(0,r0[j][1], beta*(nt-len(neighbours_list))))
	i = i + 1 
np.savetxt('testing_last_iter_not.txt', wl_array, delimiter=',')
total_l = list()
total_d = list()
good_l = list()
good_d = list()
good_wl = list()
for j in range(no_glowworms):
	l_list = leakage_read(wl_array[j])
        d_list = delay_read(wl_array[j])
	if l_list[0]+l_list[1] < dl_threshold[0]+dl_threshold[1] :
		good_l.append(l_list)
		good_d.append(d_list)
		good_wl.append(wl_array[j])
		
        total_l.append(l_list)
        total_d.append(d_list)
np.savetxt('ntesting_last_iter_not_leakage.txt', total_l, delimiter=',')
np.savetxt('ntesting_last_iter_not_delay.txt', total_d, delimiter=',')
np.savetxt('ntesting_good_not_leakage.txt', good_l, delimiter=',')
np.savetxt('ntesting_good_not_delay.txt', good_d, delimiter=',')
np.savetxt('ntesting_good_not_wl.txt', good_wl, delimiter=',')
