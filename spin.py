import numpy as np
import matplotlib.pyplot as plt
import random
import csv

spin_num = 1000 # number of spin-sites
switch_time = 10 # number of sweeps before switching to KMC
test_num = 30001

# creating neighbor list for easy access
neighbor_cable =  [[] for y in range(spin_num)] 
for i in range(1, len(neighbor_cable)-1):
    neighbor_cable[i] = [i-1, i+1]
neighbor_cable[0] = [spin_num-1, 1]
neighbor_cable[spin_num-1] = [spin_num-2, 0]

fields = ['Spin test #', 'q (spin overlap)']
with open('spin3.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(fields)
    
def hamiltonian(spin, dict, flip = False): #compute energy of site and neighbors given a spin site number
    sum = 0
    spin_value = dict[spin]
    prev_neighbor = neighbor_cable[spin][0]
    next_neighbor = neighbor_cable[spin][1]
    if flip == True:
        spin_value = -spin_value
    sum -= s[prev_neighbor]*dict[prev_neighbor]*spin_value
    sum -= s[spin]*spin_value*dict[next_neighbor]
    return sum

def flip_check(spin, dict):
    energy_unflipped = hamiltonian(spin, dict)
    energy_flipped = hamiltonian(spin, dict, flip = True)

    # print(f"Spin site {spin}:")
    # print(f"    Energy before flip is {energy_unflipped}")
    # print(f"    Energy after flip is {energy_flipped}")

    if energy_flipped < energy_unflipped:
        dict[spin] = dict[spin] * (-1)
        # print(f"    Spin site {spin} has been flipped.")
        return True
    else:
        # print(f"    Spin site {spin} has not been flipped.")
        return False
    
def glauber(t, dict):
    t += 1/spin_num
    random.seed()
    spin = random.randint(0, spin_num-1)
    flip_check(spin, dict)
    return t

def monte_carlo(t, dict, active_dict):
    random.seed()
    spin = random.choice(list(active_dict.keys()))
    prev_neighbor = neighbor_cable[spin][0]
    next_neighbor = neighbor_cable[spin][1]
    
    t += 1/len(active_dict)
    if not flip_check(spin, dict):
        active_dict.pop(spin)
    else: # spin has been flipped, neighbor sites are now active
        active_dict[spin] = active_dict[spin] * (-1)
        if spin - 1 not in active_dict.keys():
            active_dict[prev_neighbor] = dict[prev_neighbor]
        if spin + 1 not in active_dict.keys():
            active_dict[next_neighbor] = dict[next_neighbor]
    return t

def q_infinity(N):
    overlap = 0
    for spin in spin_dict_1:
        overlap += spin_dict_1[spin]*spin_dict_2[spin]
    overlap = overlap/N
    return overlap

def spin_test():
    t = 0
    while t < switch_time: # t measured in sweeps - one sweep = N spin-flip attempts
        t = glauber(t, spin_dict_1)
    active_dict_1 = {}
    for key in spin_dict_1:
        active_dict_1[key] = spin_dict_1[key]
    while active_dict_1:
        t = monte_carlo(t, spin_dict_1, active_dict_1)

    t = 0
    while t < switch_time: # t measured in sweeps - one sweep = N spin-flip attempts
        t = glauber(t, spin_dict_2)
    
    active_dict_2 = {}
    for key in spin_dict_2:
        active_dict_2[key] = spin_dict_1[key]
    while active_dict_2:
        t = monte_carlo(t, spin_dict_2, active_dict_2)
    q = q_infinity(spin_num)
    return q

q_list = []
for test in range(1, test_num):
    s = np.random.normal(0, 1, spin_num)

    spin_dict_1 = {} # list of spins
    spin_dict_2 = {}
    for i in range(0, spin_num):
        random.seed()
        spin_dict_1[i] = random.choice([-1,1])

    for key in spin_dict_1:
        spin_dict_2[key] = spin_dict_1[key]
    
    q = spin_test()
    q_list.append(q)
    # print(f"Spin test: {test}")
    # print(f"    Overlap: {spin_test()}")
    with open('spin3.csv', 'a', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow([test, q])
mean = np.mean(q_list)
std = np.std(q_list)

with open('spin3.csv', 'a', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['Mean', mean])
    csvwriter.writerow(['Standard deviation', std])
print(f"Mean overlap: {mean}")
print(f"Standard deviation: {std}")