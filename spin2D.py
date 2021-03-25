# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
import numpy as np
import matplotlib.pyplot as plt
import random
import csv

max_runs = 30000
L = 100 # length of one side
N = L**2 # number of spin-sites
J_hor = np.random.normal(0, 1, size = (L,L))
J_vert = np.random.normal(0, 1, size = (L,L))
switch_time = 1


# neighbor0, neighbor1, neighbor2, neighbor3 = left, right, above, below neighbors respectively
neighbor_cable =  [[] for y in range(N)] 
for spin in range(0, len(neighbor_cable)):
    row = spin // L
    column = spin % L

    # N - 1 = last row, last column
    if column == 0:
        neighbor0 = L*(row) + (L-1)
        neighbor1 = L*(row) + (column+1)
    elif column == L-1:
        neighbor0 = L*(row) + (column-1)
        neighbor1 = L*(row)
    else:
        neighbor0 = L*(row) + (column-1)
        neighbor1 = L*(row) + (column+1)

    if row == 0:
        neighbor2 = L*(L-1) + column
        neighbor3 = L*(row+1) + column
    elif row == L-1:
        neighbor2 = L*(row-1) + column
        neighbor3 = column
    else:
        neighbor2 = L*(row-1) + column
        neighbor3 = L*(row+1) + column
    
    neighbor_cable[spin] = [neighbor0, neighbor1, neighbor2, neighbor3] 

# initializing J couplings for all sites
# J[i] = [right edge, down edge] where edges are couplings
J =  [[] for y in range(N)]
for spin in range(0, len(neighbor_cable)):
    row = spin // L
    column = spin % L
    if column == 0:
        left_edge = J_hor[row][L-1]
    else:
        left_edge = J_hor[row][column-1]
    if row == 0:
        above_edge = J_vert[L-1][column]
    else:
        above_edge = J_vert[row-1][column]
    right_edge = J_hor[row][column]
    below_edge = J_vert[row][column]
    
    J[spin] = [left_edge, right_edge, above_edge, below_edge]

def ham2D(spin, spin_dict, flip = False): #compute energy given a spin site number
    sum = 0
    spin_value = spin_dict[spin]
    if flip == True:
        spin_value = -spin_value
    for i in range(0,4):
        nbr = neighbor_cable[spin][i]
        sum -= J[spin][i]*spin_value*spin_dict[nbr]
    return sum

def flip_check(spin, spin_dict):
    energy_unflipped = ham2D(spin, spin_dict)

    if energy_unflipped > 0:
        spin_dict[spin] = spin_dict[spin] * (-1)
        return True
    else:
        return False

def glauber(t, spin_dict):
    t += 1/N
    random.seed()
    spin = random.randint(0, N-1)
    flip_check(spin, spin_dict)
    return t

def monte_carlo(t, spin_dict, active_dict):
    random.seed()
    spin = random.choice(list(active_dict.keys()))
    t += 1/len(active_dict)
    if not flip_check(spin, spin_dict):
        active_dict.pop(spin)
    else: # spin has been flipped, neighbor sites are now active
        active_dict[spin] = active_dict[spin] * (-1)
        for nbr in neighbor_cable[spin]:
            if nbr not in active_dict.keys():
                active_dict[nbr] = spin_dict[nbr]
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
    q = q_infinity(N)
    return q
with open('spin2d.csv', 'a', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["2D spin overlap with L = 10 and 30000 runs"])

q_list = []
for run in range(1, max_runs+1):
    J_hor = np.random.normal(0, 1, size = (L,L))
    J_vert = np.random.normal(0, 1, size = (L,L))
    
    spin_dict_1 = {} # list of spins
    spin_dict_2 = {}
    for i in range(0, N):
        random.seed()
        spin_dict_1[i] = random.choice([-1,1])

    for key in spin_dict_1:
        spin_dict_2[key] = spin_dict_1[key]
    
    q = spin_test()
    q_list.append(q)
    print(f"Spin run: {run} w/ overlap: {q}")
    # with open('spin2d.csv', 'a', newline='') as csvfile:
    #     csvwriter = csv.writer(csvfile)
    #     csvwriter.writerow([run, q])
mean = np.mean(q_list)
std = np.std(q_list)

# with open('spin2d.csv', 'a', newline='') as csvfile:
#     csvwriter = csv.writer(csvfile)
#     csvwriter.writerow(['Mean', mean])
#     csvwriter.writerow(['Standard deviation', std])
print(f"Mean overlap: {mean}")
print(f"Standard deviation: {std}")