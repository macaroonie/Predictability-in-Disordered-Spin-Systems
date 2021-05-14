import numpy as np
import matplotlib.pyplot as plt
import random
import csv


def get_shape(d, dim, k, k_count):
    if dim == 1:
        return (d)
    elif dim == 2:
        if k_count == 1:
            return (k,d)
        else:
            return (d,d)
    elif dim == 3:
        if k_count == 1:
            return (k,d,d)
        elif k_count == 2:
            return (k,k,d)
        else:
            return (d,d,d)
    elif dim == 4:
        if k_count == 1:
            return (k,d,d,d)
        elif k_count == 2:
            return (k,k,d,d)
        elif k_count == 3:
            return (k,k,k,d)
        else:
            return (d,d,d,d)
    elif dim == 5:
        if k_count == 1:
            return (k,d,d,d,d)
        elif k_count == 2:
            return (k,k,d,d,d)
        elif k_count == 3:
            return (k,k,k,d,d)
        elif k_count == 4:
            return (k,k,k,k,d)
        else:
            return (d,d,d,d,d)
    elif dim == 6:
        if k_count == 1:
            return (k,d,d,d,d,d)
        elif k_count == 2:
            return (k,k,d,d,d,d)
        elif k_count == 3:
            return (k,k,k,d,d,d)
        elif k_count == 4:
            return (k,k,k,k,d,d)
        elif k_count == 5:
            return (k,k,k,k,k,d)
        else:
            return (d,d,d,d,d,d)
    elif dim == 7:
        if k_count == 1:
            return (k,d,d,d,d,d,d)
        elif k_count == 2:
            return (k,k,d,d,d,d,d)
        elif k_count == 3:
            return (k,k,k,d,d,d,d)
        elif k_count == 4:
            return (k,k,k,k,d,d,d)
        elif k_count == 5:
            return (k,k,k,k,k,d,d)
        elif k_count == 6:
            return (k,k,k,k,k,k,d)
        else:
            return (d,d,d,d,d,d,d)
    elif dim == 8:
        if k_count == 1:
            return (k,d,d,d,d,d,d,d)
        elif k_count == 2:
            return (k,k,d,d,d,d,d,d)
        elif k_count == 3:
            return (k,k,k,d,d,d,d,d)
        elif k_count == 4:
            return (k,k,k,k,d,d,d,d)
        elif k_count == 5:
            return (k,k,k,k,k,d,d,d)
        elif k_count == 6:
            return (k,k,k,k,k,k,d,d)
        elif k_count == 7:
            return (k,k,k,k,k,k,k,d)
        else:
            return (d,d,d,d,d,d,d,d)
    elif dim == 9:
        if k_count == 1:
            return (k,d,d,d,d,d,d,d,d)
        elif k_count == 2:
            return (k,k,d,d,d,d,d,d,d)
        elif k_count == 3:
            return (k,k,k,d,d,d,d,d,d)
        elif k_count == 4:
            return (k,k,k,k,d,d,d,d,d)
        elif k_count == 5:
            return (k,k,k,k,k,d,d,d,d)
        elif k_count == 6:
            return (k,k,k,k,k,k,d,d,d)
        elif k_count == 7:
            return (k,k,k,k,k,k,k,d,d)
        elif k_count == 8:
            return (k,k,k,k,k,k,k,k,d)
        else:
            return (d,d,d,d,d,d,d,d,d)
    elif dim == 10:
        if k_count == 1:
            return (k,d,d,d,d,d,d,d,d,d)
        elif k_count == 2:
            return (k,k,d,d,d,d,d,d,d,d)
        elif k_count == 3:
            return (k,k,k,d,d,d,d,d,d,d)
        elif k_count == 4:
            return (k,k,k,k,d,d,d,d,d,d)
        elif k_count == 5:
            return (k,k,k,k,k,d,d,d,d,d)
        elif k_count == 6:
            return (k,k,k,k,k,k,d,d,d,d)
        elif k_count == 7:
            return (k,k,k,k,k,k,k,d,d,d)
        elif k_count == 8:
            return (k,k,k,k,k,k,k,k,d,d)
        elif k_count == 9:
            return (k,k,k,k,k,k,k,k,k,d)
        else:
            return (d,d,d,d,d,d,d,d,d,d)

def create_neighbor_cable(side_len, dim, k=0, k_count=0):
    spin_count = side_len**(dim - k_count) * k**(k_count)
    neighbor_cable =  [[] for _ in range(spin_count)]
    
    for spin in range(0, spin_count):
        # print(f"spin = {spin}")
        side_arr = []
        vol_arr = []
        if k_count == 0:
            n = 1
        else:
            n = k_count
        d = dim - n
        tmp_spin = spin
        while d >= 0:
            face_vol = side_len**d * k**(n-1)
            vol_arr.append(face_vol)
            side = tmp_spin // face_vol
            # print(f"side = {side}")
            side_arr.append(side)
            tmp_spin = tmp_spin % face_vol
            if n > 1:
                n -= 1
            else:
                d -= 1
                
        for i in range(0, len(side_arr)):
            spin_sum = 0
            # calculate contribution of all aisles except specified aisle i
            for j in range(0, len(side_arr)):
                if i == j:
                    continue
                spin_sum += side_arr[j] * vol_arr[j]
            
            # 2 neighbors in specified aisle i
            if k_count > i:
                if side_arr[i] == 0:
                    neighbor_1 = spin_sum + (k-1) * vol_arr[i]
                    neighbor_2 = spin_sum + (side_arr[i]+1) * vol_arr[i]
                elif side_arr[i] == k-1:
                    neighbor_1 = spin_sum + (side_arr[i]-1) * vol_arr[i]
                    neighbor_2 = spin_sum # at row 0
                else:
                    neighbor_1 = spin_sum + (side_arr[i]-1) * vol_arr[i]
                    neighbor_2 = spin_sum + (side_arr[i]+1) * vol_arr[i]
            else:
                if side_arr[i] == 0:
                    neighbor_1 = spin_sum + (side_len-1) * vol_arr[i]
                    neighbor_2 = spin_sum + (side_arr[i]+1) * vol_arr[i]
                elif side_arr[i] == side_len-1:
                    neighbor_1 = spin_sum + (side_arr[i]-1) * vol_arr[i]
                    neighbor_2 = spin_sum # at row 0
                else:
                    neighbor_1 = spin_sum + (side_arr[i]-1) * vol_arr[i]
                    neighbor_2 = spin_sum + (side_arr[i]+1) * vol_arr[i]
            neighbor_cable[spin].append(neighbor_1)
            neighbor_cable[spin].append(neighbor_2)
    return neighbor_cable

def create_coupling_arr(side_len, dim, k=0, k_count=0):
    spin_count = side_len**(dim - k_count) * k**(k_count)
    J =  [[] for _ in range(spin_count)]
    normal_arr = [[] for _ in range(dim)]
    for i in range(0, dim):
        normal_arr[i] = np.random.normal(0, 1, size = get_shape(side_len, dim, k, k_count))

    for spin in range(0, spin_count):
        side_arr = []
        if k_count == 0:
            n = 1
        else:
            n = k_count
        d = dim - n
        tmp_spin = spin
        while d >= 0:
            face_vol = side_len**d * k**(n-1)
            side = tmp_spin // face_vol
            side_arr.append(side)
            tmp_spin = tmp_spin % face_vol
            if n > 1:
                n -= 1
            else:
                d -= 1

        for i in range(0, len(side_arr)):
            spin_sum = 0
            side_arr[i] -= 1
            neighbor_1 = normal_arr[i][tuple(side_arr)]
            side_arr[i] += 1
            neighbor_2 = normal_arr[i][tuple(side_arr)]
            J[spin].append(neighbor_1)
            J[spin].append(neighbor_2)
            # print(f"neighbor 1 = {neighbor_1}")
            # print(f"neighbor 2 = {neighbor_2}")
    return J

def ham(spin, spin_dict, dim): #compute energy given a spin site number
    sum = 0
    spin_value = spin_dict[spin]
    for i in range(0,2*dim):
        nbr = neighbor_cable[spin][i]
        sum -= J[spin][i]*spin_value*spin_dict[nbr]
    return sum

def can_flip(spin, spin_dict, dim):
    energy_unflipped = ham(spin, spin_dict, dim)
    if energy_unflipped > 0:
        return True
    else:
        return False

def create_active(spin_dict, dim):
    active_dict = {}
    for i in range(len(spin_dict)):
        if can_flip(i, spin_dict, dim):
            active_dict[i] = spin_dict[i]
    return active_dict

def glauber(t, spin_dict, dim):
    t += 1/len(spin_dict)
    random.seed()
    spin = random.randint(0, len(spin_dict)-1)
    if can_flip(spin, spin_dict, dim):
        spin_dict[spin] = spin_dict[spin] * (-1)
    return t

def monte_carlo(t, spin_dict, active_dict, dim):
    random.seed()
    spin = random.choice(list(active_dict.keys()))
    
    # print(f"Active list length: {len(active_dict)}")
    t += 1/len(active_dict)
    spin_dict[spin] = spin_dict[spin]*-1
    active_dict.pop(spin)
    # print(f"Active list length after POP: {len(active_dict)}")
    # if dim == 5 and len(active_dict) == 0:
    #     return t
    # print(f"    Spin site {spin} has been removed from active list")

    # if neighbor can be flipped and not in active list, add to list
    # if neighbor cannot be flipped and in active list, remove from list
    for nbr in neighbor_cable[spin]:
        if can_flip(nbr, spin_dict, dim) and nbr not in active_dict.keys():
            active_dict[nbr] = spin_dict[nbr]
            # print(f"    Neighbor spin site {nbr} added to active list.")
        elif not can_flip(nbr, spin_dict, dim) and nbr in active_dict.keys():
            active_dict.pop(nbr)
            # print(f"    Neighbor spin site {nbr} removed from active list.")
    return t

def overlap(spin_dict_1, spin_dict_2):
    total_overlap = 0
    for spin in spin_dict_1:
        total_overlap += spin_dict_1[spin]*spin_dict_2[spin]
    total_overlap = total_overlap/len(spin_dict_1)
    return total_overlap

def spin_test(spin_dict_1, spin_dict_2, dim, switch_time):
    t1 = 0
    while t1 < switch_time: # time in sweeps - 1 sweep = N spin-flip attempts
        t1 = glauber(t1, spin_dict_1, dim)
    active_dict_1 = create_active(spin_dict_1, dim)
    while active_dict_1:
        t1 = monte_carlo(t1, spin_dict_1, active_dict_1, dim)
    # print("Twin 1 has reached absorbing state.")
    t2 = 0
    while t2 < switch_time:
        t2 = glauber(t2, spin_dict_2, dim)
    active_dict_2 = create_active(spin_dict_2, dim)
    while active_dict_2:
        t2 = monte_carlo(t2, spin_dict_2, active_dict_2, dim)
    # print("Twin 2 has reached absorbing state.")
    q = overlap(spin_dict_1, spin_dict_2)
    return [q, (t1+t2)/2]

dim = 7
d = 10 # length of one side
# k = 5
k_count = 6
switch_time = 10

if k_count >= dim:
    print("k-count cannot be greater than or equal to the number of dimensions.")
    exit(1)
k_list = [4,5,6,7,8,9]
for k in k_list:
    print(f"{k} layers")
    spin_count = d**(dim - k_count) * k**(k_count)

    neighbor_cable = create_neighbor_cable(d, dim, k=k, k_count=k_count)
    q_list = []
    t_list = []
    max_runs = 30000
    for run in range(1, max_runs+1):
        J = create_coupling_arr(d, dim, k=k, k_count=k_count)
        spin_dict_1 = {} # list of spins
        spin_dict_2 = {}
        for i in range(0, spin_count):
            random.seed()
            spin_dict_1[i] = random.choice([-1,1])

        for key in spin_dict_1:
            spin_dict_2[key] = spin_dict_1[key]
        
        results = spin_test(spin_dict_1, spin_dict_2, dim, switch_time)
        q = results[0]
        t = results[1]
        q_list.append(q)
        t_list.append(t)
        id = f"[{dim},{d},{k},{k_count},{max_runs},{switch_time}]"
        print(f"{id}: Spin test: {run}\n     Overlap: {q}\n     Time : {t} sweeps")
    mean = np.mean(q_list)
    std = np.std(q_list)
    mean_time = np.mean(t_list)

    print(f"{dim}-dimensional model with side-length {d}, k = {k}, k count = {k_count}, {max_runs} runs, switch_time = {switch_time} sweeps")
    print(f"    Mean overlap: {mean}")
    print(f"    Standard deviation: {std}")
    print(f"    Mean survival time: {mean_time}")
    
    with open('7Dlayers.csv', 'a', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow([k, mean, std, mean_time, id])



