import numpy as np
import matplotlib.pyplot as plt
import random
import csv
import os

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

# side_len - length of unchanging dimensions
# k - varied dimension length
# k_count - number of varying dimensions
def create_neighbor_cable(side_len, dim, k=0, k_count=0):
    spin_count = side_len**(dim - k_count) * k**(k_count)
    neighbor_cable =  [[] for _ in range(spin_count)]
    
    for spin in range(0, spin_count):
        coord_arr = []
        coeff_arr = []
        if k_count == 0:
            n = 1 # base case
        else:
            n = k_count
        exp = dim - n
        tmp_spin = spin
        while exp >= 0:
            face_vol = side_len**exp * k**(n-1)
            coeff_arr.append(face_vol)
            coord = tmp_spin // face_vol
            coord_arr.append(coord)
            tmp_spin = tmp_spin % face_vol
            if n > 1:
                n -= 1
            else:
                exp -= 1
                
        for i in range(0, len(coord_arr)):
            # spin_omitted - value of spin excluding contribution of specified index
            spin_omitted = spin - coord_arr[i]*coeff_arr[i]
            
            # 2 neighbors in specified aisle i
            if k_count > i:
                coord_1 = (coord_arr[i]-1) % k
                coord_2 = (coord_arr[i]+1) % k
            else:
                coord_1 = (coord_arr[i]-1) % side_len
                coord_2 = (coord_arr[i]+1) % side_len
            neighbor_1 = spin_omitted + (coord_1) * coeff_arr[i]
            neighbor_2 = spin_omitted + (coord_2) * coeff_arr[i]

            neighbor_cable[spin].append(neighbor_1)
            neighbor_cable[spin].append(neighbor_2)
    return neighbor_cable

def create_coupling_arr(side_len, dim, k=0, k_count=0):
    spin_count = side_len**(dim - k_count) * k**(k_count)
    # J - array indexed by spin number
    # J[i] returns the coupling coefficients between all neighbors and spin i
    J =  [[] for _ in range(spin_count)]
    normal_arr = [[] for _ in range(dim)]
    
    # creating all coupling coefficients
    for i in range(0, dim):
        normal_arr[i] = np.random.normal(0, 1, size = get_shape(side_len, dim, k, k_count))

    # storing spin coupling coefficients
    for spin in range(0, spin_count):
        coord_arr = []
        if k_count == 0:
            n = 1
        else:
            n = k_count
        exp = dim - n
        tmp_spin = spin
        while exp >= 0:
            face_vol = side_len**exp * k**(n-1)
            side = tmp_spin // face_vol
            coord_arr.append(side)
            tmp_spin = tmp_spin % face_vol
            if n > 1:
                n -= 1
            else:
                exp -= 1

        for i in range(0, len(coord_arr)):
            coord_arr[i] -= 1
            neighbor_1 = normal_arr[i][tuple(coord_arr)] if normal_arr[i][tuple(coord_arr)] > 0 else normal_arr[i][tuple(coord_arr)] * -1
            coord_arr[i] += 1
            neighbor_2 = normal_arr[i][tuple(coord_arr)] if normal_arr[i][tuple(coord_arr)] > 0 else normal_arr[i][tuple(coord_arr)] * -1
            J[spin].append(neighbor_1)
            J[spin].append(neighbor_2)
    return J

def ham(spin, spin_dict, dim): #compute energy given a spin site number
    sum = 0
    spin_value = spin_dict[spin]
    for i in range(0,2*dim):
        nbr = neighbor_cable[spin][i]
        sum -= J[spin][i]*spin_value*spin_dict[nbr]
    return sum

def flippable(spin, spin_dict, dim):
    energy_unflipped = ham(spin, spin_dict, dim)
    if energy_unflipped > 0:
        return True
    else:
        return False

def create_active(spin_dict, dim):
    active_dict = {}
    for i in range(len(spin_dict)):
        if flippable(i, spin_dict, dim):
            active_dict[i] = spin_dict[i]
    return active_dict

def glauber(t, spin_dict, dim):
    t += 1/len(spin_dict)
    random.seed()
    spin = random.randint(0, len(spin_dict)-1)
    if flippable(spin, spin_dict, dim):
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
        if flippable(nbr, spin_dict, dim) and nbr not in active_dict.keys():
            active_dict[nbr] = spin_dict[nbr]
            # print(f"    Neighbor spin site {nbr} added to active list.")
        elif not flippable(nbr, spin_dict, dim) and nbr in active_dict.keys():
            active_dict.pop(nbr)
            # print(f"    Neighbor spin site {nbr} removed from active list.")
    return t

def overlap(spin_dict_1, spin_dict_2):
    total_overlap = 0
    for spin in spin_dict_1:
        total_overlap += spin_dict_1[spin]*spin_dict_2[spin]
    total_overlap = total_overlap/len(spin_dict_1)
    return total_overlap

def spin_test(spin_dict_1, spin_dict_2, dim, switch_time=0):
    t1 = 0
    # while t1 < switch_time: # time in sweeps - 1 sweep = N spin-flip attempts
    #     t1 = glauber(t1, spin_dict_1, dim)
    active_dict_1 = create_active(spin_dict_1, dim)
    while active_dict_1:
        t1 = monte_carlo(t1, spin_dict_1, active_dict_1, dim)
        
    t2 = 0
    # while t2 < switch_time:
    #     t2 = glauber(t2, spin_dict_2, dim)
    active_dict_2 = create_active(spin_dict_2, dim)
    while active_dict_2:
        t2 = monte_carlo(t2, spin_dict_2, active_dict_2, dim)
        
    q = overlap(spin_dict_1, spin_dict_2)
    return [q, (t1+t2)/2]

dim = 5
d = 20 # length of one side
# k = 5
k_count = 4
max_runs = 4000
# switch_time = 10

if k_count >= dim:
    print("k-count cannot be greater than or equal to the number of dimensions.")
    exit(1)
k_list = [10]  
for k in k_list:
    if k_count == 0:
        k = 0
    print(f"k = {k}, d = {d}, dim = {dim}")
    spin_count = d**(dim - k_count) * k**(k_count)

    neighbor_cable = create_neighbor_cable(d, dim, k=k, k_count=k_count)
    q_list = []
    t_list = []
    for run in range(1, max_runs+1):
        J = create_coupling_arr(d, dim, k=k, k_count=k_count)
        spin_dict_1 = {} # list of spins
        spin_dict_2 = {} 
        for i in range(0, spin_count):
            random.seed()
            spin_dict_1[i] = random.choice([-1,1])

        for key in spin_dict_1:
            spin_dict_2[key] = spin_dict_1[key]
        
        results = spin_test(spin_dict_1, spin_dict_2, dim)
        q = results[0]
        t = results[1]
        q_list.append(q)
        t_list.append(t)
        id = f"[{dim},{d},{k},{k_count},{max_runs},RF]"
        with open(os.path.join(os.getcwd(), f"Data\RF\L{d}", f'{dim}DlayersRFL{d}k{k}.csv'), 'a', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow([q, t, id])
        print(f"{id}: Spin test: {run}\n     Overlap: {q}\n     Time : {t} sweeps")
    mean = np.mean(q_list)
    std = np.std(q_list)
    mean_time = np.mean(t_list)
    time_std = np.std(t_list)

    print(f"{dim}-dimensional model with side-length {d}, k = {k}, k count = {k_count}, {max_runs} runs")
    print(f"    Mean overlap: {mean}")
    print(f"        Standard deviation: {std}")
    print(f"    Mean survival time: {mean_time}")
    print(f"        Standard deviation: {time_std}")
    
    with open(os.path.join(os.getcwd(), f"Data\RF\L{d}", f'{dim}DlayersRFL{d}.csv'), 'a', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow([k, mean, std, mean_time, time_std, id])