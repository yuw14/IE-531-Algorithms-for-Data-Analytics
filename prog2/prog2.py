# Reference: The Code Sample in Compass, prog2.pdf

import sys
import argparse
import random
import numpy as np 
import time
import math
import matplotlib.pyplot as plt  
sys.setrecursionlimit(3000)

# sort the array and pick the k-th smallest element from the sorted-array
def sort_and_select(current_array, k) :
    # sort the array
    sorted_current_array = np.sort(current_array)
    return sorted_current_array[k]

def deterministic_select(current_array, k) :

    global median_of_median_from_sorting
    global median_of_median_from_recursion
    
    if (len(current_array) <= m) :
        # just use any method to pick the k-th smallest element in the array
        # I am using the sort-and-select method here
        return sort_and_select(current_array, k)
    else : 
        # I need this array to compute the median-of-medians...
#             C
        medians_of_smaller_arrays_of_size_five = []
        
        # first, split current_array into smaller arrays with 5 elements each
        # there might be a better way than what I am doing... but this will work... 
        for i in range(0,len(current_array),5):
#             B
            smaller_array_of_size_five = []
            smaller_array_of_size_five.extend([current_array[i]])
            if ((i + 1) < len(current_array)) :
                smaller_array_of_size_five.extend([current_array[i+1]])
            if ((i + 2) < len(current_array)) :
                smaller_array_of_size_five.extend([current_array[i+2]])
            if ((i + 3) < len(current_array)) :
                smaller_array_of_size_five.extend([current_array[i+3]])
            if ((i + 4) < len(current_array)) :
                smaller_array_of_size_five.extend([current_array[i+4]])
            
            # we need each of these cases as len(smaller_array_of_size_five) can be anything between 1 and 5
            # based on len(smaller_array_of_size_five) we are computing the median of smaller_array_of_size_five for each case
            if (len(smaller_array_of_size_five) == 5) :
                medians_of_smaller_arrays_of_size_five.extend([deterministic_select(smaller_array_of_size_five,3)])
            if (len(smaller_array_of_size_five) == 4) :
                medians_of_smaller_arrays_of_size_five.extend([(deterministic_select(smaller_array_of_size_five,2)+deterministic_select(smaller_array_of_size_five,3))/2])
            if (len(smaller_array_of_size_five) == 3) :
                medians_of_smaller_arrays_of_size_five.extend([deterministic_select(smaller_array_of_size_five,2)])
            if (len(smaller_array_of_size_five) == 2) :
                medians_of_smaller_arrays_of_size_five.extend([(smaller_array_of_size_five[0]+smaller_array_of_size_five[1])/2])
            if (len(smaller_array_of_size_five) == 1) :
                medians_of_smaller_arrays_of_size_five.extend([smaller_array_of_size_five[0]])
            
        # compute the meadian of the medians_of_smaller_arrays_of_size_five array by recursion
        p = deterministic_select(medians_of_smaller_arrays_of_size_five, int(len(medians_of_smaller_arrays_of_size_five)/2))
        
        #check if we got the median-of-medians on the first go-around of the recursion
        if (len(current_array) == array_size) :
            median_of_median_from_sorting = sort_and_select(medians_of_smaller_arrays_of_size_five, 
                                                           math.ceil(len(medians_of_smaller_arrays_of_size_five)/2))
            median_of_median_from_recursion = p
            # print("n is "+str(array_size)+" m is "+str(m)+" k is "+str(k))

        
        # split the current_array into three sub-arrays: Less_than_p, Equal_to_p and Greater_than_p
        Less_than_p = []
        Equal_to_p = []
        Greater_than_p = []
        for x in current_array : 
            if (x < p) : 
                Less_than_p.extend([x])
            if (x == p) : 
                Equal_to_p.extend([x])
            if (x > p) : 
                Greater_than_p.extend([x])
                
        if (k < len(Less_than_p)) :
            return deterministic_select(Less_than_p, k)
        elif (k >= len(Less_than_p) + len(Equal_to_p)) : 
            return deterministic_select(Greater_than_p, k - len(Less_than_p) - len(Equal_to_p))
        else :
            return p


# Experimentally verifying/refuting the existence of an optimal value for m.
number_of_trials = 10
results = np.zeros((10,5),dtype = float)
index_of_n = 0
for n in range (1000, 10001, 1000):
    index_of_m = 0
    for the_real_m in [5,7,9,11,13]:
        the_amount_of_time_of_each_trial = []
        for i in range (number_of_trials) :
            the_array = [random.randint(1,100*n) for _ in range(n)]
            k = math.ceil(n/2)
            array_size = n
            m = the_real_m
            t0 = time.time()
            the_kth_smallest = deterministic_select(the_array, k)
#             print("Median-of-Medians (Computed by Sorting)     = ", median_of_median_from_sorting)
#             print("Median-of-Medians (Computed by Recursively) = ", median_of_median_from_recursion)
#             print("the_kth_smallest is "+str(the_kth_smallest))
            t1 = time.time()
#             print("runtime is "+str(t1-t0))
            the_amount_of_time_of_each_trial.append( t1-t0 )
        mean_number = np.mean(the_amount_of_time_of_each_trial)
        results[index_of_n][index_of_m] = mean_number
        # print("mean run time"+str(mean_number))
        index_of_m += 1
    index_of_n += 1

# print(results)
index_of_n = 0
for n in range (1000, 10001, 1000):
    print("When n="+str(n)+", the running time of each m is "+str(results[index_of_n]))
    index_of_n += 1


# Plot
x = [5,7,9,11,13]
y0 = np.array(results[0],dtype=float)
y1 = results[1]
y2 = results[2]
y3 = results[3]
y4 = results[4]
y5 = results[5]
y6 = results[6]
y7 = results[7]
y8 = results[8]
y9 = results[9]
plt.plot(x,y0,'o-',color = 'r')#o-:circle
plt.title("Array Size = 1000")
plt.xlabel("m")
plt.ylabel("runtime")
plt.savefig('./n=1000.png')
plt.show()

plt.plot(x,y1,'o-',color = 'r')#o-:circle
plt.title("Array Size = 2000")
plt.xlabel("m")
plt.ylabel("runtime")
plt.savefig('./n=2000.png')
plt.show()

plt.plot(x,y2,'o-',color = 'r')#o-:circle
plt.title("Array Size = 3000")
plt.xlabel("m")
plt.ylabel("runtime")
plt.savefig('./n=3000.png')
plt.show()

plt.plot(x,y3,'o-',color = 'r')#o-:circle
plt.title("Array Size = 4000")
plt.xlabel("m")
plt.ylabel("runtime")
plt.savefig('./n=4000.png')
plt.show()

plt.plot(x,y4,'o-',color = 'r')#o-:circle
plt.title("Array Size = 5000")
plt.xlabel("m")
plt.ylabel("runtime")
plt.savefig('./n=5000.png')
plt.show()

plt.plot(x,y5,'o-',color = 'r')#o-:circle
plt.title("Array Size = 6000")
plt.xlabel("m")
plt.ylabel("runtime")
plt.savefig('./n=6000.png')
plt.show()

plt.plot(x,y6,'o-',color = 'r')#o-:circle
plt.title("Array Size = 7000")
plt.xlabel("m")
plt.ylabel("runtime")
plt.savefig('./n=7000.png')
plt.show()

plt.plot(x,y7,'o-',color = 'r')#o-:circle
plt.title("Array Size = 8000")
plt.xlabel("m")
plt.ylabel("runtime")
plt.savefig('./n=8000.png')
plt.show()

plt.plot(x,y8,'o-',color = 'r')#o-:circle
plt.title("Array Size = 9000")
plt.xlabel("m")
plt.ylabel("runtime")
plt.savefig('./n=9000.png')
plt.show()

plt.plot(x,y9,'o-',color = 'r')#o-:circle
plt.title("Array Size = 10000")
plt.xlabel("m")
plt.ylabel("runtime")
plt.savefig('./n=10000.png')
plt.show()     