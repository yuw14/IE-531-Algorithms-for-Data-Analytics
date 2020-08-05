# Reference:
# Code Sample: MCMC-MH.ipynb

#!/usr/bin/env python3
# coding: utf-8

# Gibbs-Sampling procedure to compute the Probability Matrix of a Discrete-Time Markov
# Chain whose states are the d-dimensional cartesian product of the form 
# {0,1,...n-1} x {0,1,...n-1} x ... X {0,1,...n-1} (i.e. d-many products)
# 
# The target stationary distribution is expressed over the n**d many states 
#
# Written by Prof. R.S. Sreenivas for
# IE531: Algorithms for Data Analytics
#

import sys
import argparse
import random
import numpy as np 
import time
import math
import matplotlib.pyplot as plt
import itertools as it

# need this to keep the matrix print statements to 4 decimal places
np.set_printoptions(formatter={'float': lambda x: "{0:0.4f}".format(x)})

# This function computes a random n-dimensional probability vector (whose entries sum to 1)
def generate_a_random_probability_vector(n) :
    a =[]
    for i in range(n-1):
        a.extend([np.random.uniform()])
    a=np.sort(a)
    b=[a[0]]
    for i in range(1,n-1):
        b.extend([a[i]-a[i-1]])
    b.extend([1.0-a[n-2]])
    return b
    # ** WRITE THIS PART **
    # return y

# Two d-tuples x and y are Gibbs-Neighbors if they are identical, 
# (No this situation: or they differ in value at justone coordinate
def check_if_these_states_are_gibbs_neighbors(x, y, dim) :
    # x and y are dim-tuples -- we will assume they are not equal
    # count the number of coordinates where they differ
    number_of_free_coordinates = 0
    for i in range(dim):
        if(x[i]!=y[i]):
            number_of_free_coordinates+=1
        if(number_of_free_coordinates>1):
            return False
    if(number_of_free_coordinates==0 or number_of_free_coordinates>1):
        return False
    else:
        return True

    # ** WRITE THIS PART **
    # return False or True 

# Given two Gibbs-Neighbors -- that are not identical -- find the coordinate where they differ in value
# this is the "free-coordinate" for this pair
def free_coordinates_of_gibbs_neighbors(x, y, dim) :
    # we assume x and y are gibbs neighbors, then the must agree on at least (dim-1)-many coordinates
    # or, they will disagree on just one of the (dim)-many coordinates... we have to figure out which 
    # coordinate/axes is free
    for i in range(dim):
        if(x[i]!=y[i]):
            free_index = i
    # ** WRITE THIS PART **
    return free_index

# x in a dim-tuple (i.e. if dim = 2, it is a 2-tuple; if dim = 4, it is a 4-tuple) state of the Gibbs MC
# each of the dim-many variables in the dim-tuple take on values over range(n)... this function returns 
# the lexicographic_index (i.e. dictionary-index) of the state x
def get_lexicographic_index(x, n, dim) :
   # ** WRITE THIS PART ** 
    number = 0
    coordinates = []
    for i in range(dim):
        coordinates.append(x[i])
    a = dim-1
    for item in coordinates:
        if(a>=0):
            number+=item*pow(n,a)
            a=a-1
    return number

# This is an implementaton of the Gibbs-Sampling procedure
# The MC has n**dim many states; the target stationary distribution is pi
# The third_variable_is when set to True, prints the various items involved in the procedure
# (not advisable to print for large MCs)
def create_gibbs_MC(n, dim, pi, do_want_to_print) :
    if (do_want_to_print) :
        print ("Generating the Probability Matrix using Gibbs-Sampling")
        print ("Target Stationary Distribution:")
        for x in it.product(range(n), repeat = dim) :
            number = get_lexicographic_index(x, n, dim)
            print ("\u03C0", x, " = \u03C0(", number , ") = ", pi[number])
    
    # the probability matrix will be (n**dim) x (n**dim) 
    probability_matrix = [[0 for x in range(n**dim)] for y in range(n**dim)]
    
    # the state of the MC is a dim-tuple (i.e. if dim = 2, it is a 2-tuple; if dim = 4, it is a 4-tuple)
    # got this from https://stackoverflow.com/questions/7186518/function-with-varying-number-of-for-loops-python
    for x in it.product(range(n), repeat = dim) :
        # x is a dim-tuple where each variable ranges over 0,1,...,n-1
        number_of_x = get_lexicographic_index(x, n, dim)
        for y in it.product(range(n), repeat = dim) :
            # change status x to status y
            number_of_y = get_lexicographic_index(y, n, dim)
            if(number_of_x==number_of_y):
                continue
            elif(check_if_these_states_are_gibbs_neighbors(x, y, dim)):
                free_index = free_coordinates_of_gibbs_neighbors(x, y, dim)
                Denominator = 0
                element_in_denominator = list(x)
                for i in range(n):
                    element_in_denominator[free_index]=i
                    tuple_element_in_denominator=tuple(element_in_denominator)
                    number = get_lexicographic_index(tuple_element_in_denominator, n, dim)
                    Denominator += a[number]
                probability_matrix[number_of_x][number_of_y]= (1.0/dim)*a[number_of_y]/Denominator
            else:
                probability_matrix[number_of_x][number_of_y] = 0
    # when identical            
    for x in it.product(range(n), repeat = dim) :
        number_of_x = get_lexicographic_index(x, n, dim)
        for y in it.product(range(n), repeat = dim) :    
            number_of_y = get_lexicographic_index(y, n, dim)
            if(number_of_x==number_of_y):
                p = 1
                for i in range(dim) :
                    for j in range(n):
                        if(j==x[i]):
                            continue
                        new_y = list(x)
                        new_y[i]=j
                        tuple_new_y = tuple(new_y)
                        # the probability from x to new_y
                        number_of_new_y = get_lexicographic_index(tuple_new_y, n, dim)

                        p -= probability_matrix[number_of_x][number_of_new_y]

                        # free_index = i
                        # Denominator = 0
                        # element_in_denominator = x
                        # for k in range(n):
                        #     element_in_denominator[free_index]=k
                        #     number = get_lexicographic_index(element_in_denominator, n, dim)
                        #     Denominator += a[number]
                        # p -= (1.0/dim)*a[number_of_new_y]/Denominator

                probability_matrix[number_of_x][number_of_y] = p
            else:
                continue

    return probability_matrix

# Trial 1... States: {(0,0), (0,1), (1,0), (1,1)} (i.e. 4 states)
n = 2
dim = 2
a = generate_a_random_probability_vector(n**dim)
print("(Random) Target Stationary Distribution\n", a)
p = create_gibbs_MC(n, dim, a, True) 
print ("Probability Matrix:")
print (np.matrix(p))
print ("Does the Probability Matrix have the desired Stationary Distribution?", np.allclose(np.matrix(a), np.matrix(a)* np.matrix(p)))

# Trial 2... States{(0,0), (0,1),.. (0,9), (1,0), (1,1), ... (9.9)} (i.e. 100 states)
n = 10
dim = 2
a = generate_a_random_probability_vector(n**dim)
p = create_gibbs_MC(n, dim, a, False) 
print ("Does the Probability Matrix have the desired Stationary Distribution?", np.allclose(np.matrix(a), np.matrix(a)* np.matrix(p)))

# Trial 3... 1000 states 
n = 10
dim = 3
t1 = time.time()
a = generate_a_random_probability_vector(n**dim)
p = create_gibbs_MC(n, dim, a, False) 
t2 = time.time()
hours, rem = divmod(t2-t1, 3600)
minutes, seconds = divmod(rem, 60)
print ("It took ", hours, "hours, ", minutes, "minutes, ", seconds, "seconds to finish this task")
print ("Does the Probability Matrix have the desired Stationary Distribution?", np.allclose(np.matrix(a), np.matrix(a)* np.matrix(p)))

# Trial 4... 10000 states 
n = 10
dim = 4
t1 = time.time()
a = generate_a_random_probability_vector(n**dim)
p = create_gibbs_MC(n, dim, a, False) 
t2 = time.time()
hours, rem = divmod(t2-t1, 3600)
minutes, seconds = divmod(rem, 60)
print ("It took ", hours, "hours, ", minutes, "minutes, ", seconds, "seconds to finish this task")
print ("Does the Probability Matrix have the desired Stationary Distribution?", np.allclose(np.matrix(a), np.matrix(a)* np.matrix(p)))

