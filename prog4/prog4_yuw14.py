#!/usr/bin/env python3
# coding: utf-8

# Three-sided Unfair Dice into other RV distributions
#
# IE531: Algorithms for Data Analytics
# Written by Prof. R.S. Sreenivas
#
import sys
import argparse
import random
import numpy as np 
import time
import math
import matplotlib.pyplot as plt
from scipy.stats import uniform

# Some constants -- "heads" is 1 and "tails" is 0
HEADS = 1
TAILS = 0

# reading the number of trials on command line
#no_of_trials = 1000000
no_of_trials = int(sys.argv[1])+0

# assign a random value for the probability of seeing a "1" (p1), "2" (p2) or "3" (p3) for 
# the 3-sided unfair dice such that (p1, p2, p3) is uniformly distributed over the surface
# p1+p2+p3 = 1, where p1, p2, p3 are non-negative real numbers
def assign_probabilities_to_unfair_three_sided_dice() :
    # FILL CODE HERE
    # just use np.random.uniform(), go ahead!

    # this is p1
    prob_of_one = np.random.uniform() 
    # this is p2
    prob_of_two = np.random.uniform()
    # make sure p1+p2<=1, if not, select a new p2
    while((prob_of_one+prob_of_two)>1):
    	prob_of_two = np.random.uniform()
    # since p1+p2+p3=1, so p3=1-p1-p2
    prob_of_three = 1.0-prob_of_one-prob_of_two
    return prob_of_one, prob_of_two, prob_of_three

# This function simulates a single toss of the unfair 3-sided dice
def toss_of_three_sided_unfair_dice(p1, p2, p3) :
    x = np.random.uniform()
    if (x < p1) :
        return 1
    else :
    	# I modified your code here 
        if (x < p1+p2) :
            return 2
        else :
            return 3

# This function simulates a fair-coin using the unfair 3-sided dice
def simulate_fair_coin_from_unfair_three_sided_dice(p1, p2, p3) :
    # FILL Code here
    # figure out how to convert the outcome of the unfair 3-sided dice
    # into an outcome from a Fair Coin

    # toss the unfair dice 3 times
    toss1 = toss_of_three_sided_unfair_dice(p1, p2, p3)
    toss2 = toss_of_three_sided_unfair_dice(p1, p2, p3)
    toss3 = toss_of_three_sided_unfair_dice(p1, p2, p3)

    # repeat the unfair dice toss unitil they are p1, p2, p3 (no order)
    while not ( (toss1==1 and toss2==2 and toss3==3) or (toss1==1 and toss2==3 and toss3==2) or (toss1==2 and toss2==1 and toss3==3) or (toss1==2 and toss2==3 and toss3==1) or (toss1==3 and toss2==1 and toss3==2) or (toss1==3 and toss2==2 and toss3==1)):
	    toss1 = toss_of_three_sided_unfair_dice(p1, p2, p3)
	    toss2 = toss_of_three_sided_unfair_dice(p1, p2, p3)
	    toss3 = toss_of_three_sided_unfair_dice(p1, p2, p3)
        
    # if (toss1 is 1 and toss2 is 2 and toss3 is 3) or (toss1 is 1 and toss2 is 3 and toss3 is 2) or (toss1 is 2 and toss2 is 1 and toss3 is 3) we output (Fair) Heads
    # if (toss1 is 2 and toss2 is 3 and toss3 is 1) or (toss1 is 3 and toss2 is 1 and toss3 is 2) or (toss1 is 3 and toss2 is 2 and toss3 is 1) we output (Fair) Tails
    if ( (toss1==1 and toss2==2 and toss3==3) or (toss1==1 and toss2==3 and toss3==2) or (toss1==2 and toss2==1 and toss3==3) ) :
        return HEADS
    else :
        return TAILS
    # return HEADS/TAILS

# get a U.I.I.D RV by making the unfair three sided dice into a fair coin... and tossing the 
# resulting fair-coin 32 times to get discrete RV that is uniformly distributed over the 
# integers in [0, 2^{32}-1]... dividing the resulting integer by 2^{32}-1 gives us (effectively)
# a U.I.I.D. RV
def get_uiid_rvs_by_tossing_the_unfair_three_sided_dice_32_times(p1, p2, p3) :
    result = 0
    for i in range(0, 32) :
        if (simulate_fair_coin_from_unfair_three_sided_dice(p1, p2, p3) == HEADS) :
            result = result | (1 << i)
        else :
            result = result | (0 << i)
    return float(result/(pow(2,32)-1))

# plotting the histogram of the continuous RV generated by tossing the unfair three sided dice
# sufficient number of times till we get 32 fair-coin-tosses, which are then converted into a 
# number in the unit-interval

# assigning probabilities to unfair three sided dice 
p1, p2, p3 = assign_probabilities_to_unfair_three_sided_dice()

z = []
for i in range(0,no_of_trials) :
    z.extend([get_uiid_rvs_by_tossing_the_unfair_three_sided_dice_32_times(p1, p2, p3)])
plt.hist(z, bins=50)
plt.ylabel('Histogram for ' + str(no_of_trials) + ' trials');
plt.savefig("UIID_Histogram.pdf", bbox_inches='tight')

# converting (multiple) tosses of the unfair 3-sided Dice into a unit-normal distribution
# using the Box-Muller Method
a = []
for i in range(0,no_of_trials) :
    p = get_uiid_rvs_by_tossing_the_unfair_three_sided_dice_32_times(p1, p2, p3)
    q = get_uiid_rvs_by_tossing_the_unfair_three_sided_dice_32_times(p1, p2, p3)
    theta = 2*math.pi*p 
    r = np.sqrt(-2*math.log(q))
    a.extend([r*math.cos(theta)])
    a.extend([r*math.sin(theta)])
plt.hist(a, bins=50)
plt.ylabel('Histogram for ' + str(2*no_of_trials) + ' trials');
plt.savefig("Unit_Normal_Histogram.pdf", bbox_inches='tight')
