# reference: https://stackoverflow.com/questions/12545231/tower-of-hanoi-edit-k-peg-solution
# 			 Code Samples in Compass

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('n', type=int)
parser.add_argument('k', type=int)
args = parser.parse_args()

# print("n="+str(args.n))
# print("k="+str(args.k))

import sys
import math
from collections import deque
# See https://www.geeksforgeeks.org/deque-in-python/ for details on Deques
Towers = deque()
free_pegs = deque()

# Global variable that keeps track of the number of steps in our solution 
number_of_steps = 0

# It is always a good idea to set a limit on the depth-of-the-recursion-tree in Python
sys.setrecursionlimit(3000)

def initialize(n,k) :
    for i in range(k) :
        X = deque()
        if (i == 0) :
            for j in range(n) :
                X.append(j+1)
        Towers.append(X)
        
    for i in range(k):
        free_pegs.append(i)

def is_everything_legal() :
    result = True
    for i in range(k) :
        for j in range(len(Towers[i])) :
            for p in range(j,len(Towers[i])) :
                if (Towers[i][p] < Towers[i][j]) :
                    result = False
    return(result)

def move_top_disk(source, dest):
    global number_of_steps 
    number_of_steps = number_of_steps + 1
    x = Towers[source].popleft()
    Towers[dest].appendleft(x)
    if (True == is_everything_legal()) :
        y = " (Legal)"
    else :
        y = " (Illegal)"
    
    print ("Move disk " + str(x) + " from Peg " + str(source+1) + " to Peg " + str(dest+1) + y)

def move(number_of_disks, source, dest, free_pegs, medium=-1):
            
    if(1==number_of_disks):
        move_top_disk(source,dest)
    else:
        if (len(free_pegs)-2)>1:
            p = math.floor(number_of_disks/2.0)
        else:
            p = number_of_disks-1
        
        for peg in free_pegs:
            if(peg!=source and peg!=dest):
                ith = peg
                break        

        move(p,source,ith,free_pegs)
        
        index_of_ith = free_pegs.index(ith)
        free_pegs.remove(ith) 
        move(number_of_disks-p,source,dest,free_pegs)
        
        free_pegs.insert(index_of_ith,ith)
        move(p,ith,dest,free_pegs)
    
def print_peg_state(m) :
    global number_of_steps
    print ("-----------------------------")
    print ("State of Peg " + str(m+1) + " (Top to Bottom): " + str(Towers[m]))
    print ("Number of Steps = " + str(number_of_steps))
    print ("-----------------------------")
    

n = args.n
k = args.k
print(str(k)+"-Tower of Hanoi Problem, with "+str(n)+"-many disks on leftmost peg")
initialize(n,k)
print_peg_state(0)
move(n,0,k-1,free_pegs)
print_peg_state(k-1)
