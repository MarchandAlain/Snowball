# Minimum disk algorithm in N-D - Alain R. Marchand - March 2023
import os
import math
import numpy as np
import time

# Parameters
Size_params= 100, 20_000, 1.05                                  # min, max, geometric progression
Dimension=2
Shape='square'                                                            # square   disk   Gaussian   ellipse   annulus
Repeat=1                                                                     # create a new set each time 
Precision=1+1.e-12
Step=1                                                                         # only 1 or a prime number, e.g. 911

# Debugging options
Verify=True
Tests=False                                                                  # if True, perform unit tests
Test_set = None
if Test_set: Size_list, Repeat = [len(Test_set)], 1
Print_points=False

# Combinations
C_3_2=[(0,1,2), (2,0,1), (1,2,0)]                                                         # first 2 out of 3
C_3_2_plus=[(0,1,2,3), (2,0,1,3), (1,2,0,3)]                                         # first 2 out of 3 (ignore last)
C_4_2=[(0,1,2,3), (0,2,1,3), (0,3,1,2), (1,2,0,3), (1,3,0,2), (2,3,0,1)]     # first 2 out of 4
C_4_3=[(0,1,2,3), (3,0,1,2), (2,3,0,1), (1,2,3,0)]                                  # first 3 out of 4

#====================================== create_set
def create_set(size):
    """
    Draw 2D points from a Gaussian distribution
    Parameters: size, integer
    Returns: points, a list of tuples of 2 float
    """
# uniform square
    shape=(size, Dimension)
    array=np.zeros(shape, dtype=float)
    if Shape=='square':
        for d in range(Dimension):
            array[:,d]=list(np.random.random(size)*10) 

    points=[]
    for i in range(size):
        points+=[tuple(array[i,:])]

    if Print_points: print_points(points)

    return points

#====================================== stop
def stop():
    input("\nPress Enter to exit ")
    os._exit(1)

#====================================== print_points
def print_points(points):
    """
    Print points and stop
    """
    for P in points:
        for i in range(Dimension):
            print(P[i], end="\t")
        print()
    stop()

#====================================== d_squared
def d_squared(A, B):
    """
    return squared distance between points A and B
    increment global count distances
    """
    global distances
    distances+=1                                                                       # counts as 1 operation
    d=0
    for i in range(Dimension):
        d+=(A[i]-B[i])**2
    return d
    
#====================================== inc_position
def inc_position(n, position, step, size):
    """
    add n steps to position, modulo size
    """                 
    return (position+n*step) % size  

#====================================== inside
def inside(A, O, r_squared):
    """
    check whether point A is in circle (O, r_squared)
    Returns True if inside, False if outside
    """                 
    return d_squared(A, O)<=r_squared*Precision                  # allow for rounding errors                        

#====================================== obtuse
def obtuse(subset):
    """
    test whether a triangle has an obtuse angle 
    Parameters: subset, a list of 3 points, tuples of 2 float
    Returns: True if an obtuse angle is found, False otherwise
    """
    for i,j,k in C_3_2:
        (x1, y1), (x2, y2), (x3, y3) = subset[i], subset[j], subset[k]
        if (x2-x1)*(x3-x1) + (y2-y1)*(y3-y1) < 0: return True    # negative scalar product
    return False

#====================================== point_outside
def point_outside(sample, position, O, r_squared):
    """
    look for a point outside circle
    cyclic search in sample
    stop if found
    Parameters: sample, a list of 3 points, tuples of 2 float
                  position, an integer between 0 and len(sample)-1
                  O, a point, center of the circle
                  r_squared, a float, square of circle radius
    Returns: E, point if found, else None
                  new position, an integer
    """
    start=position
    while True:
        E=sample[position]
        position=inc_position(1, position, Step, len(sample))
        if not inside(E, O, r_squared):
            return E, position                                   # a point outside
        if position==start:
            return None, position                            # no point outside

#====================================== ball
def ball(A, B):
    """
    calculate ball based on 2 points
    Parameters: A, B, two points to use as a diameter
    Returns: O, a point, tuple of float
                  r_squared, a float
    """
    # calculate ball
    O=[0]*Dimension
    for i in range(Dimension):
        O[i] = (A[i]+B[i])/2                                              # middle of AB
    O=tuple(O)
    r_squared=d_squared(A, O)
    return O, r_squared

#====================================== circum_ball
def circum_ball(subset):
    """
    find center and radius of circumscribed ball
    Parameters: subset, a set of points
    Returns: O, a point, tuple of float
                  r_squared, a float
    For each point i, r² = (xo-xi)²+(yo-yi)²+... with r, xo, yo... unknown
    u = r²-xo²-yo²-...
    xi²+yi²+...-2xixo-2yiyo-...-u = 0
    solve  using linear algebra
    test (0,3), (8,3), (4,-5) --> (4,0), 25
    coeff=[[2*0, 2*3, 1], [2*8, 2*3, 1], [2*4, 2*-5, 1]]
    const=[9, 73, 41]
    """
    if len(subset)>Dimension+1 or len(subset)<3:
        print("*** Error: wrong subset size ***")
        print(subset)
        stop()

    # build matrix
##    subset=[(0,3), (8,3), (4,-5)]
    coeff, const = [], []
    for P in subset:
        a, c = [], 0
        for i in range(Dimension):
            a+=[2*P[i]]
            c+=P[i]*P[i]
        a+=[1]
        const+=[c]
        coeff+=[a]
    
    a=np.array(coeff)
    b=np.array(const)
    x = np.linalg.solve(a, b)
##    print(x)
##    stop()

    O=x[:Dimension]
    r_squared=d_squared(O, subset[0])
    return O, r_squared

#====================================== circle2
def circle2(A, B, remain):
    """
    calculate circle based on 2 points, looking for one that encloses other points
    Parameters: A, B, two points to use as a diameter
                  remain: a set of remaining points to test for inclusion in the circle
    Returns: O, a point, tuple of 2 float
                  r_squared, a float if all points are in the ball, None otherwise
                  [A, B], a list of 2 points forming a diameter of the circle
    """
    # calculate circle
    O, r_squared=ball(A, B)
    
    # test remaining point(s)
    if not point_outside(remain, 0, O, r_squared)[0]:
        return O, r_squared                                         # all points are in the circle 
        
    return None, None         
  
#====================================== circle3
def circle3(subset):
    """
    find the circumscribed circle of a subset of 3 points 
    Parameters: subset, a list of 3 points, tuples of 2 float
    Returns: (x0, y0), a point, tuple of 2 float
                  r_squared, a float, square of circle radius
    """
    (x1, y1), (x2, y2), (x3, y3) = subset
    
    if (x3-x2)*(y2-y1) == (x2-x1)*(y3-y2) :                                # check for alignment              
        return min_diameter(subset)[:2]                                      # print(subset, "are aligned")

    return circum_ball(subset)

#====================================== min_diameter
def min_diameter(subset):
    """
    try circles based on 2 points, looking for one that encloses all 3 or 4 points of the subset
    Parameter: subset, a list of 4 (or 3) points, tuples of 2 float
    Returns: O, a point, tuple of 2 float, or None if not all insode
                  r_squared, a float
                  [A, B], a list of 2 points forming a diameter of the circle
    """
    # generate combinations of two points
    if len(subset)==4:                                                       # 4 points
        combine=C_4_2
    else:
        combine=C_3_2_plus                                              # 3 points
    
    for i, j, k, l in combine:
        A, B = subset[i], subset[j]                                        # 3 or 4 points 
        remain=[subset[k]]
        if len(subset)==4:                                                    # 4 points
            remain+=[subset[l]]

        # calculate circle
        O, r_squared = circle2(A, B, remain)
        if O:
            return O, r_squared, [A, B]                              # all points are in the circle 
            
    return None, None, [None, None]                          # no solution with 2 points           
  
#====================================== min_circumscribed
def min_circumscribed(subset):
    """
    try circles based on 3 points, find the smallest one that encloses all 4 points of the subset
    Parameters: subset, a list of 4 points, tuples of 2 float
    Returns: (x0, y0), a point, tuple of 2 float
                  r_squared, a float
                  a list of 3 points defining the minimum circle
    """
    best = [None, None], float('inf'), subset[:3]                              # initialize best

    # solve if only 3 points
    if len(subset)==3:
        O, r_squared=circle3(subset) if not obtuse(subset) else circle2(subset)
        if  r_squared<best[1]:
            best=O, r_squared, subset                                               # store best

    # solve with 4 points
    else:
        # generate combinations of three points
        for i, j, k, l in C_4_3:
            combine = [subset[i], subset[j], subset[k]]
            remain = subset[l]

            # avoid obtuse triangles
            if obtuse(combine):
                continue

            # compute circle
            O, r_squared=circle3(combine)
            if not r_squared:                                                                 # aligned points
                continue
            
            # test last point in combination
            if inside(remain, O, r_squared):
                if  r_squared<best[1]:
                    best=O, r_squared, combine                               # store best

    if best[1]!=float('inf'):
        return best                                                                       # give result

    print("\n*** Anomaly: No solution found for subset ***")
    print(subset, "\n")
    os._exit(1)      

#====================================== minimum_disk
def minimum_disk (points, subset, position):
    """
    Parameters: points, a list of points, tuples of 2 float
                  subset, a list of four points to start with
                  position, in points list
    Returns: subset, a list of two or three points
    O, a point, tuple of 2 float, center of the circumscribed circle
                  r_squared, a float
    """
    global steps
    while True:
        steps+=1
        
        # try solving with 2 points
        P, s, [A, B] = min_diameter(subset)                                  
        if P:
            O, r_squared = P, s
            subset = [A, B]

        # solve with 3 points            
        else:
            O, r_squared, [A, B, C] = min_circumscribed(subset)                
            subset = [A, B, C]

        # try adding one outside point
        D, position = point_outside(points, position, O, r_squared)
        if not D:
            return subset, O, r_squared                                     # all points are inside: return
        else:
            subset += [D]

        # try completing to four with another outside point
        if len(subset) == 3:
            D, position = point_outside(points, position, O, r_squared)
            if D:
                subset += [D]

#====================================== unit_tests
def unit_tests():
    """
    test each module
    """
    global steps, distances, position
    steps=0                                                              # number of subsets examined
    distances=0                                                        # number of distances measured
    position=0                                                          # index to points
    sample=[(3,2), (3,7), (2,1), (8,2), (9,9)]
    
    assert(d_squared((0,0), (3,4))==25)                   # d_squared
    assert(distances==1)                                          # distances
    assert(inc_position(2, position, 7, 11) == 3)       # inc_position
    
    assert(inside((0,0), (3,4), 25))                              # inside
    assert(inside((3,4), (0,0), 26))                             
    assert(not inside((0,3), (4,0), 24))

    assert(obtuse(sample[:3])==True)                        # obtuse
    assert(obtuse([(0,3), (8,3), (4,-5)])==False)                        

    assert(point_outside(sample, 0, (5,5), 25)==((9,9), 0))  # point_outside & inc_position
    assert(point_outside(sample, 3, (5,5), 32)==(None, 3))

    assert(circle2((1,2), (1,4), [(1,3)])==((1,3), 1))  # circle2
    assert(circle2((1,2), (1,3), [(1,4)])==(None, None)) 

    assert(circle3(sample[:3])==((-0.5, 4.5), 18.5))     # circle3
    assert(circle3([(1,2), (1,3), (1,4)])==((1,3), 1))
    assert(circle3([(0,3), (8,3), (4,-5)])==((4,0), 25))
    assert(circle3([(3,0), (3,8), (-5,4)])==((0,4), 25))

    assert(min_diameter(sample[1:])==((5.5, 5.0), 28.25, [(2, 1), (9, 9)]))   # min_diameter (4 pts)
    assert(min_diameter(sample[:3])==((2.5, 4.0), 9.25, [(3, 7), (2, 1)]))     # min_diameter (3 pts)
    assert(min_diameter(sample[1:4])==(None, None, [None, None]))     # min_diameter (3 pts)

    assert(minimal3([(0,3), (8,3), (4,-5)])==((4.0, 0.0), 25.0, [(0, 3), (8, 3), (4, -5)])) # minimal3 (3 pts)
    assert(minimal3([(0,3), (8,3), (4,-5), (6,2)])==((4.0, 0.0), 25.0, [(0, 3), (8, 3), (4, -5)])) # minimal3 (4 pts)

    subset, O, r_squared = minimum_disk([(0,3), (8,3), (4,-5), (1,1), (6,2)], [(0,3), (8,3), (4,-5), (6,2)], 0)     # minimum_disk (4 pts) 
    assert((set(subset), O,r_squared)==(set([(0,3), (8,3), (4,-5)]), (4.0, 0.0), 25.0))                              # order is unreliable        

    print("\nTests OK")    
    os._exit(1)
    
#====================================== main program
"""
Draw sample of 2-D points and determine their minimal enclosing circle
evaluate performance for each sample
"""
if Tests: unit_tests()

if Step>1:
    print("Warning: with a Step>1 you can avoid randomization, but it is your responsibility")
    print("to verify that Step is a prime number in order to avoid avoid subsampling")
    print("Examples of prime numbers: 53, 97, 179, 271, 367, 557, 719, 911, 1777")
    print("Your Step:", Step, "\n")

# compute increasing size list
if not Test_set:
    size, max_size, coef=Size_params
    Size_list=[]
    while size<=max_size:
        Size_list+=[size]
        size=int(size*coef)
    print(Shape, "\nSet_size\tDist./N\tSteps\tTime(ms)")
else:
    Size_list=[len(Test_set)]

for index, Set_size in enumerate(Size_list):   
    for rep in range(Repeat):

        # performance stats
        steps=0                                                              # number of subsets examined
        distances=0                                                        # number of distances measured
        
        # initialize set of points
        points=create_set(Set_size)
        if Test_set: points=Test_set
        subset=points[:4]

        # start time
        t0=time.perf_counter()

        # compute minimum circle
        position=0
        subset, O, r_squared=minimum_disk(points, subset, position)

        # end time
        t1=time.perf_counter()

        # print stats
        print(Set_size, f"\t{distances/Set_size:3.3}\t{steps}\t{(t1-t0)*1000:9.2f}")

        # just a verification
        if Verify:
            if point_outside(points, 0, O, r_squared)[0]:                     
                print("\n*** Not solved ! ***\n")
                os._exit(1)

stop()

