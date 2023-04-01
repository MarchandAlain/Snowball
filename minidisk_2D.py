# Minimum disk algorithm in 2-D - Alain R. Marchand - April 2023
import os
import math
import numpy as np
from itertools import combinations
import time
import tkinter

square, disk, Gaussian, ellipse, annulus = "square", "disk", "Gaussian", "ellipse", "annulus"
# Parameters
Size_params= 2_000_000, 2_000_000, 1.05                                  # min, max, geometric progression
Dimension=2
Shape=Gaussian                                                              # square   disk   Gaussian   ellipse   annulus
Repeat=30                                                                     # create a new set each time 
Near_zero=1.e-6                                                         # for determinant
Approximate=1+1.e-12                                               # slightly more than 1
Step=1                                                                         # only 1 or a prime number, e.g. 911 or a prime close to 5N/28
Draw=False

# Debugging options
Verify=False
Tests=False                                                                  # if True, perform unit tests
Test_set = None
if Test_set: Size_list, Repeat = [len(Test_set)], 1
Print_points=False
Max_steps=1000

# drawing
point_color='red'
line_color='black'
Scale=40
offsetX=200
offsetY=100
fenetre=0
canvas=0

##====================================== create_set
def create_set(size):
    """
    Draw 2D points from a Gaussian distribution
    Parameters: size, integer
    Returns: points, a list of tuples of 2 float
    """
# uniform square
    if Shape==square:
        x=list(np.random.random(size)*10) 
        y=list(np.random.random(size)*10)         

# uniform disk
    elif Shape==disk:
        x=list(np.random.random(2*size)*10) 
        y=list(np.random.random(2*size)*10)                  # start with a square
        xy=[(a-5,b-5) for a,b in zip(x,y) if (a-5)**2+(b-5)**2<=25]
        x, y = list(zip(*xy[:size]))                                        

# Gaussian disk
    elif Shape==Gaussian:
        Mean, Sigma=0, 5
        angle=list(np.random.random(size)*2*math.pi)
        radius=list(np.random.normal(Mean, Sigma, size))
        x=[r*math.cos(a) for r, a in zip(radius, angle)]
        y=[r*math.sin(a) for r, a in zip(radius, angle)]

# Gaussian ellipse
    elif Shape==ellipse:
        Mean, Sigma, Elongation=0, 5, 3
        angle=list(np.random.random(size)*2*math.pi)
        radius=list(np.random.normal(Mean, Sigma, size))
        x=[r*math.cos(a)*Elongation for r, a in zip(radius, angle)]
        y=[r*math.sin(a) for r, a in zip(radius, angle)]

# Gaussian annulus
    elif Shape==annulus:
        Mean, Sigma=8, 1
        angle=list(np.random.random(size)*2*math.pi)
        radius=list(np.random.normal(Mean, Sigma, size))
        x=[r*math.cos(a) for r, a in zip(radius, angle)]
        y=[r*math.sin(a) for r, a in zip(radius, angle)]
    
    points=list(zip(x,y))
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

#====================================== initFenetre
def initFenetre():
    global fenetre, canvas
    fenetre=tkinter.Tk()
    fenetre.title('Minimum disk')
    canvas=tkinter.Canvas(fenetre, width=800, height=600, bg='white')
    canvas.pack()  

#====================================== to_scale
def to_scale(x, y):
    return offsetX+x*Scale, offsetY+y*Scale

#====================================== draw_points
def draw_points(points, size=2, color=point_color):
    for (u,v) in points:
        x,y=to_scale(u,v)
        canvas.create_line(x,y, x+size, y+size, fill=color, width=size)

#====================================== draw_circle
def draw_circle(O, r_squared, color=line_color):
    x, y= to_scale(*O)
    r=math.sqrt(r_squared)*Scale
    canvas.create_oval(x-r, y-r, x+r, y+r)
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
    return d_squared(A, O)<=r_squared*Approximate                  # allow for rounding errors                        

#====================================== point_outside
def point_outside(sample, position, O, r_squared):
    """
    look for a point outside circle
    cyclic search in sample
    stop if found
    Parameters: sample, a list of points, tuples of float
                  position, an integer modulo len(sample)
                  O, a point, center of the circle
                  r_squared, a float, square of radius
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
    return O, d_squared(A, O)

#====================================== circum_ball
def circum_ball(subset):
    """
    find center and radius of circumscribed ball
    assumes points are NOT in a lower D (verify elsewhere)
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
    for Dimension>2, matrix is underdetermined with only 3 points
    center must be in the plane of the 3 points: what equation ?
    """
    if len(subset)>Dimension+1:
        print("*** Error: wrong subset size ***")
        print(subset)
        stop()

    if len(subset)==2: return ball(*subset)              # only two points

    # build matrix
    coeff, const = [], []
    for P in subset:
        a, c = [], 0
        for i in range(Dimension):
            a+=[2*P[i]]
            c+=P[i]*P[i]
        a+=[1]
        const+=[c]
        coeff+=[a]

    # solve matrix
    u, v = np.array(coeff), np.array(const)
    if abs(np.linalg.det(u))<Near_zero:
        print("*** Warning: low accuracy ***")
    solved = np.linalg.solve(u, v)
    O=tuple(solved[:Dimension])
    return O, d_squared(O, subset[0])

#====================================== find_ball
def find_ball(subset, support_size):
    """
    find smallest ball based on support_size points, enclosing remaining points
    explore all combinations of points
    Parameters: subset, a set of points
    Returns: O, a point, tuple of float, None if not found
                  r_squared, a float
    """
    best=None, float('inf'), []
    if len(subset)==2: return (*ball(*subset), subset)
    if support_size>len(subset): support_size=len(subset)

    for miniset in combinations(subset, support_size):
        # compute circle
        O, r_squared = circum_ball(miniset) if support_size>2 else ball(*miniset)
        remain=list(set(subset).difference(set(miniset)))
        if not remain:
            if r_squared<best[1] or support_size==Dimension:
                best = (O, r_squared, list(miniset))
        else:
            # does it enclose all subset
            D, position= point_outside(remain, 0, O, r_squared)
            # is it the smallest
            if not D and r_squared<best[1]:
                best = (O, r_squared, list(miniset))
    if best[0]:
        return best        
    return None, 0, []        

#====================================== minimum_ball
def minimum_ball (points, subset, position):
    """
    Parameters: points, a list of points, tuples of float
                  subset, a list of Dimension+2 points to start with
                  position, in points list
    Returns: subset, a list of points
    O, a point, tuple of float, center of the circumscribed circle
                  r_squared, a float
    """
    global steps
    while True:                                                                     # loop until no more points
        steps+=1                                                                   # count steps
        if steps>Max_steps:
            print("Too many steps", O, r_squared)
            print(subset)
            print([d_squared(O,A) for A in subset])
            stop()

        # look for best circle           
        support_size=2
        while support_size<Dimension+2:
            # find_ball
            O, r_squared, base = find_ball(subset, support_size)
            if O: break
            support_size+=1

        if not O:
            print("\n*** Anomaly: No solution found for subset ***")
            print(subset)
            stop()
             
        # add point
        subset=base
        P, position = point_outside(points, position, O, r_squared)
        if not P: return subset, O, r_squared                        # all points are inside: return
        
        subset+=[P]

#====================================== unit_tests
def unit_tests():
    """
    test each module
    """
    global steps, distances, position
    if Dimension!=2: return
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

    assert(point_outside(sample, 0, (5,5), 25)==((9,9), 0))  # point_outside & inc_position
    assert(point_outside(sample, 3, (5,5), 32)==(None, 3))

    assert(ball((1,2), (1,4))==((1,3), 1))                    # ball

    assert(circum_ball(sample[:3])==((-0.5, 4.5), 18.5))     # circum_ball
    assert(circum_ball([(0,3), (8,3), (4,-5)])==((4,0), 25))
    assert(circum_ball([(3,0), (3,8), (-5,4)])==((7.105427357601002e-16, 4.0), 24.999999999999993))

    assert(find_ball(sample[1:], 2)==((5.5, 5.0), 28.25, [(2, 1), (9, 9)]))      # find_ball (4 pts)
    assert(find_ball(sample[:3], 2)==((2.5, 4.0), 9.25, [(3, 7), (2, 1)]))        # find_ball (3 pts)

    subset, O, r_squared = minimum_ball([(0,3), (8,3), (4,-5), (1,1), (6,2)], [(0,3), (8,3), (4,-5), (6,2)], 0)     # minimum_ball (4 pts) 
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

if Draw: initFenetre()

for index, Set_size in enumerate(Size_list):   
    for rep in range(Repeat):

        # performance stats
        steps=0                                                              # number of subsets examined
        distances=0                                                        # number of distances measured
        
        # initialize set of points
        points=create_set(Set_size)
        if Test_set: points=Test_set
        subset=points[:2]

        # display
        if Draw: draw_points(points)
        
        # start time
        t0=time.perf_counter()

        # compute minimum ball
        position=2
        subset, O, r_squared=minimum_ball(points, subset, position)

        # end time
        t1=time.perf_counter()

        # draw
        if Draw: draw_circle(O, r_squared)

        # print stats
        print(Set_size, f"\t{distances/Set_size:3.3}\t{steps}\t{(t1-t0)*1000:9.2f}")

        # just a verification
        if Verify:
            if point_outside(points, 0, O, r_squared)[0]:                     
                print("\n*** Not solved ! ***\n")
                stop()

stop()
