# Maximum Lk distance algorithm in N-D - Alain R. Marchand - April 2023
import os
import math
import numpy as np
import time
import tkinter

square, Gaussian, disk, annulus, ring, worst = 'square', 'Gaussian', 'disk', 'annulus', 'ring', 'worst'

# Parameters
Size_params= 200, 20_000, 1.05                                # min, max, geometric progression
Repeat=1                                                                   # create a new set each time 
Dimension=2
Shape=square
Width=.10                                                                    # for annulus and ring   
Norm=2                                                                       # fractional norms Lk
Precision=1+1.e-12
Step=1                                                                          # only 1 or a relative prime to N (e.g. 577)
Integers=False                                                              # Integers (with square) favor duplicate points
Preprocess=True                                                          # first use strides to get elongated shape

# Debugging options
Verify=False
Test=False                                                                   # if True, perform unit tests
Test_set = None                                                          # list of test points
if Test_set: Size_list, Repeat = [len(Test_set)], 1
Print_points=False

# drawing
Draw=False
if Dimension!=2: Draw=False
point_color='black'
Scale=500
offsetX=400
offsetY=300
container=0
canvas=0

#====================================== create_set
def create_set(size):
    """
    Draw N-D points from a distribution
    Parameters: size, integer
    Returns: points, a list of tuples of N float
    """
    if Shape==square:                                              # square, uniform density
        shape=(size, Dimension)
        array=np.zeros(shape, dtype=float)
        for d in range(Dimension):
            array[:,d]=list(np.random.random(size)-0.5)
            
        if Integers: array = array.astype(int)

        points=[]
        for i in range(size):
            points+=[tuple(array[i,:])]
            
    elif Shape==disk:                                              # ball, uniform density
        shape=(size*Dimension, Dimension)
        array=np.zeros(shape, dtype=float)
        for d in range(Dimension):
            array[:,d]=list(np.random.random(Dimension*size)-0.5)
        points=[]
        for i in range(2*size):
            if len(points)>=size: break
            r_squared=0
            for d in range(Dimension):
                r_squared+=array[:,d][i]**2
            if r_squared<=0.25:
                P=[]
                for d in range(Dimension):
                    P+=[array[:,d][i]]
                points+=[tuple(P)]

    elif Dimension==2:                  
        if Shape==Gaussian:                                     # Gaussian disk
            angle=list(np.random.random(size)*2*math.pi)
            radius=list(np.random.normal(0, 0.15, size))
            x=[r*math.cos(a) for r, a in zip(radius, angle)]
            y=[r*math.sin(a) for r, a in zip(radius, angle)]
            points=list(zip(x,y))

        elif Shape==annulus:                                    # Gaussian annulus
            angle=list(np.random.random(size)*2*math.pi)
            radius=list(np.random.normal(0.2, 0.2*Width, size))
            x=[r*math.cos(a) for r, a in zip(radius, angle)]
            y=[r*math.sin(a) for r, a in zip(radius, angle)]  
            points=list(zip(x,y))
            
        elif Shape==ring:                                         # ring, uniform radius
            angle=list(np.random.random(size)*2*math.pi)
            radius=list(np.random.random(size)*Width+0.4-Width/2)
            x=[r*math.cos(a) for r, a in zip(radius, angle)]
            y=[r*math.sin(a) for r, a in zip(radius, angle)]  
            points=list(zip(x,y))
        
        elif Shape==worst:
            points=worst_case(size)

    if Print_points: print_points(points)

    return points

#====================================== worst_case
def worst_case(size):
    """
    Draw 2-D points in the worst case configuration
    close to equilateral triangle
    centered on (5,5) furthest points at (4,5), (6,5)
    Parameters: size, integer
    Returns: points, a list of tuples of 2 float
    """
    dx=0.0500
    radius=0.0200
    v=math.sqrt(3)/2
    centers=[(5-v-dx, 4.5), (5, 6.05), (5+v+dx, 4.5)]
    points=[]
    for i in range(3):
        count=0
        while count<size/3 and len(points)<size-2:
            x, y = centers[i]
            x+=np.random.random(1)[0]*2*radius
            y+=np.random.random(1)[0]*2*radius
            ok=((x-centers[i][0])**2+(y-centers[i][1])**2) < radius*radius
            if ok:
                points+=[(x, y)]
                count+=1
    points+=[(4,5), (6,5)]                                # most distant points
    for i, (x, y) in enumerate(points):              # rescale
        points[i]=(x/2-2.5, y/2-2.65)
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
    for x,y in points:
        print(x, "\t", y)
    stop()

#====================================== initcontainer
def initContainer():
    """
    prepare window to draw 2-D points
    """
    global container, canvas
    container=tkinter.Tk()
    container.title('Furthest points')
    canvas=tkinter.Canvas(container, width=800, height=600, bg='white')
    canvas.pack()  

#====================================== to_scale
def to_scale(x, y):
    """
    adjust scale to window
    """
    return offsetX+x*Scale, offsetY+y*Scale

#====================================== draw_points
def draw_points(points, size=2, color=point_color):
    for (u,v) in points:
        x,y=to_scale(u,v)
        canvas.create_line(x,y, x+size, y+size, fill=color, width=size)

#====================================== distance
def distance(A, B):
    """
    return some distance between points A and B
    increment global count distances
    Norm is a fractional distance parameter
    2: Euclidean ; 1: Manhattan ; any number: fractional distance
    """
    global distances
    distances+=1                                                                       # counts as 1 operation
    d=0
    for i in range(Dimension):
        d+=(abs(A[i]-B[i]))**Norm
    return d
       
#====================================== inc_position
def inc_position(sample, position, step, size):
    """
    add one step to position, modulo size
    skip null points that have been eliminated
    """
    position=(position+step) % size
    while not sample[position]:
##        print(position, "is null")
        position=(position+step) % size
    return position

#====================================== inside
def inside(A, O, radius):
    """
    check whether point A is in ball (O, radius)
    Returns True if inside, False if outside
    """                 
    return distance(A, O)<=radius*Precision                  # allow for rounding errors                               

#====================================== point_outside
def point_outside(sample, position, step, size, O, radius):
    """
    look for a point outside ball
    cyclic search in sample
    stop if found
    Parameters: sample, a list of N points, tuples of Dimension float
                  position, an integer between 0 and size-1
                  size: total number of points, including eliminated
                  step: a parameter to browse set
                  O, a point, center of the ball
                  radius, a float, Norm power of ball radius
    Returns: E, point if found, else None
                  new position, an integer
    """
    start=position
    while True:
        E=sample[position]
        if not inside(E, O, radius):
            return E, position                                   # a point outside
        position=inc_position(sample, position, step, size)
        if position==start:
            return None, position                            # no point outside
        
#====================================== ball
def ball(A, B):
    """
    calculate ball based on 2 points
    Parameters: A, B, two points to use as a diameter
    Returns: O, a point, tuple of Dimension float
                  radius, a float, Norm power of ball radius
    """
    # calculate ball
    O=[0]*Dimension
    for i in range(Dimension):
        O[i] = (A[i]+B[i])/2                                              # middle of AB
    O=tuple(O)
    radius=distance(A, O)
    return O, radius

#====================================== unit_tests
def unit_tests():
    """
    test each module
    """
    global Dimension, steps, distances, position
    Dimension=2
    steps=0                                                              # number of subsets examined
    distances=0                                                        # number of distances measured
    position=0                                                          # index to points
    sample=[(3,2), (3,7), (2,1), (8,2), (9,9)]

    assert(distance((0,0), (3,4))==25)                        # Euclidean distance
    assert(distances==1)                                           # distances
    assert(inc_position(sample, position, 1, 5) == 1)    # inc_position
     
    assert(inside((0,0), (3,4), 25))                                 # inside
    assert(inside((3,4), (0,0), 25))                             
    assert(not inside((0,3), (4,0), 24))

    assert(point_outside(sample, 2, 1, 5, (5,5), 7)==((2,1), 2))  # point_outside & inc_position
    assert(point_outside(sample, 3, 1, 5, (5,5), 50)==(None, 3))

    assert(ball((1,2), (1,4))==((1,3), 1))                    # ball

    print("\nTests OK")    
    stop()
    
#====================================== main program
"""
Draw sample of 2-D points and determine the most distant points
evaluate performance
"""
if Norm!=2: print("\n*** Warning: this distance is not Euclidean ***\n")

if Test: unit_tests()

if Step>1:
    print("Warning: with a Step>1 you can avoid randomization, but it is your responsibility")
    print("to verify that Step is relative prime to Set_size in order to avoid avoid subsampling")
    print("suggestion: find a prime number near 5*Set_size/18")
    print("Your Step:", Step, "\n")

# compute increasing size list
if not Test_set:
    size, max_size, coef=Size_params
    Size_list=[]
    while size<=max_size:
        Size_list+=[size]
        size=int(size*coef)
    info=" width "+str(Width) if Shape in [annulus, ring] else ""
    info2=" Norm "+str(Norm)+ " preprocess" if Preprocess else " no prep."
    print(Shape+info+" dim "+str(Dimension)+info2)
    print("\nPoints\tDist./N\tSteps\tElim%\tt/pt(µs)")
else:
    Size_list=[len(Test_set)]

if Draw: initContainer()

means=[0, 0, 0, 0]                                                      # global stats 
for index, Set_size in enumerate(Size_list):
    for rep in range(Repeat):

        # performance stats
        steps=0                                                              # number of subsets examined
        distances=0                                                        # number of distances measured
        elim=0                                                               # number of points eliminated
        
        # initialize set of points
        points=create_set(Set_size)
        if Test: points=Test
        if Verify: points_copy=list(points)

        # display
        if Draw: draw_points(points)
        
        # compute ball
        position=2*Step        
        diameter=[points[0], points[Step]]
        size=Set_size
        steps=1
        prev_steps=0
        
        # start time
        t0=time.perf_counter()

##################################################
        # first stage: elongated
        if Preprocess:
            while True:
                if steps!=prev_steps:                    # avoid repeating
                    O,radius=ball(*diameter) 
                    apart=distance(*diameter)
                    prev_steps=steps
                
                # try finding one outside point            
                C, position =point_outside(points, position, Step, Set_size, diameter[1], apart)
                if not C: break
                diameter=(diameter[1],C)
                position=inc_position(points, position, Step, Set_size)       # next
                steps+=1

        # second stage: circle                
        while True:
            if steps!=prev_steps:                             # avoid repeating 
                O, radius = ball(*diameter) 
                apart=distance(*diameter)
                prev_steps=steps
            
            # try finding one outside point
            position=inc_position(points, position, Step, Set_size)       # next
            C, position =point_outside(points, position, Step, Set_size, O, radius)
            if not C: break                                       # AB is solution

            # try finding another, more distant outside point
##            print("\nfound", C)
            C_pos=position
            position=inc_position(points, position, Step, Set_size)       # next
            D, position =point_outside(points, position, Step, Set_size, C, apart)
            if not D:
##                print("eliminate", points[position])
                points[C_pos]=None                         # eliminate C
                elim+=1
                continue
            
            diameter=[C, D]
            steps+=1

##################################################

        # end time
        t1=time.perf_counter()
        rt=(t1-t0)*1000
        
        # print stats
        print(f"{Set_size}\t{distances/Set_size:3.3}\t{steps}\t{elim*100/Set_size:3.3}\t{rt*1000/Set_size:3.3}")
        means[0]+=distances/Set_size
        means[1]+=steps
        means[2]+=elim*100/Set_size
        means[3]+=rt*1000/Set_size

        # just a verification
        if Verify:
            radius=distance(*diameter)
            for A in points_copy:
                for B in points_copy:
                    d=distance(A, B)
                    if d>radius:                     
                        print("\n*** Not solved ! ***\n")
                        stop()

        # display
        if Draw:
            draw_points(diameter, size=4, color='red')
            stop()

# print global stats            
distances=means[0]/Repeat/len(Size_list)
steps=means[1]/Repeat/len(Size_list)
elim=means[2]/Repeat/len(Size_list)
rt=means[3]/Repeat/len(Size_list)
print(f"\nMean\t{distances:3.3}\t{steps:3.3}\t{elim:3.3}\t{rt:3.3}")
stop()

