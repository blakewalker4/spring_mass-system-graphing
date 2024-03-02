'''
Description: Calculates the motion of an mass-spring system without an external force given the mass, damping constant, and spring constant
             as well as the initial conditions. Uses the matplotlib library to plot the graph.
Date Completed: 3/1/24
Creator: Blake Walker
'''

from numpy import arange
from math import sqrt, exp, sin, cos
from matplotlib import pyplot as plt

print("This will create a graph of the position of a mass-spring system given the mass, damping constant, and spring constant, as well as the initial conditions")
t_max = float(input("Enter to what time (in seconds) you would like to graph: "))
m = float(input("Enter the mass (in kg): "))
b = float(input("Enter the damping coefficient (in kg/sec): "))
k = float(input("Enter the spring constant: "))
initial_position = float(input("Enter the initial displacement of the spring from equalibrium (positive is down): "))
initial_velocity = float(input("Enter the initial velocity of the spring (positive is down): "))
t = 0
t_step = t_max/10000
t_array = list(arange(0,t_max+t_step,t_step)) # creates the array of t values used in the graph
y_array = []

#solve auxillary equation
descriminant = b*b-(4*m*k)
if(descriminant > 0): # overdampened
    aux1 = (-1*b + sqrt(b*b -(4*m*k)))/(2*m)
    aux2 = (-1*b - sqrt(b*b -(4*m*k)))/(2*m)
    # solving for the constants
    c1 = (-1*aux2*initial_position + initial_velocity)/(aux1-aux2)
    c2 = (-1*aux1*initial_position + initial_velocity)/(aux2-aux1)
    print(f"The equation describing this graph is: y ≈ {round(c1,2)}e^({round(aux1,2)}t)+{round(c2,2)}e^({round(aux2,2)}t)")
    while t < t_max:
        for t in t_array:
            y1 = c1*exp(aux1*t)+c2*exp(aux2*t)
            y_array.append(y1)
            t+=0.01

if(descriminant == 0): # critically dampened
    aux = (-1*b)/(2*m)
    # solving for the constants
    c1 = initial_position                       
    c2 = initial_velocity + initial_position
    print(f"The equation describing this graph is: y ≈ e^({round(aux,2)}t)({round(c1,2)}+{round(c2,2)}t)")
    for t in t_array:
        y1 = exp(aux*t)*(c1 + c2*t)
        y_array.append(y1)
        t+=0.01

if(descriminant < 0): # underdampened
    alpha = (-1*b)/(2*m) # real component
    beta = sqrt(abs(b*b-4*m*k))/(2*m) # imaginary component
    # solving for the constants
    c1 = initial_position
    c2 = initial_velocity/beta
    print(f"The equation describing this graph is: y ≈ e^({round(alpha,2)}t)({c1}cos({round(beta,2)}t)+{round(c2,2)}sin({round(beta,2)}t))")
    for t in t_array:
        y1 = exp(alpha*t)*(c1*cos(beta*t)+c2*sin(beta*t))
        y_array.append(y1)
        t+=0.01

plt.plot(t_array, y_array)
plt.title(f"Spring-mass system with m={m}, β={b}, k={k}, y(0)={initial_position}, and y'(0)={initial_velocity}")
plt.xlabel("time")
plt.ylabel("position")
plt.axhline(0,color = 'black', linewidth=0.5)
plt.show()