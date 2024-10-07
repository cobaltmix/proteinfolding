from vpython import *
import random
import math

c = 0

# Function to initialize bonds between amino acids
def initialize_membrane(c):
    
    for i in range(180):
        if c <= 180:
            c += 1
        elif c > 180:
            c -= 1
        for j in range(c):
            sphere(pos=vec(math.sin(c/180),math.cos(c/180),), color=color.white, radius=0.1)


# Initialize membrane
initialize_membrane(c)

# Animation loop
while True:
    rate(120)
