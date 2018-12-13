from thetis import *

import os
from subprocess import call



call(["mkdir", "inputs"])
os.chdir("inputs")
call(["rm", "mesh.msh"])

f1 = open("mesh.geo","w")


outline = [[0,0],[16100,0], [16100,6000], [0,6000]]

resolution = [200 for i in range(len(outline))]

print (len(outline), resolution)

def gmsh_generator(outline,resolution):

    for i in range(len(outline)):
        f1.write('Point('+str(i+1)+') = { ' +"{}, {}, 0, {}".format(outline[i][0], outline[i][1],resolution[i]) + "}; \n")

    for i in range (len(outline)-1):
        f1.write('Line('+str(i+1)+') = { ' +"{}, {}".format(len(outline)-i, len(outline)-1-i,) + "}; \n")

    # Final connection

    f1.write('Line('+str(len(outline))+') = { ' +"{}, {}".format(1 , len(outline),) + "}; \n")
    f1.write('Line Loop(1) = {')

    for i in range (len(outline)):
        f1.write(str(i+1))
        if i < len(outline)-1:
            f1.write(", ")

    f1.write('};\n')

    f1.write('Plane Surface(6) = {1};\n')

    for i in range (len(outline)):
        f1.write('Physical Line('+str(i+1)+') = { '+"{}".format(i+1) + "}; \n")

    f1.write('Physical Surface(11) = {6};\n' )
    f1.write('Mesh.Algorithm = 6; // frontal=6, delannay=5, meshadapt=1')

    f1.close()


gmsh_generator(outline,resolution)


call(["gmsh", "mesh.geo", "-2", "mesh.msh"])


print("done")
