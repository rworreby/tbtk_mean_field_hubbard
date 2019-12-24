import numpy as np
import math

#x = np.array()
atoms = []

band_thickness = 4
band_length = 3
repetitions = int(band_thickness/4)

y_dist = 2.46
reg_dist = 1.42
half_reg_dist = reg_dist / 2.0
spec_dist = 1.23

y = 0
x = 0

whole_x = False

print(atoms)
for i in range(band_length):
    y = half_reg_dist

    for i in range(repetitions):
        atoms += {x}
        atoms += {y}
        y += reg_dist
        atoms += {x}
        atoms += {y}
        y += (2.0 * reg_dist)

    x += spec_dist
    y = 0

    for i in range(repetitions):
        atoms += {x}
        atoms += {y}
        #alt_x = x + reg_dist
        y += (2.0 * reg_dist)
        atoms += {x}
        atoms += {y}
        y += reg_dist


    y = 0

    """
    x += (reg_dist + half_reg_dist)
    y = y_dist / 2.0

    for i in range(minor_atoms):
        atoms += {x}
        atoms += {y}
        y += y_dist
    """
    x += spec_dist

for i in atoms:
    print("%.5f" % i, end=' ')
print("")

z = 0.0
f = open("gnr_8_periodic.xyz", "w+")

f.write("%d\n\n" % (len(atoms)/2))
for i in range(0, len(atoms), 2):
    f.write("C   %9.5f %9.5f %9.5f\n" % (atoms[i], atoms[i+1], z))
