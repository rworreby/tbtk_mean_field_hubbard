import numpy as np

#x = np.array()
atoms = []

band_thickness = 9 #thickness of gnr
band_length = 1

minor_atoms = int((band_thickness-1)/2)
major_atoms = int((band_thickness+1)/2)

y_dist = 2.46
reg_dist = 1.42
half_reg_dist = reg_dist / 2.0

y = 0
x = 0

whole_x = False

# print("y: %3f" % y)
# print("x: %3f" % x)

print(atoms)
for i in range(band_length):
    y = y_dist / 2.0

    for i in range(minor_atoms):
        atoms += {x}
        atoms += {y}
        y += y_dist


    x += half_reg_dist
    y = 0

    for i in range(major_atoms):
        atoms += {x}
        atoms += {y}
        alt_x = x + reg_dist
        atoms += {alt_x}
        atoms += {y}
        y += y_dist

    x += (reg_dist + half_reg_dist)
    y = y_dist / 2.0

    for i in range(minor_atoms):
        atoms += {x}
        atoms += {y}
        y += y_dist

    x += reg_dist

for i in atoms:
    print("%.5f" % i, end=' ')
print("")

z = 0.0
f = open("gnr_17_periodic.xyz", "w+")

f.write("%d\n\n" % (len(atoms)/2))
for i in range(0, len(atoms), 2):
    f.write("C   %9.5f %9.5f %9.5f\n" % (atoms[i], atoms[i+1], z))
