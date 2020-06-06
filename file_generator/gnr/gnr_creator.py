import numpy as np
import sys

def main(argv):
    atoms = []

    band_thickness = int(argv[0])
    band_length = int(argv[1])

    minor_atoms = int((band_thickness-1)/2)
    major_atoms = int((band_thickness+1)/2)

    y_dist = 2.46
    reg_dist = 1.42
    half_reg_dist = reg_dist / 2.0

    y = 0
    x = 0

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

    """
    for i in atoms:
        print("%.5f" % i, end=' ')
    print("")
    """

    z = 0.0

    file_out_name = "{}_AGNR_len_{}.xyz".format(str(band_thickness), str(band_length))
    print("generated file is named:", file_out_name)

    f = open(file_out_name, "w+")

    f.write("%d\n\n" % (len(atoms)/2))
    for i in range(0, len(atoms), 2):
        f.write("C   %9.5f %9.5f %9.5f\n" % (atoms[i], atoms[i+1], z))

if __name__ == '__main__':
    main(sys.argv[1:])
