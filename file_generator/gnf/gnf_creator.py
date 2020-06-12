import sys


def main(argv):
    atoms = []

    side_length = int(argv[0])

    max_y_atoms = side_length + 1

    y_dist = 2.46
    reg_dist = 1.42
    half_y_dist = y_dist / 2.0
    half_reg_dist = reg_dist / 2.0

    y = 0
    x = 0

    counter = 0

    num_columns = (2 * max_y_atoms - 1)
    peak_column = max_y_atoms

    for i in range(num_columns):
        if counter % 2:
            y += half_y_dist

        distance_to_mid = abs(peak_column - i - 1)

        y_counts = 0
        y_levels = peak_column - distance_to_mid

        y_adjustment = (y_levels - 1) // 2
        y += y_adjustment * y_dist
        y_counts -= y_adjustment

        for j in range(y_levels):
            atoms += {x}
            atoms += {y}
            alt_x = x + reg_dist
            atoms += {alt_x}
            atoms += {y}
            y -= y_dist
            y_counts += 1
        if counter % 2:
            y -= half_y_dist

        x = alt_x + half_reg_dist
        counter += 1
        y += y_counts * y_dist

    """
    for i in atoms:
        print("%.5f" % i, end=' ')
    print("")
    """

    z = 0.0

    file_out_name = "gnf_size_{}x{}.xyz".format(
        str(side_length),
        str(side_length)
        )
    print("generated file is named:", file_out_name)

    f = open(file_out_name, "w+")

    f.write("%d\n\n" % (len(atoms)/2))
    for i in range(0, len(atoms), 2):
        f.write("C   %9.5f %9.5f %9.5f\n" % (atoms[i], atoms[i+1], z))


if __name__ == '__main__':
    main(sys.argv[1:])
