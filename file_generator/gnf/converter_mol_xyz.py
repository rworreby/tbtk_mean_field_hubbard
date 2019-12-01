import sys

file_name = str(sys.argv[1])
f = open(file_name)

output = []

for i in range(3):
    f.readline()

num_atoms = f.readline().split()[0]
output.append(str(num_atoms + '\n' + '\n'))

for i in range(int(num_atoms)):
    line = f.readline().split()
    output.append(line[3] + '\t')
    for i in range(3):
        str_to_add = "{:9.4f} \t".format(round(float(line[i])*284.0/165.0, 4))
        output.append(str_to_add)
        if i == 2:
            output.append('\n')
f.close()

f = open('xyz_files/'+str(file_name[10:-3])+'xyz', 'w')
f.writelines(output)
f.close()
