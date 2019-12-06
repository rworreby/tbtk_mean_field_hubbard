import sys

file_name = 'eigenval_eigenvec.txt'
f = open(file_name, 'r')

lines = f.readlines()

sorted_lines = sorted(lines, key=lambda x : float(x.split()[1]), reverse=False)

f.close()

f = open('sorted_ev_ev.txt', 'w')
f.writelines(sorted_lines)
f.close()
