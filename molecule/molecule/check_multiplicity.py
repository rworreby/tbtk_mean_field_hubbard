import pandas as pd


f = open("processed_ev.txt", "r")

num_lines = sum(1 for line in open('processed_ev.txt', "r"))
print("Number of lines in file: ", num_lines)

counter_spin_up = 0
counter_spin_down = 0

for line in f:
    splitted_line = line.split(' ')
    if splitted_line[0] == 'sd':
        counter_spin_down += 1
    elif splitted_line[0] == 'su':
        counter_spin_up += 1
    else:
        print("Error occured. Spin is neither up nor down.")

print("Total spin up:   ", counter_spin_up)
print("Total spin down: ", counter_spin_down)
