import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('convergence_data.txt')

for i in range(10):
    x = df[df['run'] == i][df['iteration'] < 100]['iteration']
    y = df[df['run'] == i][df['iteration'] < 100]['error']

    plt.plot(x, y, linewidth=1)

#plt.ylim(-5,5)
#plt.xlim(-2, 99)
plt.xlabel("Number of iterations")
plt.ylabel("Total absolute error")
#plt.xticks(df['atoms'], df['size'], rotation=40)
#ax2 = ax1.twiny()
#plt.xticks(df['atoms'], df['size'], )
#plt.savefig('gnf_scaling.svg')
plt.show()
