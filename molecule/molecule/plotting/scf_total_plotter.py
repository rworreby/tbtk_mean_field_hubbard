import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

sns.set()

df = pd.read_csv('../results/scf_total_conversion.txt')

for i in range(10):
    x = df[df['run'] == i][df['iteration'] < 150]['iteration']
    y = df[df['run'] == i][df['iteration'] < 150]['error']

    plt.plot(x, y, linewidth=1)

#plt.ylim(-5,5)
#plt.xlim(-2, 99)
plt.xlabel("Number of Iterations")
plt.ylabel("Total Absolute Error")
plt.yscale('log')
#plt.xticks(df['atoms'], df['size'], rotation=40)
#ax2 = ax1.twiny()
#plt.xticks(df['atoms'], df['size'], )
plt.savefig('scf_total_convergence.pdf')
plt.show()
