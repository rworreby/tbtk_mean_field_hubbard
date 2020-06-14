import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

sns.set()

df = pd.read_csv('results/gnf_timings.txt')

plt.plot(df['atoms'], df['time'], c='black')

#plt.ylim(-5,5)
#plt.xlim(0, 99)
plt.ylabel("Execution Time in s")
plt.xlabel("Size of Graphene Nanoflake")
plt.xticks(df['atoms'], df['size'], rotation=40)
#ax2 = ax1.twiny()
#plt.xticks(df['atoms'], df['size'], )
plt.savefig('gnf_scaling.pdf')
plt.show()
