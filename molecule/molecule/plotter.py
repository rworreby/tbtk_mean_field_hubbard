import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('gnf_timings.txt')


plt.plot(df['atoms'], df['time'], c='black')

#plt.ylim(-5,5)
#plt.xlim(0, 99)
plt.ylabel("Execution time in s")
plt.xlabel("Size of graphene nanoflake")
plt.xticks(df['atoms'], df['size'], rotation=40)
#ax2 = ax1.twiny()
#plt.xticks(df['atoms'], df['size'], )
#plt.savefig('gnf_scaling.svg')
plt.show()
