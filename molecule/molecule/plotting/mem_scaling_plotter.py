import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

sns.set()

df = pd.read_csv('../results/memory_profiling.txt')

print(df)
plt.plot(df['gnf_len'], df['memory_usage'], c='black')

#plt.ylim(-5,5)
#plt.xlim(0, 99)
plt.ylabel("Memory Usage in MiB")
plt.xlabel("Size of Graphene Nanoflake")
plt.xticks(df['gnf_len'], df['gnf_size'], rotation=40)
#ax2 = ax1.twiny()
#plt.xticks(df['atoms'], df['size'], )

plt.savefig('gnf_mem_usage.pdf')
plt.show()
