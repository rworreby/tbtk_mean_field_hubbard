import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

sns.set()

df = pd.read_csv('results.txt')
plt.plot(df, c='black')

# plt.ylim(-5,5)
# plt.xlim(0, 99)

df_size = len(df)

plt.xlabel('k')
plt.ylabel('Energy')
plt.xticks(np.arange(0, df_size, (df_size-3)/2), ('-π/a', '0', 'π/a'))

# plt.show()
plt.savefig('finished_plot.pdf')
