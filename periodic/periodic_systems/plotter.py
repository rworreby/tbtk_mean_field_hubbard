import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('results.txt')
plt.plot(df, c='black')
#plt.ylim(-5,5)
#plt.xlim(0, 99)
df_size = len(df)
plt.xticks(np.arange(0,df_size, (df_size-3)/2 ), ('-π/a', '0', 'π/a'))
plt.show()
#plt.savefig('plotted_ev.png')
