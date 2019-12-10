import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('eigenval_eigenvec.txt')
plt.plot(df, c='black')
#plt.ylim(-5,5)
plt.xlim(0, 99)
plt.savefig('figures/plotted_ev.png')
plt.show()
