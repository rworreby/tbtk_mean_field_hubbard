import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('results.txt')
plt.plot(df)
plt.ylim(-6,2)
plt.show()
plt.savefig('plotted_ev.png')
