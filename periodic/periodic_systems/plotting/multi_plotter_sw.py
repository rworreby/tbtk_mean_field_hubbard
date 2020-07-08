import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

sns.set()

df1 = pd.read_csv('data_multiplots/results_left.txt')
df2 = pd.read_csv('data_multiplots/results_left_3ul.txt')
df3 = pd.read_csv('data_multiplots/results_right_sw.txt')
df4 = pd.read_csv('data_multiplots/results_right_h_sw.txt')

df_size = len(df1)
hs_df = int(df_size/2)


df1 = df1[int(len(df1)/2):]
df2 = df2[int(len(df2)/2):]
df3 = df3[int(len(df3)/2):]
df4 = df4[int(len(df4)/2):]

ax1 = plt.subplot2grid((1, 6), (0, 0), colspan=3)

ax1.plot(df1, c='black')
# plt.xlabel('k')
ax1.tick_params(labelleft=False)
plt.ylabel('Energy')
plt.xticks(np.arange(hs_df, df_size, (df_size-3)/2), ('0', 'π/a'), size='8', ha='center')
# plt.xticks(np.arange(0, df_size, (df_size-3)/2), ('π/a', '0', 'π/a'), size='8', ha='center')

ax2 = plt.subplot2grid((1, 6), (0, 3), colspan=1)
ax2.tick_params(labelleft=False)
# plt.xlabel('k')
plt.xticks(np.arange(hs_df, df_size, (df_size-3)/2), ('0', 'π/3a'), size='8', ha='right')
ax2.plot(df2, c='black')


ax3 = plt.subplot2grid((1, 6), (0, 4), colspan=1)
ax3.tick_params(labelleft=False)
ax3.plot(df3, c='black')
# plt.xlabel('k')
plt.xticks(np.arange(hs_df, df_size, (df_size-3)/2), ('0', 'π/3a'), size='8', ha='right')

ax4 = plt.subplot2grid((1, 6), (0, 5), colspan=1)
ax4.tick_params(labelleft=False)
ax4.plot(df4, c='black')
# plt.xlabel('k')
plt.xticks(np.arange(hs_df, df_size, (df_size-3)/2), ('0', 'π/a'), size='8', ha='center')

# plt.tight_layout(pad=0.8)

# plt.show()
plt.savefig('finished_plot.pdf')
