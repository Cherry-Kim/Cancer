#pip install lifelines

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np

labels = ['Group 1', 'Group 2']
groups = [1, 1, 1, 1, 1, 1,   2, 2, 2, 2, 2, 2]
events = [1, 1, 1, 1, 1, 0,   1, 1, 1, 0, 1, 1]
times  = [1, 2, 3, 4, 4.5, 5,  0.5, 0.75, 1, 1.5, 2, 3.5]

E = np.array(events, dtype=np.int32)
T = np.array(times, dtype=np.float32)

rcParams.update({'font.size': 12})
fig, ax = plt.subplots(figsize=(5,5))
styles = ['-', '--']
colors = ['r', 'g']
lw = 3

kmf = KaplanMeierFitter()
for i, label in enumerate(labels):
		ix = np.array(groups) == (i+1)
		kmf.fit(T[ix], event_observed=E[ix], label=labels[i])
		kmf.plot(ax=ax, ci_show=False, linewidth=lw, style=styles[i], c=colors[i])

ix = np.array(groups) == 2
result = logrank_test(T[ix], T[~ix], E[ix], E[~ix], alpha=.99)
pvalue = result.p_value
ax.text(3.4,0.75,'P-value=%.3f'% pvalue)

ax.set_xlabel('Time', fontsize=14)
ax.set_ylabel('Survival', fontsize=14)
ax.legend(loc='upper right')

plt.tight_layout()
plt.savefig('ex_kmplot.png', format='png', dpi=300)
plt.show()
