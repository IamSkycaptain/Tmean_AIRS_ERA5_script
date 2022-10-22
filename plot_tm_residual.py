from matplotlib import pyplot as plt
import os
import sys
import pandas as pd
import numpy as np
import matplotlib
dir_app = os.path.dirname(sys.argv[0])
matplotlib.rc('font',**{'family':'Arial','size' : 20})

airs = os.path.join(dir_app,'result\\airs_residual.csv')
efmwf = os.path.join(dir_app,'result\\ecmwf_residual.csv')
df_airs = pd.read_csv(airs)
df_ecmw = pd.read_csv(efmwf)
std_airs = np.std(df_airs['tm'])
std_ecmw = np.std(df_ecmw['tm'])
avg_airs = np.mean(df_airs['tm'])
avg_ecmw = np.mean(df_ecmw['tm'])

cl01='black'
cl02='#3498DB'

fig, ax = plt.subplots();
fig.set_size_inches(6,5)
# the histogram of the data
plt.hist(df_ecmw['tm'], 100, density=True, facecolor=cl01)
plt.hist(df_airs['tm'], 100, density=True, facecolor=cl02, alpha=0.75)

plt.grid(color='gray', linestyle=':')
plt.xlabel('Observed minus predicted (K)')
plt.ylabel('Density')
#plt.title('Residual distribution',fontweight="bold")
#plt.title('(b)')
plt.text(0.4,.22, f'ERA5-rms={std_ecmw:.1f} K')
plt.text(0.4,.17, f' AIRS-rms={std_airs:.1f} K',fontdict={'color':  cl02})
plt.tight_layout()
plt.savefig(os.path.join(dir_app, 'tm_residual.png'), dpi=300)
plt.show()