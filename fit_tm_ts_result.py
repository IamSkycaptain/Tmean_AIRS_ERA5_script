import os
import sys
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as PathEffects
import numpy as np
#sns.set_theme(style="darkgrid")
#g = sns.jointplot(x="ts", y="tm", data=df, kind="reg")
#g.fig.set_figwidth(9)
#g.fig.set_figheight(6)
#plt.sh
#matplotlib.rcParams['text.usetex'] = True
matplotlib.rc('font',**{'family':'Arial','size' : 22})
cl01='black'
cl02='#3498DB'
fig, ax = plt.subplots();
fig.set_size_inches(8,6)

dir_app = os.path.dirname(sys.argv[0])

#Data 1 Scetter ECMWF 2015-2020
data = os.path.join(dir_app, 'result\ecmwf_tmean.txt')
df = pd.read_table(data,delim_whitespace=True)
x1,y1 = df['ts'],df['tm']
ax.plot(x1, y1, '.', color=cl01, markersize=2)

xx = [280,315]
#Data 2 Scetter AIRS 2015-2020
data = os.path.join(dir_app, 'result\hdf__all_QC_0.csv')
df = pd.read_csv(data)
df = df[(df['doy']>2016) & (df['tm']>280)]
x2,y2 = df['ts'], df['tm']
ax.plot(x2, y2, '.', color=cl02, markersize=2, alpha=0.3)

#Fitting ECMWF
result1 = np.polyfit(x1,y1, 1 , full=True)
fit_fn1 = np.poly1d(result1[0])
lf1, = plt.plot(xx, fit_fn1(xx), '-',color=cl01, linewidth=3.0)
#yy = [0.27*xx[0] + 211.6, 0.27*xx[1] + 211.6] # ECMWF
#lf1 = ax.plot(xx,yy,'-k', linewidth=3.0)
plt.setp(lf1, path_effects=[
        PathEffects.withStroke(linewidth=7, foreground="w")])
#Fitting AIRS ใช้จริงมีการลดจำนวนทศนิยม       
result2 = np.polyfit(x2,y2, 1 , full=True)
ab = result2[0]
print(ab)
fit_fn2 = np.poly1d(result2[0])
lf2, = ax.plot(xx, fit_fn2(xx), '-',color=cl02, linewidth=3.0)
#yy = [0.3764*xx[0] + 177.1393, 0.3764*xx[1] + 177.1393] # AIRS
#lf2 = ax.plot(xx,yy,'-b', linewidth=3.0)
plt.setp(lf2, path_effects=[
        PathEffects.withStroke(linewidth=7, foreground="w")])
        

# Setting Histrogram
dv = make_axes_locatable(ax)
axHistx = dv.append_axes("top"  , 1.25, pad=0.2, sharex=ax)
axHisty = dv.append_axes("right", 1.25, pad=0.1, sharey=ax)
axHistx.xaxis.set_tick_params(labelbottom=False)
axHisty.yaxis.set_tick_params(labelleft=False)
axHistx.grid(color='gray', linestyle=':')
axHisty.grid(color='gray', linestyle=':')
#axHisty.set_yticks([0, 0.1, 0.2])
axHisty.set_xlabel('Density')
axHistx.set_ylabel('Density')
#axHistx.set_title('Ts-Tm Linear Regression models', fontsize="18", fontweight="bold")   
#axHistx.set_title('(a)')            
# Data 1 ECMWF    
xmax,xmin = np.max(x1),np.min(x1)
binwidth = (xmax-xmin)/100
bins = np.arange(xmin,xmax+binwidth, binwidth)
axHistx.hist(x1, bins=bins,density=True,facecolor=cl01)

ymax,ymin = np.max(y1),np.min(y1)
binwidth = (ymax-ymin)/100
bins = np.arange(ymin,ymax+binwidth, binwidth)
axHisty.hist(y1, bins=bins,density=True,facecolor=cl01,orientation='horizontal')

# Data 2 AIRS    
xmax,xmin = np.max(x2),np.min(x2)
binwidth = (xmax-xmin)/100
bins = np.arange(xmin,xmax+binwidth, binwidth)
axHistx.hist(x2, bins=bins,density=True,facecolor=cl02, alpha=0.75)

ymax,ymin = np.max(y2),np.min(y2)
binwidth = (ymax-ymin)/100
bins = np.arange(ymin,ymax+binwidth, binwidth)
axHisty.hist(y2, bins=bins,density=True,facecolor=cl02, alpha=0.75, orientation='horizontal',)
           
ax.set_ylabel('Tm (K)')
ax.set_xlabel('Ts (K)')
ax.grid(True)
ax.axis('scaled')
ax.set_ylim([275, 305])
leg = ax.legend([lf1,lf2], [r"derived by ERA5, R=0.54",r"derived by AIRS, R=0.49"], fancybox=True, loc=2)
leg.legendPatch.set_path_effects([PathEffects.withSimplePatchShadow()])


plt.tight_layout()
plt.savefig(os.path.join(dir_app, 'fit_airs_ecmwf.png'), dpi=300)
plt.show()