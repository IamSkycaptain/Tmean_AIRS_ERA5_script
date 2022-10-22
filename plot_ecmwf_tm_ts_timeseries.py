import pandas as pd
import fnmatch
import os,sys
from matplotlib import pyplot as plt
import matplotlib
dir_app = os.path.dirname(sys.argv[0])
def extract(outfile):
    #hdf__all_QC_0
    work_dir = r'E:\ECMWF\tmean'
    px = 0.12
    lat,lon = 14.36431,100.576766   #AYYA
    xmin,xmax = lon-px,lon+px
    ymin,ymax = lat-px,lat+px

    with open(outfile,'w')as outf:
        outf.write('doy,lat,lon,tm,ts,ydoy\n')
        for root,dirs,files in os.walk(work_dir,topdown=True):
            for fname in fnmatch.filter(files,'*.csv'):
                data = os.path.join(root,fname)
                df = pd.read_csv(data)
                df = df[df['hr']==6]
                nf = df[((df['lon'] >= xmin) & (df['lon'] < xmax)) & ((df['lat'] >= ymin) & (df['lat'] < ymax))] 
                if len(nf)>0:
                    print(len(nf))
                    nf = nf[['doy','lat','lon','tm','ts']]
                    nf['ydoy'] = int(fname[0:4]) + nf['doy']/366 
                    outf.write(nf.to_csv(index=False,line_terminator='').replace('doy,lat,lon,tm,ts,ydoy',''))

outfile = os.path.join(dir_app,'result\ecmwf_tm_ts_timeseries.csv')
#extract(outfile)
col01 = '#C70039'
col02 = '#754405'
df = pd.read_csv(outfile)
nf= df[df['ts']>295]
#fig,(ax1,ax2) = plt.subplots(2, sharex=True, sharey=False)
matplotlib.rc('font',**{'family':'Arial','size' : 24})
fig,ax1 = plt.subplots()
fig.set_size_inches(9,3.5)
fig.subplots_adjust(hspace=0)
#ax1.set_title('Tm and Ts from ERA5',fontweight="bold")
#ax1.set_title('(b)')
ax1.plot(nf['ydoy']-1, nf['tm'],'-s' ,color=col01, markersize=2.0, linewidth=1,label="Tm")
ax1.grid(color='gray', linestyle=':')
ax1.set_ylabel('ERA5 Tm (K)')

ax1.plot(nf['ydoy']-1, nf['ts'],'-x' ,color=col02, markersize=2.0, linewidth=1,label="Ts")
#ax2.grid(color='gray', linestyle=':')
#ax2.set_ylabel('Ts (K)')
ax1.legend(fancybox=True, loc=4, ncol=2)
plt.tight_layout()
plt.savefig(os.path.join(dir_app, 'ecmwf_tm_ts_timeseries.png'), dpi=300)
plt.show()
