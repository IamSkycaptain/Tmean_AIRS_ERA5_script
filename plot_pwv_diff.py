import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patheffects as PathEffects
from sympy import cos, init_printing, Matrix
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1 import make_axes_locatable

dir_app = os.path.dirname(sys.argv[0])
def pwv_timeseries():
    matplotlib.rc('font',**{'family':'Arial','size' : 14}) # Place this script on the top before call library
    fig,ax = plt.subplots()
    fig.set_size_inches(6.7,4)
    for root,dirs,files in os.walk(os.path.join(dir_app,r'result\pwvs'),topdown=True):
        for file in files:
            data = os.path.join(root,file)
            df = pd.read_csv(data)   
            sta = file[0:4]
            nf = df[(df['pwv_gpt3']> 10) & (df['pwv_gpt3'] < 70)]
            ax.plot(nf['ydec'],nf['pwv_gpt3'],':',label=sta,markersize=0.5, linewidth=1)
            #break
    plt.xticks([2021.041,2021.126,2021.2055,2021.2904,2021.3726,2021.4575,2021.5397,2021.6245,2021.7096,2021.7917,2021.8767,2021.9589],
                ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
    yyy = [10,70]
    ax.plot([2021.126,2021.126],yyy,'-k',linewidth=1)
    ax.plot([2021.3726,2021.3726],yyy,'-k',linewidth=1)
    ax.plot([2021.7917,2021.7917],yyy,'-k',linewidth=1)
    plt.text(2021.00,60, 'Cool',dict(fontsize=16))
    plt.text(2021.22,60, 'Hot',dict(fontsize=16))
    plt.text(2021.54,60, 'Rainy',dict(fontsize=16))
    plt.text(2021.88,60, 'Cool',dict(fontsize=16))
    
    ax.grid(color='gray', linestyle=':')
    ax.set_ylabel('PWV (mm)') # ,**{'size': 12}
    leg = plt.legend(bbox_to_anchor=(1.0,-0.08),ncol=8, handletextpad=0.1, labelspacing=0.1,columnspacing=0.1,prop={'size':10})
    #llines = leg.get_lines()
    #plt.title('GPS-PWVs in 2021', fontsize="14", fontweight="bold")   
    plt.tight_layout()
    plt.savefig(os.path.join(dir_app, 'pwv_time_series.png'), dpi=300)
    plt.show()
    
def pwv_mean_bias():
    outfile = os.path.join(dir_app,'result\pwv_compare_with_h_corr_2021.csv')
    df = pd.read_csv(outfile)
    nf = df[df['tm_model'] == 'meo']
    matplotlib.rc('font',**{'family':'Arial','size' : 14}) # Place this script on the top before call library
    fig,ax = plt.subplots()
    fig.set_size_inches(10,4)
    ax.plot(df['sta'],df['airs'] ,'-k' ,markersize=3.0, linewidth=2, label="AIRS-Tm")
    ax.plot(df['sta'],df['ecmwf'],'-+y' ,markersize=3.0, linewidth=2, label="ERA5-Tm")
    ax.plot(df['sta'],df['rata'] ,'-+b' ,markersize=3.0, linewidth=2, label="Suwantong")
    ax.plot(df['sta'],df['bevis'],'-+' ,markersize=3.0, linewidth=1, label="Bevis")
    ax.plot(df['sta'],df['mendes'],'-+',markersize=3.0, linewidth=1, label="Mendes")
    ax.plot(df['sta'],df['schue'],'-+' ,markersize=3.0, linewidth=1, label="Schueler")
    ax.plot(df['sta'],df['gpt3'] ,'-+' ,markersize=3.0, linewidth=1, label="GPT3")

    ax.grid(color='gray', linestyle=':')
    ax.set_xlabel('Station')
    ax.set_ylabel('Mean bias (mm).')
    #leg = ax.legend(fancybox=True, loc=2)
    #leg.legendPatch.set_path_effects([PathEffects.withSimplePatchShadow()])
    plt.legend(ncol=7, handletextpad=0.1, labelspacing=0.1,columnspacing=0.1,loc=4)
    old_label = 'ARNG,ATRG,AUPG,AYYA,BLMG,BTAK,CHAN,CHDN,CHMA,CHPM,CLPK,CMPN,CTAK,DBRM,DCMI,DCRI,DMSN,DPLK,DPNB,DPT9,DSSK,DYLA,ECMI,ENMA,HACH,KNKN,KNSN,LCNT,LKPT,LLPN,LMDH,LNAN,LNBP,LNSN,LPBR,MHGS,MRBR,MSOD,NKBI,NKNY,NKRM,NKSW,NRTW,PJRK,PNST,PNTG,PTLG,PYAO,RAND,RAYG,SATN,SICN,SISA,SKNK,SNCK,SOKA,SURN,TEPA,THAI,THKP,THPP,TKRI,TNPM,TNPT,TNST,TPHN,TPKT,TSKA,TSKW,TSRI,TUBN,UDON,UTOG,UTTD,VCBR,WSPG'.split(',')
    new_label = 'ARNG,,AUPG,,BLMG,,CHAN,,CHMA,,CLPK,,CTAK,,DCMI,,DMSN,,DPNB,,DSSK,,ECMI,,HACH,,KNSN,,LKPT,,LMDH,,LNBP,,LPBR,,MRBR,,NKBI,,NKRM,,NRTW,,PNST,,PTLG,,RAND,,SATN,,SISA,,SNCK,,SURN,,THAI,,THPP,,TNPM,,TNST,,TPKT,,TSKW,,TUBN,,UTOG,,VCBR,'.split(',')
    plt.xticks(old_label,new_label)
    plt.xticks(rotation=90,**{'size': 10})
    plt.tight_layout()
    plt.savefig(os.path.join(dir_app, 'pwv_mean_bias.png'), dpi=300)
    plt.show()
                
def up_tmseries(times, a, b, c, d):
    #return (a + b*times + c*np.cos(2*np.pi*(times)/365.5) + d*np.sin(2*np.pi*(times)/365.5))
    return (a + b*times + c*np.cos(2*np.pi*(times+d)/365.5))

def differenct_Pwv():    
    data = os.path.join(dir_app,'result\with_h_corr_merge.csv')
    df = pd.read_csv(data)
    df = df.sort_values(by='doy', ascending=True)
    doy = (df['doy']-2020)*365.5
    diff_pwv = df['ref_pwv']-df['rata_meo']
    coff, pcov = curve_fit(up_tmseries, doy, diff_pwv)
    print(f'a={coff[0]:.2f}, b={coff[1]:.3f}, c={coff[2]:.2f}, d={coff[3]:.2f}')
     
    matplotlib.rc('font',**{'family':'Arial','size' : 24})
    #df = df[df['sta']=='AUPG']
    fig,ax = plt.subplots()
    fig.set_size_inches(8,5.5)
    #plt.title('              The Difference of PWVs in 2021',fontsize="18", fontweight="bold")
    plt.plot(df['doy']+1,diff_pwv,'.',markersize=3.0 ,label='RMS of the different PWVs = 2.7 mm')
    plt.xticks([2021.041,2021.126,2021.2055,2021.2904,2021.3726,2021.4575,2021.5397,2021.6245,2021.7096,2021.7917,2021.8767,2021.9589],
               ['1','2','3','4','5','6','7','8','9','10','11','12'])
    #['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    plt.plot(df['doy']+1, up_tmseries(doy, *coff), '-k',linewidth=4.0)
    #plt.xticks(rotation=20)

    dv = make_axes_locatable(ax)
    axHisty = dv.append_axes("right", 1.25, pad=0.1, sharey=ax)
    axHisty.yaxis.set_tick_params(labelleft=False)
    axHisty.grid(color='gray', linestyle=':')
    ymax,ymin = np.max(diff_pwv),np.min(diff_pwv)
    binwidth = 0.5
    bins = np.arange(ymin,ymax+binwidth, binwidth)
    axHisty.hist(diff_pwv, bins=bins,density=False,facecolor='#3498DB',orientation='horizontal')
    axHisty.set_xlabel('Count')
    ann = ax.annotate("Seasonal variation with \n a maximum amplitude of 0.7 mm.",
                      xy=(2021.3, 0.0), xycoords='data',
                      xytext=(2021.5, 5), textcoords='data',
                      size=24, va="center", ha="center",
                      bbox=dict(boxstyle="round4", fc="w"),
                      arrowprops=dict(arrowstyle="simple",
                                      connectionstyle="arc3,rad=0.2",
                                      fc="r"),
                      )
                      
    ax.set_ylim([-10, 10])
    ax.set_ylabel('diff. PWVs (mm)')
    ax.set_xlabel('Month of 2021')
    ax.grid(color='gray', linestyle=':')
    #ax.legend(fancybox=True, loc=3)
    plt.tight_layout()
    plt.savefig(os.path.join(dir_app, 'plot_pwv_diff.png'), dpi=300)
    plt.show()
#pwv_timeseries()
pwv_mean_bias()
#differenct_Pwv()
    