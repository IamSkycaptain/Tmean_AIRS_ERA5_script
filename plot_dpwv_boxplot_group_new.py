from matplotlib import pyplot as plt
import matplotlib
import os
import sys
import pandas as pd
import numpy as np
#from pylab import setp
                
def setBoxColors(bp):
    bp['boxes'][0].set(facecolor='pink')
    bp['boxes'][1].set(facecolor='lightblue') 
    
    bp['medians'][0].set(color='black')
    bp['medians'][1].set(color='black') 

    
matplotlib.rc('font',**{'family':'Arial','size' : 20}) # Place this script on the top before call library
dir_app = os.path.dirname(sys.argv[0])

fig, ax = plt.subplots();
fig.set_size_inches(8,5)

data = os.path.join(dir_app, 'result\with_h_corr_merge.csv')
df = pd.read_csv(data)
dPWV = df['ref_pwv'] - df['rata_meo']
pwv_std = np.std(dPWV)
df = df[np.absolute(dPWV)<pwv_std*1.7]  
         
ref = df['ref_pwv']
airs_meo = ref-df['airs_meo']
airs_gpt = ref-df['airs_gpt3']

ecmwf_meo = ref-df['ecmwf_meo']
ecmwf_gpt = ref-df['ecmwf_gpt3']

rata_meo = ref-df['rata_meo']
rata_gpt = ref-df['rata_gpt3']

bevis_meo = ref-df['bevis_meo']
bevis_gpt = ref-df['bevis_gpt3']

mendes_meo = ref-df['mendes_meo']
mendes_gpt = ref-df['mendes_gpt3']

schue_meo = ref-df['schue_meo']
schue_gpt = ref-df['schue_gpt3']

gpt3 = ref-df['gpt3_pwv']

result = f'airs:{np.mean(airs_meo):.1f}, {np.std(airs_meo):.1f}\n'
result+= f'ecmwf:{np.mean(ecmwf_meo):.1f}, {np.std(ecmwf_meo):.1f}\n'
result+= f'rata:{np.mean(rata_meo):.1f}, {np.std(rata_meo):.1f}\n'
result+= f'bevis:{np.mean(bevis_meo):.1f}, {np.std(bevis_meo):.1f}\n'
result+= f'mendes:{np.mean(mendes_meo):.1f}, {np.std(mendes_meo):.1f}\n'
result+= f'schue:{np.mean(schue_meo):.1f}, {np.std(schue_meo):.1f}\n'
result+= f'gpt:{np.mean(gpt3):.1f}, {np.std(gpt3):.1f}\n'
print('MEO\n'+result)

result = f'airs:{np.mean(airs_gpt):.1f}, {np.std(airs_gpt):.1f}\n'
result+= f'ecmwf:{np.mean(ecmwf_gpt):.1f}, {np.std(ecmwf_gpt):.1f}\n'
result+= f'rata:{np.mean(rata_gpt):.1f}, {np.std(rata_gpt):.1f}\n'
result+= f'bevis:{np.mean(bevis_gpt):.1f}, {np.std(bevis_gpt):.1f}\n'
result+= f'mendes:{np.mean(mendes_gpt):.1f}, {np.std(mendes_gpt):.1f}\n'
result+= f'schue:{np.mean(schue_gpt):.1f}, {np.std(schue_gpt):.1f}\n'
result+= f'gpt:{np.mean(gpt3):.1f}, {np.std(gpt3):.1f}\n'

print('GPT\n'+result)


airs = pd.DataFrame({ 
    'meo': airs_meo.to_numpy(), 
    'gpt': airs_gpt.to_numpy()})

ecmwf = pd.DataFrame({ 
    'meo': ecmwf_meo.to_numpy(), 
    'gpt': ecmwf_gpt.to_numpy()})   

rata = pd.DataFrame({ 
    'meo': rata_meo.to_numpy(), 
    'gpt': rata_gpt.to_numpy()})

bevis = pd.DataFrame({ 
    'meo': bevis_meo.to_numpy(), 
    'gpt': bevis_gpt.to_numpy()})    

mendes = pd.DataFrame({ 
    'meo': mendes_meo.to_numpy(), 
    'gpt': mendes_gpt.to_numpy()})   
    
schue = pd.DataFrame({ 
    'meo': schue_meo.to_numpy(), 
    'gpt': schue_gpt.to_numpy()})   
    
gpt3 = pd.DataFrame({ 
    'GPT3': gpt3.to_numpy()})  

nt = False
green_diamond = dict(color='#cccccc',markerfacecolor='#cccccc', marker='s',markersize=3)
bp = ax.boxplot(airs, notch=nt, patch_artist=True, positions = [0.5, 1.5], widths = 0.6, flierprops=green_diamond);setBoxColors(bp)
bp = ax.boxplot(ecmwf,notch=nt, patch_artist=True, positions = [3.5, 4.5], widths = 0.6, flierprops=green_diamond);setBoxColors(bp)
bp = ax.boxplot(rata, notch=nt, patch_artist=True, positions = [6.5, 7.5], widths = 0.6, flierprops=green_diamond);setBoxColors(bp)
bp = ax.boxplot(bevis,notch=nt, patch_artist=True, positions = [9.5, 10.5], widths = 0.6, flierprops=green_diamond);setBoxColors(bp)
bp = ax.boxplot(mendes,notch=nt,patch_artist=True, positions = [12.5,13.5], widths = 0.6, flierprops=green_diamond);setBoxColors(bp)
bp = ax.boxplot(schue,notch=nt, patch_artist=True, positions = [15.5,16.5], widths = 0.6, flierprops=green_diamond);setBoxColors(bp)
bp = ax.boxplot(gpt3, notch=nt, patch_artist=True, positions = [18.5], widths = 0.6, flierprops=green_diamond);
bp['boxes'][0].set(facecolor='#CCCCCC')
bp['medians'][0].set(color='black')

ax.grid(color='gray', linestyle=':')
ax.set_ylabel('dPWV (mm)')
#ax.set_title('Estimated GPS-PWV with different Tm models',fontweight="bold")
plt.xticks([1,4,7,10,13,16,18.5], ['AIRS','ERA5','Suwantong','Bevis','Mendes','Schueler','GPT3'])
plt.xticks(rotation=20)

# draw temporary red and blue lines and use them to create a legend
hB, = ax.plot([1,1],'-', color='pink',     linewidth=7.0)
hR, = ax.plot([1,1],'-', color='lightblue',linewidth=7.0)
ax.legend((hB, hR),('Ts obtained from AWS', 'Ts obtained from GPT3'),loc=4,prop={'size':20})
hB.set_visible(False)
hR.set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(dir_app, 'dPWV_boxPlot.png'), dpi=300)
plt.show() 
