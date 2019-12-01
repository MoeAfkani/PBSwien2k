################################################################################
fermi=0 #fermi energy of bulk

################################################################################
import matplotlib.pyplot as plt
#from matplotlib.offsetbox import AnchoredText
#import matplotlib.ticker as ticker
#import pandas as pd
import matplotlib as mpl
#from matplotlib.lines import Line2D

myfont = {'family': 'Times New Roman',
        'size': 12,
}
#plt.rcParams["font.family"] = "serif"
#mpl.rc('text', usetex=True)

Ene = []
k=[]
dos=[[],[],[],[],[],[]]
dosE=[]
bandfile= open("MnTe-scf.spaghettiup_ene")
for line in bandfile:
    if line[0]=="#" or 'bandindex:' in line :
        Ene.append([])
        continue
    try:
        temp=[float(x) for x in line.split()]
    except ValueError: continue
    try:
        Ene[-1].append(temp[-1])
    except: None
    try:
        if temp[-2] not in k:
            k.append(temp[-2])
    except: None
    
dosfile1= open("MnTe-scf-dos.dos1evup")
dosfile2= open("MnTe-scf-dos.dos2evup")
for l1,l2 in zip(dosfile1,dosfile2):
    if l1[0]=="#" or l2[0]=="#":
        continue
    temp1=[float(x) for x in l1.split()]
    temp2=[float(x) for x in l2.split()]
    try:
        dosE.append(temp1[0])
        for i in range(len(dos)):
            dos[i].append(temp1[2+i]+temp2[2+i])
    except: None


    
# reduce e range

for e in Ene:
    if e==[] or abs(sum(e)/len(e)) >8 : Ene.remove(e)

Ene = [e for e in Ene if len(e)!=0]
Ene = [e for e in Ene if abs(sum(e)/len(e))<8]



for e in dosE:
    if abs(e)> 8:
        inx = dosE.index(e)
        for d in dos: d.pop(inx)
        dosE.pop(inx)


def color(e):
    inx = dosE.index(min(dosE, key=lambda x:abs(x-e)))
    try:
        tempTot = dos[1][inx] + dos[2][inx] + dos[5][inx]
        Sper = dos[1][inx]/tempTot
        Pper = dos[2][inx]/tempTot
        Dper = dos[5][inx]/tempTot
        col = [((Pper+2)/4,(Dper+2)/4,(Sper+2)/4)]
    except ZeroDivisionError: col= [(0,0,0)]
    return col

for band in Ene:
    for e in range(len(band)):
        plt.scatter(k[e],band[e],s=20 , color=color(band[e]) ,alpha=.9)
plt.show()


'''
###########################
colorbandup = [[] for _ in range(len(bandup[2]))]
for i in range(len(bandup[2])):
    for j in range(len(bandup[2][0])):
        indexenergy = ene_b.index(min(ene_b, key=lambda x:abs(x-bandup[2][i][j])))
        colorbandup[i].append((persentO_dosup[indexenergy]+persentCa_dosup[indexenergy],persentCs_dosup[indexenergy],persentN_dosup[indexenergy]+persentCa_dosup[indexenergy]))
    

###########################    
    
    
    
k_points=[bandup[1][0],bandup[1][200], bandup[1][400], bandup[1][600],bandup[1][700],bandup[1][800]]
fig, (bdw, dos, bup) = plt.subplots(1,3,figsize=(7,7), sharey='all',
                        gridspec_kw={'hspace': 0, 'wspace': 0})
E=(-3,4)
dos.set_xlim([-4.9,4.9])
dos.set_ylim(E)
dos.fill_between(tot_b_dosup,0,ene_b,facecolor='r')
dos.fill_between(tot_b_dosdn,0,ene_b,facecolor='r')
dos.axvline(x=0,c='black',ls='--')
dos.axhline(y=0,c='black',ls='--')
plt.setp((bup, bdw), 
         xticks=k_points,
         xticklabels=['W', 'L', '$\Gamma$', 'X', 'W', 'K'],
        yticks=range(E[0],E[1]+1))
plt.setp((dos), 
        yticks=range(E[0],E[1]+1))
dos.tick_params(
   axis='y',          # changes apply to the x-axis
   right='on',      # ticks along the bottom edge are off
   left='on',         # ticks along the top edge are off
   direction='in')
bdw.tick_params(
   axis='y',          # changes apply to the x-axis
   which='both',      # both major and minor ticks are affected
   right='off',      # ticks along the bottom edge are off
   left='on',         # ticks along the top edge are off
)
bup.set_xlim([0,bandup[1][-1]])
bdw.set_xlim([0,bandup[1][-1]])

for i in range(len(bandup[2])):
    for j in range(len(bandup[2][0])):
        bup.scatter(bandup[1][j], bandup[2][i][j], color=colorbandup[i][j],marker=3)

for i in range(len(banddw[2])):
    for j in range(len(banddw[2][0])):
        bdw.scatter(banddw[1][j], banddw[2][i][j], color=colorbanddw[i][j],marker=3)


for i in k_points[:-1]:  
    bup.axvline(x=i,c='black',lw=1)
    bdw.axvline(x=i,c='black',lw=1)
bup.axhline(y=0,c='black',ls='--')
bdw.axhline(y=0,c='black',ls='--')

legend_elements =[Line2D([0], [0], marker='o', color=(0,1,0), label='Cs',
                          markerfacecolor=(0,1,0), markersize=12),
		  Line2D([0], [0], marker='o', color=(1,0,1), label='Ca',
                          markerfacecolor=(1,0,1), markersize=12),
		  Line2D([0], [0], marker='o', color=(0,0,1), label='N',
                          markerfacecolor=(0,0,1), markersize=12),
		  Line2D([0], [0], marker='o', color=(1,0,0), label='O',
                          markerfacecolor=(1,0,0), markersize=12)]
bup.legend(handles=legend_elements)

bup.set_title('(c) Minority-spin')
bdw.set_title('(a) Majority-spin')
dos.set_title('(b) DOS (state/eV)')
bdw.set_ylabel('$E-E_F$ (eV)',fontsize=16)

fig.tight_layout()
plt.savefig('CsCaNObandos',format='jpeg',dpi=300)
'''