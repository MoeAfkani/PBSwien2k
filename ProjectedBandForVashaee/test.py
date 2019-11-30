import matplotlib.pyplot as plt
fermi=-0.1676706872  #fermi energy of bulk
Ene = []
k=[]
bandfile= open("graphene.spaghetti_ene")
for line in bandfile:
    if line[0]=="#" or 'bandindex:' in line :
        Ene.append([])
        continue
    temp=[float(x) for x in line.split()]
    try:
        Ene[-1].append(temp[-1])
    except: None
    try:
        if temp[-2] not in k:
            k.append(temp[-2])
    except: None
for band in Ene:
    plt.plot(k,band )
plt.show()