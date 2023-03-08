import random
#parents is a list of two admixture .ped lines, each corresponding to an individual
#the female is always the first element of the list
#out is an open file for writing the genotype
def hibrida(parents,n):
    mit=parents[0].split(" ")[1].split("_")[1]
    ind=[]
    
    nloc=(len(parents[0].split(" "))-6)/2
    for l in range(nloc):
        
        g=[]
        for p in parents:
            pspl=p.split(" ")
            g.append(random.sample([pspl[2*l+6],pspl[2*l+7]],1)[0])
        ind.append(g)
    
    string="pop_1 X_"+mit+"_"+str(n)+" 0 0 0 0"
    for x in ind:
        string+=" "+x[0]+" "+x[1]
    string+="\n"
    
    return string
def define_grupos(fich_admixture, threshold,ninds):
    #print threshold
    grupos=[[],[],[]]
    inputfile=file(fich_admixture.split("2.Q")[0]+"ped","r")
    adfile=file(fich_admixture,"r")
    a=0
    b=inputfile.readline().replace("\n","").replace("\r","")
    c=adfile.readline().replace("\n","").replace("\r","")
    while a<ninds:
        e=c.split(" ")
        if float(e[0])>threshold:
            grupos[0].append(b)
        else:
            if float(e[1]) < threshold:
                grupos[2].append(b)
        if float(e[1])>threshold:
            grupos[1].append(b)
        a+=1
        b=inputfile.readline().replace("\n","").replace("\r","")
        c=adfile.readline().replace("\n","").replace("\r","")
    return grupos
def divide_mf(grupo):
    m=[]
    f=random.sample(grupo,len(grupo)/2)
        
    for x in grupo:
        if x not in f:
            m.append(x)
    if len(grupo)%2==1:
        sl=random.sample(["m","f"],1)[0]
        if sl == "f":
            f.append(m[-1])
            m.remove(m[-1])
 
    return [f,m]

def gruposnovos(grupos):
    novosgrupos=[]
    for grup in grupos:
        a=divide_mf(grup)
        novosgrupos.append(a[0])
        novosgrupos.append(a[1])
    return novosgrupos
    
def countsalleles(ind):
    #print ind.split(" ")[1]
    indsp=ind.replace("\n","").split(" ")[6:]

    c=0
    tot=0
    for m in indsp:
        if m!="0":
            tot+=1
            if m=="1":
                c+=1
    return([c*1.0/tot,1-c*1.0/tot])
    
    
    


    





        

