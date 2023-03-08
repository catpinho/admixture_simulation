import random
def convertspedtolist(pedfile):
    a=file(pedfile,"r").read().split("\n")
    genos=[]
    for i in range(len(a)):
        if a[i]!="":
            ind=[]
            ind.append(a[i].split(" ")[1])
            ind.append(a[i].split(" ")[0])
            nl=len(a[i].split(" ")[6:])/2
            for l in range(nl):
                g=[a[i].split(" ")[2*l+6],a[i].split(" ")[2*l+7]]
                ind.append(g)
            genos.append(ind)
            
    return genos


#in the following function inds are the set of individuals on which to apply the function (result from convertpedtolist). 
#mdatadist is the distribution of missing data, contained in an external file organized as two columns (one with the number of missing loci, the other with the number of individuals missing those loci)
#maxgerror is the maximum frequency of loci with genotyping errors. Genotyping error always changes one of the alleles in the genotype, not both.

def getsallelesloc(inds):
    als=[]
    nloc=len(inds[0])-2
    for w in range(nloc):
        a=[]
        for i in inds:
            
            for al in i[2+w]:
                if al not in ["-9","0"]:
                    if al not in a:
                        a.append(al)
        als.append(a)
    return als
def throwsmissingwrong_geno_adm(inds,mdatadist,maxgerror):
    als=getsallelesloc(inds)
    nloc=len(inds[0])-2
    #print nloc
    allinds=range(len(inds))
    discardedinds=[]
    indstochoose=[]
    for x in allinds:
        if x not in discardedinds:
            indstochoose.append(x)
    ##print indstochoose
    mdata=file(mdatadist,"r").read().split("\n")
    for u in mdata:
        if u!="":
           # #print "indstochoose: "
           # #print indstochoose
            ninds=int(u.split("\t")[1])
            nmloc=int(u.split("\t")[0])
            ##print nmloc
            selinds=random.sample(indstochoose,ninds)
            ##print selinds
            for x in selinds:
            #    #print "x:" +str(x)
            #    #print "inds to choose: "+ str(indstochoose)
                i=inds[x]
            #    #print "individual: "+str(i)
            #    #print "inds before: "+str(inds)
                alloci=range(nloc)
             #   #print "alloci: "+str(len(alloci))
              
                missl=random.sample(alloci,nmloc)
                avloci=[]
                for loc2 in alloci:
                    if loc2 not in missl:
                        avloci.append(loc2)
                
             
                #print "avloci: "+str(len(avloci))
                for loc in missl:
                    i[loc+2]=["0","0"]
                 #   #print "nr of available loci "+str(len(avloci))
                    ##print nloc - nmloc
                maxnwrong=int(round(maxgerror*len(avloci),0))
                #print "maxnwrong: "+str(maxnwrong)
                wrong=random.sample(avloci,random.randint(0,maxnwrong))
                #print "wrong: "
                #print wrong
                for w in wrong:
                    alleles=als[w]
                    geno=i[w+2]
                    #print inds[x][w+2]
                    if len(alleles)==1:
                        #print "case1"

                        inds[x][w+2]=["1","2"]
                        
                    elif len(alleles)==2:
                        if geno[0]!=geno[1]:
                            #print "case2"
                            inds[x][w+2]=random.sample([["1","1"],["2","2"]],1)[0]
                           
                        elif geno[0]==geno[1]:
                            #print "case3"
                            inds[x][w+2]=alleles
                           
                    #print inds[x][w+2]
                ##print "inds after: "+str(inds)
                discardedinds.append(x)
            indstochoose=[]
            for indiv in allinds:
                if indiv not in discardedinds:
                    indstochoose.append(indiv)   
                
    return inds
