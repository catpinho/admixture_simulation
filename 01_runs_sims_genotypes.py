#this script uses functions contained in functions_sim_genotypes.py to set up a simulation of admixture between genotypes along several generations

from functions_sim_genotypes import *
import os
#NUMBER OF INDIVIDUALS OF SIMULATED POPULATION (needs to be larger than test population to avoid drift, but beware that large populations take longer to compute)
ninds=1000
#optional: allowed proportion of interspecific matings (all matings beyond this proportion will be either between pure individuals or between hybrids and any individual)
prop=0.028846
#optional: proportion of migration from pure individuals into contact zone (if applicable)
migprop=0.4
#number of different simulations to conduct
nsims=100
#proportion of species A in the real dataset
pc=0.818
#stop generation: number of generations to develop in each simulation
ngen=50

#THEN PICK APPROPRIATE CONDITIONS BY COMMENTING OUT UNDESIRED BLOCKS

#THERE ARE 2 BLOCKS RELATIVE TO THE TYPES OF MATING: RANDOM MATING + NONRANDOM MATING
#ONE OF THESE (AND ONLY ONE) NEEDS TO BE ACTIVE;

#THERE IS A 3RD BLOCK RELATIVE TO MIGRATION - COMMENT OUT OR UNCOMMENT DEPENDING ON NECESSITY

for sim in range(nsims):
    o=os.getcwd()
    #WINDOWS
#    os.system("md sim"+str(sim))
#    os.system("copy gen0* sim"+str(sim))
#    os.chdir(o+"\\sim"+str(sim))
    #UNIX
    os.system("mkdir sim"+str(sim))
    os.system("cp gen0* sim"+str(sim))
    os.chdir(o+"/sim"+str(sim))  
    for gen in range(ngen):
        out3=file("popcounts.txt","a")
        out=file("gen"+str(gen+1)+".ped","w")
        a=define_grupos("gen"+str(gen)+".2.Q",0.9,ninds)
        b=gruposnovos(a)
	#print b
        c=map(str,(map(len,b)))
        out3.write("gen"+str(gen))
        for x in c:
            out3.write("\t"+x)
        out3.write("\n")
        out3.close()
        i=0
        u=random.sample(b[0]+b[1]+b[2]+b[3]+b[4]+b[5],1)[0]
        nloc=(len(u.split(" "))-6)/2
        newlist=[]
##      #RANDOM MATING BLOCK (COMMENT OUT IF MATING IS NOT TO BE RANDOM)
##        for j in range(i,ninds):
##            reprF=random.sample(b[0]+b[2]+b[4],1)[0]
##            reprM=random.sample(b[1]+b[3]+b[5],1)[0]
##            newlist.append(hibrida([reprF,reprM],j))
        #END OF RANDOM MATING BLOCK
        
        #NONRANDOM MATING BLOCK(COMMENT OUT IF MATING IS TO BE RANDOM)
            
        #F1s 
        t=[len(b[0]),len(b[3])]
        w=[len(b[1]),len(b[2])]
        if 0 not in t:
            m1nF1=int(round(prop*sum(t),0))
            while i<m1nF1:
                a=random.sample(b[0],1)[0]
                d=random.sample(b[3],1)[0]
                newlist.append(hibrida([a,d],i))
                i+=1
                
        if 0 not in w:
            m2nF1=int(round(prop*sum(w),0))
            while i<m1nF1+m2nF1:
                be=random.sample(b[1],1)[0]
                c=random.sample(b[2],1)[0]

                newlist.append(hibrida([c,be],i))
                i+=1
        #OTHER MATINGS
        for j in range(i,ninds):
            reprF=random.sample(b[0]+b[2]+b[4],1)[0]
            
            if reprF in b[0]:
                reprM=random.sample(b[1]+b[5],1)[0]
                
            elif reprF in b[2]:
                reprM=random.sample(b[3]+b[5],1)[0]
                
            elif reprF in b[4]:
                reprM=random.sample(b[1]+b[3]+b[5],1)[0]
                
            #print "other "+str(j)+" "+str(cou)+" "+repra[4:20]+"x"+reprr[4:20]
            newlist.append(hibrida([reprF,reprM],j))

        #END OF NONRANDOM-MATING BLOCK

        #MIGRATION BLOCK (COMMENT OUT IF NO MIGRATION FROM OUTSIDE IS TO BE SIMULATED)

        rep=int(round(migprop*ninds,0))
        carb=int(round(rep*pc,0))
        boc=int(round(rep-rep*pc,0))
        indsout=random.sample(newlist,rep)
        #print len(newlist)
        #print len(indsout)
       
       
        for ind in indsout:
            if ind in newlist:
                newlist.remove(ind)
        for ind in newlist:
            if ind in indsout:
                newlist.remove(ind)
       
        for ca in range(carb):
            newlist.append("pop_1 PC_C_mgen"+str(gen)+"-"+str(ca)+" 0 0 0 0"+" 2"*nloc*2+"\n")
        for bo in range(boc):
            newlist.append("pop_1 PB_B_mgen"+str(gen)+"-"+str(bo)+" 0 0 0 0"+" 1"*nloc*2+"\n")
        #END OF MIGRATION BLOCK

        for indiv in newlist:
            out.write(indiv)
        out.close()
        out2=file("gen"+str(gen+1)+".2.Q","w")
        readfile=file("gen"+str(gen+1)+".ped","r")
        u=readfile.readline()
        while u!="":
            out2.write(str(countsalleles(u)[0])+" "+str(countsalleles(u)[1])+"\n")
            u=readfile.readline()
        out2.close()
    os.chdir(o)
    
    
    
    
    
    
    
