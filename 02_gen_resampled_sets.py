#this script requires file "errors_mdata_functions.py" in the same folder.
#It resamples a set of simulated data sets so that it contains the same number of individuals as the real data set
#it also throws missing data across individuals according to a predefined distribution contained in a separate file (see errors_mdata_functions.py for details about this file's format)
#finally, it incorporates some genotyping errors based on a maximum acceptable error (which can be based on e.g. replicate testing)
#This script was originally devised to transform in such manner a series of simulated datasets contained in a variety of folders.
#The script needs to be placed in the root folder to those data-containing folders, and no other folders should be present

#NUMBER OF INDIVIDUALS:
ni = 115
#file with missing data distribution:
md = "md_counts_cz.txt"
#maximum error rate
mer = 0.00244
#OS type (comment out which does not apply):
#ost="unix"
ost="windows"
if ost=="windows":
    sep="\\"
else:
    sep="/"
import os
from errors_mdata_functions import *
a=os.getcwd()
b=os.listdir(a)
for i in b:
    if os.path.isdir(i) == True:
        c=os.listdir(a+sep+i)
        for sim in c:
            if ".ped" in sim:
                print sim
                #converts .ped file to a list of genotypes
                t=convertspedtolist(a+sep+i+sep+sim)
                #sample ni individuals at random from the global population
                x=random.sample(t,ni)

                #throws in missing data and genotyping errors

                u=throwsmissingwrong_geno_adm(x,md,mer)

                #writes new data file; outfile names can be changed if desired below
                out=file(a+sep+i+sep+sim.replace(".ped","_"+str(ni)+"mdge.ped"),"w")
                for ind in u:
                    out.write(ind[1]+" "+ind[0]+" 0 0 0 0")
                    for g in ind[2:]:
                        out.write(" "+g[0]+" "+g[1])
                    out.write("\n")
                out.close()
                             
                         
                         
                         
