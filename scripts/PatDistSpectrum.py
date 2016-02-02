import sys
import numpy as np
from sys import *
import os
import time
import dendropy
from dendropy.calculate import treemeasure


#Ensure all of the required components were retrieved prior to running program
try:
    output = sys.argv[1]
    if "/" in str(output):
        outputPath = str(sys.argv[1]).rsplit("/",1)[0]+"/"
    else:
        outputPath = ""
    TreeFile = sys.argv[2]
    if "/" in str(TreeFile):
        TreeFileSh = str(TreeFile).rsplit("/",1)[1]
        filePath = str(sys.argv[2]).rsplit("/",1)[0]+"/"
    else:
        TreeFileSh = str(TreeFile)
        filePath = ""
    infileSp = str(TreeFileSh).rsplit(".",1)
    infilePrefix = infileSp[0]
    PatDistOutput = outputPath+infilePrefix+".IntraGenotype.PatDist.txt"
    cutoff = float(sys.argv[3])
    length = int(sys.argv[4])
    start = int(sys.argv[5])-1
except KeyboardInterrupt:
    try:
        os.remove(PatDistOutput)
    except:
        pass
except:
    print "Usage: python PatDistSpectrum.py infile idLen idStart"
    print "infile: Tree file"
    print "idLen: Number of characters for unique identifier"
    print "idStart: Starting location of unique identifier"
    sys.exit(1)

#Generate pairwise pat distances for all variants below the specified threshold 
def dendro(infile, dendroOut):
            PatDistOutput = open(dendroOut,'w')
            tree = dendropy.Tree.get_from_path(infile, "newick")
            pdm = treemeasure.PatristicDistanceMatrix(tree)
            suppress_internal_node_taxa=True
            for i, t1 in enumerate(tree.taxon_namespace):
                for t2 in tree.taxon_namespace[i+1:]:
                    if pdm(t1,t2) < float(cutoff):
                        idSp = t1.label.split(" ")
                        if idSp >1:
                            t1.label = (t1.label.replace(" ","_"))
                            t2.label =(t2.label.replace(" ","_"))
                            PatDistOutput.write("%s,%s,%s,%s" % (t1.label, t2.label, pdm(t1, t2),str(pdm.mrca(t1, t2).label) +"\n"))
                        else:
                            PatDistOutput.write("%s,%s,%s,%s" % (t1.label, t2.label, pdm(t1, t2),str(pdm.mrca(t1, t2).label) +"\n"))
                            
#Store pat dist in a dictionary for subsequent percentile calculations
def PatDistAppend(line):
    linesplit = line.split(",")
    node = linesplit[0]
    mate = linesplit[1]
    patdist = linesplit[2].rstrip("\n")
    nodeshort = node[start:length]
    mateshort = mate[start:length]
    comb = nodeshort+"__"+mateshort
    comb2 = mateshort+"__"+nodeshort
    shortdouble = nodeshort+"__"+mateshort
    if nodeshort != mateshort:      
        if not comb2 in dictInterPatDist:
            if not comb in dictInterPatDist:
                dictInterPatDist[comb]=[]
            dictInterPatDist[comb].append(float(patdist))
        else:
            dictInterPatDist[comb2].append(float(patdist))
    else:
        if not shortdouble in dictInterPatDist:
            dictInterPatDist[shortdouble]=[]
        dictInterPatDist[shortdouble].append(float(patdist))
    return dictInterPatDist
 
#Calculate pat dist percentiles   
def Spectrum():
    Percentiles={}
    for k,v in dictInterPatDist.iteritems():
        dictInterPatDist[k].sort()
        a = np.asarray(dictInterPatDist[k])
        Q0 = np.percentile(a, 0)
        Q1 = np.percentile(a, 1)
        Q5 = np.percentile(a, 5)
        Q10 = np.percentile(a, 10)
        Q20 = np.percentile(a, 20)
        Q25 = np.percentile(a, 25)
        Q30 = np.percentile(a, 30)
        Q35 = np.percentile(a, 35)
        Q40 = np.percentile(a, 40)
        Q45 = np.percentile(a, 45)
        Q50 = np.percentile(a, 50)
        Q75 = np.percentile(a, 75)
        Q90 = np.percentile(a, 90)
        Q99 = np.percentile(a, 99)
        Q100 = np.percentile(a, 100)
        if not k in Percentiles:
            Percentiles[k]=[]
        Percentiles[k].extend([Q0,Q1,Q5,Q10,Q20,Q25,Q30,Q35,Q40,Q45,Q50,Q75,Q90,Q99,Q100])
    return Percentiles
                
if __name__ == '__main__':
    
    if not os.path.isfile(PatDistOutput):
        dendro(TreeFile,PatDistOutput)
    else:
        print "Pairwise patristic distances extracted already, calculating percentiles..."
    dictInterPatDist={}
    outfile=open(outputPath+infilePrefix+'.IntraGenotype.patdistspectrum.txt','w')
    outfile.write("Samples"+"\t"+"0"+"\t"+"1"+"\t"+"5"+"\t"+"10"+"\t"+"20"+"\t"+"25"+"\t"+"30"+"\t"+"35"+"\t"+"40"+"\t"+"45"+"\t"+"50"+"\t"+"75"+"\t"+"90"+"\t"+"99"+"\t"+"100"+"\n")
    with open(PatDistOutput) as f:
        for line in f:
            PatDistAppend(line.rstrip("\n"))
        PDspectrum = Spectrum()
        for k,v in PDspectrum.iteritems():
            outfile.write(k +"\t"+str(v).replace("[","").replace("]","").replace(",","\t")+"\n")