import sys
import itertools
from ete2 import Tree
import dendropy
from dendropy.calculate import treemeasure
import os

class ClusterIdentification(object):
    
    def __init__(self):
        self.PercentileThreshold = {}
        self.dictSharedReads = {}
        self.dictClusters = {}
        self.monoFinalRes = []
        self.count = 0
        self.SerialNodes = {}
        self.t = Tree(TreeFile)
        self.nodesRemoved = []
        self.nodecheck = []
    
    ##The Split method identifies the percentile threshold for each sample from the results of PatDistSpectrum.py
    ##This threshold is determined from user input in the command line
    def Split(self,infile):
        
        Percentiles = {'0':'1','1':'2','5':'3','10':'4','20':'5','25':'6','30':'7','35':'8','40':'9','45':'10','50':'11','75':'12','90':'13','99':'14','100':'15'}
        with open(Spectrum, "r") as file1:
            for line in file1:
                if not "Samples" in line:
                    linerep = line.replace(" ","")
                    if percentile in Percentiles:
                        cutoff = Percentiles[percentile]    
                    else:
                        sys.stdout.write("Please specify the percentile as a number (0,1,5,10,20,25,30,35,40,50,75,90,100)")
                        sys.exit(1) 
                    linesp = linerep.rstrip("\n").split("\t")
                    nodes = linesp[0]
                    nodesSp = nodes.split("__")
                    if nodesSp[0] == nodesSp[1]:
                        combNode = nodesSp[0]+"__"+nodesSp[1]
                        self.PercentileThreshold[combNode] = linesp[int(cutoff)]
        
        return self.PercentileThreshold
    
    #Identifies all variants passing the threshold defined in the Split method
    def variantCollection(self):
        
        PatDistSpec = self.Split(Spectrum)
        
        with open(PatDistOutput, "r") as file2:
            for line in file2:   
                linesp=line.split(",")
                node = linesp[0]
                nodeshort=node[idStart:idLen]
                nodedouble=nodeshort+"__"+nodeshort
                mateshort = linesp[1][idStart:idLen]
                matedouble = mateshort+"__"+mateshort
                mate = linesp[1]
                patdist=(linesp[2])
                support = linesp[3].rstrip("\n")
                comb=nodeshort+"__"+mateshort
                comb2=mateshort+"__"+nodeshort
                #Store variants that are below the respective pat dist threshold defined in PatDistSpec
                if nodeshort != mateshort:
                    if float(PatDistSpec[nodedouble]) <= float(PatDistSpec[matedouble]):
                        target = float(PatDistSpec[nodedouble])
                    else:
                        target = float(PatDistSpec[matedouble])
                    if float(patdist)<=float(target):
                        if str(supportInput) == "PASS":
                            if not comb2 in self.dictSharedReads:
                                if not comb in self.dictSharedReads:
                                    self.dictSharedReads[comb]=[]
                                if not node in self.dictSharedReads[comb]:
                                    self.dictSharedReads[comb].append(node)
                                if not mate in self.dictSharedReads[comb]:
                                    self.dictSharedReads[comb].append(mate)
                            else:
                                if not node in self.dictSharedReads[comb2]:
                                    self.dictSharedReads[comb2].append(node)
                                if not mate in self.dictSharedReads[comb2]:
                                    self.dictSharedReads[comb2].append(mate)
                        elif not support == "None":
                            if float(supportInput) <= float(support):
                                if not comb2 in self.dictSharedReads:
                                    if not comb in self.dictSharedReads:
                                        self.dictSharedReads[comb]=[]
                                    if not node in self.dictSharedReads[comb]:
                                        self.dictSharedReads[comb].append(node)
                                    if not mate in self.dictSharedReads[comb]:
                                        self.dictSharedReads[comb].append(mate)
                                else:
                                    if not node in self.dictSharedReads[comb2]:
                                        self.dictSharedReads[comb2].append(node)
                                    if not mate in self.dictSharedReads[comb2]:
                                        self.dictSharedReads[comb2].append(mate)
                                        
                                        
    ##Identifying potential outliers is optional (-oR flag from the command line )
    #Based on the retrieve common ancestor function, it identifies outliers as those which contain < 3 intra-variants associated with a given sample
    def PhylyOutlierRem(self, n, node1, node2,OutlierFile,idStart,idLen):
        PhyloOutliers = {}
        ancestorList = []
        ancshort = []
        node1short = node1[idStart:idLen]
        node2short = node2[idStart:idLen]
        nodecomb = node1short+"__"+node2short
        nodecombRev = node2short+"__"+node2short
        
        if not nodecomb or not nodecombRev in self.monoFinalRes:
            if not node1 in self.nodesRemoved:
                if not node2 in self.nodesRemoved:
                    #Collect all common ancestors for each pair of variants
                    ancestor = self.t.get_common_ancestor(n)
                    for i in ancestor:
                        ancestorList.append(i.name)
                        ancestorShort = i.name[idStart:idLen]
                        if not ancestorShort in PhyloOutliers:
                            PhyloOutliers[ancestorShort] = []
                        PhyloOutliers[ancestorShort].append(1)
                    
                    #Sum the variants for each sample, if < 3, store variant as outlier
                    for k,v in PhyloOutliers.iteritems():
                        vsum = sum(v)
                        if vsum < 3:
                            for item in ancestorList:
                                if item[idStart:idLen] == k:
                                    if not item in self.nodesRemoved:
                                        if node1short in self.SerialNodes:
                                            if not node2short in self.SerialNodes[node1short]:
                                                ancestorList.remove(item)
                                                self.nodesRemoved.append(item)
                                        elif node2short in self.SerialNodes:
                                            if not node1short in self.SerialNodes[node2short]:
                                                ancestorList.remove(item)
                                                self.nodesRemoved.append(item)
                                        else:
                                            ancestorList.remove(item)
                                            self.nodesRemoved.append(item)
            for i in self.nodesRemoved:
                if not i in self.nodecheck:
                    try:
                        item = self.t.search_nodes(name=item)[0]
                        i.delete()
                        self.nodecheck.append(i)
                    except:
                        pass
        
        
        return ancestorList
                    
    #Create all combinations of intrahost sample identifiers for each respective sequential sample set
    #These results are used to assist in PhylyOutlierRem
    def intraComb(self,infile):
        
        with open(IntraFile) as f:
            for line in f:
                line = line.rstrip("\n")
                linesp = line.split(",")
                
                length = len(linesp)
                comb = int(length)
                for i in linesp:
                    self.SerialNodes[i]= []
                    for pair in itertools.combinations(linesp,2):
                        
                        for item in pair:
                            #if not item in SerialDict.keys():
                                if i != item:
                                    if not item in self.SerialNodes[i]:
                                        self.SerialNodes[i].append(item)
        
        return self.SerialNodes
    
    #First step of merging overlapping pairs of connected samples
    def ClusterKeys(self,values,node1,node2):
        
        
        if node1 in [x for v in values for x in v if type(v)==list] or node1 in values:
            if not node2 in [x for v in values for x in v if type(v)==list] or node2 in values:
                for k,v in self.dictClusters.iteritems():
                    if node1 in self.dictClusters[k]:
                        self.dictClusters[k].append(node2)
                                
        if node2 in [x for v in values for x in v if type(v)==list] or node2 in values:
            if not node1 in [x for v in values for x in v if type(v)==list] or node1 in values:
                for k,v in self.dictClusters.iteritems():
                    if node2 in self.dictClusters[k]:
                        self.dictClusters[k].append(node1)
                        
        if node1 in [x for v in values for x in v if type(v)==list] or node1 in values:
            if node2 in [x for v in values for x in v if type(v)==list] or node2 in values:
                for k,v in self.dictClusters.iteritems():
                    if node1 in self.dictClusters[k]:
                        self.dictClusters[k].append(node2)
                        self.dictClusters[k].append(node1)
                        
                    if node2 in self.dictClusters[k]:
                        self.dictClusters[k].append(node2)
                        self.dictClusters[k].append(node1)
                        
        if not node1 in [x for v in values for x in v if type(v)==list] or node1 in values:
            if not node2 in [x for v in values for x in v if type(v)==list] or node2 in values:
                self.count+=1
                if not self.count in self.dictClusters:
                    self.dictClusters[self.count]=[]
                self.dictClusters[self.count].append(node1)
                self.dictClusters[self.count].append(node2.rstrip("\n"))
                 
    #Second step of merging overlapping pairs of connected samples
    def ClusterKeys2(self,dictClusters):
        Clustvals = {}
        sysvers = str(sys.version_info[0])+"."+str(sys.version_info[1])
        if float(sysvers) == 2.7:
            ##For python 2.7
            Clustvals = {k: set(val) for k, val in self.dictClusters.items()}
        elif float(sysvers) == 2.6:
            ##For python 2.6
            Clustvals = dict((k, val) for (k, val) in self.dictClusters.items())
        merged = set()
        srt = sorted(self.dictClusters.keys())
        srt2 = srt[:]
        for key in srt:
            for k in srt2:       
                if not k == key:
                    if Clustvals[k].intersection(self.dictClusters[key]) and key not in merged:
                        merged.add(k)
                        self.dictClusters[key] = list(Clustvals[k].union(self.dictClusters[key]))
                        srt2.remove(k)
        for k in merged:    
            del self.dictClusters[k]
        try:
            if len(self.dictClusters) > 0:
                del self.dictClusters[0]
        except:
            pass
        ValLengths=[]
        ItemNumber=[]
        
        for k,v in self.dictClusters.iteritems():
            ValLengths.append(int(len(set(v))))
            for i in v:
                if not i in ItemNumber:
                    ItemNumber.append(i)
        ValLengths[:]=[]
        for k,v in self.dictClusters.iteritems():          
            vset = set(v)    
            v[:]=[]
            vset = list(vset)
            self.dictClusters[k]=str(vset)
        return self.dictClusters
    
    #Retrieve common ancestors
    def CommonAncestor(self,nodes):
        ancestors=[]
        ancestor = self.t.get_common_ancestor(nodes)
        for i in ancestor:
            ancestors.append(i.name)
        return ancestors
    
    #Prepare a list of lists for parsing       
    def NodeParser(self,nodes):
        listitem = []
        vsp = str(nodes).replace("[['","").replace("']]","").replace("['","").replace("']","").replace('"','').split("', '")
        for i in vsp:
            listitem.append(i)
        
        return listitem
        
    #Identify poly- , para-, and monophyletic pairs of variants
    def CheckMono(self,ncomb,PhyloVarRemoval,Rejects,monoFinal):
        monoResult = str(self.t.check_monophyly(values=PhyloVarRemoval, ignore_missing=True,target_attr="name"))
        monoResultSp = monoResult.split(",")
        mR = monoResultSp[1].replace("'","").replace(")","").replace(" ","")                                   
        #print monoResult
        if 'monophyletic' in mR:
            if not ncomb in monoFinal:
                monoFinal[ncomb]=[]
            monoFinal[ncomb].append(mR)
            if not ncomb in self.monoFinalRes:
                
                self.monoFinalRes.append(ncomb)
            
            return True
        elif 'paraphyletic' in mR:
            if not ncomb in monoFinal:
                monoFinal[ncomb]=[]
            monoFinal[ncomb].append(mR)
            if not ncomb in self.monoFinalRes:
                
                self.monoFinalRes.append(ncomb)
        elif not ncomb in Rejects:
            
            Rejects.append(ncomb)
    
    ##Analysis identifies all ancestors to variants passing the required percentile thresholds
    #Following the removal of outliers, it parses through every combination of these variants to determine whether the pair is polyphyletic or not
    def variantAnalysis(self):
        monoFinal = {}
        self.variantCollection()
        if outlierFlag == "TRUE":
                OutlierFile = open(outputPath+TreeShort+"."+percentile+"."+supportInput+".Outlier.txt",'w')
            
        try:
            self.intraComb(IntraFile)
        except:
            pass
        Rejects = []
        
        for k,v in self.dictSharedReads.iteritems():
            
                x=0
                count = 0
                ksp = k.split("__")
                krev = ksp[1]+"__"+ksp[0]
                FinalList = []
                clusters = self.NodeParser(self.dictSharedReads[k])
                
                                    
                for pair in itertools.combinations(clusters,2):
                    n = list(pair)
                    node1 = n[0]
                    node2 = n[1]
                    if not node1 in self.nodesRemoved:
                        if not node2 in self.nodesRemoved:
                            nodeList = []
                            node1short = str(pair)[(idStart+2):(idLen+2)]
                            pairSp = str(pair).split(",")
                            node2short = pairSp[1].replace(" ","").replace("'","")[idStart:idLen]
                            nshort = [node1short,node2short]
                            ncomb = node1short+"__"+node2short
                            ncombRev = node2short+"__"+node1short
                            if node1short != node2short:
                                if not ncomb or not ncombRev in self.monoFinalRes:
                                    if outlierFlag == "TRUE":
                                        PhyloVarRemoval = self.PhylyOutlierRem(n,node1,node2, OutlierFile,idStart,idLen)
                                    else:
                                        PhyloVarRemoval = self.CommonAncestor(n)
                                    
                                    for i in PhyloVarRemoval:                    
                                        node = i[idStart:idLen]
                                        if not node in nodeList:
                                            nodeList.append(node)
                            if not node1 in self.nodesRemoved:
                                if not node2 in self.nodesRemoved:            
                                    if len(nodeList) == 2:
                                        if not ncomb or not ncombRev in self.monoFinalRes:
                                            
                                            if self.CheckMono(ncomb,PhyloVarRemoval,Rejects,monoFinal):
                                                break
                                    elif len(nodeList)>2:
                                        monoPos = 0
                                        lengthNode = len(nodeList)
                                        flag = 0
                                        nodeCheck = 0
                                        nodeRemoval = []
                                        
                                        for i in set(nodeList):
                                            if not i in nshort:
                                                if not i in self.SerialNodes:
                                                    nodeRemoval.append(i)
                                                    flag = 1
                                        if flag == 0:
                                            if len(PhyloVarRemoval)>1:
                                                if not ncomb or not ncombRev in self.monoFinalRes:
                                                    if self.CheckMono(ncomb,PhyloVarRemoval,Rejects,monoFinal):
                                                        break           
                                        else:
                                            for item in nodeRemoval:
                                                nodeRemovalShort = item[idStart:idLen]
                                                if not nodeRemovalShort+"__"+node1short in self.dictSharedReads.keys():
                                                    if not node1short+"__"+nodeRemovalShort in self.dictSharedReads.keys():
                                                        if not node2short+"__"+nodeRemovalShort in self.dictSharedReads.keys():
                                                            if not nodeRemovalShort+"__"+node2short in self.dictSharedReads.keys():
                                                                for i in PhyloVarRemoval:
                                                                    nodeShort = i[idStart:idLen]
                                                                    if nodeShort in nodeRemoval:
                                                                        PhyloVarRemoval.remove(i)
                                            if len(PhyloVarRemoval) > 1:
                                                if not ncomb in self.monoFinalRes:
                                                    if self.CheckMono(ncomb,PhyloVarRemoval,Rejects,monoFinal):
                                                        break
                                            flag = 0
        if outlierFlag == "TRUE":
            for i in set(self.nodesRemoved):
                OutlierFile.write('%s\n' % i)
                           
        self.dictClusters ={}
        self.count = 0
        for i in self.monoFinalRes:
            if not 'polyphyletic' in i:
                if not self.count in self.dictClusters:
                    self.dictClusters[self.count]=[]        
                values = self.dictClusters.values()
                
                isp = i.split("__")
                node1 = isp[0]
                node2 = isp[1].split("\t")[0]
                ClusterKeys = self.ClusterKeys(values,node1,node2)
        try:
            FinalClustering = self.ClusterKeys2(ClusterKeys)
        except:
            print "WARNING: Patristic Distance Data files may be empty..."
            sys.exit(1)
        for k,v in monoFinal.iteritems():
            print k +"\t"+str(v)
        print "Clusters that are polyphyletic: "+str(Rejects)
        return FinalClustering

#Organizes filepathways and ensures files can be located
Spectrum = sys.argv[1]
PatDistOutput = sys.argv[2]
percentile = sys.argv[3]
supportInput = sys.argv[4].upper()
IntraFile = sys.argv[5]
TreeFile = sys.argv[6]
idLen = int(sys.argv[7])
idStart = int(sys.argv[8])-1
outlierFlag = sys.argv[9].upper()
output = sys.argv[10]

if "/" in output:
    outputPath = str(sys.argv[10]).rsplit("/",1)[0]+"/"
else:
    outputPath = ""
if "/" in str(TreeFile):
    filePath = str(sys.argv[6]).rsplit("/",1)[0]+"/"
    TreeShort = str(sys.argv[6]).rsplit("/",1)[1]
else:
    filePath = ""
    TreeShort = str(TreeFile).rsplit(".", 1)[0]
if os.path.isfile(TreeFile):
    TreeFile = TreeFile
    filePath = str(sys.argv[6]).rsplit("/",1)[0]+"/"
elif os.path.isfile(filePath+TreeShort):
    TreeFile = filePath+TreeShort  
elif os.path.isfile(outputPath+TreeShort):
    TreeFile = outputPath+TreeShort
else:
    print "WARNING: Input tree file could not be located..."
    sys.exit(1)

if IntraFile != 'False':
    if not os.path.isfile(IntraFile):
        
        print "WARNING: Intrahost data could not be located..."
        sys.exit(1)
    
if __name__ == '__main__':  
    
    ClusterAnalysis = ClusterIdentification()
    
    Clusters = ClusterAnalysis.variantAnalysis()

    ClusterOutputFile = open(output,'w')
    
    for k,v in Clusters.iteritems():
        ClusterOutputFile.write(str(k) +"\t"+ v+"\n")