from ete2 import PhyloTree
from ete2 import Tree
import dendropy
from dendropy.calculate import treemeasure
import sys
import itertools
import os

class SangerTreeParsing(object):
    def __init__(self):
        self.monoFinal = []
        self.monoPairs = []
        self.Clusters = []
        self.dictClusters = {}
        self.SerialNodes = {}
    #If tree support is not required, this retrieves variants with pat distances below the specified threshold
    def TreeParsePatDist(self,newickFile,cutoff):
        tree = dendropy.Tree.get_from_path(newickFile, "newick")
        pdm = treemeasure.PatristicDistanceMatrix(tree)
        pdfile = open(outputPath+TreeShort+"."+str(cutoff)+"."+str(support)+".patdist.txt",'w')
        PatDist = []
        suppress_internal_node_taxa=True
        for i, t1 in enumerate(tree.taxon_namespace):
            for t2 in tree.taxon_namespace[i+1:]:
                if pdm(t1,t2) < float(cutoff):
                    pdfile.write("%s,%s,%s,%s" % (t1.label, t2.label, pdm(t1, t2),str(pdm.mrca(t1, t2).label) )+"\n")
                    pdlist = []
                    idSp = t1.label.split(" ")
                    #this accounts for underscores on UIDs as they are split with dendropy
                    if idSp >1:
                        pdlist.append(t1.label.replace(" ","_"))
                        pdlist.append(t2.label.replace(" ","_"))
                        PatDist.append(pdlist)
                    else:
                        pdlist.append(t1.label)
                        pdlist.append(t2.label)
                        PatDist.append(pdlist)
        return PatDist
    #This retrieves variants with pat distances below the specified pat dist and tree support thresholds
    def TreeParsePatDistSupport(self,newickFile,cutoff,support):
        tree = dendropy.Tree.get_from_path(newickFile, "newick")
        pdm = treemeasure.PatristicDistanceMatrix(tree)
        pdfile = open(outputPath+TreeShort+"."+str(cutoff)+"."+str(support)+".patdist.txt",'w')
        PatDist = []
        suppress_internal_node_taxa=True
        for i, t1 in enumerate(tree.taxon_namespace):
            for t2 in tree.taxon_namespace[i+1:]:
                if pdm(t1,t2) < float(cutoff):
                    if not (str(pdm.mrca(t1, t2).label)) == "None":
                        if float(str(pdm.mrca(t1, t2).label)) > float(support):
                            pdfile.write("%s,%s,%s,%s" % (t1.label, t2.label, pdm(t1, t2),str(pdm.mrca(t1, t2).label) )+"\n")
                            pdlist = []
                            idSp = t1.label.split(" ")
                            if idSp >1:
                                #print t1.label
                                pdlist.append(t1.label.replace(" ","_"))
                                pdlist.append(t2.label.replace(" ","_"))
                                PatDist.append(pdlist)
                            else:
                                pdlist.append(t1.label)
                                pdlist.append(t2.label)
                                PatDist.append(pdlist)
                    elif pdm(t1,t2) < float(cutoff):
                        pdfile.write("%s,%s,%s,%s" % (t1.label, t2.label, pdm(t1, t2),str(pdm.mrca(t1, t2).label) )+"\n")
                        pdlist = []
                        idSp = t1.label.split(" ")
                        if idSp >1:
                            pdlist.append(t1.label.replace(" ","_"))
                            pdlist.append(t2.label.replace(" ","_"))
                            PatDist.append(pdlist)
                        else:
                            pdlist.append(t1.label)
                            pdlist.append(t2.label)
                            PatDist.append(pdlist)
        return PatDist

    #Cluster samples passing specified thresholds
    def Cluster(self,ResultList):
        self.dictClusters = {}
        count = 0
        for lst in ResultList:
            for pair in itertools.combinations(lst,2):
                node1 = pair[0]
                node2 = pair[1]
                if not count in self.dictClusters:
                    self.dictClusters[count]=[]
                values = self.dictClusters.values()
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
                        count+=1
                        if not count in self.dictClusters:
                            self.dictClusters[count]=[]
                        self.dictClusters[count].append(node1)
                        self.dictClusters[count].append(node2.rstrip("\n"))
        Clustvals = {}
        Clustvals = {k: set(val) for k, val in self.dictClusters.items()}
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
        del self.dictClusters[0]
        ValLengths=[]
        ItemNumber=[]
        for k,v in self.dictClusters.iteritems():
            ValLengths.append(int(len(set(v))))
            for i in v:
                if not i in ItemNumber:
                    ItemNumber.append(i)
        #print ItemNumber
        #print "The number samples that clustered is " +"\t"+str(len(ItemNumber))
        #print "The number of clusters is "+"\t"+str(len(self.dictClusters))
        #print "The smallest cluster size is " +"\t"+str(sorted(ValLengths)[0])
        #print "The largest cluster size is " +"\t"+str(sorted(ValLengths)[-1])
        ValLengths[:]=[]
        self.Clusters = []
        for k,v in self.dictClusters.iteritems():
            self.Clusters.append(str(set(v)).replace("set(","").replace(")",""))
        return self.Clusters
    
    
    def CheckMonophyly(self,PDlist):
        t = Tree(filePath+TreeFile)
        monoShort=[]
        x=0
        for item in PDlist:
            cluster=[]
            pairL=[]
            flag = 0
            y=0
            clusterRaw=str(item).replace("[","").replace("]","").replace('"',"").replace(" ","").replace("'",'').split(',')
            for i in clusterRaw:
                if not i in cluster:
                    cluster.append(i)
            monoResult = str(t.check_monophyly(values=cluster, ignore_missing=True,target_attr="name"))
            #Identify poly- , para-, and monophyletic relationships for clusters
            if 'monophyletic' in monoResult:
                for pair in itertools.combinations(cluster,2):
                    m = list(pair)
                    #combM = str(m).replace("', '",",").replace("[","").replace("]","").replace("'","")
                    if not m in self.monoPairs:
                        self.monoPairs.append(m)
                    
            elif 'paraphyletic' in monoResult:
                for pair in itertools.combinations(cluster,2):
                    m = list(pair)
                    #combM = str(m).replace("', '",",").replace("[","").replace("]","").replace("'","")
                    if not m in self.monoPairs:
                        self.monoPairs.append(m)
            else:
                cluster2 = []
                for pair in itertools.combinations(cluster,2):
                    n = list(pair)  
                    monoResult = str(t.check_monophyly(values=n, ignore_missing=True,target_attr="name"))
                    if not 'polyphyletic' in monoResult:
                        if not n in cluster:
                            cluster.append(n)
                            for i in n:
                                if i in cluster:
                                    cluster.remove(i)
                        if not n in self.monoPairs:
                            self.monoPairs.append(n)
                #Breaksdown large clusters to identify poly- , para-, and monophyletic sub-clusters
                while y < 2:
                    y+=1
                    for pair in itertools.combinations(cluster,2):
                        pairL = list(pair)
                        cluster2 = []
                        if type(pairL[0]) is list:
                            if type(pairL[1]) is list:
                                cluster2 = pairL[0]+pairL[1]
                            else:
                                for i in pairL[0]:
                                    if not i in cluster2:
                                        cluster2.append(i)
                                if not pairL[1] in cluster2:
                                    cluster2.append(pairL[1])
                        elif type(pairL[1]) is list:
                            for i in pairL[1]:
                                if not i in cluster2:
                                    cluster2.append(i)
                            if not type(pairL[0]) is list:
                                if not pairL[0] in cluster2:
                                    cluster2.append(pairL[0])
                        else:
                            if not pairL[0] in cluster2:
                                cluster2.append(pairL[0])
                            if not pairL[1] in cluster2:
                                cluster2.append(pairL[1])
                        monoResult = str(t.check_monophyly(values=cluster2, ignore_missing=True,target_attr="name"))
                        if not 'polyphyletic' in monoResult:
                            x+=1
                            if not cluster2 in cluster:
                                cluster.append(cluster2)
                            for item in cluster2:
                                for i in item:
                                    if i in cluster:
                                        cluster.remove(i)
                            if not cluster2 in self.monoPairs:
                                self.monoPairs.append(cluster2)
        return self.monoPairs
#Ensure all of the required files and settings were input correctly
try:
    TreeFile = sys.argv[1]
    cutoff = sys.argv[2]
    cutoff = float(str(cutoff).replace("'","").replace('"',""))
    support = sys.argv[3]
    output = sys.argv[4]
    support = support.upper()
    support = str(support).replace("'","").replace('"',"")
    
    
except IndexError:
    print "Usage: python ete.ConsensusSequencingClusterIDFin5.py newickFile PDThreshold TSThreshold output"
    print "newickFile: any newick tree file"
    print "PDThreshold: Patristic Distance Threshold: ex '0.03', default '0.3'"
    print "TSThreshold: Bootstrap / Tree Support Threshold: ex '0.95' or 'pass' if you want to identify your own"
    print "Output: output filename and location"
    sys.exit(1)
    
    
if "/" in output:
    outputPath = str(sys.argv[4]).rsplit("/",1)[0]+"/"
else:
    outputPath = ""
if "/" in str(TreeFile):
    filePath = str(sys.argv[1]).rsplit("/",1)[0]+"/"
    TreeShort = str(sys.argv[1]).rsplit("/",1)[1]
else:
    filePath = ""
    TreeShort = str(TreeFile).rsplit(".", 1)[0]
if os.path.isfile(TreeFile):
    TreeFile = TreeFile
    filePath = ""
elif os.path.isfile(filePath+TreeFile):
    TreeFile = filePath+TreeFile  
elif os.path.isfile(outputPath+TreeFile):
    TreeFile = outputPath+TreeFile
else:
    print "WARNING: Input tree file could be located..."
    sys.exit(1)

if __name__ == '__main__':
    Cluster1 = SangerTreeParsing()
    
    if support != "PASS":
        PD = Cluster1.TreeParsePatDistSupport(TreeFile,cutoff,support)
    else:
        PD = Cluster1.TreeParsePatDist(TreeFile,cutoff)
    
    Cluster2 = Cluster1.Cluster(PD)
    Cluster3 = Cluster1.CheckMonophyly(Cluster2)
    Cluster4 = Cluster1.Cluster(Cluster3)
    x=0
    ClusterOutputFile = open(output,'w')
    for i in Cluster4:
        x+=1
        ClusterOutputFile.write(str(x)+"\t"+i+"\n")