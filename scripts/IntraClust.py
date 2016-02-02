#!/usr/bin/env python

from argparse import ArgumentParser
import sys
import subprocess as sp
import subprocess
import sys
import os
import datetime
current = datetime.datetime.now()
TimeStamp = current.strftime("%Y-%m-%d-%H.%M")

flag = 0


parser = ArgumentParser()

parser.add_argument("-i", "--input", dest="infile",required=True,
                    help="Input fasta or newick file", metavar="infile")
parser.add_argument("-t1", "--patristic distance", dest="patdist", default=True,required=True,
                    help="Patristic Distance (Consensus sequencing (ex 0.03))"
                    "Percentile (Deep sequencing dataset (ex 25)) Threshold", metavar="threshold")
parser.add_argument("-t2", "--tree support", dest="support", default=True,required=True,
                    help="Tree / Bootstrap support threshold", metavar="threshold")
parser.add_argument("-t3", "--genotype pat dist threshold", dest="genotype", default=0.3,
                    help="Optional: If multiple genotypes are included, this threshold prevents inter-genotype comparisons (default is 0.3)", metavar="threshold")
parser.add_argument("-a", "--analysis type", dest="analysis", default=True,required=True,
                    help="Type of analysis: Consensus/Sanger (ex C) or Deep Sequencing (ex DS)", metavar="analysis")
parser.add_argument("-path", "--pathway", dest="pathway", default=True,
                    help="Optional: Required only when a newick file is needed, enter the FastTree file pathway (ex: ../Documents/)", metavar="path")
parser.add_argument("-idS", "--idStart", dest="start", default=True,nargs='?', type= int,
                    help="Required for deep seq analysis, specify the starting position of the uniqe identifier in your reads for each sample (ex: 1)", metavar="start")
parser.add_argument("-idL", "--idLen", dest="len", default=True,nargs='?', type=int,
                    help="Required for deep seq analysis, specify the length of the uniqe identifier in your reads for each sample (ex: 3)", metavar="len")
parser.add_argument("-oR", "--outlierRem", dest="outlier", default=True,nargs='?',type = str,
                    help="Required for deep seq analysis if outlier removal is not needed (ex '-oR F'). Outliers are removed by default (ex '-oR T' ) ", metavar="outlier")
parser.add_argument("-o", "--output", dest="output", required = True,
                    help="Required: this specifies the output filename location", metavar="output")


args = parser.parse_args()

def run_FastTree(path, infile,outfile):
    cmd = "%sFastTree -gtr -nt %s >%s" % (FastTreePathway, args.infile,outputPath+inputfilePrefix+".newick")
    print "Running the following command:\n"+ cmd
    logfile.write("Running the following command:\n"+ cmd+"\n")
    os.system(cmd)
    return
def run_ConsensusClusterIdentifier(i):
    cmd = "python ete.ConsensusSequencingClusterID.py %s %s %s %s " % (i, PatDistCutoff, TreeSupport,output )
    print "Running the following command:\n"+ cmd
    logfile.write("Running the following command:\n"+ cmd+"\n")
    os.system(cmd)
    return
def run_PatDistSpectrum(i):
    cmd = "python PatDistSpectrum.py %s %s %s %s %s" % ( output, outputPath+TreeFile, genotypeCutoff, idLen, idStart )
    print "Running the following command:\n"+ cmd
    logfile.write("Running the following command:\n"+ cmd+"\n")
    os.system(cmd)
    return
def run_ETE(Spectrum):
    cmd = "python ete.DeepSequencingClusterID.py %s %s %s %s %s %s %s %s %s %s %s" % (outputPath+Spectrum, outputPath+PatDistOutput, PatDistCutoff, TreeSupport, filePath+IntraFile, filePath+FastaFile, outputPath+TreeFile, idLen, idStart, outlierFlag, output)
    print "Running the following command:\n"+ cmd
    logfile.write("Running the following command:\n"+ cmd+"\n")
    os.system(cmd)
    return

if "/" in str(args.infile):
    inputfile = str(args.infile).rsplit("/",1)[1]
    filePath = str(args.infile).rsplit("/",1)[0]+"/"
else:
    inputfile = str(args.infile)
    filePath = ""
output = str(args.output)
if "/" in str(args.output):
    outputPath = str(args.output).rsplit("/",1)[0]+"/"
else:
    outputPath = ""

inputfileSp = inputfile.rsplit(".",1)
inputfilePrefix = inputfileSp[0]
inputfileSuffixU = inputfileSp[-1].upper()
inputfileSuffix = inputfileSp[-1]
PatDistCutoff = args.patdist.upper()
TreeSupport = args.support.upper()
genotypeCutoff = args.genotype
IntraFile = 'Intrahost.txt'

if args.analysis == "C":
    
    if not args.len is True:
        print "For consensus sequencing, identifier information is not required"
        flag = 1
    if not args.start is True:
        print "For consensus sequencing, identifier information is not required"
        flag = 1
    if not args.outlier is True:
        print "For consensus sequencing, outlier removal is not possible using this approach"
        flag = 1
    try:
        if isinstance( int(args.patdist), int ) == True:
            if PatDistCutoff != "PASS":
                flag = 1
                print "For consensus sequencing, a patristic distance is required ..."
            
    except:
        pass
    if args.pathway is True:
        print "Please enter the pathway to locate FastTree (ex -path ../../) ..."
        flag = 1
    
    if flag == 1:
        sys.exit(1)
        flag = 0
elif args.analysis == "DS":
    print inputfilePrefix+'.IntraGenotype.patdistspectrum.txt'
    if args.len is True:
        print "For deep sequencing, identifier information is required, see help (-h)"
        flag = 1
    if args.start is True:
        if not flag == 1:
            print "For deep sequencing, identifier information is required, see help (-h)"
            flag = 1
    try:
        int(PatDistCutoff)
        
            
    except:
        if PatDistCutoff != "PASS":
            print "For deep sequencing, please enter one of the following percentiles:\n"
            print "0,1,5,10,20,25,30,35,40,45,50,75,90,100"
            flag = 1
    if flag == 1:
        sys.exit(1)
        flag = 0
if os.path.isfile(outputPath+'IntraClustLogfile.txt'):
    logfile = open(outputPath+'IntraClustLogfile.txt','a')
    logfile.write("\n")
    logfile.write("Analysis initiated at:"+ TimeStamp+"\n")
else:
    logfile = open(outputPath+'IntraClustLogfile.txt','w')
    logfile.write("\n")
    logfile.write("Analysis initiated at:"+ TimeStamp+"\n")
if __name__ == '__main__':
    if inputfileSuffixU == "FASTA" or inputfileSuffixU == "FAS" or inputfileSuffixU == "FA":
        
            
        if not os.path.isfile(outputPath+inputfilePrefix+".newick"):
            print "Detected FASTA file, constructing newick tree with FastTree..."
            if str(args.pathway[-1]) != "/":
                args.pathway = str(args.pathway)+"/"
            if os.path.isfile(args.pathway+"FastTree"):
                if os.path.isfile(args.infile):
                    FastTreePathway=args.pathway
                    run_FastTree(FastTreePathway,args.infile,outputPath+inputfilePrefix+".newick")
                else:
                    print "Input file not found. Please provide the correct file pathway to the input fasta / newick file"
                    sys.exit(1)
                    
            else:
                print "Please provide the correct file pathway to locate FastTree, or enter '' if not needed (ex: -path ../../)"
                sys.exit(1)
        else:
            print "Phylogenetic tree already completed for this file...."
            
                
        if args.analysis == "C" or args.analysis == "CONSENSUS":
            PatDistCutoff = args.patdist
            TreeSupport = args.support
            if os.path.isfile(outputPath+inputfilePrefix+".newick"):
                if not os.path.isfile(output):
                    run_ConsensusClusterIdentifier(outputPath+inputfilePrefix+".newick")
                else:
                    print "Analysis already completed for this file...."
            else:
                print "Input file not found. Please provide the correct file pathway to the input fasta / newick file"
                sys.exit(1)
        elif args.analysis == "DS" or args.analysis == "DEEP SEQUENCING":
            idLen = args.len
            idStart = args.start
            try:
                idLen = args.len
                idStart = args.start
                PatDistOutput = inputfilePrefix+".IntraGenotype.PatDist.txt"
                FastaFile = filePath+inputfilePrefix+"."+inputfileSuffix
                TreeFile = inputfilePrefix+".newick"
                outlierFlag = args.outlier
                outfile = args.output
            except:
                cmdError = "python IntraClust.py -h"
                os.system(cmdError)
                sys.exit(1)
            
            try:
                if not os.path.isfile(outputPath+inputfilePrefix+'.IntraGenotype.patdistspectrum.txt'):
                    run_PatDistSpectrum(outputPath+inputfilePrefix+".IntraGenotype.PatDist.txt")
            except KeyboardInterrupt:
                try:
                    os.remove(outputPath+inputfilePrefix+'.IntraGenotype.patdistspectrum.txt')
                except:
                    pass
            except:
                sys.exit(1)
            try:
                
                if not os.path.isfile(output):
                    print output
                    run_ETE(inputfilePrefix+'.IntraGenotype.patdistspectrum.txt')
                else:
                    print "Cluster analysis already completed for this file...."
            except KeyboardInterrupt:
                try:
                    os.remove(output)
                except:
                    pass
            except:
                sys.exit(1)
    elif inputfileSuffixU == "NEWICK" or inputfileSuffixU == "NWK" or inputfileSuffixU == "TREE" or inputfileSuffixU == "NW":
        print "Detected a newick phylogenetic tree, proceeding to cluster identification..."
        
        if args.analysis == "C" or args.analysis == "CONSENSUS":
            try:
                
                if not os.path.isfile(output):
                    run_ConsensusClusterIdentifier(outputPath+args.infile)
                else:
                    print "Cluster analysis already completed for this file...."
            except KeyboardInterrupt:
                try:
                    os.remove(output)
                except:
                    pass
            except:
                cmdError = "python IntraClust.py -h"
                os.system(cmdError)
                sys.exit(1)
        elif args.analysis == "DS" or analysis == "DEEP SEQUENCING":
            try:
                idLen = args.len
                idStart = args.start
                PatDistOutput = inputfilePrefix+".IntraGenotype.PatDist.txt"
                IntraFile = 'Intrahost.txt'
                TreeFile = inputfilePrefix+"."+inputfileSuffix
                FastaFile = inputfilePrefix+".fasta"
                outlierFlag = args.outlier
                outfile = args.output
            except:
                cmdError = "python IntraClust.py -h"
                os.system(cmdError)
                sys.exit(1)
            try:
                if not os.path.isfile(outputPath+inputfilePrefix+'.IntraGenotype.patdistspectrum.txt'):
                    run_PatDistSpectrum(outputPath+inputfilePrefix+".IntraGenotype.PatDist.txt")
                if not os.path.isfile(output):
                    run_ETE(inputfilePrefix+'.IntraGenotype.patdistspectrum.txt')
                else:
                    print "Cluster analysis already completed for this file...."
                
                
            except IndexError:
                print "If you generated your own phylogenetic tree, please ensure that the FASTA file and the tree filenames only differ by their file extensions...."
                sys.exit(1)
    else:
        print "A fasta file (ending in '.fasta', '.fa', or '.fas') or a newick tree file (ending in '.newick, '.tree', or '.nw') is required..."
    current2 = datetime.datetime.now()
    TimeStampFin = current2.strftime("%Y-%m-%d-%H:%M")
    logfile.write("Analysis completed at:"+ TimeStampFin+"\n")