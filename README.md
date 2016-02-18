#Intra-Clust

##Description
These python scripts are used for identifying phylogenetic clusters through the analysis of intra-host variants derived from deep sequencing data.  A Sanger sequencing cluster analysis pipeline is also provided for comparison. Clusters are identified directly from a fasta file or newick phylogenetic tree input.
For deep sequencing, clusters are identified using patristic distance thresholds determined by the intra-host diversity from each individual's viral population.
Whereas for Sanger sequencing analyses, clusters are traditionally defined using both patristic distance and tree support thresholds specified by the user.
After a preliminary run, the output files can be used to help determine which threshold is most appropriate for a given dataset.
This allows species with different mutation rates to be analyzed accordingly. This program was tested using a HCV deep sequencing dataset.

##Requirements
This package was implemented in python 2.7 and 2.6. For python version incompatibilities, please contact me.
Python packages required:

numpy
http://www.scipy.org/scipylib/download.html

dendropy
https://pythonhosted.org/DendroPy/
install dendropy: sudo pip install -U dendropy

ETE Toolkit
http://etetoolkit.org/

FastTree is optional, however, the scripts included here are set up to run it if you choose to start from a fasta file input.
http://www.microbesonline.org/fasttree/#Install