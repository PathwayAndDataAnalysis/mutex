This project aims to discover driver alterations in genomic datasets by detecting mutual exclusion (mutex) patterns in alterations of groups of genes with a common downstream effect. Below are example mutex groups from TCGA glioblastoma dataset.

![http://resources.mutex.googlecode.com/hg/gbm_tcga_groups.svg](http://resources.mutex.googlecode.com/hg/gbm_tcga_groups.svg)

Here, we show member genes of a mutex group in a compound node. Label of the compound node displays the coverage of gene alterations in the group as percentage. When we merge the detected mutex groups, we get the below aberrant signaling network composed of potential cancer-driver genes.

![http://resources.mutex.googlecode.com/hg/gbm_tcga.svg](http://resources.mutex.googlecode.com/hg/gbm_tcga.svg)

And below is the oncoprint for the group TRPV4-TP53-CDKN2A-MDM2-MYH1.

```
TRPV4  MMMM............................................................................................................
TP53   ....MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMD............................................................................
CDKN2A ...D.....................DDDDDDDDDDDMMDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD....................
MDM2   ........................M...................................................................AAAAAAA.............  
MYH1   ...................................................................................................MM...........
```

In the oncoprint, each column represents a patient sample. A dot (.) means that the gene is not altered in that patient sample. M: mutated, A: copy number gained, D: copy number lost. In the analysis (also in the oncoprint), we use the copy number alterations only if they are also verified with gene expression.

The directed signaling network we use for finding candidate genes with a common downstream is compiled using databases Pathway Commons, SPIKE, and SignaLink. In the graphs above, blue solid edges represent post-translational modifications, and green dashed edges represent transcriptional regulations.



Please see UserGuide to download, customize, and run Mutex code.

A publication about analysis of TCGA datasets using Mutex is still in progress. A [preprint](http://dx.doi.org/10.1101/009878) is available in the meantime.