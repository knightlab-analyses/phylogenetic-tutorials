
# Kernel Methods

Kernel methods aim to quantify dissimilarity between microbial communities.  By defining a distance between two microbial communities, it is possible to perform classification and cluster communities obtain a high level understanding of how microbial communities are impacted by specific environmental variables.

One of the most commonly used distance metrics used in microbial ecology is [Unifrac](http://aem.asm.org/content/71/12/8228.short).  Unifrac incorporates information from prior knowledge of the evolutionary history of the microbes, in addition to the observed microbial abundances.


```python
%load_ext rpy2.ipython
```

    The rpy2.ipython extension is already loaded. To reload it, use:
      %reload_ext rpy2.ipython



```python
%%R -o OTUTable -o disturbance_frequency

library(ape)
library(phytools)
library(magrittr)
library(nlme)
set.seed(1)
tree <- rtree(10)
n <- 25
disturbance_frequency <- rexp(n) %>% sort %>% log

Q <- diag(9) %>% cbind(rep(0,9),.) %>% rbind(.,rep(0,10)) 
Q <- Q+t(Q)-diag(10)
# Q <- Q
RNAcopyNumber <- sim.history(tree, Q, anc = '3')


abundances <- function(dst,RNAcopy){
  m <- length(RNAcopy)
  muTot <- 1e4   ##a mean of 10,000 sequence counts per sample.
  logmu <- 3 * dst * log(RNAcopy) #model to yield linear changes in log-ratios
  muRel <- exp(logmu) / sum(exp(logmu))  #mean relative abundances 
  mu = muRel * muTot
  size = 1
  N <- rnbinom(m, size, mu=mu)
}

OTUTable <- sapply(as.list(disturbance_frequency),
                   FUN=function(dst, c) abundances(dst, c),
                   c=as.numeric(RNAcopyNumber$states)) %>% matrix(., ncol=n, byrow=F)
OTUTable <- as.data.frame(OTUTable)
rownames(OTUTable) <- names(RNAcopyNumber$states)
write.tree(tree, file='tree.nwk', tree.name=TRUE)
```


    Some columns (or rows) of Q don't sum to 0.0. Fixing.
    Done simulation(s).




```python
from skbio import TreeNode
from skbio.diversity import beta_diversity
tree = TreeNode.read('tree.nwk')
table = OTUTable.astype(int).T
wu_dm = beta_diversity("weighted_unifrac", table.values, table.index, tree=tree,
                        otu_ids=table.columns)
```


```python
wu_dm
```




![svg](output_4_0.svg)




```python
import pandas as pd
samp_md = pd.DataFrame(disturbance_frequency, index=table.index, columns=['Disturbance_Frequency'])
```


```python
from skbio.stats.ordination import pcoa

pc = pcoa(wu_dm)
pc.plot(samp_md, 'Disturbance_Frequency',
        axis_labels=('PC 1', 'PC 2', 'PC 3'),
        title='Disturbance_Frequency', cmap='jet', s=50)
```

    /Users/mortonjt/miniconda3/envs/phylogenetic-tutorials/lib/python3.5/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:111: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.4781015679113227 and the largest is 21.026077067875075.
      RuntimeWarning





![png](output_6_1.png)




```python

```
