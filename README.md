## Analyses and data associated with manuscript: [_Germline DNA methylation in six species of reef corals: patterns and potential roles in response to environmental change_ ](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.13414)

---

###Description of repository directory structure

_Each subdirectory has its own readme file describing the contents._

---   

`data/` - Raw data, including publically available fasta data. These are housed in separate subdirectories for each species.

Coral species subdirectories include:

  `data/Ahya` - *Acropora hyacinthus* 
  
  `data/Amil` - *Acropora millepora* 
  
  `data/Apalm` - *Acropora palmata* 
  
  `data/Past` - *Porites astreoides* 
  
  `data/Pdam` - *Pocillopora damicornis* 
  
  `data/Spist` - *Stylophora pistillata* 
  


--- 

`ipynb/` - jupyter notebooks.


---

`analyses/` - Includes subdirectories for each coral species analyzed, which include output from analyses performed in jupyter (formerly IPython) notebooks. Also includes a `scripts/` directory where R scripts used for graphics and statistical analyses are housed, as well as scripts called on from jupyter notebooks. 

---



### Sofware versions originally used in this analyses (on Mac OS X v10.10.3) include: 

* Python: 2.7.9  
* IPython: 3.1.0
* R: 3.1.3  
* NCBI Blast: 2.2.3 


---

###General notes on workflow

The workflow for each species starts with jupyter notebooks ending with the suffix "_blast_anno.ipynb", which performs a blast annotation of the transcriptome and provides instructions for further annotation with GOSlim terms. Next, CpG O/E analysis is carried out in jupyter notebooks with the suffix "_CpG_ratio.ipynb". The remainder of analyses are conducted in R using the following scripts: 

- CpG_Density.R: Plots and compares mixture models for CpG O/E data.
- CpG_GOslim.R: Deals with reading in a tab delimited file containing CpGo/e and GOSlim information. The script performs Fisher's exact tests and plots various types of figures. Note that some files are derived from analyses created in CpG_Density.R, so that script should be run prior to and alongside this one.
- CpG_deg.R: Evaluates differentiallly expressed genes in Acropora hyacinthus, A. millepora, and A. palmata.
- Expression.R: Analyzes gene expression vs. CpG O/E.

---

We are actively trying to improve this realizing that we are likely missing dependancies, etc. Any suggestions and feedback is welcome. 

