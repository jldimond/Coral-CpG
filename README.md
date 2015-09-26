## Juptyer notebooks and other analyses associated with manuscript: _Germline DNA methylation in six species of reef corals: patterns and potential roles in response to environmental change_

---

###Description of repository directory structure

_Each subdirectory has its own readme file describing the contents._

---   

`data/` - Raw data, including publically available fasta data. These are housed in separate subdirectories for each species.

--- 

`ipynb/` - jupyter notebooks.


---

`analyses/` - Includes subdirectories for each coral species analyzed, which include output from analyses performed in jupyter (formerly IPython) notebooks. Also includes a `scripts/` directory where R scripts used for graphics and statistical analyses are housed, as well as scripts called on from jupyter notebooks. Coral species subdirectories include:

  `analyses/Ahya` - *Acropora hyacinthus* 
  
  `analyses/Amil` - *Acropora millepora* 
  
  `analyses/Apalm` - *Acropora palmata* 
  
  `analyses/Past` - *Porites astreoides* 
  
  `analyses/Pdam` - *Pocillopora damicornis* 
  
  `analyses/Spist` - *Stylophora pistillata* 
  
---



### Sofware versions originally used in this analyses (on Mac OS X v10.10.3) include: 

* Python: 2.7.9  
* IPython: 3.1.0
* R: 3.1.3  
* NCBI Blast: 2.2.3 

---

###Instructions for full annotation and analysis (interactive execution) notebook

1) **Before you get started**

To execute the jupyter (Ipython, `.ipynb`) notebooks in their entirety you will need:   

* IPython - [install instructions](http://ipython.org/install.html)    
* NBCI Blast -  [install instructions](http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)  
* R - [install instructions](http://www.r-project.org/)  
* SQLShare-pythonclient (optional; current workflow uses manual upload and download) [install instructions](https://github.com/uwescience/sqlshare-pythonclient)

---

In addition you will need a local copy of the UniProt/SwissProt Blast Database. 
If you do not already have this database you can create it once you install NCBI Blast. Create a blast database first download a fasta file from <http://www.uniprot.org/downloads> that of Reviewed (Swiss-Prot) then run make blastdb commands.
Below is an example code if you wanted to create the database in a subdirectory named `blastdb`. This will result in files > 300 MB.

`$ mkdir blastdb`

`$ cd blastdb`

`$ curl -O ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz`

`$ gunzip uniprot_sprot.fasta.gz`

`$ makeblastdb -in uniprot_sprot.fasta -dbtype prot -out uniprot_sprot`

This will generate a Protein database that you can use to blast sequences. 

---

2) Download the repository zip file to a local directory and uncompress. 

3) Launch IPython from the repository primary directory. 
For example, using Terminal on MacOSX.


```
$ cd /Desktop/Coral-CpG-ratio-MS
$ ipython notebook

```
This will launch IPython in your web browser.  


4) Open notebook by clicking on any `.ipynb` files. This will open a new tab in your browser.

5) Execute code as written or modify to your liking. To execute cell type *shift-enter*. Relative pathnames assume you are within the Coral-CpG-ratio-MS directory.


---

**Workflow** 

The workflow for each species starts with jupyter notebooks ending with the suffix "_blast_anno.ipynb", which performs a blast annotation of the transcriptome and provides instructions for further annotation with GOSlim terms. Next, CpG O/E analysis is carried out in jupyter notebooks with the suffix "_CpG_ratio.ipynb". The remainder of analyses are conducted in R using the following scripts: 

- CpG_Density.R: Plots and compares mixture models for CpG O/E data.
- CpG_GOslim.R: Deals with reading in a tab delimited file containing CpGo/e and GOSlim information. The script performs Fisher's exact tests and plots various types of figures. Note that some files are derived from analyses created in CpG_Density.R, so that script should be run prior to and alongside this one.
- CpG_deg.R: Evaluates differentiallly expressed genes in Acropora hyacinthus, A. millepora, and A. palmata.
- Expression.R: Analyzes gene expression vs. CpG O/E.

---

We are actively trying to improve this realizing that we are likely missing dependancies, etc. Any suggestions and feedback is welcome. 

