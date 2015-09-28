## jupyter (IPython) notebooks
These notebooks detail analysis procedures and can be rendered natively within GitHub by clicking on them. Notebook files start with the species prefix: 


  `Ahya` - *Acropora hyacinthus* 
  
  `Amil` - *Acropora millepora* 
  
  `Apalm` - *Acropora palmata* 
  
  `Past` - *Porites astreoides* 
  
  `Pdam` - *Pocillopora damicornis* 
  
  `Spist` - *Stylophora pistillata* 
  

Notebooks end with the following suffixes:
  

`_blast_anno.ipynb` - Performs transcriptome annotation.

`_CpG_ratio.ipynb` - Analyzes CpG O/E.

`_zoox_removal.ipynb` - Identifies putative *Symbiodinium* sequences for removal (Pdam and Spist only).

`_expression.ipynb` - Joins expression count data with CpG O/E data (only for Ahya, Amil, and Apalm).

`_exp_CpG_ratio.ipynb` - Only for Amil; performs CpG O/E calculation on Amil transcriptome using GenBank accession numbers for contig IDs. This was necessary to join this data with expression data.

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

5) Execute code as written or modify to your liking. To execute cell type *shift-enter*. Relative pathnames assume you are within the Coral-CpG directory.


