This directory includes subdirectories for each coral species analyzed, where output from analyses performed in jupyter (formerly IPython) notebooks is directed. Coral species subdirectories include:

  `Ahya` - *Acropora hyacinthus* 
  
  `Amil` - *Acropora millepora* 
  
  `Apalm` - *Acropora palmata* 
  
  `Past` - *Porites astreoides* 
  
  `Pdam` - *Pocillopora damicornis* 
  
  `Spist` - *Stylophora pistillata* 
  

Each coral species subdirectory includes files with the following suffixes:

`_blastx_uniprot.tab` - Uniprot/Swissprot annotation

`_cpg_GOslim.tab` - File with contig ID, CpG O/E, and GOSlim Biological Process term

`Ahya`, `Amil`, and `Apalm` subdirectories aso include the following:

`_exp_CpG2` - File with contig ID, expression counts for each sample, and CpG O/E in the last column


The directory also includes a `scripts/` directory where R scripts used for graphics and statistical analyses are housed, as well as scripts called on from jupyter notebooks. 