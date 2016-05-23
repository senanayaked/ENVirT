# ENVirT
Repository for the Matlab implementation of ENVirT : simultaneous estimation of metavirome richness and average genome length

ENVirT is an algorithm which provides estimations of the viral richness of a metagenomic sample without prior assumptions of 
the average genome length being required as an input. 

## Input Specification
To run the algorithm on a dataset, use the following input specification: 
* ``` spectrum_file``` : the file containing the contig spectrum that needs to be analyzed
* ```n_reads``` : the number of reads in the contig spectrum. This needs to be explicitly given since trimming the contig spectrum
  may change the number of reads calculated from it. 
* ```avg_read_length``` : the average length of a single read in the sample. 
* ```overlap_length``` : The overlap length (in base pairs) considered for the sample.
* ```trim_length``` : trim length for the contig spectrum. this can help with reducing the running time.
* ```M_min``` : lower bound for the richness estimate
* ```M_max``` : upper bound for the richness estimate
* ```p_res``` : resolution for the p value
* ```p_max``` : maximum for the p value
* ```L_res``` : resolution for average genome length estimates
* ```L_min``` : lower bound for average genome length estimates
* ```L_max``` : upper bound for average genome length estimates
* ```n_windows```: number of niching windows
* ```result_folder``` : path of the results folder
* ```persist_intermediates``` : True for a file with intermediate results, False to skip writing intermediate results to 
  files. recommended to skip writing for faster calculations. 
  
## Running a Sample

A sample run would require adding ENVIRT.m to the working directory path. For demonstration purposes, consider the following
folder structure

+ Working_Directory
|-contig_spectrum.txt
|-  ENVIRT.m
|-  ENVIRT_CORE.m

Suppose the ```contig_spectrum.txt``` contains the contig spectrum of 10000 reads. An example scenario can be seen as below.

  ```ENVIRT('contig_spectrum.txt',10000,100,35,50,1,15000,'',5,'',10000,310000,29,'./Result_Set',true) ```
  
Note that some inputs are given empty strings, since they contain default values. However, the spectrum file, number of 
reads, the average read length, and the overlap length cannot be left empty. 

For further information, refer to the original publication of ENVirT, by D. Jayasundara et al. 
