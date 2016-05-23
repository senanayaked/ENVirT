# ENVirT
Repository for the Matlab implementation of ENVirT : simultaneous estimation of metavirome richness and average genome length

ENVirT is an algorithm which provides simultaneous estimates of the viral richness and average genome length of a metagenomic sample without prior assumptions of the average genome length. 

## Input Specification
To run the algorithm on a dataset, use the following input specification: 
* spectrum_file : File containing the contig spectrum that needs to be analyzed [required]
* n_reads       : Number of reads in the contig spectrum. [required]
* avg_read_length : Average length of a single read in the sample (bp). [required]
* overlap_length : Overlap length considered when generating the contig spectrum (bp). [required]
* trim_length : Number of contig spectrum elements to consider in the algorithm. This can help with reducing the running time.
* M_min : Lower bound of the search space for richness (M)
* M_max : Upper bound of the search space for richness (M)
* d_res : Resolution for the distribution specific parameter (d)
* d_max : Upper bound of the search space for the distribution specific parameter (d)
* L_res : Initial resolution of average genome length (L)
* L_min : Lower bound of the search space for average genome length (L)
* L_max : Upper bound of the search space for average genome length (L)
* n_windows : Mumber of niching windows along the L dimension
* result_folder : Path of the results folder
* persist_intermediates : True for a file with intermediate results, False to skip writing intermediate results to files. Recommended to skip writing for faster calculations. 
  
## Running a Sample

A sample run would require adding ENVIRT.m to the working directory path. For demonstration purposes, consider the following
folder structure

+ Working_Directory
|-contig_spectrum.txt
|-  ENVIRT.m
|-  ENVIRT_CORE.m

Suppose the ``contig_spectrum.txt`` contains the contig spectrum of 10000 reads. An example scenario can be given as below.

ENVIRT('contig_spectrum.txt',10000,100,35,50,1,15000,'',5,'',10000,310000,29,'./Result_Set',true)

## Outputs

At the end of the execution of the function 'ENVirT', a new folder called 'Result_Set' will be created as follows.

+ Working_Directory
|-  contig_spectrum.txt
|-  ENVIRT.m
|-  ENVIRT_CORE.m
|+  Result_Set
  |-  ga_results_1.txt
  |-  ga_results_2.txt
  |-  ga_results_3.txt
  |-  ga_results_4.txt
  |-  model_results.txt
  |-  results.txt
  |-  results_final.txt

The final results, with the best fitting model, will be recorded in the 'results_final.txt' file as follows. 

 location: contig_spectrum.txt
 richness: 4879
 average genome length: 297733
 d: 2.48900
 a: lognormal
 residual error: 1.697270e-04
 run time: 8121.02s
 evenness: 0.6569

In ases where the estimates of each model is needed, you may read the results in 'model_results.txt'. And intermediate results are written in 'results.txt'. The optimization steps for each model are written in ga_results_<model>.txt files. 
(for the models, 1- Power Law, 2- Exponential, 3- Logarithmic, 4-Lognormal).

For further information, refer to the original publication of ENVirT, by D. Jayasundara et al. 
