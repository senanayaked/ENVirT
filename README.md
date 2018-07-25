# ENVirT
Repository for the Matlab implementation of ENVirT : simultaneous estimation of metavirome richness and average genome length

ENVirT is an algorithm which provides simultaneous estimates of the viral richness and average genome length of a metagenomic sample without prior assumptions of the average genome length. 

## Running the Algorithm

To run the algorithm on a dataset, open ENVIRT.m using a matlab distribution and execute the file without any input parameters. 
`` >> ENVIRT()``

You will be prompted to select a contig spectrum file. 

The next input window will have the following inputs to be completed. 

### Config Spectrum Configuration
* Number of reads in the contig spectrum. [required]
* Average length of a single read in the sample (bp). [required]
* Overlap length considered when generating the contig spectrum (bp). [required]
* Number of contig spectrum elements to consider in the algorithm. This can help with reducing the running time.

The next window will prompt you to input the search configuration for the algorithm. 

### Search Configuration 
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

We have included a sample dataset in the ``Contig_Spectrum_Example`` directory. The running instructions are as below. 

* Open matlab and execute ENVIRT.m as follows

`` >> ENVIRT(); ``

* Select the ``AB_1_NG_10000_GL_25_NR_10000_RL_100_S_0_SPECTRUM`` inside the ``Contig_Spectrum_Example`` directory 
 
* Complete the _Contig Spectrum Configuration Prompt_ with the following inputs 
	- No. of reads : 10000
	- Average read length : 100 
	- Overlap Length : 35
	- Trim Length : 50 
* Leave the _search configuration_ in their default values

* Select result save location

* Select if you want to save the files that are intermediate results

Now the script will be executed. 


## Outputs

At the end of the execution of the function 'ENVirT', a new folder called 'Result_Set' will be created as follows.

+  Result_Set
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
