# Generalised Telegraph Model for scRNA-seq transcriptional bursting. 

Transcription is a discrete process with stochastic properties, resulting in the expression of RNA in bursts. The classical telegram model (CTM) often used to describe transcriptional bursting is a two- state model where promoters switch between on and off states given exponentially distributed dwell times. It was first implemented to describe transcriptional dynamics using scRNA-seq by Larsson et al. (2019) [1]. The model captures part of the diversity of transcription regulation.
However, the initiation of transcription is a multi-step process that doesn’t always follow exponentially distributed dwell times.Thus, we explore the limitations of the CTM to describe certain gene regulation patterns and compare its performance to a generalized telegram model (GTM) that allows arbitrary dwell-time distributions. The GTM was first implemented by Luo et al. 2022 [2] in Matlab, we implemented a python version for better comparison and increased scalability. 
While both models estimate similar transcriptional bursting parameters for most genes. Some known bursting genes, certain transcription factors regulated by extracellular signaling (i.e. Fos, Atf4) are better described by the GTM. This highlights the importance of a more comprehensive transcriptional dynamic model to capture more of the diversity of transcriptional regulation mechanisms.

* Example script (python [script_name] [input_CSV_file] [output_directory])
> python .\generate_gtm_example.py .\..\data\inputs\example.csv .\..\data\results\total_allelic\  

The data used here is publicly available data originating by the Sandberg lab[1,3]

Coming soon : 
- more detailed tutorial
- more detail example dataset
- benchmarking of the GTM
- script to simulate bursting data
- visualisation functions


Citations : 

[1] Songhao Luo, Zhenquan Zhang, Zihao Wang, Xiyan Yang, Xiaoxuan Chen, Tianshou Zhou, Jiajun Zhang, Inferring transcriptional bursting kinetics from single-cell snapshot data using a generalized telegraph model. bioRxiv 2022.07.17.500373; doi: https://doi.org/10.1101/2022.07.17.500373 
[2] Larsson AJM, Johnsson P, Hagemann-Jensen M, Hartmanis L, Faridani OR, Reinius B, Segerstolpe Å, Rivera CM, Ren B, Sandberg R. Genomic encoding of transcriptional burst kinetics. Nature. 2019 Jan;565(7738):251-254. doi: 10.1038/s41586-018-0836-1. Epub 2019 Jan 2. PMID: 30602787; PMCID: PMC7610481. 
[3] Larsson AJM, Ziegenhain C, Hagemann-Jensen M, Reinius B, Jacob T, Dalessandri T, Hendriks GJ, Kasper M, Sandberg R. Transcriptional bursts explain autosomal random monoallelic expression and affect allelic imbalance. PLoS Comput Biol. 2021 Mar 9;17(3):e1008772. doi: 10.1371/journal.pcbi.1008772. PMID: 33690599; PMCID: PMC7978379. 
