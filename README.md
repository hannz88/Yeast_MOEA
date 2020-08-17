# Yeast MOEA
![R Version](https://img.shields.io/badge/R-3.6.1-brightgreen)
![C++ Version](https://img.shields.io/badge/C++-17-brightgreen)

The project contains the codes I've used to perform biclustering analysis on Synthetic Physical Interaction (SPI) from Thorpe Lab (Ólafsson and Thorpe, 2015).
To explain it simply, conventional clustering is looking for groups of rows with similar patterns across columns, biclustering is looking for subset of rows under a subset columns (Padilha and Campello, 2017).
The project aim was to detect dark proteome using SPI. This is done by finding biclusters of genes that have similar growth patterns when the gene products interact with certain proteins. 

Due to confidentiality, the original dataset is not included here. 
The project starts off by performing PCA as an exploratory analysis and then using genetic algorithm (GA) to optimise it. The codes for them are in `PCA.R` and `Genalg_PCA.R`, respectively. 
The input required is a csv of float/ double, where the rows consists of values for the genes readout and the columns are proteins.

The project then uses NSGA2 to optimise 3 objectives that have been documented to detect different types of biclusters (Dale, Zhao and Obafemi-Ajayi, 2019). R package "ecr" is used to perform NSGA2 (Bossek, 2017).
Since the project uses a "white-box" approach, some of the subfunctions was tweaked before being used. The NSGA2 subfunctions are in `NSGA2_subfunctions.R` . Import them when running `NSGA2.R`.
The 3 objectives are used as fitness functions for NSGA2. As pure R was not time-efficient enough, the 3 objectives were written in C++ and imported as library into R. The codes for the objectives are in `yeastfitness` folder.
To install yeastfitness as a library in R, run the following command in the terminal in the directory where yeastfitness is stored:
`R CMD INSTALL yeastfitness`

`NSGA2.R` contains codes for NSGA2. Users could change the number iterations and populations here. for them are in `PCA.R` and `Genalg_PCA.R`, respectively. The RData output should contains the genome of the optimal solutions.

`NSGA2_Postprocessing.R` contains codes for:
- obtaining the Pareto optimal set genome.
- combining solutions that are highly overlapped using Jaccard coefficient.
- plotting the progression of the 3 objectives over 10 iterations (a small number due to visibility issue).

`GO_enrichment.R` contains codes for performing gene overrepresentation test. Here, user needs to provide the original dataset that was used for `NSGA2.R` as well as the output from `NSGA2.R` . The codes initially look through the
Pareto optimal set and find the biclusters that are significantly enriched using PANTHER API (Mi et al., 2019). The output is a list of index of which bicluster in the input that are significantly enriched. The user can then either look 
them up on PANTHER (or other GO enrichment site) or programmatically through a local IDLE. In the first case, there are also codes that will write out the candidate gene list and the background genes list into separate txt files. Lastly,
the file also include codes to visually inspect the bicluster of interest via heatmap.

## References:
- Bossek, J., 2017, July. ecr 2.0: a modular framework for evolutionary computation in R. In Proceedings of the Genetic and Evolutionary Computation Conference Companion (pp. 1187-1193).
- Dale, J., Zhao, J. and Obafemi-Ajayi, T., 2019, July. Multi-objective Optimization Approach to find Biclusters in Gene Expression Data. In 2019 IEEE Conference on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB) (pp. 1-8). IEEE.
- Mi, H., Muruganujan, A., Huang, X., Ebert, D., Mills, C., Guo, X. and Thomas, P.D., 2019. Protocol Update for large-scale genome and gene function analysis with the PANTHER classification system (v. 14.0). Nature protocols, 14(3), pp.703-721.
- Ólafsson, G. and Thorpe, P.H., 2015. Synthetic physical interactions map kinetochore regulators and regions sensitive to constitutive Cdc14 localization. Proceedings of the National Academy of Sciences, 112(33), pp.10413-10418.
- Padilha, V.A. and Campello, R.J., 2017. A systematic comparative evaluation of biclustering techniques. BMC bioinformatics, 18(1), p.55.
