# AICg
R functions for computing the AICg for the phylogenomic quartet models of R package MSCquartets

The file AICgRFunctionsAnalysis.R applies the NANUQ algorithm of Allman et. al. (2019) to infer a hybridization network from a collection of gene trees, under the level-1 network multispecies coalescent (NMSC) model.

MSCquartets uses the hypothesis tests of Mitchell et. al. (2019) to infer either a star tree, a resolved tree or a level-1 network for each set of four taxa. An N-taxon tree or level-1 network is inferred from the inferred 4-taxon trees and level-1 networks.

The function quartetTableResolvedAICg infers one of the same three models for each set of four taxa, but uses the AICg instead of hypothesis tests.

Users must input a file of topological or metric gene trees in Newick format. Users must also specify:
1) a level for a Holm-Bonferroni-like adjustment to control false discoveries of 4-taxon level-1 networks,
2) a level for the AICg weight of the star tree model, at or below which the star tree model is not selected and
3) a quantity to specify the penalty of the star tree model.
