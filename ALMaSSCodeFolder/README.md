This folder contains the ALMaSS files that I have made changes to during my thesis.  
Here follows a description of what I have added to each file!  

In the file-duo `GeneticMaterial.cpp` and `GeneticMaterial.h`:  
* A new class `GeneticMaterialNEW` containing 2 genomes and a paternal and maternal lineage ID, with functions:
  - A constructor function `GeneticMaterialNew()`; <br>
Calls `InitializeGenome();`.
  - `InitializeGenome()`; <br>
    Randomly generates two genomes for the vole from a file of allele frequency input, by producing a discrete distribution with the probabilites provided by the allele frequency intput. <br>
    If an empty string is passed to the `InitializeGenome()`; the vole keeps the genome it already has. (For keeping the inherited genes from parents, if the vole is an offspring.)
  - `PrintGenomes()`;
    Prints genenomes in cout, mainly for testing.
  - `EraseGenomes()`;
    Erases genomes, maninly for testing.
  - `InitializeTestingGenomesFemale()`;
    A function that gives all females initialized a genome og just 0's and one of just 1's. For testing.
  - `InitializeTestingGenomesMale()`;
    A function that gives all males initialized a genome og just 2's and one of just 3's. For testing. These two testing genome functions practical for testing inheritance dynamics, for instance for testing the function `MakeGameteSimple()`;
  - `MakeGameteSimple()`;
    This function makes a gamete from an instance of `GeneticMaterialNew`. It takes the two genomes contained in the `GeneticMaterialNew` (Genome0 and Genome1) and for each locus, it is randomly decided if the allele from Genome0 or Genome1 is going to be included into the gamete. This means, that there is is no linkage in inheritance. This function is applied to the vole mother's `GeneticMaterialNew` and the vole father's `GeneticMaterialNew` when offspring are being produced.

In the file-duo `VolePopulationManager.cpp` and `VolePopulationManager.h`:  
* New functions for the class PopulationManager:
  - `Four_QuadrantBasedGeneticOutput()`;
    This function divides the landscape map into 4 different quadrants with buffer zones between them and extensively records the genetics present and other information about the voles in them. It records:
    * Birth postition and current position of the voles (for inference of movement).
    * The full genome of each vole as well as its position in the map.
    * The allele frequencies present in each quadrant.
    * FST values between all pairs of quadrants. 
    * FIS values of all quadrants.
    * Genome-wide, population-wide average and observed and expected heterozygosity rate.
  - `Nine_QuadrantBasedGeneticOutput()`;
    This function divides the landscape map into 9 (3*3) different quadrants WITHOUT buffer zones between them and extensively records the genetics present and other information about the voles in them. It records:
    * Birth postition and current position of the voles (for inference of movement).
    * The full genome of each vole as well as its position in the map.
    * The allele frequencies present in each quadrant.
    * FST values between all pairs of quadrants. 
    * FIS values of all quadrants.
    * Genome-wide, population-wide average and observed and expected heterozygosity rate.
  - `FourQuadrantsPopulationSizeProbe()`;
    This functions counts how many voles are present in each quadrant of four at sampling time, divided by males and females. Also counts how many voles are excluded by being in the buffer zone.
  - `NineQuadrantsPopulationSizeProbe()`;
  This functions counts how many voles are present in each of nine quadrant at sampling time, divided by males and females.
  - `MakeAlleleInputFile()`;
    This function is to be called very early in the population initiation (before creating voles) to generate allele frequencies for each chromosome and locus. This funciton takes two sets of two alpha parameters. Half of the chromosomes' allele frequencies are generated according to the first set of alpha parameters, and the second half is generated according to the second set of alpha parameters. The alpha-parameters are used in making a beta distribution from which an allele frequency is drawn. In this project, we use di-allelic loci, so the allele frequency of the other allele is then given by the 1-first allele frequency. 
  - `CompareGenomes()`;
    This function takes the full genome of each vole recorded by `Four_QuadrantBasedGeneticOutput()`; og `Nine_QuadrantBasedGeneticOutput()`; samples 500 random voles, and compares their genomes. For each locus, a score is given based on similarity between the two voles.
* If they share both alleles: 2 points
* If they share one allele: 1 point
* If they share no allele: 0 points.
The sum is then divided by the max possible score (256) to get a normalized similarity score between 0 and 1. Outputs a file of comma seperates values of the position of the voles, the distance between them, and their similarity score.
  - `GenerationCountOutput()`;
    Outputs the generation count for the voles.
  - `LineagesOutput()`;
    Outputs the paternal and maternal lineages of the voles at a provided sampling rate (here: 1).

* Added to the PopulationManager functions `PopulationManager::CreateObjects` and `PopulationManager::CreateObjectsInit`:
  - A vole that is being born from parents keeps the genome from its parents' gametes, so it calls the construction function `NAME` with an empty string for the filename of the allele frequency input.
  - A vole under initialization calls the function with a string that contains the allele frequency file and will get its own genome based on allele frequency inputs.

In the file-duo `Vole_all.cpp` and `Vole_all.h`:  
* Each vole `class struct_Vole_Adult` was given the attributes:
  - `yborn` (y coordinate of birth)
  - `xborn` (x coordinate of birth)
  - `new_Genes`, an instance of `GeneticMaterialNEW`
  - `new_Mates_Genes`, an instance of `GeneticMaterialNEW`
  - `GenerationCount`
  - `mitochondrial lineage`
  - `ychromosome lineage`

* Added to the `Vole_Female::st_Lactating()` function (which is used for giving birth):
  - Voles record their birthplace into yborn and xborn.
  - Voles inherit a generation count +1 from their mother.
  - The mother vole recombines her genes and the father vole's genes with the function `MakeGameteSimple()`.
  - The mother passes on her mitochondrial lineage and the father's ychromosome lineage to the offspring.

In the file-duo `Vole_toletoc.cpp` and `Vole_toletoc.h`:  
* The only modification here is that when a vole looks up the movement quality of a large road, it will in 10% of lookup cases draw the value "1" (the vole can then walk onto the road) and in 90% of the cases draw the value "-1" (the vole cannot walk onto to road). 
