This folder contains the ALMaSS files that I have made changes to during my thesis.  
Here follows a description of what I have added to each file!   
_____________________________________________
In the file-duo `GeneticMaterial.cpp` and `GeneticMaterial.h`:  
* A new class `GeneticMaterialNEW` containing 2 genomes and a paternal and maternal lineage ID, with functions:
  - A constructor function `GeneticMaterialNew()`; (line 66-74 in the cpp-file) <br>
Calls `InitializeGenome();`.
  - `InitializeGenome()`; (line 85-132 in the ccp-file) <br>
    Randomly generates two genomes for the vole from a file of allele frequency input, by producing a discrete distribution with the probabilites provided by the allele frequency intput. <br>
    If an empty string is passed to the `InitializeGenome()`; the vole keeps the genome it already has. (For keeping the inherited genes from parents, if the vole is an offspring.)
  - `PrintGenomes()`; (line 146-163 in the cpp-file)
    Prints genenomes in cout (the terminal), mainly for testing.
  - `EraseGenomes()`; (line 76-79 in the cpp-file)
    Erases genomes, mainly for testing.
  - `InitializeTestingGenomesFemale()`; (line 140-144 in the cpp file)
    A function that gives all females initialized a genome og just 0's and one of just 1's. For testing.
  - `InitializeTestingGenomesMale()`; (line 134-137 in the cpp-file)
    A function that gives all males initialized a genome og just 2's and one of just 3's. For testing. These two testing genome functions practical for testing inheritance dynamics, for instance for testing the function `MakeGameteSimple()`;
  - `MakeGameteSimple()`; (line 166-182 in the cpp-file)
    This function makes a gamete from an instance of `GeneticMaterialNew`. It takes the two genomes contained in the `GeneticMaterialNew` (Genome0 and Genome1) and for each locus, it is randomly decided if the allele from Genome0 or Genome1 is going to be included into the gamete. This means, that there is is no linkage in inheritance. This function is applied to the vole mother's `GeneticMaterialNew` and the vole father's `GeneticMaterialNew` when offspring are being produced.
___________________________________________
In the file-duo `VolePopulationManager.cpp` and `VolePopulationManager.h`:  
* New functions for the class PopulationManager:
  - `Four_QuadrantBasedGeneticOutput()`;(line 3551 to 4024 in the cpp-file)
    This function divides the landscape map into 4 different quadrants with buffer zones between them and extensively records the genetics present and other information about the voles in them. It records:
    * Birth postition and current position of the voles (for inference of movement).
    * The full genome of each vole as well as its position in the map.
    * The allele frequencies present in each quadrant.
    * FST values between all pairs of quadrants. 
    * FIS values of all quadrants.
    * Genome-wide, population-wide average and observed and expected heterozygosity rate.
  - `Nine_QuadrantBasedGeneticOutput()`; (line 4029-4463 in the ccp-file)
    This function divides the landscape map into 9 (3*3) different quadrants WITHOUT buffer zones between them and extensively records the genetics present and other information about the voles in them. It records:
    * Birth postition and current position of the voles (for inference of movement).
    * The full genome of each vole as well as its position in the map.
    * The allele frequencies present in each quadrant.
    * FST values between all pairs of quadrants. 
    * FIS values of all quadrants.
    * Genome-wide, population-wide average and observed and expected heterozygosity rate.
  - `FourQuadrantsPopulationSizeProbe()`; (line 4652-4782 in cpp-file)
    This functions counts how many voles are present in each quadrant of four at sampling time, divided by males and females. Also counts how many voles are excluded by being in the buffer zone.
  - `NineQuadrantsPopulationSizeProbe()`; (line 4783-4902 in cpp-file)
  This functions counts how many voles are present in each of nine quadrant at sampling time, divided by males and females.
  - `MakeAlleleInputFile()`; (line 4469-4540 in the cpp-file)
    This function is to be called very early in the population initiation (before creating voles) to generate allele frequencies for each chromosome and locus. This funciton takes two sets of two alpha parameters. Half of the chromosomes' allele frequencies are generated according to the first set of alpha parameters, and the second half is generated according to the second set of alpha parameters. The alpha-parameters are used in making a beta distribution from which an allele frequency is drawn. In this project, we use di-allelic loci, so the allele frequency of the other allele is then given by the 1-first allele frequency. 
  - `CompareGenomes()`; (line in 4543-4650 cpp-file)
    This function takes the full genome of each vole recorded by `Four_QuadrantBasedGeneticOutput()`; og `Nine_QuadrantBasedGeneticOutput()`; samples 500 random voles, and compares their genomes. For each locus, a score is given based on similarity between the two voles.
    * If they share both alleles: 2 points
    * If they share one allele: 1 point
    * If they share no allele: 0 points. <br>
    The sum is then divided by the max possible score (256) to get a normalized similarity score between 0 and 1. Outputs a file   of   comma seperates values of the position of the voles, the distance between them, and their similarity score.
  - `GenerationCountOutput()`; (line 4904-4953 in the cpp-file)
    Outputs the generation count for the voles.
  - `LineagesOutput()`; (line 4955-5008 in the ccp-file)
    Outputs the paternal and maternal lineages of the voles at a provided sampling rate (here: 1).
  - `Four_AssignQuadrant();` ( line 5010-5032 in the cpp-file)
    Is a function that gets called by `Four_QuadrantBasedGeneticOutput();` and FourQuadrantsPopulationSizeProbe();`       for putting each vole into a quadrant (0-3) in the 2*2 quadrant grid with buffer zones.
  - `Nine_AssignQuadrant();` (line 5034-5062 in the cpp-file)
    Is a function that gets called by `Nine_QuadrantBasedGeneticOutput();` and NineQuadrantsPopulationSizeProbe();`       for putting each vole into a quadrant (0-3) in the 2*2 quadrant grid with buffer zones.

* Added to the PopulationManager functions `PopulationManager::CreateObjects()` (line 2750-2876 in the cpp-file) and `PopulationManager::CreateObjectsInit()`(line 2883-2931 in the cpp-file):
  - A vole that is being born from parents keeps the genome from its parents' gametes, so it calls the construction function (`GeneticMaterialNew()`;) with an empty string for the filename of the allele frequency input.
  - A vole under initialization calls the function with a string that contains the allele frequency file and will get its own genome based on allele frequency inputs.
__________________________________________________________
In the file-duo `Vole_all.cpp` and `Vole_all.h`:  
* Each vole `class struct_Vole_Adult` (line 122-152 in the header file) was given the attributes:
  - `new_Genes`, an instance of `GeneticMaterialNEW`
  - `new_Mates_Genes`, an instance of `GeneticMaterialNEW`
  - `yborn` (y coordinate of birth)
  - `xborn` (x coordinate of birth)
  - `GenerationCount`
  - `mitochondrial lineage`
  - `ychromosome lineage`
  
* Added to the `Vole_Female::st_Lactating()` function (which is used for giving birth) (line 2067-2210 in the cpp-file):
  - Voles record their birthplace into yborn and xborn.
  - Voles inherit a generation count +1 from their mother.
  - The mother vole recombines her genes and the father vole's genes with the function `MakeGameteSimple()`.
  - The mother passes on her mitochondrial lineage and the father's ychromosome lineage to the offspring.
* Numerous small adjustments to enable inheritance in many different constructor functions.
___________________________________________________________
In the file-duo `Vole_toletoc.cpp` and `Vole_toletoc.h`: <br>
In this file-duo, only two minor changes were:
* In the function `vole_tole_move_quality` (line 5-348 in the cpp-file) that when a vole looks up the movement quality of a large road, it will in 10% of lookup cases draw the value "1" (the vole can then walk onto the road) and in 90% of the cases draw the value "-1" (the vole cannot walk onto to road).
* In the `vole_tole_assess_barrier` function (line 818- 966 in the cpp-file), the boolean barrier status of the landscape category "large road" was set as "false" to enable movement across large roads.
