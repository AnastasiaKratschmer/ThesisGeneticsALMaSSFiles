This folder contains the ALMaSS files that I have made changes to during my thesis.  
Here follows a description of what I have added to each file!  

In the file-duo `GeneticMaterial.cpp` and `GeneticMaterial.h`:  
* A new class "GeneticMaterialNEW" with functions:
  - A constructor function GeneticMaterialNew();
  - InitializeGenome();
  - PrintGenomes();
  - EraseGenomes();
  - InitializeTestingGenomesFemale();
  - InitalizeTestingGenomesMale();
  - MakeGameteSimple();

In the file-duo `VolePopulationManager.cpp` and `VolePopulationManager.h`:  
* New functions for the class PopulationManager:
  - Four_QuadrantBasedGeneticOutput();
  - Nine_QuadrantBasedGeneticOutput();
  - FourQuadrantsPopulationSizeProbe();
  - NineQuadrantsPopulationSizeProbe();
  - MakeAlleleInputFile();
  - CompareGenomes();
  - GenerationCountOutput();
  - LineagesOutput()

In the file-duo `Vole_all.cpp` and `Vole_all.h`:  
* Added to the Vole_Female::st_Lactating();-function (which is used for giving birth):
  - Voles record their birthplace.
  - Voles inherit a generation count +1 from their mother.
  - The mother vole takes her genes and and recombines them, and takes the father vole's genes and recombines those with the function MakeGameteSimple();
  - The mother also passes on her mitochondrial lineage and the fathers ychromosome lineage to the offspring.

In the file-duo `Vole_toletoc.cpp` and `Vole_toletoc.h`:  
* The only modification here is that when a vole looks up 
