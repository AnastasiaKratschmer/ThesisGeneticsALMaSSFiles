This folder contains the ALMaSS files that I have made changes to during my thesis.  
Here follows a description of what I have added to each file!  

In the file-duo `GeneticMaterial.cpp` and `GeneticMaterial.h`:  
* A new class `GeneticMaterialNEW` containing 2 genomes and a paternal and maternal lineage ID, with functions:
  - A constructor function `GeneticMaterialNew()`; <br>
Calls `InitializeGenome();`.
  - `InitializeGenome()`; <br>
    Randomly generates two genomes for the vole from a file of allele frequency input. <br>
    If an empty string is passed to the `InitializeGenome()`; the vole keeps the genome it already has. (For keeping the inherited genes from parents, if the vole is an offspring.
  - `PrintGenomes()`;
  - `EraseGenomes()`;
  - `InitializeTestingGenomesFemale()`;
  - `InitializeTestingGenomesMale()`;
  - `MakeGameteSimple()`;

In the file-duo `VolePopulationManager.cpp` and `VolePopulationManager.h`:  
* New functions for the class PopulationManager:
  - `Four_QuadrantBasedGeneticOutput()`;
  - `Nine_QuadrantBasedGeneticOutput()`;
  - `FourQuadrantsPopulationSizeProbe()`;
  - `NineQuadrantsPopulationSizeProbe()`;
  - `MakeAlleleInputFile()`;
  - `CompareGenomes()`;
  - `GenerationCountOutput()`;
  - `LineagesOutput()`;

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
* The only modification here is that when a vole looks up ...
