/* 
*******************************************************************************************************
Copyright (c) 2011, Christopher John Topping, Aarhus University
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided
that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the
following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and 
the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
********************************************************************************************************
*/

/**
\file
\brief
<B>GeneticMaterial.h This file contains the headers for the genetic material classes</B> \n
*/
/**
\file 
 by Chris J. Topping \n
 Version of 1st Nov 2000 \n
 All rights reserved. \n
 \n
 With additions as noted in: \n
 Doxygen formatted comments in May 2008 \n
*/
//---------------------------------------------------------------------------
//CfgStr chromoFilePath1("CHROMO_FILE_PATH1", CFG_CUSTOM, "");

extern std::string filename_AlleleInput;
extern int ChromosomeCounting;
extern int LocusCounting;
extern int DiffBaseCounting;

#ifndef GeneticMaterialH
#define GeneticMaterialH


#ifdef __LINUX
#include "ALMaSSDefines.h"

#else
#include "../ALMaSSDefines.h"
//#include "../Vole/vole_all.h"

#endif
//---------------------------------------------------------------------------
//------------------------------------------------------------------------------

class Landscape;

//------------------------------------------------------------------------------

/**
\brief
Class to handle statistics and constructs based on allele frequencies
*/
class AlleleFreq
{
  protected:
    int AlleleNumber[32][16];
    float AlleleFrequency[32][16];
    float HE[32];
    float HO[32];
    int NoAlleles[32];
  public:
    AlleleFreq();
    int SupplyAN(int loc, int al) { return AlleleNumber[loc][al];}
/*
    void SetAF(int loc, int al, int val) { AlleleFrequency[loc][al]=val;}
    void IncAN(int loc, int al) { AlleleNumber[loc][al]++;}
    void IncAF(int loc, int al) { AlleleFrequency[loc][al]++;}
    void IncHO(int loc) { HO[loc]++;}
    float SupplyAF(int loc, int al) { return AlleleFrequency[loc][al];}
    float SupplyHE(int loc) { return HE[loc];}
    float SupplyHO(int loc) { return HO[loc];}
    int SupplyNoAlleles(int loc) { return NoAlleles[loc];}
    void CalcNoAlleles();
    void CalcAF();
    void CalcHE();
    void CalcHO(int si);
    void Flush();
*/
};
//------------------------------------------------------------------------------

class GeneticMaterialNEW {
public:
    vector<vector<int>> Genome0; //Two genome versions= Diploid organism.
    vector<vector<int>> Genome1;
    int MitochondrialLine;
    int YChromoLine;
    GeneticMaterialNEW(const string& filename_AlleleInput, bool new_Genes_empty); //constructor with input
    void InitializeGenomes(const string& filename, vector<vector<int>>& Genome0, vector<vector<int>>& Genome1); // reads in allele freqs from file and generates a genome for the vole
    void PrintGenomes();
    void EraseGenomes();
    void InitializeTestingGenomesFemale(); //will put all initalized female genomes to 2 (genome0) and 3 (genome1) for testing purposes
    void InitializeTestingGenomesMale(); //will put all initalized male genomes to 0 (genome0) and 1 (genome1) for testing purposes
    vector<vector<int>> MakeGameteSimple() const;
    GeneticMaterialNEW() { //default constructor
      Genome0 = vector<vector<int>>(ChromosomeCounting, vector<int>(LocusCounting, -1));
      Genome1 = vector<vector<int>>(ChromosomeCounting, vector<int>(LocusCounting, -1));
    }
  };
/**
\brief
Class for the genetic material optionally carried by animals in ALMaSS
*/
class GeneticMaterial
{
 protected:
   uint32 Chromosome[6];
 public:
   GeneticMaterial();
   void ReadFrequencies();
   void SetAllele(int pos, uint32 value, int Chromosome);
   uint32 GetAllele(int pos, int Chromosome);
   void PrintChromosome(char *C, int Chromosome);
   int HomozygosityCount();
   int HeterozygosityCount();
   void Recombine(GeneticMaterial* Gen21, GeneticMaterial* Gene2);
   void Initiation(AlleleFreq* Al);
   float ScoreReproduction();
   float ScoreHQThreshold();
	void SetGeneticFlag();
	void SetDirectFlag();
	void UnsetGeneticFlag();
	void UnsetDirectFlag();
	uint32 GetGeneticFlag();
	uint32 GetDirectFlag();
   void Mutation_1(); // random mutation
   void Mutation_1ab(); // random mutation only a&b
   void Mutation_2(); // mutation one step either way
   void Mutation_3(); // mutation from one state to another only (a->b) (b->c)}
   void Mutation_4(); // mutation from one state to another only (a->b) (b->a) for first locus
};

//---------------------------------------------------------------------------

class AlleleFreq1616
{
  protected:
    int AlleleNumber[16][16];
    float AlleleFrequency[16][16];
    float HE[16];
    float HO[16];
    int NoAlleles[16];
  public:
    AlleleFreq1616();
	int SupplyAN(int loc, int al) { return AlleleNumber[loc][al];}
/*
    void SetAF(int loc, int al, int val) { AlleleFrequency[loc][al]=val;}
    void IncAN(int loc, int al) { AlleleNumber[loc][al]++;}
    void IncAF(int loc, int al) { AlleleFrequency[loc][al]++;}
    void IncHO(int loc) { HO[loc]++;}
    float SupplyAF(int loc, int al) { return AlleleFrequency[loc][al];}
    float SupplyHE(int loc) { return HE[loc];}
    float SupplyHO(int loc) { return HO[loc];}
    int SupplyNoAlleles(int loc) { return NoAlleles[loc];}
    void CalcNoAlleles();
    void CalcAF();
    void CalcHE();
    void CalcHO(int si);
    void Flush();
*/
};
//------------------------------------------------------------------------------
class GeneticMaterial1616
{
 protected:
	uint32 Chromosome[4]; // = 32 loci (in 2 chromosomes) & 4 bits each
 public:
	GeneticMaterial1616();
	void SetAllele(unsigned int locus, uint32 value, unsigned int Chromo);
	uint32 GetAllele( unsigned int locus, unsigned int Chromo );
	void PrintChromosome(char *C, unsigned int Chromosome);
	void SetGeneticFlag();
	void SetDirectFlag();
	void UnsetGeneticFlag();
	void UnsetDirectFlag();
	uint32 GetGeneticFlag();
	uint32 GetDirectFlag();
	int HomozygosityCount();
	int HeterozygosityCount();
	void Recombine(GeneticMaterial1616* Gene1, GeneticMaterial1616* Gene2);
	void Initiation(AlleleFreq1616* Al);
	void Mutation_1(); // random mutation
	void Mutation_2(); // 'next in line' mutation + (16 becomes 0)
	void Mutation_3(); // 'next in line' mutation +/- (16 becomes 0 and -1 becomes 15)
};
//---------------------------------------------------------------------------


class AlleleFreq256_16
{
  protected:
    int AlleleNumber[256][16];
    float AlleleFrequency[16][256];
    float HE[16];
    float HO[16];
    int NoAlleles[256];
  public:
    AlleleFreq256_16();
	int SupplyAN(int loc, int al) { return AlleleNumber[al][loc];}
};
//------------------------------------------------------------------------------
class GeneticMaterial256_16
{
 protected:
	unsigned char Chromosome[32]; // = 32 loci (in 2 chromosomes) & 8 bits each
 public:
	GeneticMaterial256_16();
	void SetAllele(unsigned int locus, uint32 value, unsigned int Chromo);
	uint32 GetAllele( unsigned int locus, unsigned int Chromo );
	void Mutation_3(); // 'next in line' mutation +/- 
	void SetGeneticFlag();
	void SetDirectFlag();
	int HomozygosityCount() { return 0; } // To add if needed
	int HeterozygosityCount() { return 0; } // To add if needed
	void UnsetGeneticFlag();
	void UnsetDirectFlag();
	uint32 GetGeneticFlag();
	uint32 GetDirectFlag();
	void PrintGenes();
	void Recombine(GeneticMaterial256_16* Gene1, GeneticMaterial256_16* Gene2);
	void Initiation(AlleleFreq256_16* Al);
/*
	void Mutation_1(); // random mutation
	void Mutation_2(); // 'next in line' mutation + (16 becomes 0)
*/
};
//---------------------------------------------------------------------------

#endif
