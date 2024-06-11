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
<B>GeneticMaterial.cpp This file contains the source for the genetic material classes</B> \n
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
//

//#include <stdio.h>
#include "../Landscape/ls.h"
#include "../Vole/GeneticMaterial.h"
//#include "../Vole/random_device.hpp"
#include "../ALMaSS_all/BatchALMaSS/ALMaSS_Random.h"
#include <chrono>
#include <random>
#include "../Vole/vole_all.h"

extern CfgBool cfg_Fixed_random_sequence;
extern CfgInt cfg_FixedRandomSeed;


CfgStr filename_AlleleInputCfg("ALLELE_FREQ_FILE", CFG_CUSTOM, "AlleleFreqs.txt");

CfgInt ChromosomeCountCfg("CHROMO_COUNT", CFG_CUSTOM, 4);
CfgInt LocusCountCfg("LOCUS_COUNT", CFG_CUSTOM, 32); //apparently does not actually read from the config, just reads from here..
CfgInt DiffBaseCountCfg("DIFF_BASES", CFG_CUSTOM, 2);

int ChromosomeCounting = static_cast<int>(ChromosomeCountCfg.value());
int LocusCounting = static_cast<int>(LocusCountCfg.value());
int DiffBaseCounting = static_cast<int>(DiffBaseCountCfg.value());

string filename_AlleleInput= filename_AlleleInputCfg.value();
//---------------------------------------------------------------------------
GeneticMaterialNEW::GeneticMaterialNEW(const string& filename_AlleleInput, bool new_Genes_empty) {
  Genome0 = vector<vector<int>>(ChromosomeCounting, vector<int>(LocusCounting, -1)); //init the genes empty for filling out (if the genes were already empty before being passed to this function)
  Genome1 = vector<vector<int>>(ChromosomeCounting, vector<int>(LocusCounting, -1));
  if (new_Genes_empty) { //the genes were empty before the constructor func was called
    if (!filename_AlleleInput.empty()) { //if a filename was specified (if we desire to generate a genoms)
      InitializeGenomes(filename_AlleleInput,Genome0,Genome1); //call the init function on the two newly made empty genes, to be filled out from the file
    }
  }
}

void GeneticMaterialNEW::EraseGenomes(){ //just set everything to -1 
  Genome0 = vector<vector<int>>(ChromosomeCounting, vector<int>(LocusCounting, -1));
  Genome1 = vector<vector<int>>(ChromosomeCounting, vector<int>(LocusCounting, -1));
};


// Read in chromosome and locus-combination specific probabilites and generate alleles randomly from that into 2 genomes (0 and 1).
// The CSV file needs the correct amount of rows (ChromosomeCounting*LocusCounting) and correct amount of columns (DiffBasesCount).
// This function could be made less sensitive to needing the correct amount of rows and cols by making a nested for loop instead of while(getline)...
void GeneticMaterialNEW::InitializeGenomes(const string& filename_AlleleInput, vector<vector<int>>& Genome0, vector<vector<int>>& Genome1) {
    const string filename=filename_AlleleInput;
    ifstream MyReadFile(filename); //check if file can be opened
      if (!MyReadFile.is_open()) {
        cout << "Error: Unable to open file " << filename << "\n";
      }

    string myText; //string for holding each line (which is a serios of locus specific allele frequencies)
    int i = 0; //keep track of iterations
    while (getline(MyReadFile, myText)) { //iterate through lines which all describe a set of chomo and locus specific probabilites of alleles
        istringstream ss(myText); //make a string stream of the line (each locus)
        vector<float> Floats; // vector to hold extracted floats (allele frequencies)
        string Token; //to hold one float at a time while iterating through the line (the locus)
        while (getline(ss, Token, ',')) { //iterating through the line (the locus)
            try {
                Floats.push_back(stof(Token)); //put each found float into the Floats vector
            } catch (...) {
                cerr << "Error: Failed to convert token to float: " << Token << endl;
            }
        }
        
        if (Floats.size() == 2) { //control if the previous steps went well
          float float1 = Floats[0]; //extract single floats
          float float2 = Floats[1];
          vector<float> Weights = {float1, float2}; //making floats into weights for discrete distribution :)

          std::ostringstream weigts_corr_syn; //for streaming the Weights vector and converting it to have " " between the floats for the distribution syntax
          for (auto it = Weights.begin(); it != Weights.end(); ++it) { //making syntax of weights for distribution
            weigts_corr_syn << *it;
            if (it != Weights.end() - 1) { // Add space if it's not the last element
              weigts_corr_syn << " ";
            }
          }
          probability_distribution p1("DISCRETE", weigts_corr_syn.str()); //now let's make a discrete dist with the weights from the file (locus specific)
          auto seed = std::chrono::steady_clock::now().time_since_epoch().count(); //get a random seed
          int result0=p1.Geti(); // draw 0 or 1 based from the p1 distribution
          int result1=p1.Geti(); // draw 0 or 1 based from the p1 distribution

          int curr_ch = i / LocusCounting; // indexing
          int curr_loc = i % LocusCounting;
          Genome0[curr_ch][curr_loc] = result0; //Put generated value into genomes
          Genome1[curr_ch][curr_loc] = result1;
        }else {
          cerr << "Error: Incorrect number of floats in line: " << myText << endl;
      }
        i++; //increase nr of iterations to keep track of indexing    
    };
}
//makes a 2 genomes for a male, one with only 0, one with only 1, for testing purposes
void GeneticMaterialNEW::InitializeTestingGenomesMale() {
  Genome0 = vector<vector<int>>(ChromosomeCounting, vector<int>(LocusCounting, 0));
  Genome1 = vector<vector<int>>(ChromosomeCounting, vector<int>(LocusCounting, 1));      
}

//makes a 2 genomes for a female, one with only 2, one with only 3, for testing purposes
void GeneticMaterialNEW::InitializeTestingGenomesFemale() {
  Genome0 = vector<vector<int>>(ChromosomeCounting, vector<int>(LocusCounting, 2));
  Genome1 = vector<vector<int>>(ChromosomeCounting, vector<int>(LocusCounting, 3));
      
}

//For printing your genomes, probably mostly for debugging :) 
void GeneticMaterialNEW::PrintGenomes(){
    cout << " \nGenome0:" << endl;
    for (const auto& row : Genome0) { //go trough chromos
        for (int value : row) { //go through loci
            cout << value << " ";
        }
        cout << endl;
    }
    
    cout << "\n Genome1:" << endl;
    for (const auto& row : Genome1) {
        for (int value : row) {
            cout << value << " ";
        }
        cout << endl;
    }
}

//this function is resposible for making the two gametes that go into an offspring during reproduction (makes 1 at a time....)
vector<vector<int>> GeneticMaterialNEW::MakeGameteSimple() const {
    vector<vector<int>> Gamete(ChromosomeCounting, vector<int>(LocusCounting, -1)); //init "empty" gamete
    for (int i = 0; i < ChromosomeCounting; i++) { // iterate in a chromosome wise fashion
        for (int j = 0; j < LocusCounting; j++) { // through each locus as well
            int WhichGenome;
            //randomly drawing which chromo to get genes from
            WhichGenome = rand() % 2; //get either 0 or 1 (by generating radnom nr, div by 2, find remainder (either 0 or 1))
            //extract the from genome0 or genome1 the desired base.
            if (WhichGenome == 0) {
                Gamete[i][j] = Genome0[i][j];
            } else if (WhichGenome == 1) {
                Gamete[i][j] = Genome1[i][j];
            }
        }
    }
    return Gamete;
}
//--------------------------------- THE END of my stuff (anastasia)
double MutationChance;
unsigned char g_MaxAllele;

AlleleFreq::AlleleFreq( ) {
  FILE * FreqFile = fopen("GeneticFrequencies.txt", "r" );
  int data;
  if ( !FreqFile ) {
    g_msg->Warn( "GeneticFrequencies File missing", "" );
    exit( 0 );
  }
  for ( int i = 0; i < 16; i++ ) {
    for ( int j = 0; j < 4; j++ ) {
      fscanf( FreqFile, "%d", & data );
      AlleleNumber[ i ] [ j ] = data;
    }
  }
  for ( int i = 16; i < 32; i++ ) {
    for ( int j = 0; j < 16; j++ ) {
      fscanf( FreqFile, "%d", & data );
      AlleleNumber[ i ] [ j ] = data;
    }
  }
  fclose( FreqFile );
}

//---------------------------------------------------------------------------

/*

void AlleleFreq::Flush() {
  for ( int loc = 0; loc < 16; loc++ ) {
    HO[ loc ] = 0;
    HE[ loc ] = 0;
    NoAlleles[ loc ] = 0;
    for ( int al = 0; al < 4; al++ ) {
      AlleleFrequency[ loc ] [ al ] = 0;
      AlleleNumber[ loc ] [ al ] = 0;
    }
  }
  for ( int loc = 16; loc < 32; loc++ ) {
    HO[ loc ] = 0;
    HE[ loc ] = 0;
    NoAlleles[ loc ] = 0;
    for ( int al = 0; al < 16; al++ ) {
      AlleleFrequency[ loc ] [ al ] = 0;
      AlleleNumber[ loc ] [ al ] = 0;
    }
  }
}

//---------------------------------------------------------------------------


void AlleleFreq::CalcNoAlleles() {
  // Counts the number of alleles existing for each locus
  for ( int loc = 0; loc < 16; loc++ ) {
    for ( int al = 0; al < 4; al++ ) {
      if ( AlleleNumber[ loc ] [ al ] > 0 ) NoAlleles[ loc ] ++;
    }
  }
  for ( int loc = 16; loc < 32; loc++ ) {
    for ( int al = 0; al < 16; al++ ) {
      if ( AlleleNumber[ loc ] [ al ] > 0 ) NoAlleles[ loc ] ++;
    }
  }
}

//---------------------------------------------------------------------------


void AlleleFreq::CalcHE() {
  for ( int loc = 0; loc < 16; loc++ ) {
    float NAlleles = 0;
    // How many alleles for each locus
    for ( int al = 0; al < 4; al++ ) {
      NAlleles += ( float )AlleleNumber[ loc ] [ al ];
    }
    // Calculated the proportion
    for ( int al = 0; al < 4; al++ ) {
      AlleleFrequency[ loc ] [ al ] = ( float )AlleleNumber[ loc ] [ al ] / NAlleles;
    }
    //HE=2ab+2bc+2bd+2ac+2ad+2cd
    HE[ loc ] = ( 2 * AlleleFrequency[ loc ] [ 0 ] * AlleleFrequency[ loc ] [ 1 ] )
         + ( 2 * AlleleFrequency[ loc ] [ 1 ] * AlleleFrequency[ loc ] [ 2 ] )
         + ( 2 * AlleleFrequency[ loc ] [ 1 ] * AlleleFrequency[ loc ] [ 3 ] )
         + ( 2 * AlleleFrequency[ loc ] [ 0 ] * AlleleFrequency[ loc ] [ 2 ] )
         + ( 2 * AlleleFrequency[ loc ] [ 0 ] * AlleleFrequency[ loc ] [ 3 ] )
         + ( 2 * AlleleFrequency[ loc ] [ 2 ] * AlleleFrequency[ loc ] [ 3 ] );
  }
}

//---------------------------------------------------------------------------


void AlleleFreq::CalcHO( int si ) {
  for ( int loc = 0; loc < 16; loc++ ) {
    HO[ loc ] /= ( float )si;
  }
}

//---------------------------------------------------------------------------


void AlleleFreq::CalcAF() {
  for ( int loc = 0; loc < 16; loc++ ) {
    int NoAlleles = 0;
    for ( int al = 0; al < 4; al++ ) NoAlleles += AlleleNumber[ loc ] [ al ];
    for ( int al = 0; al < 4; al++ ) {
      AlleleFrequency[ loc ] [ al ] = AlleleNumber[ loc ] [ al ] / ( float )NoAlleles;
    }
  }
}

*/

//---------------------------------------------------------------------------

void GeneticMaterial::SetGeneticFlag() {
	SetAllele(0,1,0);
}
//---------------------------------------------------------------------------
void GeneticMaterial::SetDirectFlag() {
	SetAllele(0,1,1);
}
//---------------------------------------------------------------------------

void GeneticMaterial::UnsetGeneticFlag() {
	SetAllele(0,0,0);
}
//---------------------------------------------------------------------------
void GeneticMaterial::UnsetDirectFlag() {
	SetAllele(0,0,1);
}
//---------------------------------------------------------------------------

uint32 GeneticMaterial::GetGeneticFlag() {
	return GetAllele(0,0);
}
//---------------------------------------------------------------------------
uint32 GeneticMaterial::GetDirectFlag() {
	return GetAllele(0,1);
}
//---------------------------------------------------------------------------

void GeneticMaterial::SetAllele( int locus, uint32 value, int Chromo ) {
  Chromo*=3; // now 0 or 3
  uint32 mask;
  if (locus<16) {
  // Get the right chromosome
  // Create the mask
  // Does it twice because 32 bits coding for 16 loci
  mask = 0x03 << locus;
  mask = mask << locus;
  // just to make make sure it is 0-3
  value = value & 0x03;
  // create the value mask
  value = value << locus;
  value = value << locus;
  // clear the locus
  Chromosome[ Chromo ] &= ~mask;
  // write the value
  Chromosome[ Chromo ] |= value;
  } else {
    Chromo++; // now 1 or 4
    locus-=16;
    if (locus>=8) {
      Chromo++;
      locus-=8;
    }
    mask = 0x0F << (locus*4);
    value = value & 0x0f; // make sure there was no extra stuff added!
    // create the value mask
    value = value << (locus*4);
    Chromosome[ Chromo ] &= ~mask;
    // write the value
    Chromosome[ Chromo ] |= value;
    }
}

//---------------------------------------------------------------------------

uint32 GeneticMaterial::GetAllele( int locus, int Chromo ) {
  uint32 value;
  // Get the right chromosome
  // if Chromo==0 then 0-2, else 3-5
  Chromo *= 3; // 0 or 3
  if ( locus < 16 ) {
    // Shift it so the locus is in the last two bits
    // Does it twice because 32 bis coding for 16 loci
    value = Chromosome[ Chromo ] >> locus;
    value = ( value >> locus ) & 0x03;
  } else {
    Chromo++; // 1 or 4
    locus -= 16; // Now 0 to 16
    if ( locus >= 8 ) {
      Chromo++; // 2 or 5
      locus -= 8;
    }
    value = Chromosome[ Chromo ] >> ( locus * 4 );
    value = value & 0x0f;
  }
  return value;
}

//---------------------------------------------------------------------------

void GeneticMaterial::PrintChromosome( char * C, int Chromo ) {
  for ( int i = 0; i < 16; i++ ) {
    uint32 allele = GetAllele( i, Chromo );
    switch ( allele ) {
      case 0:
        C[ i ] = 'a';
      break;
      case 1:
        C[ i ] = 'b';
      break;
      case 2:
        C[ i ] = 'c';
      break;
      case 3:
        C[ i ] = 'd';
      break;
      case 4:
        C[ i ] = 'e';
      break;
      case 5:
        C[ i ] = 'f';
      break;
      case 6:
        C[ i ] = 'g';
      break;
      case 7:
        C[ i ] = 'h';
      break;
      case 8:
        C[ i ] = 'i';
      break;
      case 9:
        C[ i ] = 'j';
      break;
      case 10:
        C[ i ] = 'k';
      break;
      case 11:
        C[ i ] = 'l';
      break;
      case 12:
        C[ i ] = 'm';
      break;
      case 13:
        C[ i ] = 'n';
      break;
      case 14:
        C[ i ] = 'o';
      break;
      case 15:
        C[ i ] = 'p';
      break;
    }
  }
  C[ 16 ] = 0;
}

//---------------------------------------------------------------------------

int GeneticMaterial::HomozygosityCount() {
  // OK OK there is an easy way to do this by calling HeterozygosityCount and
  // subtracting this from 32, but just is case that little bit of saved time is useful:
  int homozyg=0;
  for ( int i = 0; i < 32; i++ ) {
    if ( GetAllele( i, 0 ) == GetAllele( i, 1 ) ) homozyg++;
  }
  return homozyg;
}
//---------------------------------------------------------------------------

int GeneticMaterial::HeterozygosityCount() {
  int heterozyg = 0;
  for ( int i = 0; i < 32; i++ ) {
    if ( GetAllele( i, 0 ) != GetAllele( i, 1 ) ) heterozyg++;
  }
  return heterozyg;
}

//---------------------------------------------------------------------------

void GeneticMaterial::Recombine( GeneticMaterial * Gene1, GeneticMaterial * Gene2 ) {
  for ( int i = 0; i < 32; i++ ) {
    // For each locus
    // Choose which chromosome for each parent
    int g0 = g_random_fnc( 2 );
    int g1 = g_random_fnc( 2 );
    // get the two alleles
    uint32 a0 = Gene1->GetAllele( i, g0 );
    uint32 a1 = Gene2->GetAllele( i, g1 );
    //  put a0 into chromo0 & a1 to chromo1 & vice versa
    SetAllele( i, a0, 0 );
    SetAllele( i, a1, 1 );
  }
}

//---------------------------------------------------------------------------

GeneticMaterial::GeneticMaterial() {
  // ensure zeros in all loci
  for ( int i = 0; i < 6; i++ ) Chromosome[ i ] = 0;
}

//---------------------------------------------------------------------------

/**
The method called to intialise genes on initiation of the simulation. \n
Gene frequencies are based on an external text file input read in on construction.
*/
void GeneticMaterial::Initiation( AlleleFreq * Al ) {
  uint32 value;
  for ( int l = 0; l < 32; l++ ) {
    //if ( l < 16 ) c = 0; else if ( l < 24 ) c = 1; else c = 2;

    int chance = g_random_fnc( 1000 );
    uint32 index = 0;
    while ( chance > Al->SupplyAN( l, index ) ) {
      index++;
    }
    value = index;
    // set the value
    SetAllele( l, value, 0 );
    chance = g_random_fnc( 1000 );
    index = 0;
    while ( chance > Al->SupplyAN( l, index ) ) {
      index++;
    }
    value = index;
    // set the value
    SetAllele( l, value, 1 );
  }
}

//---------------------------------------------------------------------------
/**
This function can be used to alter reproductive effects based on genetic codes. These are only used in population genetic research.
*/
float GeneticMaterial::ScoreReproduction() {
  return 1.0;
  /* OLD CODE OUTDATED uint32 allele0a = GetAllele(0,0); // loci 0 uint32 allele1a = GetAllele(1,0); // loci 1
  uint32 allele0b = GetAllele(0,1); // loci 0 uint32 allele1b = GetAllele(1,1); // loci 1
  // Initial rules are that locus 0 and 1 are ac or bd then OK // likewise they may be ca or db
  // any other combination is bad // Same for loci 1 bool IsOK0=false; switch(allele0a) { case 0:   //a
  if (allele1a==2) IsOK0=true; break; case 1:   //b if (allele1a==3) IsOK0=true; break; case 2:   //c
  if (allele1a==0) IsOK0=true; break; case 3:   //d if (allele1a==1) IsOK0=true; break; default:
  FILE* errfile=fopen("GeneticErrorFile.Txt","w"); fprintf(errfile,"Unknown Allele Number\n"); fclose(errfile); exit(10);
  break; } bool IsOK1=false; switch(allele0b) { case 0: if (allele1b==2) IsOK1=true; break; case 1:
  if (allele1b==3) IsOK1=true; break; case 2: if (allele1b==0) IsOK1=true; break; case 3: if (allele1b==1) IsOK1=true; break;
  default: FILE* errfile=fopen("GeneticErrorFile.Txt","w"); fprintf(errfile,"Unknown Allele Number\n"); fclose(errfile);
  exit(11); break; } // determine the effect of the genetics // In the simple case it is good or bad
  if (IsOK0 && IsOK1) return 1.0; else return 0.05; */
}

//---------------------------------------------------------------------------
/**
This function can be used to alter fitness based on associated genetic codes. These are only used in population genetic research, e.g. to create hybrid zones.
*/
float GeneticMaterial::ScoreHQThreshold() {
  return 1.0;
  // OLD CODE OUTDATED
  /* Ditte's Simulation Version uint32 allele0a = GetAllele(0,0); // loci 0 uint32 allele1a = GetAllele(1,0); // loci 1
  uint32 allele0b = GetAllele(0,1); // loci 0 uint32 allele1b = GetAllele(1,1); // loci 1
  // Initial rules are that if 0a and 1a are 0 & 2 or 2 & 0 then OK   (a,c)(c,a)
  // likewise they may be 1 & 3 or 3 & 1  (b,d) (d,b) // any other combination is bad // Same for loci 1 bool IsOK0=false;
  switch(allele0a) { case 0: if (allele1a==2) IsOK0=true; break; case 1: if (allele1a==3) IsOK0=true; break; case 2:
  if (allele1a==0) IsOK0=true; break; case 3: if (allele1a==1) IsOK0=true; break; default: assert(NULL); break; }
  bool IsOK1=false; switch(allele0b) { case 0: if (allele1b==2) IsOK1=true; break; case 1: if (allele1b==3) IsOK1=true; break;
  case 2: if (allele1b==0) IsOK1=true; break; case 3: if (allele1b==1) IsOK1=true; break; default: assert(NULL); break; }
  // determine the effect of the genetics // In the simple case it is good or bad
  if (IsOK0 && IsOK1) return 1.0; else return 0.9;

  // Lar Bach's version uint32 allele2a = GetAllele(2,0); // loci 0 uint32 allele2b = GetAllele(2,1); // loci 0
  float result= -0.5; switch (allele2a) { case 0: case 2: case 3: break; default: result+=0.5; } switch (allele2b) { case 0:
  case 2: case 3: break; default: result+=0.5; break; } return result;
  // returns -0.5 if homozygous aa, 0.5 if homozygous bb & het=0

  */
}

//---------------------------------------------------------------------------

/**
random allele choice
*/
// Each locus is checked and if mutation chance equals 1 the allele and chromosome
//is chosen at random

void GeneticMaterial::Mutation_1()
{
  for ( int i = 0; i < 16; i++ ) {
    if ( g_rand_uni_fnc() < MutationChance ) // one chance in Mutation Chance
    {
      SetAllele( i, g_random_fnc( 4 ), g_random_fnc( 2 ) );
    }
  }
  for ( int i = 16; i < 32; i++ ) {
    if ( g_rand_uni_fnc() < MutationChance ) // one chance in Mutation Chance
    {
      SetAllele( i, g_random_fnc( 16 ), g_random_fnc( 2 ) );
    }
  }
}

//---------------------------------------------------------------------------

/**
random allele choice a & b only
*/
void GeneticMaterial::Mutation_1ab() // random allele choice
{
  // Only used when all loci have only two alleles!
  for ( int i = 0; i < 32; i++ ) {
    if ( g_rand_uni_fnc() < MutationChance ) // one chance in Mutation Chance
    {
      SetAllele( i, g_random_fnc( 2 ), g_random_fnc( 2 ) );
    }
  }
}

//---------------------------------------------------------------------------

/**
Move one allele +/-
*/
void GeneticMaterial::Mutation_2() // move one allele +/-
{
  for ( int i = 0; i < 16; i++ ) {
    if ( g_rand_uni_fnc() < MutationChance ) // one chance in Mutation Chance
    {
      int strand = g_random_fnc( 2 );
      int allele = GetAllele( i, strand );
      if ( g_random_fnc( 2 ) == 1 ) allele++; else allele--;
      if ( allele == -1 ) allele = 3; else if ( allele == 4 ) allele = 0;
      SetAllele( i, allele, strand );
    }
  }
  for ( int i = 16; i < 32; i++ ) {
    if ( g_rand_uni_fnc() < MutationChance ) // one chance in Mutation Chance
    {
      int strand = g_random_fnc( 2 );
      int allele = GetAllele( i, strand );
      if ( g_random_fnc( 2 ) == 1 ) allele++; else allele--;
      if ( allele == -1 ) allele = 15; else if ( allele == 16 ) allele = 0;
      SetAllele( i, allele, strand );
    }
  }
}

//---------------------------------------------------------------------------

/**
switch a<->b & c<->d
*/
void GeneticMaterial::Mutation_3()
{
  // NB Only works for the first 16 loci
  for ( int i = 0; i < 16; i++ ) {
    if ( g_rand_uni_fnc() < MutationChance ) // one chance in Mutation Chance
    {
      int strand = g_random_fnc( 2 );
      int allele = GetAllele( i, strand );
      switch ( allele ) {
        case 0:
          allele = 1;
        break;
        case 1:
          allele = 0;
        break;
        case 2:
          allele = 3;
        break;
        case 3:
          allele = 2;
        break;
      }
      SetAllele( i, allele, strand );
    }
  }
}
//---------------------------------------------------------------------------

void GeneticMaterial::Mutation_4()
{
	/** Specially mutation of only the first locus with two options 0/1 */
	if (g_rand_uni_fnc() < MutationChance) // one chance in Mutation Chance
	{
		int strand = g_random_fnc(2);
		int allele = GetAllele(0, strand);
		allele = (allele + 1) && 1; // only 1 or 0 is allowed
	}
}

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
// GeneticMaterial1616
//---------------------------------------------------------------------------

// const int MutChance = 2000;

AlleleFreq1616::AlleleFreq1616( ) {
  FILE * FreqFile = fopen("GeneticFrequencies_1616_Mut.txt", "r" );
  int data;
  if ( !FreqFile ) {
    g_msg->Warn( "GeneticFrequencies_1616_Mut.txt File missing - AllelFreq1616", "" );
    exit( 0 );
  }
  for ( int i = 0; i < 16; i++ ) {
    for ( int j = 0; j < 16; j++ ) {
      fscanf( FreqFile, "%d", & data );
      AlleleNumber[ i ] [ j ] = data;
    }
  }

  fclose( FreqFile );
}
//---------------------------------------------------------------------------
AlleleFreq256_16::AlleleFreq256_16( ) {
  FILE * FreqFile = fopen("GeneticFrequencies_256_16_Mut.txt", "r" );
  int data;
  if ( !FreqFile ) {
    g_msg->Warn( "GeneticFrequencies_256_16_Mut.txt File missing - AllelFreq1616", "" );
    exit( 0 );
  }
  for ( int i = 0; i < 256; i++ ) {
    for ( int j = 0; j < 16; j++ ) {
      fscanf( FreqFile, "%d", & data );
      AlleleNumber[ i ] [ j ] = data;
    }
  }

  fclose( FreqFile );
}
//---------------------------------------------------------------------------
GeneticMaterial1616::GeneticMaterial1616() {
  // ensure zeros in all loci
  for ( int i = 0; i < 4; i++ ) Chromosome[ i ] = 0;
}

//---------------------------------------------------------------------------

uint32 GeneticMaterial1616::GetAllele( unsigned int locus, unsigned int Chromo ) {
	// This is for 32 bit machines, 64bit is easier
	// locus must be 0 to 15
	//Chromo must be either 0 or 1
	// These debug tests below are costly so turn off in release code
	#ifdef __GENDEBUG
		if (Chromo>1) {
			g_msg->Warn( "Chromo > 1 in GeneticMaterial1616 - get allele", NULL );
			exit( 0 );
		}
		if (locus>15) {
			g_msg->Warn( "locus > 15 in GeneticMaterial1616 - get allele", NULL );
			exit( 0 );
		}
	#endif

    uint32 segment=((Chromo<<1) | ((locus & 0x08)>>3)); // locates on which segment in the comosome we are in
	uint32 allele=0x0F & (Chromosome[segment]>>((locus & 0x07)<<2)); //Locates the allele and shifts it down to the first postition to read
	return allele;
}
//---------------------------------------------------------------------------

void GeneticMaterial1616::SetAllele( unsigned int locus, uint32 value, unsigned int Chromo ) {
	// This is for 32 bit machines, 64bit is easier
	// locus must be 0 to 15
	//Chromo must be either 0 or 1
	// These debug tests below are costly so turn off in release code
	#ifdef __GENDEBUG
		if (Chromo>1) {
			g_msg->Warn( "Chromo > 1 in GeneticMaterial1616 - set allele", NULL );
			exit( 0 );
		}
		if (locus>15) {
			g_msg->Warn( "locus > 15 in GeneticMaterial1616 - set allele", NULL );
			exit( 0 );
		}
	#endif

    uint32 segment=((Chromo<<1) | ((locus & 0x08)>>3));
	uint32 mask = 0x0F; //0000 1111F
	// Need to shift the mask over the correct allele
	mask=mask<<((locus&7)<<2);
	value = value & 0x0f; // make sure there was no extra stuff added!
    // create the value mask
    value = value << ((locus&7)<<2);
    Chromosome[ segment ] &= ~mask; // get rid of the current info
    Chromosome[ segment ] |= value; // write the value
}
//---------------------------------------------------------------------------

void GeneticMaterial1616::PrintChromosome( char * C, unsigned int Chromo ) {
  for ( int i = 0; i < 16; i++ ) {
    uint32 allele = GetAllele( i, Chromo );
    switch ( allele ) {
      case 0:
        C[ i ] = 'a';
      break;
      case 1:
        C[ i ] = 'b';
      break;
      case 2:
        C[ i ] = 'c';
      break;
      case 3:
        C[ i ] = 'd';
      break;
      case 4:
        C[ i ] = 'e';
      break;
      case 5:
        C[ i ] = 'f';
      break;
      case 6:
        C[ i ] = 'g';
      break;
      case 7:
        C[ i ] = 'h';
      break;
      case 8:
        C[ i ] = 'i';
      break;
      case 9:
        C[ i ] = 'j';
      break;
      case 10:
        C[ i ] = 'k';
      break;
      case 11:
        C[ i ] = 'l';
      break;
      case 12:
        C[ i ] = 'm';
      break;
      case 13:
        C[ i ] = 'n';
      break;
      case 14:
        C[ i ] = 'o';
      break;
      case 15:
        C[ i ] = 'p';
      break;
    }
  }
  C[ 16 ] = 0;
}
//---------------------------------------------------------------------------

int GeneticMaterial1616::HomozygosityCount() {
  // OK OK there is an easy way to do this by calling HeterozygosityCount and
  // subtracting this from 32, but just is case that little bit of saved time is useful:
  int homozyg=0;
  for ( int i = 0; i < 16; i++ ) {
    if ( GetAllele( i, 0 ) == GetAllele( i, 1 ) ) homozyg++;
  }
  return homozyg;
}
//---------------------------------------------------------------------------

int GeneticMaterial1616::HeterozygosityCount() {
  int heterozyg = 0;
  for ( int i = 0; i < 16; i++ ) {
    if ( GetAllele( i, 0 ) != GetAllele( i, 1 ) ) heterozyg++;
  }
  return heterozyg;
}
//---------------------------------------------------------------------------------------

void GeneticMaterial1616::Recombine(GeneticMaterial1616 * Gene1, GeneticMaterial1616 * Gene2 ) {
	// Is called with hers and his genes
  for ( int i = 0; i < 16; i++ ) {
    // For each locus
    // Choose which chromosome for each parent
    int g0 = g_random_fnc(2);
    int g1 = g_random_fnc(2);
    // get the two alleles
    uint32 a0 = Gene1->GetAllele( i, g0 );
    uint32 a1 = Gene2->GetAllele( i, g1 );
    //  put a0 into chromo0 & a1 to chromo1 & vice versa
    SetAllele( i, a0, 0 );
    SetAllele( i, a1, 1 );
  }
}

//---------------------------------------------------------------------------
void GeneticMaterial1616::SetGeneticFlag() {
	SetAllele(0,1,0);
}
//---------------------------------------------------------------------------
void GeneticMaterial1616::SetDirectFlag() {
	SetAllele(0,1,1);
}
//---------------------------------------------------------------------------

void GeneticMaterial1616::UnsetGeneticFlag() {
	SetAllele(0,0,0);
}
//---------------------------------------------------------------------------
void GeneticMaterial1616::UnsetDirectFlag() {
	SetAllele(0,0,1);
}
//---------------------------------------------------------------------------

uint32 GeneticMaterial1616::GetGeneticFlag() {
	return GetAllele(0,0);
}
//---------------------------------------------------------------------------
uint32 GeneticMaterial1616::GetDirectFlag() {
	return GetAllele(0,1);
}
//---------------------------------------------------------------------------

/**
The method called to intialise genes on initiation of the simulation. \n
Gene frequencies are based on an external text file input read in on construction.
*/
void GeneticMaterial1616::Initiation( AlleleFreq1616 * Al ) {
  uint32 value; //, c;
  for ( int l = 0; l < 16; l++ ) {
   // if ( l < 16 ) c = 0; else c = 2;

    int chance = g_random_fnc( 1000 );
    uint32 index = 0;
    while ( chance >= Al->SupplyAN( l, index ) ) {
      index++;
    }
    value = index;
    // set the value
    SetAllele( l, value, 0 );

    chance = g_random_fnc( 1000 );
    index = 0;
    while ( chance >= Al->SupplyAN( l, index ) ) {
      index++;
    }
    value = index;
    // set the value
    SetAllele( l, value, 1 );
  }
}
//---------------------------------------------------------------------------

/**
random allele choice
*/
// Each locus is checked and if mutation chance equals 1 the allel and chromosome
//is chosen at random

void GeneticMaterial1616::Mutation_1()
{
	if (MutationChance != 0)
	{
  for ( int i = 0; i < 16; i++ )
  {
    if ( g_rand_uni_fnc() < MutationChance ) // one chance in Mutation Chance
    {
      SetAllele( i, g_random_fnc( 16 ), g_random_fnc( 2 ) );
	}
  }
	}
}

//---------------------------------------------------------------------------

/**
Move one allele + and 16 becomes 0
*/
void GeneticMaterial1616::Mutation_2() // move one allele +
{
	if (MutationChance != 0){
  for ( int i = 0; i < 16; i++ ) {
    if (g_rand_uni_fnc() < MutationChance)  // one chance in Mutation Chance for the locus
    {
      unsigned int strand = g_random_fnc( 2 ); // kromosom 0 or 1
      unsigned int allele = GetAllele( i, strand );
      allele++;
      allele&=0x0f; //and with 1111 - allels curls up 15 -> 0
      SetAllele( i, allele, strand );
    }
  }
}
}

//---------------------------------------------------------------------------

/**
Move one allele one up or down
*/
void GeneticMaterial1616::Mutation_3() // move one allele + or -
{
	if (MutationChance != 0){
  for ( int i = 0; i < 16; i++ ) {
    if ( g_rand_uni_fnc() < MutationChance) // one chance in Mutation Chance for the locus
    {
      unsigned strand = g_random_fnc( 2 ); // kromosom 0 or 1
      unsigned allele = GetAllele( i, strand );
      //if (strand == 1) allele++; else allele--;
	  if ( g_random_fnc( 2 ) == 1 ) allele++; else allele--;

	  //For mutations less than 0 and more than 15 the mutation should result in 1 or 14

	  if (allele > 0x0f) // Test to see if mutation to negative or above 15
	  {
		  allele |= 0x1C;
		  if(allele > 0x01C) { allele &= 0x01;}
		  else allele >>=1;
		  //allele|= 0xF0; // (240) Sets neg to 0, and 16 to 15 -> to avoid curling up
		  //allele = (~allele);
	  }

      SetAllele( i, allele, strand );
	}
  }
}
}

//---------------------------------------------------------------------------
GeneticMaterial256_16::GeneticMaterial256_16() {
  // ensure zeros in all loci
  for ( int i = 0; i < 32; i++ ) Chromosome[ i ] = 0;
}

//---------------------------------------------------------------------------

void GeneticMaterial256_16::SetAllele( unsigned int locus, uint32 value, unsigned int Chromo ) {
	if (Chromo==1) locus +=16;
	Chromosome[ locus ] = (unsigned char) value;
}
//---------------------------------------------------------------------------

uint32 GeneticMaterial256_16::GetAllele( unsigned int locus, unsigned int Chromo ) {
	if (Chromo==1) locus +=16;
	return (uint32) Chromosome[ locus ];
}
//---------------------------------------------------------------------------

/**
Move one allele one up or down
*/
void GeneticMaterial256_16::Mutation_3() // move one allele + or -
{
	if (MutationChance != 0){
		for ( int i = 0; i < 16; i++ ) {
			if ( g_rand_uni_fnc() < MutationChance ) // one chance in Mutation Chance for the locus
			{
				unsigned strand = g_random_fnc( 2 ); // kromosom 0 or 1
				int allele = GetAllele( i, strand );
				if ( g_random_fnc( 2 ) == 1 ) allele++; else allele--;
				//For mutations less than 0 and more than 256 the mutation should result in 1 or g_MaxAllele-1
				if (allele > g_MaxAllele) allele-=2;
				else if (allele < 0) allele = 1;
				SetAllele( i, (uint32) allele, strand );
			}
		}
	}
}
//---------------------------------------------------------------------------

void GeneticMaterial256_16::SetGeneticFlag() {
	SetAllele(0,1,0);
}
//---------------------------------------------------------------------------
void GeneticMaterial256_16::SetDirectFlag() {
	SetAllele(0,1,1);
}
//---------------------------------------------------------------------------

void GeneticMaterial256_16::UnsetGeneticFlag() {
	SetAllele(0,0,0);
}
//---------------------------------------------------------------------------
void GeneticMaterial256_16::UnsetDirectFlag() {
	SetAllele(0,0,1);
}
//---------------------------------------------------------------------------

uint32 GeneticMaterial256_16::GetGeneticFlag() {
	return GetAllele(0,0);
}
//---------------------------------------------------------------------------
uint32 GeneticMaterial256_16::GetDirectFlag() {
	return GetAllele(0,1);
}
//---------------------------------------------------------------------------

void GeneticMaterial256_16::Recombine(GeneticMaterial256_16 * Gene1, GeneticMaterial256_16 * Gene2 ) {
	// Is called with hers and his genes
  for ( int i = 0; i < 16; i++ ) {
    // For each locus
    // Choose which chromosome for each parent
    int g0 = g_random_fnc(2);
    int g1 = g_random_fnc(2);
    // get the two alleles
    uint32 a0 = Gene1->GetAllele( i, g0 );
    uint32 a1 = Gene2->GetAllele( i, g1 );
    //  put a0 into chromo0 & a1 to chromo1 & vice versa
    SetAllele( i, a0, 0 );
    SetAllele( i, a1, 1 );
  }
}
//---------------------------------------------------------------------------

/**
The method called to intialise genes on initiation of the simulation. \n
Gene frequencies are based on an external text file input read in on construction.
*/
void GeneticMaterial256_16::Initiation( AlleleFreq256_16 * Al ) {
  uint32 value; //, c;
  for ( int l = 0; l < 16; l++ ) {
    int chance = g_random_fnc( 1000 );
    uint32 index = 0;
    while ( chance >= Al->SupplyAN( l, index ) ) {
      index++;
    }
    value = index;
    // set the value
    SetAllele( l, value, 0 );
    chance = g_random_fnc( 1000 );
    index = 0;
    while ( chance >= Al->SupplyAN( l, index ) ) {
      index++;
    }
    value = index;
    // set the value
    SetAllele( l, value, 1 );
  }
}
//---------------------------------------------------------------------------
