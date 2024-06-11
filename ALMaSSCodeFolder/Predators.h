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
/** \file Predators.h
\brief <B>The header file for all predator lifestages and population manager classes</B>
*/
/**  \file Predators.h
Version of  28 January 2001. \n
By Chris J. Topping

*/

//---------------------------------------------------------------------------
#ifndef PredatorsH
#define PredatorsH
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

class TPredator;
class TPredator_Population_Manager;

//------------------------------------------------------------------------------
/**
Used for the population manager's list of predators
*/
//typedef vector<TPredator*> TListOfPredators;
//---------------------------------------------------------------------------

/**
Predators like other ALMaSS animals work using a state/transition concept. These are the predator behavioural states.
*/
typedef enum
{
      tops_InitialState=0,
      tops_Hunting,
      tops_Dispersal,
      tops_Movement
} TTypeOfPredatorState;

/**
\brief
Used for creation of a new predator object
*/
class struct_Predator
{
 public:
  int x;
  int y;
  int species;
  Landscape* L;
  TPredator_Population_Manager * PM;
};

/**
\brief
The base class for predators encompsassing all their general behaviours
*/
/**
Defines simple predators that are really nothing more than moving mortality probabilities of different sizes and strengths.
Breeding is once a year if enough prey are eaten, and death occurs if a starvation criteria is reached. The last individual cannot die so
the population can never go extinct.
If predators do not eat for a definable number of hunts then they may disperse, otherwise they move only locally.
Different types of predators can be defined in the same simulation by defining the Weasel and Owl classes using
configuration variables to create e.g. specialist or generalist predators, with different movement patterns, hunting efficiency and numerical responses.
*/
class TPredator : public TAnimal
{
   /*
   A predator must have some simple functionality:

   Search Area:  where it looks for prey
   HomeRange: containing all its current search areas
   Kill Efficiency: How good it is at finding and killing prey
   Movement: How much it moves around in its home range
   Dispersal: How often and how it changes its home range
   */

   // Inherits m_Location_x, m_Location_y, m_OurLandscape from TAnimal
   // NB All areas are squares of size length X length

protected:
   TTypeOfPredatorState CurrentPState;
   unsigned SpeciesID;
   unsigned m_DispersalMax;
   unsigned m_SearchArea;// length of side of a square area
   int m_kills_this_season;
   int PreyResponse1;
   int PreyResponse2;
   int m_Search_x;
   int m_Search_y;
   int SimW;
   int SimH;
   unsigned m_FailureCount;
   unsigned m_NoFailuresBeforeDispersal;
   unsigned m_HomeRange; // TL corner = m_Location_x,m_Location_y
   bool m_HaveTerritory;
   int m_KillEfficiency; // 0-1000 = 0-100%
   Vole_Population_Manager* m_Prey;
   TPredator_Population_Manager* m_OurPopulationManager;
   vector<Vole_Base*>* CurrentPrey;
public:
   TPredator(Vole_Population_Manager* ThePrey, int p_x, int p_y,
                           Landscape* p_L, TPredator_Population_Manager* p_PPM);
   ~TPredator();
   bool OverlapMyTerritory(unsigned x, unsigned y);
   virtual void st_Dispersal();
   virtual void st_Movement();
   virtual int st_Hunting();
   virtual void BeginStep       (void) {}
   virtual void Step            (void) {}
   virtual void EndStep         (void) {}

   unsigned SupplySpeciesID() {return SpeciesID;}
   int SupplyKill() {return m_kills_this_season;}
   bool SupplyTerr() {return m_HaveTerritory;}
   int SupplyKillEff() {return m_KillEfficiency;}
   int SupplyHomeRange() {return m_HomeRange;}
};

/**
\brief
The class to handle all predator population related matters
*/
class TPredator_Population_Manager : public Population_Manager
{
public:
// Methods
   virtual void Run(int);
   bool InOtherTerritory(unsigned sp, int p_x, int p_y, TPredator* p_Pred);
   TPredator_Population_Manager(Landscape* L,Vole_Population_Manager* VPM);
   virtual ~TPredator_Population_Manager (void);
   void CreateObjects(int ob_type, TAnimal *pvo,/*void* null ,*/
                                         struct_Predator* data,int number);
   void inc_inds(unsigned list) {m_no_individuals[list]++;}
   void dec_inds(unsigned list) {m_no_individuals[list]--;}
   unsigned supply_no_inds(unsigned list) {return m_no_individuals[list];}

protected:
// Attributes
   Vole_Population_Manager* m_Prey;
   unsigned NoPredatorTypes;
   unsigned m_no_individuals[2];

// Methods
   virtual bool StepFinished();
   virtual void DoFirst(){}
   virtual void DoBefore(){}
   virtual void DoAfter(){}
   virtual void DoLast(){}
   void CloseTheReallyBigOutputProbe() {}; // This will always be done by the main population manager
   void PredSampleFile();
   void PredAutumnSample();
   void PredSpringSample();
   void PredSpringAutumnSample();

};

/**
\brief
The Weasel class is one of two current implementations of TPredator.
*/
/**
It is configurable via config parameters and in other than name and default configuration it is identical to the Owl class
*/
class Weasel : public TPredator
{
public:
//   virtual void st_Movement()  {};
//   virtual int st_Hunting()   {};
//   virtual void st_Dispersal() {};
   virtual void BeginStep       (void);
   virtual void Step            (void);
   virtual void EndStep         (void) {}
   Weasel(Vole_Population_Manager* ThePrey, int p_x, int p_y, Landscape * p_L,
                                          TPredator_Population_Manager * p_PPM);
   virtual ~Weasel();
};

/**
\brief
The Owl class is one of two current implementations of TPredator.
*/
/**
It is configurable via config parameters and in other than name and default configuration it is identical to the Weasel class
*/
class Owl : public TPredator
{
 public:
//   virtual void st_Movement()  {};
//   virtual void st_Hunting()   {};
//   virtual void st_Dispersal() {};
   virtual void BeginStep       (void);
   virtual void Step            (void);
   virtual void EndStep         (void) {}
   Owl(Vole_Population_Manager* ThePrey, int p_x, int p_y, Landscape * p_L,
                                          TPredator_Population_Manager * p_PPM);
   virtual ~Owl();
};

//---------------------------------------------------------------------------
#endif
