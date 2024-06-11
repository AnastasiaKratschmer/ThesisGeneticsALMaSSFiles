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
/** \file Predators.cpp
\brief <B>The main source code for all predator lifestage and population manager classes</B>
*/
/**  \file Predators.cpp
Version of  28 January 2001. \n
By Chris J. Topping \n \n
*/

//---------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include<vector>
#include "../Landscape/ls.h"
#include "../BatchALMaSS/PopulationManager.h"
#include "../Vole/GeneticMaterial.h"
#include "../Vole/vole_all.h"
#include "../Vole/Predators.h"
#include "../BatchALMaSS/BinaryMapBase.h"
#include "../BatchALMaSS/MovementMap.h"
#include "../Vole/VolePopulationManager.h"

//---------------------------------------------------------------------------

#define WEASEL 0
#define OWL 1

static CfgBool cfg_PredMortalityDataUsed("PRED_MORTALITY_DATA_USED",CFG_CUSTOM, false);
static CfgBool cfg_PredSampleDataUsed("PRED_SAMPLE_DATA_USED",CFG_CUSTOM, false);

static CfgInt  cfg_pred_first_sample_day("PRED_SAMPLE_FILE_DAY_ONE", CFG_CUSTOM,109);
static CfgInt  cfg_pred_second_sample_day( "PRED_SAMPLE_FILE_DAY_TWO", CFG_CUSTOM,287 );

static CfgInt  cfg_weasel_breed_threshold( "WEASEL_BT", CFG_CUSTOM,5 );
static CfgInt  cfg_owl_breed_threshold( "OWL_BT", CFG_CUSTOM,5000000 );
static CfgInt  cfg_weasel_death_threshold( "WEASEL_DT", CFG_CUSTOM,10 );
static CfgInt  cfg_owl_death_threshold( "OWL_DT", CFG_CUSTOM,-1 );
static CfgInt  cfg_weasel_breed_day( "WEASEL_BD", CFG_CUSTOM,115 );
static CfgInt  cfg_owl_breed_day( "OWL_BD", CFG_CUSTOM,120 );
static CfgInt cfg_weasel_kill_efficiency("WEASEL_KILL_EFF",CFG_CUSTOM,250);
static CfgInt cfg_weasel_home_range("WEASEL_HOME_RANGE",CFG_CUSTOM,250);
static CfgInt cfg_weasel_search_area("WEASEL_SEARCH_AREA",CFG_CUSTOM,100);
static CfgInt cfg_weasel_NoFailuresBeforeDispersal("WEASEL_FAILURES",CFG_CUSTOM,10);
static CfgInt cfg_weasel_DispersalMax("WEASEL_DISPERSAL_MAX",CFG_CUSTOM,200);
static CfgInt cfg_weasel_StartingNo("WEASEL_START_NO",CFG_CUSTOM,0);
static CfgInt cfg_owl_StartingNo("OWL_START_NO",CFG_CUSTOM,0);

// Constants for the predator species that need to be fast
int weasel_breed_threshold;
int owl_breed_threshold;
int weasel_death_threshold;
int owl_death_threshold;
int weasel_breed_day;
int owl_breed_day;
int owl_StartingNo;
int weasel_StartingNo;

//---------------------------------------------------------------------------

TPredator_Population_Manager::~TPredator_Population_Manager (void)
{
   // Should all be done by the Population_Manager destructor
}
//---------------------------------------------------------------------------

TPredator_Population_Manager::TPredator_Population_Manager(Landscape* L,
                                               Vole_Population_Manager* VPM)
 : Population_Manager(L , 2)
{
  
  m_ListNames[WEASEL] = "Weasel";
  m_ListNames[OWL] = "Owl";
  m_ListNameLength = 2;
  m_SimulationName = "Vole Predators";

  // Constants for the predator species
  weasel_breed_threshold=cfg_weasel_breed_threshold.value();
  owl_breed_threshold=cfg_owl_breed_threshold.value();
  weasel_death_threshold=cfg_weasel_death_threshold.value();
  owl_death_threshold=cfg_owl_death_threshold.value();
  weasel_breed_day=cfg_weasel_breed_day.value();
  owl_breed_day=cfg_owl_breed_day.value();
  weasel_StartingNo=cfg_weasel_StartingNo.value();
  owl_StartingNo=cfg_owl_StartingNo.value();
  // Remember the prey
  m_Prey=VPM;
  m_no_individuals[WEASEL]=0;
  m_no_individuals[OWL]=0;
  // Create some weasels and owls
  struct_Predator* sp;
  sp = new struct_Predator;
  sp->PM = this;
  sp->L = m_TheLandscape;
  m_population_type = TOP_Predators;
  for (int i=0; i<weasel_StartingNo; i++)
  {
    sp->x = g_random_fnc(SimW);
    sp->y = g_random_fnc(SimH);
    CreateObjects(0,NULL,sp,1); // 0 = weasel
  }
  for (int i=0; i<owl_StartingNo; i++)
  {
    sp->x = g_random_fnc(SimW);
    sp->y = g_random_fnc(SimH);
    CreateObjects(1,NULL,sp,1); // 1 = Owl
  }
  delete sp;
  ReallyBigOutputPrb=0;
  RipleysOutputPrb=0;
}

//---------------------------------------------------------------------------
void TPredator_Population_Manager::CreateObjects(int ob_type,
           TAnimal * ,struct_Predator * data,int number)
{
   Weasel*  new_Weasel;
   Owl*  new_Owl;
   for (int i=0; i<number; i++)
   {
    if (ob_type == WEASEL)   // Weasel
    {
       new_Weasel = new Weasel(m_Prey,data->x, data->y, data->L, data->PM);
       PushIndividual(ob_type, new_Weasel);
       inc_inds(WEASEL);
    }
    if (ob_type == OWL)  // Owl
    {
       new_Owl = new Owl(m_Prey,data->x, data->y, data->L, data->PM);
        PushIndividual(ob_type, new_Owl);
       inc_inds(OWL);
    }
   }
}

//---------------------------------------------------------------------------

void TPredator_Population_Manager::PredSampleFile(){

	int today = m_TheLandscape->SupplyDayInYear();
	int year = m_TheLandscape->SupplyYearNumber();
	int month = m_TheLandscape->SupplyMonth();
	int dayinMo = m_TheLandscape->SupplyDayInMonth();

	if ((cfg_weasel_StartingNo.value())&& (cfg_owl_StartingNo.value()> 0)){

		FILE* pref=fopen("PredSampleData.txt","a");
		if (pref == NULL) {
			m_TheLandscape->Warn("TPredator_Population_Manager::PredSampleFile","Could Not Open PredSampleData.txt File");
			exit(0);
		}

		Weasel* Wp;
		unsigned int SizeWeasel = (unsigned int) SupplyListSize(0);
		Owl* Op;
		unsigned int SizeOwl = (unsigned int) SupplyListSize(1);
		unsigned int TotalSize = SizeWeasel+SizeOwl;

		int w = m_TheLandscape->SupplySimAreaWidth();
		int h = m_TheLandscape->SupplySimAreaHeight();
		fprintf(pref,"%d\t %d\t %d\t %d\t %u\t %u\n", 0,w ,0, h, SizeWeasel, TotalSize);

		for (unsigned i=0; i<SizeWeasel; i++){

			Wp=dynamic_cast<Weasel*>(SupplyAnimalPtr(0, i));
            int Wx{0}, Wy{0};
            SupplyLocXY(0, i, Wx, Wy);

			int Wpoly = m_TheLandscape->SupplyPolyRef(Wx, Wy);
			int Wele = m_TheLandscape->SupplyElementType(Wpoly);
			int Wveg = m_TheLandscape->SupplyVegType(Wpoly);

			unsigned Wspecies = Wp->SupplySpeciesID();
			int kill = Wp->SupplyKill();
			int KillEff = Wp->SupplyKillEff();
			int Terr = Wp->SupplyTerr();
			int HomeR = Wp->SupplyHomeRange();

			fprintf(pref,"%d\t",year);
			fprintf(pref,"%d\t",month);
			fprintf(pref,"%d\t",dayinMo);
			fprintf(pref,"%d\t",today);

			fprintf(pref,"%d\t",Wspecies);
			fprintf(pref,"%d\t",kill);
			fprintf(pref,"%d\t",KillEff);
			fprintf(pref,"%d\t",Terr);
			fprintf(pref,"%d\t",HomeR);

			fprintf(pref,"%d\t",Wx);
			fprintf(pref,"%d\t",Wy);
			fprintf(pref,"%d\t",Wpoly);
			fprintf(pref,"%d\t",Wele);
			fprintf(pref,"%d\t",Wveg);
			fprintf(pref,"\n");
			}

		fprintf(pref,"%d\t %d\t %d\t %d\t %u\t %u\n", 0,w ,0, h, SizeOwl, TotalSize);

		for (unsigned i = 0; i < SizeOwl; i++){
			Op=dynamic_cast<Owl*>(SupplyAnimalPtr(1, i));
			int Ox{0}, Oy{0};
            SupplyLocXY(1, i, Ox, Oy);
			int Opoly = m_TheLandscape->SupplyPolyRef(Ox, Oy);
			int Oele = m_TheLandscape->SupplyElementType(Opoly);
			int Oveg = m_TheLandscape->SupplyVegType(Opoly);

			int species = Op->SupplySpeciesID();
			int kill = Op->SupplyKill();
			int KillEff = Op->SupplyKillEff();
			int Terr = Op->SupplyTerr();
			int HomeR = Op->SupplyHomeRange();

			fprintf(pref,"%d\t",year);
			fprintf(pref,"%d\t",month);
			fprintf(pref,"%d\t",dayinMo);
			fprintf(pref,"%d\t",today);

			fprintf(pref,"%d\t",species);
			fprintf(pref,"%d\t",kill);
			fprintf(pref,"%d\t",KillEff);
			fprintf(pref,"%d\t",Terr);
			fprintf(pref,"%d\t",HomeR);

			fprintf(pref,"%d\t",Ox);
			fprintf(pref,"%d\t",Oy);
			fprintf(pref,"%d\t",Opoly);
			fprintf(pref,"%d\t",Oele);
			fprintf(pref,"%d\t",Oveg);
			fprintf(pref,"\n");
		}
		fclose(pref);
	}

	else if (cfg_weasel_StartingNo.value()> 0){

		FILE* pref=fopen("PredSampleData.txt","a");
		if (pref == NULL) {
			m_TheLandscape->Warn("TPredator_Population_Manager::PredSampleFile","Could Not Open PredSampleData.txt File");
			exit(0);
		}

		Weasel* Wp;
		unsigned int SizeWeasel = (unsigned int) SupplyListSize(0);
		unsigned int TotalSize = (unsigned int) SizeWeasel;

		int w = m_TheLandscape->SupplySimAreaWidth();
		int h = m_TheLandscape->SupplySimAreaHeight();
		fprintf(pref,"%d\t %d\t %d\t %d\t %u\t %u\n", 0,w ,0, h, SizeWeasel, TotalSize);

		for (unsigned i = 0; i < SizeWeasel; i++){
			Wp=dynamic_cast<Weasel*>(SupplyAnimalPtr(0, i));
			int Wx{0}, Wy{0};
            SupplyLocXY(0, i, Wx, Wy);
			int Wpoly = m_TheLandscape->SupplyPolyRef(Wx, Wy);
			int Wele = m_TheLandscape->SupplyElementType(Wpoly);
			int Wveg = m_TheLandscape->SupplyVegType(Wpoly);

			int species = Wp->SupplySpeciesID();
			int kill = Wp->SupplyKill();
			int KillEff = Wp->SupplyKillEff();
			int Terr = Wp->SupplyTerr();
			int HomeR = Wp->SupplyHomeRange();

			fprintf(pref,"%d\t",year);
			fprintf(pref,"%d\t",month);
			fprintf(pref,"%d\t",dayinMo);
			fprintf(pref,"%d\t",today);

			fprintf(pref,"%d\t",species);
			fprintf(pref,"%d\t",kill);
			fprintf(pref,"%d\t",KillEff);
			fprintf(pref,"%d\t",Terr);
			fprintf(pref,"%d\t",HomeR);

			fprintf(pref,"%d\t",Wx);
			fprintf(pref,"%d\t",Wy);
			fprintf(pref,"%d\t",Wpoly);
			fprintf(pref,"%d\t",Wele);
			fprintf(pref,"%d\t",Wveg);
			fprintf(pref,"\n");
		}
		fclose(pref);
	}

	else if (cfg_owl_StartingNo.value()> 0) {

		FILE* pref=fopen("PredSampleData.txt","a");
		if (pref == NULL) {
			m_TheLandscape->Warn("TPredator_Population_Manager::PredSampleFile","Could Not Open PredSampleData.txt File");
			exit(0);
		}

		  Owl* Op;
		  unsigned int SizeOwl = (unsigned int) SupplyListSize(1);
		  unsigned int TotalSize = SizeOwl;
		  int w = m_TheLandscape->SupplySimAreaWidth();
		  int h = m_TheLandscape->SupplySimAreaHeight();
		  fprintf(pref,"%d\t %d\t %d\t %d\t %u\t %u\n", 0,w ,0, h, SizeOwl, TotalSize);

		  for (unsigned i = 0; i < SizeOwl; i++){
			  Op=dynamic_cast<Owl*>(SupplyAnimalPtr(1, i));
              int Wx{0}, Wy{0};
              SupplyLocXY(1, i, Wx, Wy);

			  int Wpoly = m_TheLandscape->SupplyPolyRef(Wx, Wy);
			  int Wele = m_TheLandscape->SupplyElementType(Wpoly);
			  int Wveg = m_TheLandscape->SupplyVegType(Wpoly);

			  int Wspecies = Op->SupplySpeciesID();
			  int Wkill = Op->SupplyKill();
			  int WKillEff = Op->SupplyKillEff();
			  int Terr = Op->SupplyTerr();
			  int HomeR = Op->SupplyHomeRange();

			  fprintf(pref,"%d\t",year);
			  fprintf(pref,"%d\t",month);
			  fprintf(pref,"%d\t",dayinMo);
			  fprintf(pref,"%d\t",today);

			  fprintf(pref,"%d\t",Wspecies);
			  fprintf(pref,"%d\t",Wkill);
			  fprintf(pref,"%d\t",WKillEff);
			  fprintf(pref,"%d\t",Terr);
			  fprintf(pref,"%d\t",HomeR);

			  fprintf(pref,"%d\t",Wx);
			  fprintf(pref,"%d\t",Wy);
			  fprintf(pref,"%d\t",Wpoly);
			  fprintf(pref,"%d\t",Wele);
			  fprintf(pref,"%d\t",Wveg);
			  fprintf(pref,"\n");
		  }
		  fclose(pref);
	}
}
//-------------------------------------------------------------------------------------------------------
void TPredator_Population_Manager::PredAutumnSample()
{
		int today = m_TheLandscape->SupplyDayInYear();
		int year = m_TheLandscape->SupplyYearNumber();

		if ((cfg_weasel_StartingNo.value()> 0)&& (cfg_owl_StartingNo.value()> 0)){
			FILE* pref=fopen("PredAutumnOutput.txt","a");
			if (pref == NULL) {
				m_TheLandscape->Warn("TPredator_Population_Manager::PredAutumnSample","Could Not Open PredAutumnOutput.txt File");
				exit(0);
			}

			unsigned int SizeWeasel = (unsigned int) SupplyListSize(0);
			unsigned int SizeOwl = (unsigned int) SupplyListSize(1);
			unsigned int TotalSize = SizeWeasel+SizeOwl;

			fprintf(pref,"%d\t", year);
			fprintf(pref,"%d\t", today);
			fprintf(pref,"%d\t", SizeWeasel);
			fprintf(pref,"%d\t", SizeOwl);
			fprintf(pref,"%d\t", TotalSize);
			fprintf(pref,"\n");

			fclose(pref);
		}

		  else if (cfg_weasel_StartingNo.value()> 0){
		  unsigned int SizeWeasel = (unsigned int) SupplyListSize(0);

		  FILE* pref = fopen("PredAutumnOutput.txt","a");
		  if (pref == NULL) {
			  m_TheLandscape->Warn("TPredator_Population_Manager::PredAutumnSample","Could Not Open PredAutumnOutput.txt File");
			  exit(0);
		  }
		  	fprintf(pref,"%d\t", year);
			fprintf(pref,"%d\t", today);
			fprintf(pref,"%d\t", SizeWeasel);
			fprintf(pref,"%d\t", 0);
			fprintf(pref,"%d\t", SizeWeasel);
			fprintf(pref,"\n");

			fclose(pref);
		  }

		  else if (cfg_owl_StartingNo.value()> 0) {
		  unsigned int SizeOwl = (unsigned int) SupplyListSize(1);

		  FILE* pref = fopen("PredAutumnOutput.txt","a");
		  if (pref == NULL) {
			  m_TheLandscape->Warn("TPredator_Population_Manager::PredAutumnSample","Could Not Open PredAutumnOutput.txt File");
			  exit(0);
		  }
		  fprintf(pref,"%d\t", year);
		  fprintf(pref,"%d\t", today);
		  fprintf(pref,"%d\t", 0);
		  fprintf(pref,"%d\t", SizeOwl);
		  fprintf(pref,"%d\t", SizeOwl);
		  fprintf(pref,"\n");

		  fclose(pref);
		  }
}
//-------------------------------------------------------------------------------------------------------
void TPredator_Population_Manager::PredSpringSample()
{
		int today = m_TheLandscape->SupplyDayInYear();
		int year = m_TheLandscape->SupplyYearNumber();

		if ((cfg_weasel_StartingNo.value())&& (cfg_owl_StartingNo.value()> 0)){
			unsigned int SizeWeasel = (unsigned int) SupplyListSize(0);
			unsigned int SizeOwl = (unsigned int) SupplyListSize(1);
			unsigned int TotalSize = SizeWeasel+SizeOwl;

			FILE* pref=fopen("PredSpringOutput.txt","a");
			if (pref == NULL) {
				m_TheLandscape->Warn("TPredator_Population_Manager::PredSpringSample","Could Not Open PredSpringOutput.txt File");
				exit(0);
			}
			fprintf(pref,"%d\t", year);
			fprintf(pref,"%d\t", today);
			fprintf(pref,"%d\t", SizeWeasel);
			fprintf(pref,"%d\t", SizeOwl);
			fprintf(pref,"%d\t", TotalSize);
			fprintf(pref,"\n");

			fclose(pref);
		}

		  else if (cfg_weasel_StartingNo.value()> 0){
		  unsigned int SizeWeasel = (unsigned int) SupplyListSize(0);

		  FILE* pref=fopen("PredSpringOutput.txt","a");
		  if (pref == NULL) {
			  m_TheLandscape->Warn("TPredator_Population_Manager::PredSpringSample","Could Not Open PredSpringOutput.txt File");
			  exit(0);
		  }
		  	fprintf(pref,"%d\t", year);
			fprintf(pref,"%d\t", today);
			fprintf(pref,"%d\t", SizeWeasel);
			fprintf(pref,"%d\t", 0);
			fprintf(pref,"%d\t", SizeWeasel);
			fprintf(pref,"\n");

			fclose(pref);
		  }

		  else if (cfg_owl_StartingNo.value()> 0) {
		  unsigned int SizeOwl = (unsigned int) SupplyListSize(1);

		  FILE* pref=fopen("PredSpringOutput.txt","a");
		  if (pref == NULL) {
			  m_TheLandscape->Warn("TPredator_Population_Manager::PredSpringSample","Could Not Open PredSpringOutput.txt File");
			  exit(0);
		  }
		  fprintf(pref,"%d\t", year);
		  fprintf(pref,"%d\t", today);
		  fprintf(pref,"%d\t", 0);
		  fprintf(pref,"%d\t", SizeOwl);
		  fprintf(pref,"%d\t", SizeOwl);
		  fprintf(pref,"\n");

		  fclose(pref);
		  }
}

//-------------------------------------------------------------------------------------------------------
void TPredator_Population_Manager::PredSpringAutumnSample()
{
		int today = m_TheLandscape->SupplyDayInYear();
		int year = m_TheLandscape->SupplyYearNumber();

		if ((cfg_weasel_StartingNo.value())&& (cfg_owl_StartingNo.value()> 0)){
			unsigned int SizeWeasel = (unsigned int) SupplyListSize(0);
			unsigned int SizeOwl = (unsigned int) SupplyListSize(1);
			unsigned int TotalSize = SizeWeasel+SizeOwl;

			FILE* pref=fopen("PredSpringAutumnOutput.txt","a");
			if (pref == NULL) {
			m_TheLandscape->Warn("TPredator_Population_Manager::PredSpringAutumnSample","Could Not Open PredSpringAutumnOutput.txt File");
			exit(0);
			}

			fprintf(pref,"%d\t", year);
			fprintf(pref,"%d\t", today);
			fprintf(pref,"%d\t", SizeWeasel);
			fprintf(pref,"%d\t", SizeOwl);
			fprintf(pref,"%d\t", TotalSize);
			fprintf(pref,"\n");

			fclose(pref);
		}

		  else if (cfg_weasel_StartingNo.value()> 0){
		  unsigned int SizeWeasel = (unsigned int) SupplyListSize(0);

		  FILE* pref=fopen("PredSpringAutumnOutput.txt","a");
		  if (pref == NULL) {
		  m_TheLandscape->Warn("TPredator_Population_Manager::PredSpringAutumnSample","Could Not Open PredSpringAutumnOutput.txt File");
		  exit(0);
		  }
		  	fprintf(pref,"%d\t", year);
			fprintf(pref,"%d\t", today);
			fprintf(pref,"%d\t", SizeWeasel);
			fprintf(pref,"%d\t", 0);
			fprintf(pref,"%d\t", SizeWeasel);
			fprintf(pref,"\n");

			fclose(pref);
		}
		  else if (cfg_owl_StartingNo.value()> 0) {
		  unsigned int SizeOwl = (unsigned int) SupplyListSize(1);

		  FILE* pref=fopen("PredSpringAutumnOutput.txt","a");
		  if (pref == NULL) {
		  m_TheLandscape->Warn("TPredator_Population_Manager::PredSpringAutumnSample","Could Not Open PredSpringAutumnOutput.txt File");
		  exit(0);
		  }
		  fprintf(pref,"%d\t", year);
		  fprintf(pref,"%d\t", today);
		  fprintf(pref,"%d\t", 0);
		  fprintf(pref,"%d\t", SizeOwl);
		  fprintf(pref,"%d\t", SizeOwl);
		  fprintf(pref,"\n");

		  fclose(pref);
		  }
}



//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------

void TPredator_Population_Manager::Run(int)
{
  DoFirst();
 // begin step actions ...
 // set all stepdone to false.... is this really necessary??
  int today = m_TheLandscape->SupplyDayInYear();
  int Y=m_TheLandscape->SupplyYearNumber();

  if (cfg_PredSampleDataUsed.value() ){
	  if (today==cfg_pred_first_sample_day.value() || today ==cfg_pred_second_sample_day.value())
	  {
	  PredSampleFile();
	  PredSpringAutumnSample();
	  }
	  if (today == cfg_pred_first_sample_day.value())
	  {
		  PredSpringSample();
	  }
	  if (today == cfg_pred_second_sample_day.value())
  	  {
		  PredAutumnSample();
	  }
  }

  if (today==364 ) {
	  FILE* pref = fopen("PredProbe.txt","a");
	  if (!pref) {
		  m_TheLandscape->Warn("Predator_Population_Manager Destructor","Could Not Open PredProbe.txt File");
		  exit(0);
	  }
	  if ((SupplyListSize(0)> 0)&& (SupplyListSize(1)> 0)){
	  unsigned int PredNo0 = (unsigned int) SupplyListSize(0);
	  unsigned int PredNo1 = (unsigned int) SupplyListSize(1);
	  fprintf(pref,"%d\t%u\t%u\t%u\t%u\n", Y, 0, PredNo0, 1, PredNo1);
	  }
	  else if (SupplyListSize(0)> 0){
	  unsigned int PredNo0 = (unsigned int) SupplyListSize(0);
	  fprintf(pref,"%d\t%u\t%u\n", Y, 0, PredNo0);
	  }
	  else if (SupplyListSize(1)> 0) {
	  unsigned int PredNo1 = (unsigned int) SupplyListSize(1);
	  fprintf(pref,"%d\t%u\t%u\n", Y, 1, PredNo1);
	  }
	  fclose(pref);
  }
// cal
 for (unsigned listindex=0; listindex<SupplyListIndexSize(); listindex++)
 {
   for (unsigned j=0; j< SupplyListSize(listindex); j++)
   {
       SupplyAnimalPtr(listindex, j)->SetStepDone(false);
   }
 }
 // call the begin-step-method of all objects
 for (unsigned listindex=0; listindex<SupplyListIndexSize();listindex++)
 {
   for (unsigned j=0; j< SupplyListSize(listindex); j++)
       SupplyAnimalPtr(listindex, j)->BeginStep();
 }
  DoBefore();
 // call the step-method of all objects
 do
 {
   for (unsigned listindex=0; listindex<SupplyListIndexSize();listindex++)
   {
     for (unsigned j=0; j< SupplyListSize(listindex); j++)
     {
         SupplyAnimalPtr(listindex, j)->Step();
     }
   } // for listindex
 } while (!StepFinished());
 DoAfter();
 // call the end-step-method of all objects
 for (unsigned listindex=0; listindex<SupplyListIndexSize();listindex++)
 {
   for (unsigned j=0; j< SupplyListSize(listindex); j++)
   {
       SupplyAnimalPtr(listindex, j)->EndStep();
   }
 }
 // ----------------
 // end of this step actions

 // For each animal list
  for (unsigned listindex=0; listindex<SupplyListIndexSize();listindex++)
 {
   // Must check each object in the list for m_CurrentStateNo==-1
   int TAend=(int) SupplyListSize(listindex)-1;
   for (int j=TAend; j>=0; j--)  // Search backwards is more efficicent
   {
     if (SupplyAnimalPtr(listindex, j)->GetCurrentStateNo() == -1) // code for kill it
     {
       delete SupplyAnimalPtr(listindex,j); // this line will become redundant if we use smart pointers instead of new
       RemoveFromList(listindex, j);
     }
   }
 }
 DoLast();
}
//---------------------------------------------------------------------------

// returns true if and only if all objects have finished the current step
bool TPredator_Population_Manager::StepFinished (void)
{
   for (unsigned listindex=0; listindex<SupplyListIndexSize();listindex++)
   {
     for (unsigned j=0; j< SupplyListSize(listindex); j++)
     {
       if (SupplyAnimalPtr(listindex,j)->GetStepDone() == false)
       {
         return false;
       }
     }
   }
 return true;
}
//---------------------------------------------------------------------------

bool TPredator_Population_Manager::InOtherTerritory(unsigned sp, int x, int y,
          TPredator* p_Pred)
{
  // Go through sp species and see if a territory at x,y will overlap with
  // theirs
  for (unsigned i=0; i< SupplyListSize(sp); i++)
  {
     TPredator* APredator=(TPredator *) SupplyAnimalPtr(sp,i);
     if (APredator->OverlapMyTerritory(x,y))
     {
        if (APredator!=p_Pred) return true;
     }
  }
  return false;
}
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
//       TPREDATOR CODE
//---------------------------------------------------------------------------

/**
Tpredator constructor
*/
TPredator::TPredator(Vole_Population_Manager* ThePrey, int p_x, int p_y,
     Landscape* p_L, TPredator_Population_Manager* p_PPM) : TAnimal(p_x,p_y)
{
    CurrentPState=tops_InitialState;
    m_KillEfficiency=0;
    m_kills_this_season=0;
    m_FailureCount=0;
    m_NoFailuresBeforeDispersal=5; // Default
    m_OurPopulationManager=p_PPM;
    m_Prey = ThePrey;
    m_Search_x=m_Location_x;
    m_Search_y=m_Location_y;
    m_HomeRange=0;
    m_HaveTerritory=false;
    m_SearchArea=0;
    m_DispersalMax=0;
    SimH=m_OurLandscape->SupplySimAreaHeight();
    SimW=m_OurLandscape->SupplySimAreaWidth();
    PreyResponse1=0;
    PreyResponse2=0;
}
//---------------------------------------------------------------------------

TPredator::~TPredator()
{

}

//---------------------------------------------------------------------------

int TPredator::st_Hunting()
{
   unsigned kills=0;
   // count days since last kill
   // Takes the Search_x, Search_y, SearchArea. Applies KillEfficiency to
   // all voles defined by this square
   CurrentPrey=m_Prey->SupplyVoleList(m_Search_x,m_Search_y,m_SearchArea);
   //int s=CurrentPrey->size();  // **CJT** to help with debug
   for (unsigned i=0; i<CurrentPrey->size(); i++)
   {
     if (g_random_fnc(1000)<m_KillEfficiency)
     {
      (*CurrentPrey)[i]->OnKilled();
      kills++;
     }
   }
   // Must tidy up here because m_Prey cannot know when to do it
   CurrentPrey->clear();
   delete CurrentPrey;
   // record the kills
   m_kills_this_season+=kills;
   return kills;
}


//---------------------------------------------------------------------------

void TPredator::st_Movement()
{
   // Can relocate search_x & search_y to be up somewhere in the homerange
   // but must have all the square in the HomeRange
   // 1. Drift a bit
   m_Location_x+=g_random_fnc(3)-1;
   m_Location_y+=g_random_fnc(3)-1;
   m_Location_x=(SimW+m_Location_x)%SimW;
   m_Location_y=(SimH+m_Location_y)%SimH;
   m_OurPopulationManager->UpdateGuardMap(m_Location_x,m_Location_y,m_guard_cell_x,m_guard_cell_y);
   // 2. determine search area
   int max_dist=m_HomeRange-m_SearchArea;
   m_Search_x=(m_Location_x+g_random_fnc(max_dist))%SimW;
   m_Search_y=(m_Location_y+g_random_fnc(max_dist))%SimH;
}

//---------------------------------------------------------------------------

bool TPredator::OverlapMyTerritory(unsigned x, unsigned y)
{
   // ensure we can't go negative
   x+=SimW;
   y+=SimH;
   unsigned mx = m_Location_x+SimW;
   unsigned my = m_Location_y+SimH;
   // most likely that it is not in so test for false
   if (x<mx-m_HomeRange) return false;
    else if (x>=mx+m_HomeRange) return false;
      else if (y<my-m_HomeRange) return false;
        else if (y>=my+m_HomeRange) return false;
         else return true;
}

//---------------------------------------------------------------------------

void TPredator::st_Dispersal()
{

  /** 
  * Moves the home range to an area where it does not overlap with a conspecific
  */
  bool found=false;
  unsigned Count=0;
  while ((!found)&&(Count<100))
  {
    // Simple random walk
    Count++;
    m_Location_x=((m_Location_x+(g_random_fnc(2*m_DispersalMax)-m_DispersalMax)))%SimW;
    m_Location_y=((m_Location_y+(g_random_fnc(2*m_DispersalMax)-m_DispersalMax)))%SimH;
    m_OurPopulationManager->UpdateGuardMap(m_Location_x,m_Location_y,m_guard_cell_x,m_guard_cell_y);

    if (!m_OurPopulationManager->InOtherTerritory(SpeciesID,m_Location_x,
                                                                 m_Location_y,this))
    {
      m_HaveTerritory=true;
      found=true;
    }
  }
}


//---------------------------------------------------------------------------
//       WEASEL CODE
//---------------------------------------------------------------------------


Weasel::Weasel(Vole_Population_Manager* ThePrey, int p_x, int p_y,
                           Landscape * p_L, TPredator_Population_Manager* p_PPM)
                                          : TPredator(ThePrey,p_x,p_y,p_L,p_PPM)
{
    SpeciesID=0;
    m_KillEfficiency=cfg_weasel_kill_efficiency.value();
    m_HomeRange=cfg_weasel_home_range.value();
    m_SearchArea=cfg_weasel_search_area.value();
    m_DispersalMax=cfg_weasel_DispersalMax.value();
    m_NoFailuresBeforeDispersal=cfg_weasel_NoFailuresBeforeDispersal.value();
    PreyResponse1=1;
    PreyResponse2=1;
	SpeciesID=0;
}
//---------------------------------------------------------------------------


Weasel::~Weasel()
{
    //Nothing to do
}
//---------------------------------------------------------------------------

void Weasel::BeginStep()
{
  int day= m_OurLandscape->SupplyDayInYear();
  if (day==weasel_breed_day)
  {
    int noToMake=m_kills_this_season/weasel_breed_threshold;
    for (int k=0; k<noToMake; k++)
    {
      // make a new weasel
      struct_Predator* sp;
      sp = new struct_Predator;
      sp->PM = m_OurPopulationManager;
      sp->L = m_OurLandscape;
      sp->x = m_Location_x;
      sp->y = m_Location_y;
      m_OurPopulationManager->CreateObjects(0,NULL,sp,1); // 0 = weasel
      delete sp;
    }
  }
  else if (day==364)
  {
    if (m_kills_this_season<weasel_death_threshold)
    {
      if (m_OurPopulationManager->supply_no_inds(WEASEL)>1)
      {
        m_OurPopulationManager->dec_inds(WEASEL);
        KillThis();
      }
    }
    m_kills_this_season=0; // reset the count
  }
}
//---------------------------------------------------------------------------

void Weasel::Step()
{
  if (m_StepDone || m_CurrentStateNo == -1) return;
  switch (CurrentPState)
  {
   case tops_InitialState: // Initial state
    CurrentPState=tops_Dispersal;
    m_HaveTerritory=false;
    break;
   case tops_Hunting:
    if (st_Hunting()<PreyResponse1) CurrentPState=tops_Movement;
    m_StepDone=true;
    break;
   case tops_Dispersal:
    st_Dispersal();
    if (m_HaveTerritory) CurrentPState=tops_Hunting;
    m_StepDone=true;
    break;
   case tops_Movement:
    st_Movement();
    if (st_Hunting()<PreyResponse2)   // alter this figure to increase functional response
     m_FailureCount++;
    else m_FailureCount=0;
    if (m_FailureCount>m_NoFailuresBeforeDispersal)
    {
     m_HaveTerritory=false;
     CurrentPState=tops_Dispersal;
    }
    else CurrentPState=tops_Hunting;
    m_StepDone=true;
    break;
   default:
	   m_OurLandscape->Warn("Weasel::Step()","unknown state - default");
    exit(1);
   }
}
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
//      OWL CODE
//---------------------------------------------------------------------------


Owl::Owl(Vole_Population_Manager* ThePrey, int p_x, int p_y,
                           Landscape * p_L, TPredator_Population_Manager* p_PPM)
                                          : TPredator(ThePrey,p_x,p_y,p_L,p_PPM)
{
    SpeciesID=1;
    m_KillEfficiency=250;  // 25%
    m_HomeRange=800;
    m_SearchArea=50;
    m_DispersalMax=1000;
    m_NoFailuresBeforeDispersal=200;
    PreyResponse1=1;
    PreyResponse2=1;
	SpeciesID=1;
}
//---------------------------------------------------------------------------

void Owl::BeginStep()
{
  int day= m_OurLandscape->SupplyDayInYear();
  if (day==owl_breed_day)
  {
    int noToMake=m_kills_this_season/owl_breed_threshold;
    for (int k=0; k<noToMake; k++)
    {
      // make a new owl
      struct_Predator* sp;
      sp = new struct_Predator;
      sp->PM = m_OurPopulationManager;
      sp->L = m_OurLandscape;
      sp->x = m_Location_x;
      sp->y = m_Location_y;
      m_OurPopulationManager->CreateObjects(OWL,NULL,sp,1); // 1 = owl
      delete sp;
    }
  }
  else if (day==364)
  {
    if (m_kills_this_season<owl_death_threshold)
    {
      if (m_OurPopulationManager->supply_no_inds(OWL)>1)
      {
        m_OurPopulationManager->dec_inds(OWL);
        KillThis();
      }
    }
    m_kills_this_season=0; // reset the count
  }
}
//---------------------------------------------------------------------------

void Owl::Step()
{
  if (m_StepDone || m_CurrentStateNo == -1) return;
  switch (CurrentPState)
  {
   case tops_InitialState: // Initial state
    CurrentPState=tops_Dispersal;
    m_HaveTerritory=false;
    break;
   case tops_Hunting:
    if (st_Hunting()<PreyResponse1) CurrentPState=tops_Movement;
    m_StepDone=true;
    break;
   case tops_Dispersal:
    st_Dispersal();
    CurrentPState=tops_Hunting;
    m_StepDone=true;
    break;
   case tops_Movement:
    st_Movement();
    if (st_Hunting()<PreyResponse2) // alter this figure to increase functional response
     m_FailureCount++;
    else m_FailureCount=0;
    if (m_FailureCount>m_NoFailuresBeforeDispersal)
    {
      m_HaveTerritory=false;
      CurrentPState=tops_Dispersal;
    }
    else CurrentPState=tops_Hunting;
    m_StepDone=true;
    break;
   default:
	   m_OurLandscape->Warn("Owl::Step()","unknown state - default");
    exit(1);
   }
}
//---------------------------------------------------------------------------

Owl::~Owl()
{
    //Nothing to do
}
//---------------------------------------------------------------------------


