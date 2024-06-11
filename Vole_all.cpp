/** \file
* \brief <B>This file contains the code for all vole lifestage classes</B> \n
*
* by Chris J. Topping \n
* Version of 28th Jan 2001 \n
* \n
* With additions as noted in: \n
* April 2006 \n
* November 2007 \n
* Doxygen formatted comments in May 2008 \n
*/

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

//---------------------------------------------------------------------------
//
#define _CRT_SECURE_NO_DEPRECATE

#include <iostream>
#include <fstream>
#include "../Landscape/ls.h"
#include "../BatchALMaSS/PopulationManager.h"
#include "../BatchALMaSS/BinaryMapBase.h"
#include "../BatchALMaSS/MovementMap.h"
#include "../Vole/vole_all.h"
#include "../Vole/VolePopulationManager.h"
#include "../Vole/GeneticMaterial.h"
#include "../Vole/Predators.h"
#include "../Vole/Vole_toletoc.h"

extern CfgBool cfg_RecordVoleMort;
extern CfgFloat l_pest_daily_mort;
extern string filename_AlleleInput;

#ifdef __VOLEPESTICIDEON
int Vole_Female::m_EndoCrineDisruptionGestationLength = 5; // Default value - will be altered by Vole_Female constructor
#endif

/** \brief Input parameter for daily extra mortality chance while dispersing */
CfgFloat cfg_extradispmort("VOLEDISPMORT", CFG_CUSTOM, 0.055);
/** \brief A parameter determining the density dependence effect. Zero will have no density dependence effect, 1.0 will have infinite impact */
CfgFloat cfg_VoleResourceRegrowth("VOLE_RESOURCEREGROWTH", CFG_CUSTOM, 0.5);

static CfgFloat cfg_volepcidebiodegredrate("VOLE_PCIDE_BIODEGREDATIONRATE", CFG_CUSTOM, 0.0); // 1.0 = no degradation, 0 = 100%
static CfgFloat cfg_InfanticideRangeRelToTerRange("VOLE_IINFANTICIDERANGE_RELTO_TERRANGE", CFG_CUSTOM, 0.5);
static CfgBool cfg_VoleMortalityDataUsed("VOLE_MORTALITY_DATA_USED", CFG_CUSTOM, false);
static CfgInt cfg_VoleEndoCrineDisruptionGestationLength("VOLE_ENDOCRINEDISPRUPTORGESTATION", CFG_CUSTOM, 21);
static CfgFloat cfg_InfanticideProbability("INFANTI_PROBA", CFG_CUSTOM, 0.01);
/** Scales to the average minimum territory quality acceptable */
static CfgFloat cfg_vole_habqualscaler("VOLE_HABQUALSCALER", CFG_CUSTOM, 2.1);
/** Minimum territory size for females */
CfgInt cfg_MinFemaleTerritorySize("VOLE_MINFEMALETERRITORYSIZE", CFG_CUSTOM, 8);
/** Maximum territory size for females */
static CfgInt cfg_MaxFemaleTerrSize("VOLE_MAXFEMALETERRITORYSIZE", CFG_CUSTOM, 8);
/** Minimum territory size for males */
CfgInt cfg_MinMaleTerritorySize("VOLE_MINMALETERRITORYSIZE", CFG_CUSTOM, 9);
/** Maximum territory size for males  */
static CfgInt cfg_MaxMaleTerrSize("VOLE_MAXMALETERRITORYSIZE", CFG_CUSTOM, 23);
/** Average physiological ifespan in months above 1 yr */
static CfgInt cfg_vole_LifeMonths("VOLE_LIFEMONTHS", CFG_CUSTOM, 3);

/** The max number of starvation days before death */
static CfgInt cfg_MaxStarvationDays("VOLE_MAXSTARVATIONDAYS", CFG_CUSTOM, 29);
extern double MoveToLessFavourable;
/** array of reciprocals to use instead of divding */
extern double g_speedy_Divides[2001];
/** Toxicological variable for specific pesticide effects */
extern CfgInt cfg_productapplicstartyear;
/** Toxicological variable for specific pesticide effects */
static CfgFloat cfg_PesticideAccumulationThreshold("VOLE_PESTICDEACCUMULATIONTHRESHOLD", CFG_CUSTOM, 20.0); //mg/kg
static CfgFloat cfg_PesticideAccumulationThresholdModelink2("VOLE_PESTICDEACCUMULATIONTHRESHOLD_MODELINKTWO", CFG_CUSTOM, 20.0); //mg/kg
static CfgFloat cfg_PesticideAccumulationThreshold2("VOLE_PESTICDEACCUMULATIONTHRESHOLDTWO", CFG_CUSTOM, 44.1); //mg/kg
static CfgFloat cfg_PesticideFemaleMaturityDelay("VOLE_PCIDE_FEMALEMATURITYDELAY", CFG_CUSTOM, 0.0192);
static CfgFloat cfg_PesticideLitterSizeReduction("VOLE_PCIDE_LITTERSIZEREDUCTION", CFG_CUSTOM, 0.0038);
static CfgFloat cfg_PesticideWeaningReduction("VOLE_PCIDE_WEANINGREDUCTION", CFG_CUSTOM, 0.0017);
CfgBool cfg_ResistanceDominant("VOLE_RESISTANCEDOMINANT", CFG_CUSTOM, false);

extern CfgInt cfg_productapplicendyear;
extern CfgBool cfg_ReallyBigOutputMonthly_used;
extern CfgBool cfg_CfgRipleysOutputUsed;
extern CfgBool cfg_AorOutput_used;
extern CfgBool l_pest_enable_pesticide_engine;
extern CfgBool cfg_pest_residue_or_rate; // true if residue, otherwise it is rate per m2
/** Minimum territory size for males from Erlinge et al 1990 */
CfgInt cfg_VoleDDepConst("VOLE_DDEPCONST", CFG_CUSTOM, 4);
CfgInt cfg_MinReproAgeM("VOLE_MINREPROAGEM", CFG_CUSTOM, 30);
CfgInt cfg_MinReproAgeF("VOLE_MINREPROAGEF", CFG_CUSTOM, 23);

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Constants
//---------------------------------------------------------------------------

/** 1st september female voles cannot mature after this */
const int FemNoMature = September;
/** Probability of maturing in each month of age - will do this if 5 months old unless it is winter */
//const double MaturationChance [5] = {0.020,0.100,0.200,0.500,1.000};
/** 21 days gestation */
const unsigned TheGestationPeriod = 21;
/** Age at weaning */
const unsigned WeanedAge = 13;
// Growth curve data for weight of voles
const int WeanedWeight = 5;
/** The maximum weight of adult males */
const double MaxWeightM = 60; // is 'double' to make MaleTerritoryRangeSlope a double
/** The maximum weight of adult females */
const double MaxWeightF = 55; // is 'double' to make FemaleTerritoryRangeSlope a double
/** How long it takes to get to max weight at max rate of growth */
const double DaysAtMaxGrowth = 90.0;
/** g per day during growing for males */
const double growthperdayM = (MaxWeightM - WeanedWeight) / DaysAtMaxGrowth;
/** g per day during growing for females */
const double growthperdayF = (MaxWeightF - WeanedWeight) / DaysAtMaxGrowth;
/** The day it is assumed that grass stops growing */
const int GrowStopDate = August;
/** The youngest a male can be reproductive at.*/
int MinReproAgeM;
/**The youngest a female can be reproductive at.*/
int MinReproAgeF;
/**The smallest a male can be reproductive at.*/
const unsigned MinReproWeightM = 40; //The smallest a male can be reproductive at.
/**The smallest a female can be reproductive at.*/
const int MinReproWeightF = 20; //The smallest a female can be reproductive at.
/** Age difference required before eviction by older male */
int g_sigAgeDiff = 30;
/** This is only here in case the resource requirement ought to change per month - it is not currently used. \n
Need 13 below because SupplyMonth() returns 1-12 */
const double FemaleResourceReq[13] =
{
	0,
	1.0,
	1.0,
	1.0,
	1.0,
	1.0,
	1.0,
	1.0,
	1.0,
	1.0,
	1.0,
	1.0,
	1.0
};

/** Max distance between fix distance in Mols radiotracking data was 48m */
const unsigned FemaleMovement = 50; // is normally 50
//const unsigned FemaleMovement = 100; // for walking DOUBLE
//const unsigned FemaleMovement = 200; // for walking FOUR TIMES
//const unsigned FemaleMovement = 400; // for walking EIGHT TIMES
/** arbitrary minimum movement of 5 - because we need one */
const unsigned MinFemaleMovement = 5; //is normally 5
//const unsigned MinFemaleMovement = 10; //for walking DOUBLE
//const unsigned MinFemaleMovement = 20; //for walking FOUR TIMES
//const unsigned MinFemaleMovement = 40; //for walking EIGHT TIMES
const unsigned VoleMoveInterval = 30; // will have max movement after 120 days, is normally 30
//const unsigned VoleMoveInterval = 60; // for walking DOUBLE
//const unsigned VoleMoveInterval = 120; // for walking FOUR TIMES
//const unsigned VoleMoveInterval = 240; // for walking EIGHT TIMES
/** Max distance between fix distance in Mols radiotracking data was 110m. For the male the distance that he can move increases with size */
const unsigned MaleMovement[4] = {10, 40, 70, 110}; // is normally 10,40,70,110
//const unsigned MaleMovement[4] = {20, 80, 140, 220}; // FOR WALKING DOUBLE
//const unsigned MaleMovement[4] = {40, 160, 280, 440}; // FOR WALKING FOUR TIMES
//const unsigned MaleMovement[4] = {80, 320, 560, 880}; // FOR WALKING EIGHT TIMES
/** arbitrary minimum movement of 5 - because we need one */
const unsigned MinMaleMovement = 5; //is normally 5
//const unsigned MinMaleMovement = 10; //for walking DOUBLE
//const unsigned MinMaleMovement = 20; //for walking FOUR TIMES
//const unsigned MinMaleMovement = 40; //for walking EIGHT TIMES
/** 1/10 deduced from Erlinge et al, 1983 */
int g_MaleReproductFinish = 0;
/** Data from M.arvalis (Heise & Lippke 1997) on infanticide attempts success rate with age of litter */
const double InfanticideChanceByAge[9] = {0.97, 0.86, 0.75, 0.64, 0.54, 0.43, 0.32, 0.21, 0.11};

// Mortality Data TODO - make these config variables
/** Mortality per day as a background mortality encompassing all things not directly modelled. This may include specialist predators if they are not simulated - if they are then reduce this value */
double g_DailyMortChance = 0.003;
double g_DailyMortChanceMaleTerr = 0.003;
/** Farm operation Mortality */
const double VoleHarvestMort = 0.20;
/** Farm operation Mortality */
const double VoleStriglingMort = 0.50;
/** Farm operation Mortality */
const double VoleSoilCultivationMort = 0.75;
/** Farm operation Mortality */
const double VolePigGrazingMort = 0.25;
/** Vole herbicide diect mortality */
const double VoleHerbicicideMort = 0.0;
/** Vole insecticide diect mortality */
const double VoleInsecticideMort = 0.0;

double g_extradispmort;
double g_NoFemalesMove = 0.01;
//---------------------------------------------------------------------------
//                          Vole_Base
//---------------------------------------------------------------------------

unsigned int Vole_Base::m_MaxMaleTerritorySize = 0;
unsigned int Vole_Base::m_MaxFemaleTerritorySize = 0;
unsigned int Vole_Base::m_MinMaleTerritorySize = 0;
unsigned int Vole_Base::m_MinFemaleTerritorySize = 0;
double Vole_Base::m_MinFVoleHabQual = 0;
double Vole_Base::m_MinJMVoleHabQual = 0;
double Vole_Base::m_MinMVoleHabQual = 0;
double Vole_Base::m_MaleTerritoryRangeSlope = 0;
double Vole_Base::m_FemaleTerritoryRangeSlope = 0;
bool Vole_Base::m_BreedingSeason = false;
#ifdef __VoleStarvationDays
	int Vole_Base::m_MaxStarvationDays = 21;
#endif

/**
\brief
Constructor for Vole_Base
*/
//Vole_Base::Vole_Base(struct_Vole_Adult* a_AVoleStruct_ptr, const string& filename_AlleleInput) 
  //  : TAnimal(a_AVoleStruct_ptr->x, a_AVoleStruct_ptr->y) { 
  //  Init(a_AVoleStruct_ptr, filename_AlleleInput); 
//}


//Vole_Base::Vole_Base(struct_Vole_Adult* a_AVoleStruct_ptr, const string& filename_AlleleInput) 
  //  : TAnimal(a_AVoleStruct_ptr->x, a_AVoleStruct_ptr->y), new_Genes(filename_AlleleInput), new_Mates_Genes("") { 
    //Init(a_AVoleStruct_ptr, filename_AlleleInput); 
//}

//Vole_Base::Vole_Base(struct_Vole_Adult* a_AVoleStruct_ptr, const string& filename_AlleleInput) 
 //   : TAnimal(a_AVoleStruct_ptr->x, a_AVoleStruct_ptr->y),new_Genes(filename_AlleleInput), new_Mates_Genes("") { 
   // Init(a_AVoleStruct_ptr, filename_AlleleInput); 
//}

//Vole_Base::Vole_Base(struct_Vole_Adult* a_AVoleStruct_ptr, const string& filename_AlleleInput) 
  //  : TAnimal(a_AVoleStruct_ptr->x, a_AVoleStruct_ptr->y) { 
    //Init(a_AVoleStruct_ptr, filename_AlleleInput); 
//}

Vole_Base::Vole_Base(struct_Vole_Adult* a_AVoleStruct_ptr, const string& filename_AlleleInput, int i) //THE MOST RECENTLY USED ONE!!
    : TAnimal(a_AVoleStruct_ptr->x, a_AVoleStruct_ptr->y), struct_Vole_Adult(*a_AVoleStruct_ptr) {
	//cout << "\n I am now in hte Vole_Base  constructor and will print positions: ";
	//cout << a_AVoleStruct_ptr->x << "," << a_AVoleStruct_ptr->y;
	//cout << "\n I am in the Vole base constructor,and we are gonna print hte genes!!!!";
	//a_AVoleStruct_ptr->new_Genes.PrintGenomes();
	//a_AVoleStruct_ptr->new_Mates_Genes.PrintGenomes();
	bool new_Genes_empty = (a_AVoleStruct_ptr->new_Genes.Genome0[0][0] == -1);
	if (new_Genes_empty){
		//cout << "\n new genes were empty AF!";
		a_AVoleStruct_ptr->new_Genes = GeneticMaterialNEW(filename_AlleleInput, new_Genes_empty);
    	a_AVoleStruct_ptr->new_Mates_Genes = GeneticMaterialNEW("", new_Genes_empty);
		a_AVoleStruct_ptr->new_Genes.MitochondrialLine=i;
		a_AVoleStruct_ptr->new_Genes.YChromoLine=i;
		//cout << "\n mito was " <<  a_AVoleStruct_ptr->new_Genes.MitochondrialLine;
		//cout << "\n ychrom was " <<  a_AVoleStruct_ptr->new_Genes.YChromoLine;
	}
    //Initialize new_Genes and new_Mates_Genes

    Init(a_AVoleStruct_ptr, filename_AlleleInput); 
}

//Vole_Base::Vole_Base(struct_Vole_Adult* a_AVoleStruct_ptr, const string& filename_AlleleInput, bool isOffspring) 
 //   : TAnimal(a_AVoleStruct_ptr->x, a_AVoleStruct_ptr->y), struct_Vole_Adult(*a_AVoleStruct_ptr) {
    // Check if a valid pointer to genes is provided
   // if (isOffspring=true) {
        // Initialize new_Genes and new_Mates_Genes using the pointer to genes
     //   new_Genes = a_AVoleStruct_ptr->new_Genes;
       // new_Mates_Genes = new_Mates_Genes = GeneticMaterialNEW("");
    //} else {
        // If the pointer is null, initialize the genes from the provided filename
      //  new_Genes = GeneticMaterialNEW(filename_AlleleInput);
        //new_Mates_Genes = GeneticMaterialNEW("");
    //}
    //Init(a_AVoleStruct_ptr, filename_AlleleInput); 
//}

/* Vole_Base::Vole_Base(struct_Vole_Adult* a_AVoleStruct_ptr, const string& filename_AlleleInput) 
    : TAnimal(a_AVoleStruct_ptr->x, a_AVoleStruct_ptr->y), struct_Vole_Adult(*a_AVoleStruct_ptr) {
    // Initialize new_Genes and new_Mates_Genes by copying genes from the pointer
    if (!a_AVoleStruct_ptr->new_Genes.Genome0.empty() && !a_AVoleStruct_ptr->new_Genes.Genome1.empty()) {
        new_Genes = a_AVoleStruct_ptr->new_Genes;
        new_Mates_Genes = GeneticMaterialNEW("");
    } else {
        // If the pointer is null or doesn't contain genetic information, initialize new_Genes and new_Mates_Genes using the provided filename
        new_Genes = GeneticMaterialNEW(filename_AlleleInput);
        new_Mates_Genes = GeneticMaterialNEW("");
    }
    Init(a_AVoleStruct_ptr, filename_AlleleInput); 
}
*/


//---------------------------------------------------------------------------

void Vole_Base::ReInit(struct_Vole_Adult* a_AVoleStruct_ptr, const string& filename_AlleleInput, int i) { Init(a_AVoleStruct_ptr, filename_AlleleInput); }
//---------------------------------------------------------------------------

void Vole_Base::Init(struct_Vole_Adult* a_AVoleStruct_ptr, const string& filename_AlleleInput, int i) {
	//cout << "\n I am now in volebase init! positions are ";
	//a_AVoleStruct_ptr->new_Genes.PrintGenomes();
	//a_AVoleStruct_ptr->new_Mates_Genes.PrintGenomes();
	//cout << a_AVoleStruct_ptr->x << "," << a_AVoleStruct_ptr->y;
	m_Location_x = a_AVoleStruct_ptr->x;
	m_Location_y = a_AVoleStruct_ptr->y;
	m_OurLandscape = a_AVoleStruct_ptr->L;
	m_CurrentStateNo = 0;
	m_OurPopulation = a_AVoleStruct_ptr->VPM;
	m_Mature = false;
	m_Age = 0;
	m_Weight = 0.5;
	m_NoOfYoungTotal = 0;
	m_DispVector = -1;
	m_Have_Territory = false;
	m_TerrRange = 0; // The size of the territory
	m_MinTerrRange = 0;
	m_Reserves = 3;
	m_LifeSpan = 365 + 60 + static_cast<int>(g_rand_uni_fnc() * (30 * cfg_vole_LifeMonths.value())); //
	CurrentVState = tovs_InitialState;
	m_intrappos.m_inAtrap = false;
	m_BornLastYear = false;

	m_MyGenes = a_AVoleStruct_ptr->Genes;
	//new_Genes = new_Genes(filename_AlleleInput);
	//new_Genes = a_AVoleStruct_ptr->new_Genes;
	//new_Mates_Genes = GeneticMaterialNEW("");
	new_Genes= a_AVoleStruct_ptr->new_Genes;
	new_Mates_Genes= a_AVoleStruct_ptr->new_Mates_Genes;
	//cout << "\n In Vole Base:: Init we just set hte genes to be:";
	//new_Genes.PrintGenomes();
	//new_Mates_Genes.PrintGenomes();

	
	m_SimH = m_OurPopulation->SupplySimH();
	m_SimW = m_OurPopulation->SupplySimW();
	IDNo = m_OurPopulation->IDNumber++;

	const unsigned FatherIdNo = a_AVoleStruct_ptr->FatherId;
	Set_FatherId(FatherIdNo);
	const unsigned MotherIdNo = a_AVoleStruct_ptr->MotherId;
	Set_MotherId(MotherIdNo);
	Set_XBorn(m_Location_x);
	Set_YBorn(m_Location_y);
	Set_PolyRefBorn(m_Location_x, m_Location_y);
	Set_VegBorn(m_Location_x, m_Location_y);
	Set_ElemBorn(m_Location_x, m_Location_y);
	Set_BirthYear(a_AVoleStruct_ptr->BirthYear);

	m_MinMaleTerritorySize = cfg_MinMaleTerritorySize.value();
	m_MaxMaleTerritorySize = cfg_MaxMaleTerrSize.value();
	if (m_MaxMaleTerritorySize < m_MinMaleTerritorySize) m_MaxMaleTerritorySize = m_MinMaleTerritorySize; // Corrects a bug caused by incorrect input settings
	m_MinFemaleTerritorySize = cfg_MinFemaleTerritorySize.value();
	m_MaxFemaleTerritorySize = cfg_MaxFemaleTerrSize.value();
	if (m_MaxFemaleTerritorySize < m_MinFemaleTerritorySize) m_MaxFemaleTerritorySize = m_MinFemaleTerritorySize; // Corrects a bug caused by incorrect input settings
	m_MinFVoleHabQual = cfg_vole_habqualscaler.value() * (4 * m_MinFemaleTerritorySize * m_MinFemaleTerritorySize);
	m_MinMVoleHabQual = cfg_vole_habqualscaler.value() * (4 * m_MinMaleTerritorySize * m_MinMaleTerritorySize);
	m_MinJMVoleHabQual = cfg_vole_habqualscaler.value() * (4 * m_MinFemaleTerritorySize * m_MinFemaleTerritorySize); // Yes this should be female
	m_MaleTerritoryRangeSlope = (m_MaxMaleTerritorySize - m_MinMaleTerritorySize) / (MaxWeightM - MinReproWeightM);
	m_FemaleTerritoryRangeSlope = (m_MaxFemaleTerritorySize - m_MinFemaleTerritorySize) / (MaxWeightF - MinReproWeightF);
#ifdef __VoleStarvationDays
		m_StarvationDays = 0;
		m_MaxStarvationDays = cfg_MaxStarvationDays.value();
#endif
#ifdef __VOLEPESTICIDEON
		ClearPesticidLoad();
		m_pesticideBioDegradeRate = cfg_volepcidebiodegredrate.value();
		SetPesticideInfluenced1(false);
		SetPesticideInfluenced2(false);
		m_pesticideInfluenced3 = 0.0;
#endif
	m_fertile = true;
}

//---------------------------------------------------------------------------

Vole_Base::~Vole_Base() {
}

//---------------------------------------------------------------------------
/**
\brief
Duplicates a vole
*/
/** Method used to duplicate a vole - most commonly used for experimental manipulation of populations e.g. return rate experiments */
void Vole_Base::CopyMyself(VoleObject a_vtype) {
	const auto av = new struct_Vole_Adult;
	av->VPM = m_OurPopulation;
	av->L = m_OurLandscape;
	av->m_flag = true; // Used to signal pesticide effect to CreateObjects
	av->PolyRefBorn = m_OurLandscape->SupplyPolyRef(m_Location_x, m_Location_y);
	av->ElemBorn = m_OurLandscape->SupplyElementType(av->PolyRefBorn);
	av->VegBorn = m_OurLandscape->SupplyVegType(m_Location_x, m_Location_y);
	av->x = (m_Location_x + 1) % m_OurPopulation->SupplySimW();
	av->y = m_Location_y % m_OurPopulation->SupplySimH();
	// Do the genetics
	av->Genes.Recombine(&m_MyGenes, &m_MyGenes);
	m_OurPopulation->CreateObjects(a_vtype, this, av, 1, filename_AlleleInput, true);
	// object will be destroyed by death state
	// but must let Dad know anyway
	delete av;
}

//-----------------------------------------------------------------------------
/**
\brief
All voles age at the end of the day
*/
void Vole_Base::EndStep() {
	m_Age++; // Once a day increment m_Age
}

//---------------------------------------------------------------------------
/**
\brief
All voles end here on death
*/
/**
Called when a vole dies. Just removes itself from the map and sets a flag to destroy the object in the endStep
*/
void Vole_Base::st_Dying() {
	FreeLocation();
	m_CurrentStateNo = -1;
}

//---------------------------------------------------------------------------

/**
\brief
Do a mortality test
*/
/**
Takes both physiological lifespan and background mortality into account to determine whether the vole should die - repeated calls increase the risk of dying
*/
bool Vole_Base::MortalityTest() {
	// returns true if the vole should die
	if (g_rand_uni_fnc() < g_DailyMortChance) { return true; }
	// returns true if the vole should die
	return false;
}

//---------------------------------------------------------------------------

/** Returns the sum of the qualities in the area covered by the territory (even if the vole does not have one) \n\n
* Uses the algorithm for fast searching of square arrays with wrap around co-ordinates.\n
* parameters = x,y starting coordinates top left\n
* range = extent of the space to search\n
* bottom left coordinate is therefore x+range, y+range \n
* 22/09/2000 \n
* \n
* For each polygon it gets the quality and multiplies by area of that polygon in the territory.
* This is summed for all polygons in the territoryto get overall quality. The value is then divided by
* the number of voles present less a threshold value.\n
*/
double Vole_Base::CalculateCarryingCapacity(int p_x, int p_y, int a_ddep) const {
	int NoPolygons = 0;
	int PolyRefData[500][2];
	// First convert centre co-rdinates to square co-ordinates
	double quality = 0;
	int x = p_x - m_MinTerrRange;
	if (x < 0) x += m_SimW;
	int y = p_y - m_MinTerrRange;
	if (y < 0) y += m_SimH;
	const int range_x = m_MinTerrRange + m_MinTerrRange;
	const int range_y = m_MinTerrRange + m_MinTerrRange;
	// Stage 1 make a list of polygons
	// create the extent variables
	const int xextent0 = x + range_x;
	const int yextent0 = y + range_y;
	const int xextent1 = x + range_x - m_SimW;
	const int yextent1 = y + range_y - m_SimH;
	// Create the looping variables needed
	int Dfinx = xextent0;
	int Dfiny = yextent0;
	int Afinx = 0; // unless the finx values for A-C are changed
	int Bfinx = 0; // the value of zero will stop the A-C loops from executing
	int Cfinx = 0;
	int Afiny = 0; // this one assigned 0 to stop compiler complaints
	// Now create the loop values;
	if (x + range_x <= m_SimW)
	{
		// Type B & D (overlap bottom, no overlap)
		if (yextent0 > m_SimH)
		{
			// Type B (overlap bottom only)
			Dfiny = m_SimH;
			Bfinx = x + range_x;
		}
	}
	else
	{
		// Type A & C overlap left edge or bottom & left
		if (yextent0 > m_SimH)
		{
			// Type C overlap bottom and left
			Afinx = xextent1;
			Afiny = m_SimH;
			Bfinx = m_SimW;
			Cfinx = xextent1;
			Dfinx = m_SimW;
			Dfiny = m_SimH;
		}
		else
		{
			// Type A overlap left edge
			Afinx = xextent1;
			Afiny = yextent0;
			Dfinx = m_SimW;
		}
	}
	// the default is:
	// Type D no overlap

	// A Loop
	for (int i = 0; i < Afinx; i++)
	{
		for (int j = y; j < Afiny; j++)
		{
			// Get the polyref
			const int PRef = m_OurLandscape->SupplyPolyRefIndex(i, j);
			// check if we have had this one already
			///*
			bool found = false;
			for (int k = 0; k < NoPolygons; k++)
			{
				if (PolyRefData[k][0] == PRef)
				{
					PolyRefData[k][1]++;
					found = true;
					break;
				}
			}
			if (!found)
			//*/
			{
				// don't have this one so get the height & type to the PolyDatas arrays
				PolyRefData[NoPolygons][0] = PRef;
				PolyRefData[NoPolygons][1] = 1;
				NoPolygons++;
			}
		}
	}
	// B Loop
	for (int i = x; i < Bfinx; i++)
	{
		for (int j = 0; j < yextent1; j++)
		{
			// Get the polyref
			const int PRef = m_OurLandscape->SupplyPolyRefIndex(i, j);
			// check if we have had this one already
			/*
				int index=-1;
				    for (int k=0; k<NoPolygons; k++)
			      {
			        if (PolyRefData[k][0]==PRef)
			        {
			          index=k;
			          break;
			        }
			      }
			      if (index!=-1)
			      {
			        PolyRefData[index][1]++;
			      }
				  else
			*/
			///*
			bool found = false;
			for (int k = 0; k < NoPolygons; k++)
			{
				if (PolyRefData[k][0] == PRef)
				{
					PolyRefData[k][1]++;
					found = true;
					break;
				}
			}
			if (!found)
			//*/
			{
				// don't have this one so get the height & type to the PolyDatas arrays
				PolyRefData[NoPolygons][0] = PRef;
				PolyRefData[NoPolygons][1] = 1;
				NoPolygons++;
			}
		}
	}
	// C Loop
	for (int i = 0; i < Cfinx; i++)
	{
		for (int j = 0; j < yextent1; j++)
		{
			// Get the polyref
			const int PRef = m_OurLandscape->SupplyPolyRefIndex(i, j);
			// check if we have had this one already
			/*
				int index=-1;
				    for (int k=0; k<NoPolygons; k++)
			      {
			        if (PolyRefData[k][0]==PRef)
			        {
			          index=k;
			          break;
			        }
			      }
			      if (index!=-1)
			      {
			        PolyRefData[index][1]++;
			      }
				  else
			*/
			///*
			bool found = false;
			for (int k = 0; k < NoPolygons; k++)
			{
				if (PolyRefData[k][0] == PRef)
				{
					PolyRefData[k][1]++;
					found = true;
					break;
				}
			}
			if (!found)
			//*/
			{
				// Don't have this one so get the height & type to the PolyDatas arrays
				PolyRefData[NoPolygons][0] = PRef;
				PolyRefData[NoPolygons][1] = 1;
				NoPolygons++;
			}
		}
	}
	// D Loop
	int last = -99999;
	int k_index = -1;
	for (int i = x; i < Dfinx; i++)
	{
		for (int j = y; j < Dfiny; j++)
		{
			// Get the polyref
			const int PRef = m_OurLandscape->SupplyPolyRefIndex(i, j);
			// check if we have had this one already
			if (last == PRef) { PolyRefData[k_index][1]++; }
			else
			{
				bool found = false;
				for (int k = 0; k < NoPolygons; k++)
				{
					if (PolyRefData[k][0] == PRef)
					{
						PolyRefData[k][1]++;
						found = true;
						last = PRef;
						k_index = k;
						break;
					}
				}
				if (!found)
				{
					// don't have this one so get the height & type to the PolyDatas arrays
					PolyRefData[NoPolygons][0] = PRef;
					PolyRefData[NoPolygons][1] = 1;
					NoPolygons++;
				}
			}
		}
	}
	// End of search algorithm
	/* for each polygon get the quality and multiply by amount of that polygon in the territory sum this to get overall quality */
	for (int i = 0; i < NoPolygons; i++) { quality += PolyRefData[i][1] * m_OurPopulation->GetHabitatQuality(PolyRefData[i][0]); }
	//    return quality;

	const int Voles = m_OurPopulation->SupplyHowManyVoles(m_Location_x, m_Location_y, m_TerrRange);
	if (Voles > a_ddep) quality /= Voles - a_ddep;
	// return the total quality
	return quality;
}

//---------------------------------------------------------------------------

/**
Do an assessment of the quality of the voles territory, and move.
\brief
*/
/**
This function does the same as the other CalcuateCarryingCapacity
     except that it alters the two values stand_x, stand_y to set them
     to be in the best polygon in the territory
\n
\n
Returns the mean of the qualities in the area covered by the territory (even if the vole does not have one) \n\n
 Uses the algorithm for fast searching of square arrays with wrap around co-ordinates.
\n   parameters = x,y starting coordinates top left
\n   range = extent of the space to search
\n   bottom left coordinate is therefore x+range, y+range
\n   22/09/2000
\n
\n For each polygon it gets the quality and multiplies by area of that polygon in the territory. This is summed for all polygons in the territory to get overall quality.
*/
//double Vole_Base::CalculateCarryingCapacity(int p_x, int p_y,
//                                              int & p_stand_x, int & p_stand_y)
//{
//  int NoPolygons = 0;
//  int PolyRefData[200][2];
//  // First convert centre co-rdinates to square co-ordinates
//  double quality = 0;
//  int x=p_x - m_MinTerrRange;
//  if (x<0) x+=m_SimW;
//  int y=p_y - m_MinTerrRange;
//  if (y<0) y+=m_SimH;
//  int range_x=m_MinTerrRange+m_MinTerrRange;
//  int range_y=m_MinTerrRange+m_MinTerrRange;
//  // Stage 1 make a list of polygons
//  // create the extent variables
//  int xextent0 = x+range_x;
//  int yextent0 = y+range_y;
//  int xextent1 = (x+range_x)-m_SimW;
//  int yextent1 = (y+range_y)-m_SimH;
//  // Create the looping variables needed
//  int Dfinx=xextent0;
//  int Dfiny=yextent0;
//  int Afinx=0;  // unless the finx values for A-C are changed
//  int Bfinx=0;  // the value of zero will stop the A-C loops from executing
//  int Cfinx=0;
//  int Afiny=0; // this one is only assigned as 0 to stop compiler complaints
//  // Now create the loop values;
//  if (x+range_x<=m_SimW)
//  {
//    // Type B & D (overlap bottom, no overlap)
//    if (yextent0>m_SimH)
//    {
//      // Type B (overlap bottom only)
//      Dfiny=m_SimH;
//      Bfinx=x+range_x;
//    }
//  }
//  else
//  {
//    // Type A & C overlap left edge or bottom & left
//    if (yextent0>m_SimH)
//    {
//      // Type C overlap bottom and left
//      Afinx=xextent1;
//      Afiny=m_OurPopulation->m_SimH;
//      Bfinx=m_OurPopulation->m_SimW;
//      Cfinx=xextent1;
//      Dfinx=m_OurPopulation->m_SimW;
//      Dfiny=m_OurPopulation->m_SimH;
//    }
//    else
//    {
//      // Type A overlap left edge
//      Afinx=xextent1;
//      Afiny=yextent0;
//      Dfinx=m_SimW;
//    }
//  }
//  // the default is:
//  // Type D no overlap
//
//
//  // A Loop
//  for (int i=0; i<Afinx; i++)
//  {
//    for (int j=y; j<Afiny; j++)
//    {
//      // Get the polyref
//      int PRef = m_OurLandscape->SupplyPolyRefIndex(i,j);
//      // check if we have had this one already
//      int index=-1;
//      for (int k=0; k<NoPolygons; k++)
//      {
//        if (PolyRefData[k][0]==PRef)
//        {
//          index=k;
//          break;
//        }
//      }
//      if (index!=-1)
//      {
//        PolyRefData[index][1]++;
//      }
//      else
//      {
//        // don't have this one so get the height & type to the PolyDatas arrays
//        PolyRefData[NoPolygons][0]=PRef;
//        PolyRefData[NoPolygons][1]=1;
//        NoPolygons++;
//      }
//    }
//  }
//  // B Loop
//  for (int i=x; i<Bfinx; i++)
//  {
//    for (int j=0; j<yextent1; j++)
//    {
//      // Get the polyref
//      int PRef = m_OurLandscape->SupplyPolyRefIndex(i,j);
//      // check if we have had this one already
//      int index=-1;
//      for (int k=0; k<NoPolygons; k++)
//      {
//        if (PolyRefData[k][0]==PRef)
//        {
//          index=k;
//          break;
//        }
//      }
//      if (index!=-1)
//      {
//        PolyRefData[index][1]++;
//      }
//      else
//      {
//        // don't have this one so get the height & type to the PolyDatas arrays
//        PolyRefData[NoPolygons][0]=PRef;
//        PolyRefData[NoPolygons][1]=1;
//        NoPolygons++;
//      }
//    }
//  }
//  // C Loop
//  for (int i=0; i<Cfinx; i++)
//  {
//    for (int j=0; j<yextent1; j++)
//    {
//      // Get the polyref
//      int PRef = m_OurLandscape->SupplyPolyRefIndex(i,j);
//      // check if we have had this one already
//      int index=-1;
//      for (int k=0; k<NoPolygons; k++)
//      {
//        if (PolyRefData[k][0]==PRef)
//        {
//          index=k;
//          break;
//        }
//      }
//      if (index!=-1)
//      {
//        PolyRefData[index][1]++;
//      }
//      else
//      {
//        // Don't have this one so get the height & type to the PolyDatas arrays
//        PolyRefData[NoPolygons][0]=PRef;
//        PolyRefData[NoPolygons][1]=1;
//        NoPolygons++;
//      }
//    }
//  }
//  // D Loop
//  for (int i=x; i<Dfinx; i++)
//  {
//    for (int j=y; j<Dfiny; j++)
//    {
//      // Get the polyref
//      int PRef = m_OurLandscape->SupplyPolyRefIndex(i,j);
//      // check if we have had this one already
//      int index=-1;
//      for (int k=0; k<NoPolygons; k++)
//      {
//        if (PolyRefData[k][0]==PRef)
//        {
//          index=k;
//          break;
//        }
//      }
//      if (index!=-1)
//      {
//        PolyRefData[index][1]++;
//      }
//      else
//      {
//        // don't have this one so get the height & type to the PolyDatas arrays
//        PolyRefData[NoPolygons][0]=PRef;
//        PolyRefData[NoPolygons][1]=1;
//        NoPolygons++;
//      }
//    }
//  }
//  // End of search algorithm
//    // for each polygon get the quality and multiply by amount of
//    // that polygon in the territory
//    // sum this to get overall quality
//    double bestqual=-2;
//    unsigned found=0;
//    for ( int i=0; i<NoPolygons; i++)
//    {
//		double qual = m_OurPopulation->GetHabitatQuality(PolyRefData[i][0]);
//       if (bestqual<qual)
//       {
//         bestqual=qual;
//         found = i;
//       }
//       quality+=PolyRefData[i][1]*qual;
//    }
//    // now find a location in polygon i
//  // D Loop
//  for (int i=x; i<Dfinx; i++)
//  {
//    for (int j=y; j<Dfiny; j++)
//    {
//      // Get the polyref
//       if (m_OurLandscape->SupplyPolyRefIndex(i,j)==PolyRefData[found][0])
//      {
//        p_stand_x = p_x;
//        p_stand_y = p_y;
//        return quality;
//      }
//    }
//  }
//  // A Loop
//  for (int i=0; i<Afinx; i++)
//  {
//    for (int j=y; j<Afiny; j++)
//    {
//      // Get the polyref
//      if (m_OurLandscape->SupplyPolyRefIndex(i,j)==PolyRefData[found][0])
//      {
//        p_stand_x = p_x;
//        p_stand_y = p_y;
//        return quality;
//      }
//    }
//  }
//  // B Loop
//  for (int i=x; i<Bfinx; i++)
//  {
//    for (int j=0; j<yextent1; j++)
//    {
//      // Get the polyref
//      if (m_OurLandscape->SupplyPolyRefIndex(i,j)==PolyRefData[found][0])
//      {
//        p_stand_x = p_x;
//        p_stand_y = p_y;
//        return quality;
//      }
//    }
//  }
//  // C Loop
//  for (int i=0; i<Cfinx; i++)
//  {
//    for (int j=0; j<yextent1; j++)
//    {
//      if (m_OurLandscape->SupplyPolyRefIndex(i,j)==PolyRefData[found][0])
//      {
//        p_stand_x = p_x;
//        p_stand_y = p_y;
//        return quality;
//      }
//    }
//  }
//  m_OurLandscape->Warn("Vole_Base::CalculateCarryingCapacity - error",NULL);
//  exit(1);
//}
//---------------------------------------------------------------------------

/**
\brief
Movement
*/
/**
This will alter m_Location_x & m_Location_y. \n
It will give a rather directed movement towards p_Vector. \n
Generally the vole will stay in the best habitat, but occasional mis-steps occur.\n
p_Vector gives the preferred direction (0-7), p_Distance is the number of steps.
*/
void Vole_Base::MoveTo(int p_Vector, int p_Distance, int p_iterations) {
	//cout << "\n \n \n I am in MoveTo now!";
	//cout << "\n current position is " <<  m_Location_x << " , " << m_Location_y << "\t";
	bool RightBefore;
	bool RightAfter;
	if (m_Location_x >=5050){
		RightBefore=true;
	}else{
		RightBefore=false;
	}
	int m_Location_x_orig=m_Location_x;
	const int offset = p_Distance * p_iterations;
	if (m_Location_x - offset < 0 || m_Location_x + offset >= m_SimW || m_Location_y - offset < 0 || m_Location_y + offset >= m_SimH)
	{
		//cout<< "\n ALERT A vole wants to leave!!";
		// Need correct coords
		// Make sure that the coords can't become -ve
		int vx = m_Location_x + m_SimW;
		int vy = m_Location_y + m_SimH;
		do
		{
			DoWalkingCorrect(p_Distance, p_Vector, vx, vy);
			p_Distance = 1;
		}
		while (GetLocation(vx % m_SimW, vy % m_SimH) && p_iterations-- > 0);
		// alter the voles location (& correct coords)
		FreeLocation();
		m_Location_x = vx % m_SimW;
		m_Location_y = vy % m_SimH;
		//if (vx >= 2*m_SimW ||vx <= m_SimW  ){
		//	cout << "\n ALERT B: the non normalized width is: " << vx;
		//}
		//if (vy >= 2*m_SimH ||vy <= m_SimH  ){
		//	cout << "\n ALERT C : the non normalized height is: " << vy;
		//}
		if (m_Location_x >=5050){
			RightAfter=true;
		}else{
			RightAfter=false;
		}
		//if (RightBefore!=RightAfter){
		//	cout << "ALERT D, something moved across the middle!!!";
		//	cout << "\n Before: " << m_Location_x_orig << ", now: " << m_Location_x;

		//}
		//cout << "\n new position (after trying to leave) is " <<  m_Location_x << " , " << m_Location_y << "\t";
		m_OurPopulation->UpdateGuardMap(m_Location_x, m_Location_y, m_guard_cell_x, m_guard_cell_y);
		SetLocation();
	}
	else
	{
		// Don't need correct coords
		int vx = m_Location_x;
		int vy = m_Location_y;
		do
		{
			DoWalking(p_Distance, p_Vector, vx, vy);
			p_Distance = 1;
		}
		while (GetLocation(vx, vy) && p_iterations-- > 0);
		// alter the voles location (& correct coords)
		FreeLocation();
		m_Location_x = vx;
		m_Location_y = vy;
		if (m_Location_x >=5050){
			RightAfter=true;
		}else{
			RightAfter=false;
		}
		//if (RightBefore!=RightAfter){
			//cout << "\n ALERT D, something moved across the middle!!!";
			//cout << "\n Before: " << m_Location_x_orig << ", now: " << m_Location_x;
		//}
		//cout << "\n new position (after not trying to leave) is " <<  m_Location_x << " , " << m_Location_y << "\t";
		m_OurPopulation->UpdateGuardMap(m_Location_x, m_Location_y, m_guard_cell_x, m_guard_cell_y);
		SetLocation();
	}
}

//---------------------------------------------------------------------------

/**
\brief
Walking
*/
/**
This method does the actual stepping - there is no look ahead here, so steps are taken one at a time based on the habitat type and vector given.
*/
void Vole_Base::DoWalking(int p_Distance, int& p_Vector, int& vx, int& vy) const {
	int t[5], q[5];
	///*
	t[0] = p_Vector;
	t[1] = p_Vector + 1 & 0x07;
	t[2] = p_Vector + 7 & 0x07;
	t[3] = p_Vector + 2 & 0x07;
	t[4] = p_Vector + 6 & 0x07;
	//*/
	for (int i = 0; i < p_Distance; i++)
	{
		// test the squares at Vector, Vector+1+2, Vector-1-2
		// They have either a quality or are inaccessible (water,buildings)
		// if can go to one of these squares then pick the best one
		// if all or some are equal then take a random pick.
		// if all are 'bad' then add or subtract one from vector and try again

		/*	t[0] = p_Vector;
		    t[1] = (p_Vector+1) & 0x07;
		    t[2] = (p_Vector+7) & 0x07;
		    t[3] = (p_Vector+2) & 0x07;
		    t[4] = (p_Vector+6) & 0x07; */

		q[0] = MoveQuality(vx + g_vector_x[t[0]], vy + g_vector_y[t[0]]);
		q[1] = MoveQuality(vx + g_vector_x[t[1]], vy + g_vector_y[t[1]]);
		q[2] = MoveQuality(vx + g_vector_x[t[2]], vy + g_vector_y[t[2]]);
		q[3] = MoveQuality(vx + g_vector_x[t[3]], vy + g_vector_y[t[3]]);
		q[4] = MoveQuality(vx + g_vector_x[t[4]], vy + g_vector_y[t[4]]);
		if (g_rand_uni_fnc() < MoveToLessFavourable)
		{
			// allow a mistake once in a while
			for (int j = 1; j < 5; j++) q[j] = -1;
		}
		// Now pick the best of these
		int best = q[0];
		int score = 1;
		for (int ii = 1; ii < 5; ii++)
		{
			if (q[ii] > best)
			{
				best = q[ii];
				score = 1;
			}
			else if (q[ii] == best) score++;
		}
		if (best == -1)
		{
			// can't go anywhere so change the vector
			if (g_random_fnc(2)) ++p_Vector;
			else --p_Vector;
			p_Vector &= 0x07;
			///*
			t[0] = p_Vector;
			t[1] = p_Vector + 1 & 0x07;
			t[2] = p_Vector + 7 & 0x07;
			t[3] = p_Vector + 2 & 0x07;
			t[4] = p_Vector + 6 & 0x07;
			//*/
		}
		else
		{
			// Can go to one of score squares
			const int scored = g_random_fnc(score); // pick one
			int loop = 0;
			for (int ii = 0; ii < 5; ii++)
			{
				if (best == q[ii]) loop++; // count the squares with 'best' quality
				if (loop > scored)
				{
					loop = ii; // go to i-square
					break;
				}
			}
			// change co-ordinates
			vx += g_vector_x[t[loop]];
			vy += g_vector_y[t[loop]];
		}
	}
}

//------------------------------------------------------------------------------
/**
\brief
Walking where there is a danger of stepping off the world.
*/
/**
This method does the actual stepping - there is no look ahead here, so steps are taken one at a time based on the habitat type and vector given.
\n This version corrects coords for wrap around. This is slower so is only called when necessary.
*/
void Vole_Base::DoWalkingCorrect(int p_Distance, int& p_Vector, int& vx, int& vy) const {
	//cout << "\n ALERT I am now in DoWalkingCorrect!!";
	int t[5], q[5];
	///*
	t[0] = p_Vector;
	t[1] = p_Vector + 1 & 0x07;
	t[2] = p_Vector + 7 & 0x07;
	t[3] = p_Vector + 2 & 0x07;
	t[4] = p_Vector + 6 & 0x07;
	//*/
	for (int i = 0; i < p_Distance; i++)
	{
		// test the squares at Vector, Vector+1+2, Vector-1-2
		// They have either a quality or are inaccessible (water,buildings)
		// if can go to one of these squares then pick the best one
		// if all or some are equal then take a random pick.
		// if all are 'bad' then add or subtract one from vector and try again
		/* t[0] = p_Vector;
	     t[1] = (p_Vector+1) & 0x07;
	     t[2] = (p_Vector+7) & 0x07;
	     t[3] = (p_Vector+2) & 0x07;
	     t[4] = (p_Vector+6) & 0x07; */
		//if (vx + g_vector_x[t[0]])> 11000 || vx + g_vector_x[t[0]]<0){

		//}
		
		q[0] = MoveQuality((vx + g_vector_x[t[0]]) % m_SimW, (vy + g_vector_y[t[0]]) % m_SimH);
		q[1] = MoveQuality((vx + g_vector_x[t[1]]) % m_SimW, (vy + g_vector_y[t[1]]) % m_SimH);
		q[2] = MoveQuality((vx + g_vector_x[t[2]]) % m_SimW, (vy + g_vector_y[t[2]]) % m_SimH);
		q[3] = MoveQuality((vx + g_vector_x[t[3]]) % m_SimW, (vy + g_vector_y[t[3]]) % m_SimH);
		q[4] = MoveQuality((vx + g_vector_x[t[4]]) % m_SimW, (vy + g_vector_y[t[4]]) % m_SimH);
		if (g_rand_uni_fnc() < MoveToLessFavourable)
		{
			// allow a mistake once in a while
			for (int j = 1; j < 5; j++) q[j] = -1;
		}
		//cout << "\n The found MoveQualities were: ";
		// Now pick the best of these
		int best = q[0];
		int score = 1;
		for (int ii = 1; ii < 5; ii++)
		{
			//cout << q[ii] << ", ";
			if (q[ii] > best)
			{
				best = q[ii];
				score = 1;
			}
			else
				if (q[ii] == best) score++;
		}
		if (best == -1)
		{
			//cout << "\n The best movement quality was -1.";
			// can't go anywhere so change the vector
			if (g_random_fnc(2)) ++p_Vector;
			else --p_Vector;
			p_Vector &= 0x07;
			// /*
			t[0] = p_Vector;
			t[1] = p_Vector + 1 & 0x07;
			t[2] = p_Vector + 7 & 0x07;
			t[3] = p_Vector + 2 & 0x07;
			t[4] = p_Vector + 6 & 0x07; //*/
		}
		else
		{
			// Can go to one of score squares
			const int scored = g_random_fnc(score); // pick one
			int loop = 0;
			for (int ii = 0; ii < 5; ii++)
			{
				if (best == q[ii]) loop++; // count the squares with 'best' quality
				if (loop > scored)
				{
					loop = ii; // go to i-square
					break;
				}
			}
			// change co-ordinates
			vx += g_vector_x[t[loop]];
			vy += g_vector_y[t[loop]];
		}
	}
}

//------------------------------------------------------------------------------

/**
\brief
Dispersal - directed movement.
*/
/**
This works like MoveTo above, but with a rather more directed movement aimed at moving further from the start point.
*/
void Vole_Base::Escape(int p_Vector, int p_Distance) {
	// This will alter m_Location_x & m_Location_y
	// it will give a very directed movement towards p_Vector

	// p_Vector gives the preferred direction (0-7)
	// p_Distance is the number of steps
	int vx = m_Location_x + m_SimW;
	int vy = m_Location_y + m_SimH;
	//    int counter = 0;
	for (int i = 0; i < p_Distance; i++)
	{
		// test the squares at Vector, Vector+1+2, Vector-1-2
		// They have either a quality or are inaccessible (water,buildings)
		// if can go to one of these squares then pick the best one
		// if all or some are equal then take a random pick.
		// if all are 'bad' then add or subtract one from p_Vector and try again
		//      counter ++;
		int t[3], q[3];
		t[0] = p_Vector;
		t[1] = p_Vector + 1 & 0x07;
		t[2] = p_Vector + 7 & 0x07;
		q[0] = MoveQuality((vx + g_vector_x[t[0]]) % m_SimW, (vy + g_vector_y[t[0]]) % m_SimH);
		q[1] = MoveQuality((vx + g_vector_x[t[1]]) % m_SimW, (vy + g_vector_y[t[1]]) % m_SimH);
		q[2] = MoveQuality((vx + g_vector_x[t[2]]) % m_SimW, (vy + g_vector_y[t[2]]) % m_SimH);
		// Now pick the best of these
		int noscore = 0;
		for (int ii = 0; ii < 3; ii++) { if (q[ii] == -1) { noscore++; } }
		if (noscore == 3)
		{
			// can't go anywhere so change the vector
			if (g_random_fnc(2)) p_Vector++;
			else p_Vector--;
		}
		else
		{
			// Can go to at least one of score squares
			// try the middle first
			int loop = 0;
			if (q[0] == -1)
			{
				// otherwise randomly try one side or other
				// it must be possible to go to one of these
				const int which = g_random_fnc(2);
				if (which == 1)
				{
					if (q[1] != -1) loop = 1;
					else loop = 2;
				}
				else
				{
					if (q[2] != -1) loop = 2;
					else loop = 1;
				}
			}
			// change co-ordinates
			vx += g_vector_x[t[loop]];
			vy += g_vector_y[t[loop]];
		}
	}
	// alter the voles location
	FreeLocation();
	m_Location_x = vx % m_SimW;
	m_Location_y = vy % m_SimH;
	m_OurPopulation->UpdateGuardMap(m_Location_x, m_Location_y, m_guard_cell_x, m_guard_cell_y);
	SetLocation();
}

//---------------------------------------------------------------------------

/**
\brief
Test a location for quality while moving.
*/
/**
Returns the quality of a patch of habitat at p_x,p_y. \n
Can't walk through another vole though - so test location first.
*/
int Vole_Base::MoveQuality(int p_x, int p_y) const {
	if (m_OurPopulation->m_VoleMap->GetMapValue(p_x, p_y) != nullptr)
		return -1;
	// Nobody there so carry on

	return vole_tole_move_quality(m_OurLandscape, p_x, p_y);
}

//---------------------------------------------------------------------------
//                           Vole_JuvenileFemale
//---------------------------------------------------------------------------

/**
\brief
Vole_JuvenileFemale constructor
*/
Vole_JuvenileFemale::Vole_JuvenileFemale(struct_Vole_Adult* p_aVoleStruct,const string& filename_AlleleInput, int i): Vole_Base(p_aVoleStruct,filename_AlleleInput, i) { 
	//cout << "\n I am now in hte Vole_JuvenileFemale constructor";
	//cout << "\n these are the genes after going to Juvenile Female constructor";
	//p_aVoleStruct->new_Genes.PrintGenomes();
	//p_aVoleStruct->new_Mates_Genes.PrintGenomes();
	m_Sex = false;
	m_TerrRange = m_MinFemaleTerritorySize;
	m_MinTerrRange = m_MinFemaleTerritorySize;
	m_BornLastYear = false;
	m_Weight = WeanedWeight;
	m_Age = 14;
#ifdef __VOLEPESTICIDEON
  m_maturitydelay = 0.0;
#endif
}

//---------------------------------------------------------------------------

void Vole_JuvenileFemale::ReInit(struct_Vole_Adult* p_aVoleStruct, const string& filename_AlleleInput, int i) {
	Init(p_aVoleStruct, filename_AlleleInput, i);
	m_Sex = false;
	m_TerrRange = m_MinFemaleTerritorySize;
	m_MinTerrRange = m_MinFemaleTerritorySize;
	m_BornLastYear = false;
	m_Weight = WeanedWeight;
	m_Age = 14;
#ifdef __VOLEPESTICIDEON
	m_maturitydelay = 0.0;
#endif
}

//---------------------------------------------------------------------------

Vole_JuvenileFemale::~Vole_JuvenileFemale() {
	// Nothing to do
}

//---------------------------------------------------------------------------

/**
\brief
Female vole BeginStep
*/
/**
The BeginStep is one of the three timestep divisions. This is called once for each vole before Step and EndStep. \n
The main function here is to remove voles that die before they take up CPU resources in the Step code.
\n Can also be used to check for pesticide accumulation levels in pesticide simulation
*/
void Vole_JuvenileFemale::BeginStep() {
	if (MortalityTest())
	{
		CurrentVState = tovs_FDying;
		m_StepDone = true;
		if (cfg_RecordVoleMort.value()) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_FBck, 1);
	}
	else
	{
		if (--m_LifeSpan < 1)
		{
			if (cfg_RecordVoleMort.value()) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_FLife, 1);
			CurrentVState = tovs_FDying;
			m_StepDone = true;
		}
	}
#ifdef __VOLEPESTICIDEON
	PesticideIngestion();
	ActOnPesticideDose();
#endif
}

//---------------------------------------------------------------------------

/**
\brief
JuvenileFemale vole Step
*/
/**
The Step is one of the three timestep divisions. This is called repeatedly after BeginStep and before EndStep, until all voles report that they are done with Step. \n
\n Most of the behaviours are controlled by moving voles between behavioural states in Step (for other models this is also done in BeginStep and EndStep).
\n When a vole is done for the day it will signal this by setting m_StepDone==true. NB that a call to one behaviour may trigger a call to another behaviour on the next call to step inside the same timestep. In this way a daily cycle of activity can be undertaken (i.e. do reproduction and explore)
*/
void Vole_JuvenileFemale::Step() {
	if (m_StepDone) return;
	switch (CurrentVState)
	{
	case 0: // Initial state
		CurrentVState = tovs_FEvaluateExplore;
		break;
	case tovs_FEvaluateExplore: // Evaluate & Explore
		switch (st_Evaluate_n_Explore())
		{
		case 1: // Dead from dispersal mortality
			CurrentVState = tovs_FDying; // Die
			break;
		case 2: // Starved
			if (cfg_RecordVoleMort.value()) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_FStarve, 1);
			CurrentVState = tovs_FDying; // Die
			break;
		default:
			; // Do nothing
		}
#ifdef __VOLEPESTICIDEON
		if (m_Age >= MinReproAgeF+ m_maturitydelay ) st_BecomeSubAdult();
#else
		if (m_Age == MinReproAgeF) st_BecomeSubAdult();
#endif
		m_StepDone = true;
		break;
	case tovs_FDying: // Die
		FreeLocation();
		m_StepDone = true;
		break;
	default:
		m_OurLandscape->Warn("Vole_Female::Step - unknown state", nullptr);
		exit(9);
	}
}

//---------------------------------------------------------------------------

/**
\brief
Female vole EndStep
*/
/**
The EndStep one of the three timestep divisions. This is called once for each vole after BeginStep and Step. \n
The main function here is to remove voles that have died during step and otherwise to grow if not at max weight.
It also checks if the vole was killed due to human management and determines the potential territory size.
*/
void Vole_JuvenileFemale::EndStep() {
#ifdef __VoleStarvationDays
  if ( m_StarvationDays > m_MaxStarvationDays)
  {
	  st_Dying();
	  if (cfg_RecordVoleMort.value() ) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_FStarve,1);
  }
#endif
	CheckManagement();
	if (CurrentVState == tovs_FDying) { st_Dying(); }

	else
	{
		m_Age++;
		/*              FEMALE GROWTH NOTES
	
		 Female vole grows until 20g.
		 After that she will only grow if he has matured
		 Growth continues up to 55g.
	
		 Reproduction cannot occur below 20g or 20 days
	
		 Growth only occurs between 1 March and 1st August
	
		 taken from Hanson L, 1977, Oikos 29.
		*/
		const int today = m_OurLandscape->SupplyDayInYear();
		if (today < GrowStopDate && today >
			m_OurPopulation->SupplyGrowthStartDate())
		{
			if (m_Weight < 20) { m_Weight += growthperdayF; }
			else if (m_Mature == true && m_Weight < MaxWeightF)
			{
				m_Weight += growthperdayF;
				m_TerrRange = static_cast<int>(m_MinFemaleTerritorySize + m_FemaleTerritoryRangeSlope * (m_Weight - MinReproWeightF));
			}
		}
		else if (today == 1) m_BornLastYear = true; // must be true if alive on 1st Jan.
	}
}

//---------------------------------------------------------------------------

/**
\brief
External event handler.
*/
/**
This method evaluates external events and chooses a suitable response (in this case a probability of dying because other effects will be taken up by the evaluate and explore state.
*/
bool Vole_JuvenileFemale::OnFarmEvent(FarmToDo event) {
	switch (event)
	{
//Here begins soilcultivation mortality
	case sleep_all_day: break;
	case autumn_harrow:
	case autumn_or_spring_plough:
	case autumn_plough:
	case autumn_roll:
	case autumn_sow:
	case autumn_sow_with_ferti:
	case burn_straw_stubble:
	case burn_top:
	case flammebehandling:
	case deep_ploughing:
	case heavy_cultivator_aggregate:
	case hilling_up:
	case preseeding_cultivator:
	case preseeding_cultivator_sow:
	case row_cultivation:
	case spring_harrow:
	case spring_plough:
	case spring_roll:
	case spring_sow:
	case spring_sow_with_ferti:
	case strigling_hill:
	case strigling_sow:
	case stubble_cultivator_heavy:
	case stubble_harrowing:
	case stubble_plough:
	case summer_harrow:
	case summer_plough:
	case summer_sow:
	case winter_harrow:
	case winter_plough:
	case bed_forming:
	case shallow_harrow:
	case shredding:
		if (g_rand_uni_fnc() < VoleSoilCultivationMort)
			CurrentVState = tovs_FDying;
		break;
//Here begins harvest mortality
	case swathing:
	case straw_removal:
	case straw_chopping:
	case straw_covering:
	case fiber_covering:
	case fiber_removal:
	case mow:
	case hay_bailing:
	case hay_turning:
	case harvest:
	case cut_to_hay:
	case cut_to_silage:
	case cut_weeds:
	case bulb_harvest:
	case flower_cutting:
	case green_harvest:
	case harvest_bushfruit:
	case harvestshoots:
		if (g_rand_uni_fnc() < VoleHarvestMort)
			CurrentVState = tovs_FDying;
		break;
//Here begins grazing mortality
	case pigs_out:
		if (g_rand_uni_fnc() < VolePigGrazingMort)
			CurrentVState = tovs_FDying;
		break;
//Here begins strigling mortality
	case strigling:
		if (g_rand_uni_fnc() < VoleStriglingMort)
			CurrentVState = tovs_FDying;
		break;
// Here begins insecticide mortality
	case syninsecticide_treat: /*\todo I added it because it was so for male, was there a reason it wasn't a mortality for female? */
	case insecticide_treat:
		if (g_rand_uni_fnc() < VoleInsecticideMort)
			CurrentVState = tovs_FDying;
		break;
// Here begins herbicide mortality
	case herbicide_treat: /*\todo I added it because it was so for male, was there a reason it wasn't a mortality for female? */
		if (g_rand_uni_fnc() < VoleHerbicicideMort)
			CurrentVState = tovs_FDying;
		break;
//Here begins non-mortality events
	case biocide:
	case fa_ammoniumsulphate:
	case fa_boron:
	case fa_calcium:
	case fa_cu:
	case fa_greenmanure:
	case fa_k:
	case fa_manganesesulphate:
	case fa_manure:
	case fa_n:
	case fa_nk:
	case fa_npk:
	case fa_npks:
	case fa_p:
	case fa_pk:
	case fa_pks:
	case fa_rsm:
	case fa_sk:
	case fa_sludge:
	case fa_slurry:
	case fp_ammoniumsulphate:
	case fp_boron:
	case fp_calcium:
	case fp_cu:
	case fp_greenmanure:
	case fp_k:
	case fp_liquidNH3:
	case fp_manganesesulphate:
	case fp_manure:
	case fp_n:
	case fp_nc:
	case fp_nk:
	case fp_npk:
	case fp_npks:
	case fp_ns:
	case fp_p:
	case fp_pk:
	case fp_pks:
	case fp_rsm:
	case fp_sk:
	case fp_sludge:
	case fp_slurry:
	case fungicide_treat:
	case glyphosate:
	case molluscicide:
	case org_fungicide:
	case org_herbicide:
	case org_insecticide:
	case last_treatment:
	case pheromone:
	case trial_control:
	case trial_insecticidetreat:
	case trial_toxiccontrol:
	case product_treat:
	case cattle_out:
	case cattle_out_low:
	case water:
	case manual_weeding:
	case pruning:
	case suckering:
	case growth_regulator:
		break;
	default:
		g_msg->Warn(WARN_FILE, "Vole_Female::OnFarmEvent(): Unknown event type:",
		            m_OurLandscape->EventtypeToString(event));
		exit(1);
	}
	if (CurrentVState == tovs_FDying)
	{
		if (cfg_RecordVoleMort.value()) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_FFarm, 1);
		return true;
	}
	return false;
}

//----------------------------------------------------------------------------

/**
\brief
Death from external entity.
*/
/**
External event has caused death - probably eaten by a explicitly modelled predator
*/
void Vole_JuvenileFemale::OnKilled() {
	if (cfg_RecordVoleMort.value()) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_FPred, 1);
	CurrentVState = tovs_FDying;
}

//---------------------------------------------------------------------------

void Vole_JuvenileFemale::st_BecomeSubAdult() {
	//cout << "\n I am in Vole_JuvenileFemale::stBecomeSubAdult() and the locations are: ";
	//cout << m_Location_x << " , " << m_Location_y;
	const auto av = new struct_Vole_Adult;
	av->VPM = m_OurPopulation;
	av->L = m_OurLandscape;
	av->m_flag = true; // Used to signal pesticide effect to CreateObjects
	// Create the new voles (50% chance of male/female)
	av->PolyRefBorn = m_OurLandscape->SupplyPolyRef(m_Location_x, m_Location_y);
	av->ElemBorn = m_OurLandscape->SupplyElementType(m_Location_x, m_Location_y);
	av->VegBorn = m_OurLandscape->SupplyVegType(m_Location_x, m_Location_y);
	av->FatherId = m_FatherId;
	av->MotherId = m_MotherId;
	av->xborn=xborn;
	av->yborn=yborn;
	av->x = m_Location_x;
	av->y = m_Location_y;
	av->weight = m_Weight;
	av-> GenerationCount=GenerationCount;
	// Do the genetics
	av->Genes = m_MyGenes;
	av->new_Genes=new_Genes;
	av->new_Mates_Genes=new_Mates_Genes;
	av->age = m_Age;
	m_OurPopulation->CreateObjects(vob_Female, this, av, 1, filename_AlleleInput,true);
	// Remove the current object
	CurrentVState = tovs_FDying;
	FreeLocation();
	m_CurrentStateNo = -1;
	m_StepDone = true;
}

//---------------------------------------------------------------------------

/**
\brief
Main territory evaluation behaviour.
*/
/**
Evaluates the quality of her habitat and does some limited exploration in the surrounding area to see if she can improve it by moving.
*/
int Vole_JuvenileFemale::st_Evaluate_n_Explore() {
#ifdef __USINGTRAPLINES
	CheckTraps();
#endif
	const double Qual = CalculateCarryingCapacity(m_Location_x, m_Location_y, cfg_VoleDDepConst.value());
	// This function determines whether he must leave and whether he does
	// an explore. He will do this if he is not in optimal conditions

	if (Qual < m_MinFVoleHabQual)
	{
		// check an area MinFemaleMovement to FemaleMovement metres away
		m_DispVector = g_random_fnc(8);
		// Not quality dependent dispersal, very directed
		return Dispersal(-1, g_random_fnc(FemaleMovement) + 1);
	}
	return Dispersal(Qual, 1 + g_random_fnc(FemaleMovement) + 1);
	// returns 1 if it has died of extra dispersal mortality or starvation
}

//---------------------------------------------------------------------------

/**
\brief
Location map function
*/
inline void Vole_JuvenileFemale::SetLocation() { m_OurPopulation->m_VoleMap->SetMapValue(m_Location_x, m_Location_y, this); };
/**
\brief
Location map function
*/
inline void Vole_JuvenileFemale::FreeLocation() { m_OurPopulation->m_VoleMap->ClearMapValue(m_Location_x, m_Location_y); };
/**
\brief
Location map function
*/
inline bool Vole_JuvenileFemale::GetLocation(int px, int py) {
	if (m_OurPopulation->m_VoleMap->GetMapValue(px, py) != nullptr) return true;
	return false;
};
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//                           Vole_Female
//---------------------------------------------------------------------------
/**
\brief
Female dispersal
*/
/**
Checks p_Distance away to see if it can find a territory in the next FHabQualThreshold category or with an improved quality of 1.1* p_OldQual if already in optimal habitat
\n This entails some risk though, so there is a fixed 2.5% increase in the mortality chance when it does this.
*/
int Vole_JuvenileFemale::Dispersal(double p_OldQual, int p_Distance) {
	// Returns 1 for die, 0 for carry on
	if (m_DispVector == -1) m_DispVector = g_random_fnc(8); // Choose direction 0-7
	// Go that far in that direction (assuming it is possible to do that)
	const int oldx = m_Location_x;
	const int oldy = m_Location_y;
	MoveTo(m_DispVector, p_Distance, 10);
	//  Now we are there so what is the new quality
	// 1. Get the carrying capacity
	const double CC = CalculateCarryingCapacity(m_Location_x, m_Location_y, cfg_VoleDDepConst.value());
	if (CC <= p_OldQual)
	{
		// Don't want to move
		FreeLocation();
		m_Location_x = oldx;
		m_Location_y = oldy;
		m_OurPopulation->UpdateGuardMap(m_Location_x, m_Location_y, m_guard_cell_x, m_guard_cell_y);
		SetLocation();
	}
	else
	{
		if (g_rand_uni_fnc() < g_extradispmort)
		{
			if (cfg_RecordVoleMort.value()) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_FDisp, 1);
			return 1;
		}
	}
	return 0;
}

//---------------------------------------------------------------------------

/**
\brief
Vole_Female constructor
*/
Vole_Female::Vole_Female(struct_Vole_Adult* p_aVoleStruct, const string& filename_AlleleInput, int i) : Vole_JuvenileFemale(p_aVoleStruct, filename_AlleleInput, i) {
	//cout << "\n I am now in hte Vole_Female constructor and will print genes:";
	//p_aVoleStruct->new_Genes.PrintGenomes();
	//p_aVoleStruct->new_Mates_Genes.PrintGenomes();

	m_NoOfYoung = 0;
	m_Pregnant = false;
	m_TerrRange = m_MinFemaleTerritorySize;
	m_MinTerrRange = m_MinFemaleTerritorySize;
	m_MatesIdNo = -1;
	m_Weight = WeanedWeight;
	m_Age = MinReproAgeF + 1;
	m_DaysUntilBirth = 0;
#ifdef __VOLEPESTICIDEON
	Vole_Female::m_EndoCrineDisruptionGestationLength = cfg_VoleEndoCrineDisruptionGestationLength.value();
	m_pesticideloadindex = 0;
	for (int i = 0; i<21; i++) m_pesticideloadarray[i] = 0;
#endif
}

//---------------------------------------------------------------------------
void Vole_Female::ReInit(struct_Vole_Adult* p_aVoleStruct, const string& filename_AlleleInput, int i) {
	Init(p_aVoleStruct, filename_AlleleInput, i);
	m_NoOfYoung = 0;
	m_Pregnant = false;
	m_TerrRange = m_MinFemaleTerritorySize;
	m_MinTerrRange = m_MinFemaleTerritorySize;
	m_MatesIdNo = -1;
	m_Weight = WeanedWeight;
	m_Age = MinReproAgeF + 1;
	m_DaysUntilBirth = 0;
#ifdef __VOLEPESTICIDEON
	Vole_Female::m_EndoCrineDisruptionGestationLength = cfg_VoleEndoCrineDisruptionGestationLength.value();
	m_pesticideloadindex = 0;
	for (int i = 0; i<21; i++) m_pesticideloadarray[i] = 0;
#endif
}

//---------------------------------------------------------------------------

Vole_Female::~Vole_Female() {
	// Nothing to do
}

//---------------------------------------------------------------------------

/**
\brief
Female vole Step
*/
/**
The Step is one of the three timestep divisions. This is called repeatedly after BeginStep and before EndStep, until all voles report that they are done with Step. \n
\n Most of the behaviours are controlled by moving voles between behavioural states in Step (for other models this is also done in BeginStep and EndStep).
\n When a vole is done for the day it will signal this by setting m_StepDone==true. NB that a call to one behaviour may trigger a call to another behaviour on the next call to step inside the same timestep. In this way a daily cycle of activity can be undertaken (i.e. do reproduction and explore)
*/
void Vole_Female::Step() {
	if (m_StepDone) return;
	switch (CurrentVState)
	{
	case 0: // Initial state
		CurrentVState = tovs_FEvaluateExplore;
		break;
	case tovs_FEvaluateExplore: // Evaluate & Explore
		switch (st_Evaluate_n_Explore())
		{
		case 1:
			CurrentVState = tovs_FDying; // Die
			break;
		case 2:
			if (cfg_RecordVoleMort.value()) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_FStarve, 1);
			CurrentVState = tovs_FDying; // Die
			break;
		default:
			CurrentVState = tovs_ReproBehaviour; // repro behaviour
		}
		m_StepDone = true;
		break;
	case tovs_ReproBehaviour: // ReproBehaviour
		CurrentVState = st_ReproBehaviour();
		break;
	case tovs_Lactating: // Lactating
		CurrentVState = st_Lactating();
		m_StepDone = true;
		break;
	case tovs_GiveBirth:
		CurrentVState = st_GiveBirth();
		m_StepDone = true; // wait until next day
		break;
	case tovs_FMaturation: // Maturation
		CurrentVState = st_BecomeReproductive(); // only returns Eval n Explore
		break;
	case tovs_Mating: // Mating
		CurrentVState = st_Mating(); // only returns Eval n Explore
		break;
	case tovs_UpdateGestation: // Update Gestation
		CurrentVState = st_UpdateGestation();
		break;
	case tovs_SpecialExplore: // Special Explore
		st_Special_Explore();
		CurrentVState = tovs_ReproBehaviour; // Repro behaviour
		m_StepDone = true;
		break;
	case tovs_FDying: // Die
		FreeLocation();
		m_StepDone = true;
		break;
	default:
		m_OurLandscape->Warn("Vole_Female::Step - unknown state", nullptr);
		exit(9);
	}
}

//---------------------------------------------------------------------------

/**
\brief
Reproductive switch
*/
/**
* This state is simply a switch determing what behaviour to exhibit. In the case of pesticide exposure to endocrine disruptors
* there may be an impact of mating with a steriel male leading to a false pregnancy. In this case m_DaysUntilBirth will be
* > 0, but she will not be pregnant. This will have the effect of delaying mating until the counter is zero.
*/
TTypeOfVoleState Vole_Female::st_ReproBehaviour() const {
	if (!m_Mature) return tovs_FMaturation; // Go to BecomeReproductive;
	if (m_NoOfYoung > 0) return tovs_Lactating; //Go to lactating
	if (!m_Pregnant && m_DaysUntilBirth == 0) return tovs_Mating; // Go to Mating;
	return tovs_UpdateGestation;
	// Go to UpdateGestation;
}

//---------------------------------------------------------------------------

/**
\brief
Gestation control.
*/
/**
Decreases the number of days until birth by 1, checks for special pesticide effects for vinclozolin like pesticide simulations
*/
TTypeOfVoleState Vole_Female::st_UpdateGestation() {
	m_DaysUntilBirth--;
	if (m_DaysUntilBirth == 0) return tovs_GiveBirth;
	return tovs_FEvaluateExplore;
}

//---------------------------------------------------------------------------

/**
\brief
Female vole maturation control.
*/
/**
Tests to see if the female should mature
*/
TTypeOfVoleState Vole_Female::st_BecomeReproductive() {
	if (m_BreedingSeason && m_Have_Territory) { if (m_Age >= MinReproAgeF) { Setm_Mature(); } }
	return tovs_FEvaluateExplore;
}

//---------------------------------------------------------------------------

/**
\brief
Litter production.
*/
/**
Produces a litter and records the information if necessary
*/
TTypeOfVoleState Vole_Female::st_GiveBirth() {
	// TODO must locate the nest position if this is considered important - currently it is set as the centre of the territory
	if (m_Have_Territory && m_Pregnant)
	{
		// Must produce young
		int yr = 0;
		if (m_BornLastYear) yr = 1;
		m_NoOfYoung = m_OurPopulation->ReproTable[yr][m_OurLandscape->SupplyDayInYear() / 30];
		m_Pregnant = false;
#ifdef __VOLEPESTICIDEON
		if ((m_OurLandscape->SupplyPesticideType() == ttop_ModelinkPesticide) || (m_OurLandscape->SupplyPesticideType() == ttop_ModelinkPesticide21TWA))
		{
			double dose = 0.0;
			double survive = 1.0;
			if (m_OurLandscape->SupplyPesticideType() == ttop_ModelinkPesticide21TWA)
			{
				for (int i = 0; i < 21; i++) dose += m_pesticideloadarray[i];
				dose /= 21;
				dose -= cfg_PesticideAccumulationThresholdModelink2.value();
			}
			else
			{
				dose = m_pesticideInfluenced3 - cfg_PesticideAccumulationThreshold.value();
				m_pesticideInfluenced3 = 0;
			}
			if (dose > 0) // this has the maximum dose received
			{
				// Calc dose response in terms of chance of death for individual offspring
				survive = 1 - (cfg_PesticideLitterSizeReduction.value() * dose);
			}
			// This applies the survival probability to each individual - why do we not do this to the litter as a whole
			// - this is because the probability may be too small to get a % response from e.g. 4 young
			int ny = 0;
			for (int y = 0; y < m_NoOfYoung; y++) if (g_rand_uni() < survive)
			{
				ny++;
				m_OurPopulation->AddToNotImpacted();
			}
			else m_OurPopulation->AddToImpacted();
			m_NoOfYoung = ny;
		}
		if (m_NoOfYoung > 0) m_OurPopulation->incLittersProduced(); else m_OurPopulation->incLittersLost();
#endif
		// Must age young
		m_YoungAge = 0;
		// Update the YoungProduced Today Counter
		m_OurPopulation->AddToYoung(m_NoOfYoung);
		m_NoOfYoungTotal += m_NoOfYoung;
		return tovs_ReproBehaviour; //NB if m_NoOfYoung is zero due to pesticide effects this will quietly be dealt with by the repro switch
	}
	m_NoOfYoung = 0;
	m_Pregnant = false;
	return tovs_FEvaluateExplore;
}

//---------------------------------------------------------------------------
/**
\brief
Lactation.
*/
/**
Once the litter reaches weaning age then individual voles are created and set in motion by this method. Pesticide effects may be specified here too.
*/
TTypeOfVoleState Vole_Female::st_Lactating() {
	if (m_YoungAge++ >= WeanedAge)
	{
		auto av = new struct_Vole_Adult;
		av->VPM = m_OurPopulation;
		av->L = m_OurLandscape;
		av->m_flag = true; // Used to signal pesticide effect to CreateObjects
#ifdef __VOLEPESTICIDEON
	  double mdelay = 0;
	  if ((m_OurLandscape->SupplyPesticideType() == ttop_ModelinkPesticide)|| (m_OurLandscape->SupplyPesticideType() == ttop_ModelinkPesticide21TWA))
	  {
		  double dose = 0.0;
		  double survive = 1.0;
		  if (m_OurLandscape->SupplyPesticideType() == ttop_ModelinkPesticide21TWA)
		  {
			  for (int i=0; i<21; i++) dose+=m_pesticideloadarray[i];
			  dose /=21;
			  dose -= cfg_PesticideAccumulationThresholdModelink2.value();
		  } else
		  {
			dose = m_pesticideInfluenced3-cfg_PesticideAccumulationThreshold.value();
			m_pesticideInfluenced3=0;
		  }
		  if ( dose > 0 ) // this has the maximum dose received
		  {
			  // Calc dose response in terms of chance of death for individual offspring
			  survive = 1-(cfg_PesticideWeaningReduction.value() * dose);
			  mdelay = cfg_PesticideFemaleMaturityDelay.value() * dose;
		  }
		  // This applies the survival probability to each individual - why do we not do this to the litter as a whole
		  // - this is because the probability may be too small to get a % response from e.g. 4 young
		  int ny=0;
		  for (int y=0; y<m_NoOfYoung; y++) if (g_rand_uni()<survive) ny++;
		  m_NoOfYoung = ny;
	  }
#endif
		const unsigned MotherIdNo = SupplyIDNo();
		const unsigned FatherIdNo = SupplyMateId();

		// Create the new voles (50% chance of male/female)
		for (int i = 0; i < m_NoOfYoung; i++)
		{
			av->PolyRefBorn = m_OurLandscape->SupplyPolyRef(m_Location_x, m_Location_y);
			av->ElemBorn = m_OurLandscape->SupplyElementType(m_Location_x, m_Location_y);
			av->VegBorn = m_OurLandscape->SupplyVegType(m_Location_x, m_Location_y);
			av->FatherId = FatherIdNo;
			av->MotherId = MotherIdNo;
	
			av->x = (m_Location_x + g_random_fnc(10)) % m_OurPopulation->SupplySimW();
			av->y = (m_Location_y + g_random_fnc(10)) % m_OurPopulation->SupplySimH();

			av->xborn=av->x; // fill into the pointer the x coordinate where the animal was born
			av->yborn=av->y;// fill into the pointer the x coordinate where the animal was born

			// Do the genetics (old)
			av->Genes.Recombine(&m_MyGenes, &m_MatesGenes); //old genetics
			av->Genes.Mutation_4();

			//Doing the new genetics and putting them to the new pointer
			av->new_Genes.Genome0= this->new_Genes.MakeGameteSimple(); //make gamete from mom's genes
			av->new_Genes.Genome1= this->new_Mates_Genes.MakeGameteSimple(); //make gamete from dads' genes
			av->new_Genes.MitochondrialLine = this -> new_Genes.MitochondrialLine; //pass on mitochrondrial/maternal line from mom
			av -> new_Genes.YChromoLine= this -> new_Mates_Genes.YChromoLine; //pass on ychromosome/paternal line from mom dad
			
			//Getting generationcount intp pointer
			int GenCountHolder=this->GenerationCount;
			av->GenerationCount = GenCountHolder + 1;
#ifdef __Mutation
			//av->Genes256_16.Mutation_3(); // This is only one of a number of mutation types we can use
#endif
			if (g_random_fnc(2) == 1) // Even sex ratio
			{
				// Males
#ifdef __VOLEPESTICIDEON
			if (m_OurLandscape->SupplyPesticideType() == ttop_Vinclozolin)
			{
			// If her mate was pesticide influenced then ensure all male offspring are also
			// by setting the gene on if we are useing Vinclozolin
			// If either the male was directly affected, or the female directly during gestation
			if ( m_pesticideInfluenced1 )
			{
				av->m_dflag = true;
				av->m_gflag = true;
				m_OurPopulation->AddToImpacted();
			}
			// Otherwise if the male was genetically affected
			else if (MatesGenes.GetGeneticFlag())
			{
				av->m_gflag = true;
				av->m_dflag=false;
				m_OurPopulation->AddToGeneticImpacted();
			}
			else
			{
				av->m_gflag=false;
				av->m_dflag=false;
				m_OurPopulation->AddToNotImpacted();
			}
		}
		else if (m_OurLandscape->SupplyPesticideType() == ttop_ReproductiveEffects )
		{
			// If either the male was directly affected, or the female directly during gestation
			if ( m_pesticideInfluenced1 )
			{
				av->m_dflag = true;
				av->m_gflag = false;
				m_OurPopulation->AddToImpacted();
			}
			else
			{
				av->m_gflag=false;
				av->m_dflag=false;
				m_OurPopulation->AddToNotImpacted();
			}
		}
		else
		{
			av->m_gflag=false;
			av->m_dflag=false;
		}
#endif
				m_OurPopulation->CreateObjects(vob_JuvenileMale, this, av, 1, filename_AlleleInput, true);
			}
			else
			{
#ifdef __VOLEPESTICIDEON
			if ((m_OurLandscape->SupplyPesticideType() == ttop_ModelinkPesticide)|| (m_OurLandscape->SupplyPesticideType() == ttop_ModelinkPesticide21TWA))
			{
				av->misc_use = mdelay;
			}
			else av->misc_use = 0.0;
#endif
				m_OurPopulation->CreateObjects(vob_JuvenileFemale, this, av, 1, filename_AlleleInput, true);
			}
		}
		m_OurPopulation->AddToJuvs(m_NoOfYoung);
		m_Reserves = 0; // all reserves used up
		m_NoOfYoung = 0; // No more young to feed
		delete av; // clean up
		return tovs_SpecialExplore; //  CurrentVState=14 Special Explore
	}
	return tovs_Lactating;
	// carry on feeding young
}

//---------------------------------------------------------------------------

/**
\brief
Main territory evaluation behaviour.
*/
/**
Evaluates the quality of her habitat and does some limited exploration in the surrounding area to see if she can improve it by moving.
*/
int Vole_Female::st_Evaluate_n_Explore() {
#ifdef __USINGTRAPLINES
    CheckTraps();
#endif
	/*
		unsigned MovementIndex;
	    if (m_Age>90) MovementIndex=3;
	     else if (m_Age>60) MovementIndex=2;
	       else if (m_Age>30) MovementIndex=1;
	         else MovementIndex=0;
	*/
	const double Qual = CalculateCarryingCapacity(m_Location_x, m_Location_y, cfg_VoleDDepConst.value());
	// This function determines whether he must leave and whether he does
	// an explore. Will do this if he is not in optimal conditions
	if (Qual <= m_MinFVoleHabQual)
	{
#ifdef __VoleStarvationDays
	   ++m_StarvationDays;
#endif
		m_Have_Territory = false;
		// check an area MinFemaleMovement to FemaleMovement metres away
		m_DispVector = g_random_fnc(8);
		// Not quality dependent dispersal, very directed
		return Dispersal(-1, g_random_fnc(FemaleMovement) + MinFemaleMovement);
	}
	m_Have_Territory = true;
#ifdef __VoleStarvationDays
	   m_StarvationDays = 0;
#endif
	return Dispersal(Qual, g_random_fnc(FemaleMovement) + MinFemaleMovement); // returns 1 if it has died of extra dispersal mortality or starvation
}

//---------------------------------------------------------------------------

/**
\brief
Post weaning territory expansion.
*/
/**
Called after weaning of a litter when the female needs to expand her territory again. This does not use the Dispersal function so avoids extra mortality at this point.
*/
int Vole_Female::st_Special_Explore() {
	//The vole must now find a good territory
	// She will explore in all 8 directions from the nest and pick the best
	// from these and here present position
	unsigned coords[9][2];
	double qualities[9];
	// Check the quality here
	coords[8][0] = m_Location_x;
	coords[8][1] = m_Location_y;
	if (!m_OurPopulation->SupplyOlderFemales(m_Location_x, m_Location_y, m_Age, m_TerrRange))
		qualities[8] = -1;
	else { qualities[8] = CalculateCarryingCapacity(m_Location_x, m_Location_y, cfg_VoleDDepConst.value()); }
	// Now check the quality in all 8 directions
	const unsigned steps = 2 * FemaleMovement + MinFemaleMovement;
	for (unsigned i = 0; i < 8; i++)
	{
		MoveTo(i, steps, 10); // changes m_Location_x & y
		coords[i][0] = m_Location_x;
		coords[i][1] = m_Location_y;
		if (!m_OurPopulation->SupplyOlderFemales(m_Location_x, m_Location_y, m_Age,
		                                         m_TerrRange))
			qualities[i] = -1;
		else { qualities[i] = CalculateCarryingCapacity(m_Location_x, m_Location_y, cfg_VoleDDepConst.value()); }
	}
	// find the best out of the nine qualities
	double best = qualities[8]; // here
	unsigned found = 8; // default stay here
	for (unsigned i = 0; i < 8; i++)
	{
		if (qualities[i] > best)
		{
			best = qualities[i];
			found = i;
		}
		if (best > 1) { if (m_Mature) { m_Have_Territory = true; } }
		else { m_Have_Territory = false; }
		FreeLocation();
		m_Location_x = coords[found][0];
		m_Location_y = coords[found][1];
		m_OurPopulation->UpdateGuardMap(m_Location_x, m_Location_y, m_guard_cell_x, m_guard_cell_y);
		SetLocation();
	}
	return 0;
}

//---------------------------------------------------------------------------

/**
\brief
Female mating.
*/
/**
The mating state is where genes are passed and any genetic effects from the male need to be evaluated (e.g. pesticide induced infertility)
*/
TTypeOfVoleState Vole_Female::st_Mating() {
	//if it is the breeding season
	if (m_BreedingSeason)
	{
		// Find a male if there is one

		Vole_Male* Mate = m_OurPopulation->FindClosestMale(m_Location_x, m_Location_y, m_MaxMaleTerritorySize); // Changed from a constant of 20 to her Territory range **//TD 090519
		// Get his genes if we found a mate
		if (Mate)
		{
			m_MatesGenes = Mate->SupplyGenes();
			m_MatesIdNo = Mate->SupplyIDNo();
			new_Mates_Genes= Mate->new_Genes;
#ifdef __VOLEPESTICIDEON
		if (!Mate->GetFertile())
		{
			m_Pregnant = false;
		}
		else m_Pregnant = true;
#else
			m_Pregnant = true;
#endif
			m_DaysUntilBirth = TheGestationPeriod;
		}
	}
	return tovs_FEvaluateExplore;
}

//---------------------------------------------------------------------------

/**
\brief
Determines whether an infanticide attempt will succeed.
*/
/**
     A mate has attempted infanticide.\n
	 If they are 9 days old there is no chance of mortality
     otherwise it is proportional to age. \n
     Data from M.arvalis (Heise & Lippke 1997)
*/
void Vole_Female::OnInfanticideAttempt() {
	//try this if she has kids
	if (m_NoOfYoung)
	{
		if (m_YoungAge < 9)
		{
			const double mortchance = g_rand_uni_fnc();
			if (mortchance <= InfanticideChanceByAge[m_YoungAge])
			{
				m_OurPopulation->AddToNoYoungInfanticideCount(m_NoOfYoung, m_YoungAge);

				// Young are killed
				m_NoOfYoung = 0; // will cause the vole to be mated
				m_Pregnant = false;
			}
		}
	}
}

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//                           Vole_Male
//---------------------------------------------------------------------------

/**
\brief
Vole_Male constructor.
*/
/**
- just calls Init()
*/
Vole_Male::Vole_Male(struct_Vole_Adult* p_aVoleStruct, const string& filename_AlleleInput, int i) : Vole_JuvenileMale(p_aVoleStruct, filename_AlleleInput, i) {
	//cout << "\n I am now in hte Vole_Male  constructor and will print positions: ";
	//cout << p_aVoleStruct->x,p_aVoleStruct->y;
	//p_aVoleStruct->new_Genes.PrintGenomes();
	//p_aVoleStruct->new_Mates_Genes.PrintGenomes();
	m_Age = p_aVoleStruct->age;
	m_Weight = p_aVoleStruct->weight;
	m_BirthYear = p_aVoleStruct->BirthYear;
	m_Sex = true;
	m_TerrRange = m_MinMaleTerritorySize;
	m_MinTerrRange = m_MinMaleTerritorySize;
}

//---------------------------------------------------------------------------

void Vole_Male::ReInit(struct_Vole_Adult* p_aVoleStruct,const string& filename_AlleleInput, int i) {
	Init(p_aVoleStruct, filename_AlleleInput, i);
	m_Age = p_aVoleStruct->age;
	m_Weight = p_aVoleStruct->weight;
	m_BirthYear = p_aVoleStruct->BirthYear;
	m_Sex = true;
	m_TerrRange = m_MinMaleTerritorySize;
	m_MinTerrRange = m_MinMaleTerritorySize;
}

//---------------------------------------------------------------------------

Vole_Male::~Vole_Male() {
	// Nothing to do
}

//---------------------------------------------------------------------------

/**
\brief
Male vole Step
*/
/**
The Step is one of the three timestep divisions. This is called repeatedly after BeginStep and before EndStep, until all voles report that they are done with Step. \n
\n Most of the behaviours are controlled by moving voles between behavioural states in Step (for other models this is also done in BeginStep and EndStep).
\n When a vole is done for the day it will signal this by setting m_StepDone==true. NB that a call to one behaviour may trigger a call to another behaviour on the next call to step inside the same timestep. In this way a daily cycle of activity can be undertaken.
*/
void Vole_Male::Step() {
	if (m_StepDone || m_CurrentStateNo == -1) return;
	switch (CurrentVState)
	{
	case tovs_InitialState: //initial state
		CurrentVState = tovs_MEvaluateExplore;
		break;
	case tovs_MMaturation: // Maturation
		if (st_Maturation()) Setm_Mature();
		CurrentVState = tovs_MEvaluateExplore; // Eval and Explore
		break;
	case tovs_MEvaluateExplore: // Eval&Explore
		CurrentVState = st_Eval_n_Explore();
		m_StepDone = true;
		break;
	case tovs_Infanticide: // Infanticide
		st_Infanticide();
		CurrentVState = tovs_MEvaluateExplore;
		break;
	case tovs_MDying:
		FreeLocation();
		m_StepDone = true;
		break;
	default:
		m_OurLandscape->Warn("Vole_Male::Step - unknown return error", nullptr);
		exit(1);
	}
}

//---------------------------------------------------------------------------

/**
\brief
Male vole EndStep
*/
/**
The EndStep is one of the three timestep divisions. This is called once for each vole after BeginStep and Step. \n
The main function here is to remove voles that have died during step and otherwise to grow if not at max weight.
It also checks if the vole was killed due to human management and determines the potential territory size.
*/
void Vole_Male::EndStep() {
#ifdef __VoleStarvationDays
  if (m_Have_Territory) m_StarvationDays = 0;
  if ( m_StarvationDays > m_MaxStarvationDays)
  {
	  st_Dying();
	  if (cfg_RecordVoleMort.value() ) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_MStarve,1);
  }
  else
#endif
	{
		CheckManagement();
		if (CurrentVState == tovs_MDying) { st_Dying(); }
		else
		{
			m_Age++;
			/*              MALE GROWTH NOTES

			 Male vole grows until 20g. After that he will only grow if he has matured
			 Growth continues up to 60g.

			 Reproduction cannot occur below 40g or 40 days

			 Growth only occurs between 1 March and 1st August
			 taken from Hanson L, 1977, Oikos 29.
			*/
			const int today = m_OurLandscape->SupplyDayInYear();
			if (today < GrowStopDate && today > m_OurPopulation->SupplyGrowthStartDate())
			{
				if (m_Weight < 20) { m_Weight += growthperdayM; }
				else if (m_Mature == true && m_Weight < MaxWeightM) { m_Weight += growthperdayM; }
			}
			else if (today == 1)
			{
				m_BornLastYear = true;
				m_Mature = true;
			} // must be true if alive on 1st Jan.
			DetermineTerritorySize();
		}
	}
}

//---------------------------------------------------------------------------

/**
\brief
Male vole maturation control.
*/
/**
Decide whether to become mature or not
*/
int Vole_Male::st_Maturation(void) const {
	if (m_OurLandscape->SupplyDayInYear() <
		m_OurPopulation->SupplyGrowthStartDate())
		return 0; // don't mature
	return 1;
	// matures
}

//---------------------------------------------------------------------------

/**
\brief
Male vole infanticide behaviour.
*/
/**
 Will only enter here if have taken over a new area a reasonable
 distance from the old one and will only commit infanticide in 10% of his attempts in his minimum territory range. \n

Tells the population manager to send an infanticide message to
all females in the territory.
*/
void Vole_Male::st_Infanticide(void) const {
	if (g_rand_uni_fnc() < cfg_InfanticideProbability.value())
	{
		// int SearchInf =(int) m_TerrRange*cfg_InfanticideRangeRelToTerRange.value();
		m_OurPopulation->SendMessage(tovm_Infanticide, m_Location_x, m_Location_y, m_MinTerrRange, false, m_Weight);
	}
}

//---------------------------------------------------------------------------

/**
\brief
Male vole main territory assessment behaviour.
*/
/**
Evaluates the quality of his habitat and does some limited exploration in the surrounding area to see if she can improve it by moving.
*/
TTypeOfVoleState Vole_Male::st_Eval_n_Explore(void) {
	VoleDispersalReturns retval = vdisp_CarryOn;

	unsigned MovementIndex;
	if (m_Age > 90) MovementIndex = 3;
	else if (m_Age > 60) MovementIndex = 2;
	else if (m_Age > 30) MovementIndex = 1;
	else MovementIndex = 0;
	const int OldMales = m_OurPopulation->SupplyInOlderTerr(m_Location_x, m_Location_y, m_Age, m_TerrRange); // returns -1 if a male has p_x,p_y in his territory and is older than p_Age else returns the number of females present
	const int NoFemales = OldMales;
	const double Qual = CalculateCarryingCapacity(m_Location_x, m_Location_y, 10000);
	//if (OldMales>1) Qual/=OldMales; // Share resources with the females, but no-one else
	if (m_BreedingSeason)
	{
		if (m_Have_Territory) // Must be mature
		{
#ifdef __USINGTRAPLINES
			CheckTraps();
#endif
			if (!m_Mature)
			{
				return tovs_MMaturation; // Special case, immature male in breeding season with territory will mature
			}
			// If he his mature then he must move if in an older males territory
			bool NoMates = false;
			if (NoFemales < 1)
			{
				if (g_rand_uni_fnc() < g_NoFemalesMove)
				{
					const Vole_Female* af = m_OurPopulation->FindClosestFemale(m_Location_x, m_Location_y, MaleMovement[MovementIndex]);
					if (af != nullptr)
					{
						const AnimalPosition ap = af->SupplyPosition();
						int tx = (m_SimW + ap.m_x + (g_random_fnc(2) - 1)) % m_SimW;
						int ty = (m_SimH + ap.m_y + (g_random_fnc(2) - 1)) % m_SimH;
						while (GetLocation(tx, ty) && MoveQuality(tx, ty) < 1)
						{
							tx = (m_SimW + tx + (g_random_fnc(2) - 1)) % m_SimW;
							ty = (m_SimH + ty + (g_random_fnc(2) - 1)) % m_SimH;
						}
						FreeLocation();
						m_Location_x = tx;
						m_Location_y = ty;
						m_OurPopulation->UpdateGuardMap(m_Location_x, m_Location_y, m_guard_cell_x, m_guard_cell_y);
						SetLocation();
					}
					else NoMates = true;
				}
			}
			else if ((OldMales < 1 || Qual < m_MinFVoleHabQual || NoMates) && !m_BornLastYear) // In larger males territory, must move or are no females
			{
				m_Have_Territory = false;
				m_DispVector = g_random_fnc(8);
				retval = Dispersal(Qual, MaleMovement[MovementIndex]);
			}
		}
		else
		{
			// No territory & Breeding season
#ifdef __USINGTRAPLINES
			CheckTraps();
#endif
			// Not quality dependent dispersal, very directed
			retval = Dispersal(-1, MaleMovement[MovementIndex]);
		}
	}
	else // Not breeding season
	{
#ifdef __USINGTRAPLINES
		CheckTraps();
#endif
		if (OldMales < 1 || Qual < m_MinFVoleHabQual)
		{
#ifdef __VoleStarvationDays
				m_StarvationDays++;
#endif
			retval = Dispersal(-1, 1 + g_random_fnc(MaleMovement[MovementIndex]));
		}
	}
	// retval now holds the answer to the returning state
	switch (retval)
	{
	case vdisp_Mature:
		Setm_Mature();
		return tovs_Infanticide;
	case vdisp_Infanticide:
		return tovs_Infanticide;
	case vdisp_Die:
		return tovs_MDying;
	case vdisp_CarryOn:
	default:
		return CurrentVState;
	}
}

//---------------------------------------------------------------------------

/**
\brief
Male vole dispersal behaviour.
*/
/**
Works like female dispersal - but a return code of 3 will trigger infanticide \n
Checks p_Distance away to see if it can find a territory in the next MHabQualThreshold category or with \n
an improved quality of 1.1* p_OldQual if already in optimal habitat\n
This entails some risk though, so there is a fixed increase in the mortality chance when it does this.
*/
VoleDispersalReturns Vole_Male::Dispersal(double p_OldQual, int p_Distance) {
	// int retval = 0;
	// p_OldQuatells whether dispersal is conditional on quality or not
	// p_OldQual is set to old habitat quality or -1
	// p_Distance is the p_Distance used by the move function
	// aim is to move in a directed way traversing the landscape using the best
	// habitats
	// remember the old co-ordinates in case infanticide needs to be triggered
	const int oldx = m_Location_x;
	const int oldy = m_Location_y;
#ifdef __VOLE_SMALL_LANDSCAPE
	m_DispVector = g_random_fnc(8); // Choose direction 0-7
#else
	if (m_DispVector == -1) m_DispVector = random(8); // Choose direction 0-7
#endif
	// Go that far in that direction (assuming it is possible to do that)
	MoveTo(m_DispVector, p_Distance, 10);
	//  Now we are there so what is the new quality
	// 1. Get the carrying capacity
	int OldMales = 1;
	int NoFemales = 1;
	if (m_BreedingSeason)
	{
		OldMales = m_OurPopulation->SupplyInOlderTerr(m_Location_x, m_Location_y, m_Age, m_TerrRange);
		NoFemales = OldMales; // If can be 0 or >0 if OldMales does not return 1
	}
	const double CC = CalculateCarryingCapacity(m_Location_x, m_Location_y, cfg_VoleDDepConst.value());
	if (p_OldQual == -1)
	{
		// Do an extra mortality test if we are forced to disperse
		if (g_rand_uni_fnc() < g_extradispmort)
		{
			if (cfg_RecordVoleMort.value()) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_MDisp, 1);
			return vdisp_Die;
		}
	}

	if (CC < p_OldQual || OldMales < 1 || NoFemales < 1)
	{
		// Don't want to move
		FreeLocation();
		m_Location_x = oldx;
		m_Location_y = oldy;
		m_OurPopulation->UpdateGuardMap(m_Location_x, m_Location_y, m_guard_cell_x, m_guard_cell_y);
		SetLocation();
		return vdisp_CarryOn;
	}
	if (CC > m_MinMVoleHabQual)
	{
		m_Have_Territory = true;
		if (m_BreedingSeason)
		{
			if (!m_Mature)
			{
				return vdisp_Mature; // Will also trigger infanticde
			}
			if (abs(oldx - m_Location_x) > 2 * m_TerrRange || abs(oldy - m_Location_y) > 2 * m_TerrRange)
			{
				// Have we moved more than _TerrRange from old home
				return vdisp_Infanticide; // This will result in an infanticide attempt if females in the new territory have young and the male has moved outside his old terr range.
			}
		}
	}
	return vdisp_CarryOn;
}

//---------------------------------------------------------------------------
/**
\brief
Do a mortality test
*/
/**
Takes both physiological lifespan and background mortality into account to determine whether the vole should die - repeated calls increase the risk of dying
*/
bool Vole_Male::MortalityTest() {
	// returns true if the vole should die
	if (m_Have_Territory) { if (g_rand_uni_fnc() < g_DailyMortChanceMaleTerr) return true; }
	else if (g_rand_uni_fnc() < g_DailyMortChance) return true;
	return false;
}

//---------------------------------------------------------------------------
/**
\brief
Calculates the territory size needed for a vole of his weight
*/
inline void Vole_Male::DetermineTerritorySize() {
	if (m_Weight < MinReproWeightM) return;
	m_TerrRange = static_cast<int>(m_MinMaleTerritorySize + m_MaleTerritoryRangeSlope * (m_Weight - MinReproWeightM));
}

//---------------------------------------------------------------------------

/**
\brief
Currently not used.
*/
/**
* Check whether our location is of sufficient quality to allow use to feed, ie. is it over or under min
* suitable quality (m_MinMVoleHabQual). \n
*/
inline bool Vole_Male::CanFeed() const {
	if (CalculateCarryingCapacity(m_Location_x, m_Location_y, cfg_VoleDDepConst.value()) >= m_MinMVoleHabQual) return true;
	return false;
}

//---------------------------------------------------------------------------

void Vole_Base::CheckTraps() {
	if (m_OurLandscape->SupplyYearNumber() < 10) return; // Don't sample until the population has stabilised
	const int r = m_TerrRange;
	const int vx = m_Location_x;
	const int vy = m_Location_y;
	for (int rad = 0; rad < r; rad++)
	{
		// Two Loops
		for (int x = vx - rad; x <= vx + rad; x++)
		{
			int y = vy + rad;
			if (m_OurPopulation->IsTrap(x, y))
			{
				if (g_random_fnc(1) == 0)
				{
					m_intrappos.m_x = x;
					m_intrappos.m_y = y;
					m_intrappos.m_EleType = m_OurLandscape->SupplyElementType(x, y);
					m_intrappos.m_VegType = m_OurLandscape->SupplyVegType(x, y);
					m_intrappos.m_inAtrap = true;
					return;
				}
			}
			y = vy - rad;
			if (m_OurPopulation->IsTrap(x, y))
			{
				if (g_random_fnc(1) == 0)
				{
					m_intrappos.m_x = x;
					m_intrappos.m_y = y;
					m_intrappos.m_EleType = m_OurLandscape->SupplyElementType(x, y);
					m_intrappos.m_VegType = m_OurLandscape->SupplyVegType(x, y);
					m_intrappos.m_inAtrap = true;
					return;
				}
			}
		}
		for (int y = vy - (rad - 1); y <= vy + (rad - 1); y++)
		{
			int x = vx + rad;
			if (m_OurPopulation->IsTrap(x, y))
			{
				if (g_random_fnc(1) == 0)
				{
					m_intrappos.m_x = x;
					m_intrappos.m_y = y;
					m_intrappos.m_EleType = m_OurLandscape->SupplyElementType(x, y);
					m_intrappos.m_VegType = m_OurLandscape->SupplyVegType(x, y);
					m_intrappos.m_inAtrap = true;
					return;
				}
			}
			x = vx - rad;
			if (m_OurPopulation->IsTrap(x, y))
			{
				if (g_random_fnc(1) == 0)
				{
					m_intrappos.m_x = x;
					m_intrappos.m_y = y;
					m_intrappos.m_EleType = m_OurLandscape->SupplyElementType(x, y);
					m_intrappos.m_VegType = m_OurLandscape->SupplyVegType(x, y);
					m_intrappos.m_inAtrap = true;
					return;
				}
			}
		}
	}
}

/*
*****************************************************************************************
******************************** Juvenile Male Code *************************************
*****************************************************************************************
*/

/**
\brief
Extra movement on weaning.
*/
/**
Cause the vole to do some exploration on maturity, just to move him away from the litter centre. \n
After this first day he will go into the normal Eval_n_Explore.
*/
void Vole_JuvenileMale::st_JuvenileExplore(void) {
	// Just do some movement at first
	MoveTo(g_random_fnc(8), MinMaleMovement, 20);
}

//---------------------------------------------------------------------------
/**
\brief
JuvenileMale vole exernal event handler.
*/
/**
This method evaluates external events and chooses a suitable response (in this case a probability of dying because other effects will be taken up by the evaluate and explore state. \n
This method is inherited by the adult male class
*/
bool Vole_JuvenileMale::OnFarmEvent(FarmToDo event) {
	switch (event)
	{
	case sleep_all_day:
		break;
//Here begins soilcultivation mortality
	case autumn_harrow:
	case autumn_or_spring_plough:
	case autumn_plough:
	case autumn_roll:
	case autumn_sow:
	case autumn_sow_with_ferti:
	case burn_straw_stubble:
	case burn_top:
	case flammebehandling:
	case deep_ploughing:
	case heavy_cultivator_aggregate:
	case hilling_up:
	case preseeding_cultivator:
	case preseeding_cultivator_sow:
	case row_cultivation:
	case spring_harrow:
	case spring_plough:
	case spring_roll:
	case spring_sow:
	case spring_sow_with_ferti:
	case strigling_hill:
	case strigling_sow:
	case stubble_cultivator_heavy:
	case stubble_harrowing:
	case stubble_plough:
	case summer_harrow:
	case summer_plough:
	case summer_sow:
	case winter_harrow:
	case winter_plough:
	case bed_forming:
	case shallow_harrow:
	case shredding:
		if (g_rand_uni_fnc() < VoleSoilCultivationMort)
			CurrentVState = tovs_MDying;
		break;
		//Here begins harvest mortality
	case swathing:
	case straw_removal:
	case straw_chopping:
	case straw_covering:
	case fiber_covering:
	case fiber_removal:
	case mow:
	case hay_bailing:
	case hay_turning:
	case harvest:
	case cut_to_hay:
	case cut_to_silage:
	case cut_weeds:
	case bulb_harvest:
	case flower_cutting:
	case green_harvest:
	case harvest_bushfruit:
	case harvestshoots:
		if (g_rand_uni_fnc() < VoleHarvestMort)
			CurrentVState = tovs_MDying;
		break;
		//Here begins grazing mortality
	case pigs_out:
		if (g_rand_uni_fnc() < VolePigGrazingMort)
			CurrentVState = tovs_MDying;
		break;
		//Here begins strigling mortality
	case strigling:
		if (g_rand_uni_fnc() < VoleStriglingMort)
			CurrentVState = tovs_MDying;
		break;
		// Here begins insecticide mortality
	case syninsecticide_treat: /*\todo I added it because it was so for male, was there a reason it wasn't a mortality for female? */
	case insecticide_treat:
		if (g_rand_uni_fnc() < VoleInsecticideMort)
			CurrentVState = tovs_MDying;
		break;
		// Here begins herbicide mortality
	case herbicide_treat: /*\todo I added it because it was so for male, was there a reason it wasn't a mortality for female? */
		if (g_rand_uni_fnc() < VoleHerbicicideMort)
			CurrentVState = tovs_MDying;
		break;
		//Here begins non-mortality events
	case biocide:
	case fa_ammoniumsulphate:
	case fa_boron:
	case fa_calcium:
	case fa_cu:
	case fa_greenmanure:
	case fa_k:
	case fa_manganesesulphate:
	case fa_manure:
	case fa_n:
	case fa_nk:
	case fa_npk:
	case fa_npks:
	case fa_p:
	case fa_pk:
	case fa_pks:
	case fa_rsm:
	case fa_sk:
	case fa_sludge:
	case fa_slurry:
	case fp_ammoniumsulphate:
	case fp_boron:
	case fp_calcium:
	case fp_cu:
	case fp_greenmanure:
	case fp_k:
	case fp_liquidNH3:
	case fp_manganesesulphate:
	case fp_manure:
	case fp_n:
	case fp_nc:
	case fp_nk:
	case fp_npk:
	case fp_npks:
	case fp_ns:
	case fp_p:
	case fp_pk:
	case fp_pks:
	case fp_rsm:
	case fp_sk:
	case fp_sludge:
	case fp_slurry:
	case fungicide_treat:
	case glyphosate:
	case molluscicide:
	case org_fungicide:
	case org_herbicide:
	case org_insecticide:
	case last_treatment:
	case pheromone:
	case trial_control:
	case trial_insecticidetreat:
	case trial_toxiccontrol:
	case product_treat:
	case cattle_out:
	case cattle_out_low:
	case water:
	case manual_weeding:
	case pruning:
	case suckering:
	case growth_regulator:
		break;
	default:
		g_msg->Warn(WARN_FILE, "Vole_JuvenileMale::OnFarmEvent(): Unknown event type:",
		            m_OurLandscape->EventtypeToString(event));
		exit(1);
	}
	if (CurrentVState == tovs_MDying)
	{
		if (cfg_RecordVoleMort.value()) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_MFarm, 1);
		return true;
	}
	return false;
}

//----------------------------------------------------------------------------

/**
\brief
Map location function
*/
inline void Vole_JuvenileMale::SetLocation() { m_OurPopulation->m_VoleMap->SetMapValue(m_Location_x, m_Location_y, this); };
//----------------------------------------------------------------------------

/**
\brief
Map location function
*/
inline void Vole_JuvenileMale::FreeLocation() { m_OurPopulation->m_VoleMap->ClearMapValue(m_Location_x, m_Location_y); };
//----------------------------------------------------------------------------

/**
\brief
Map location function
*/
inline bool Vole_JuvenileMale::GetLocation(int px, int py) {
	if (m_OurPopulation->m_VoleMap->GetMapValue(px, py) != nullptr) return true;
	return false;
};
//---------------------------------------------------------------------------

/**
\brief
Juvenile Male vole BeginStep
*/
/**
The BeginStep is one of the three timestep divisions. This is called once for each vole before Step and EndStep. \n
The main function here is to remove voles that die before they take up CPU resources in the Step code.
*/
void Vole_JuvenileMale::BeginStep() {
	if (MortalityTest())
	{
		if (cfg_RecordVoleMort.value()) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_MBck, 1);
		CurrentVState = tovs_MDying;
		m_StepDone = true;
	}
	else
	{
		if (--m_LifeSpan < 1)
		{
			if (cfg_RecordVoleMort.value()) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_MLife, 1);
			CurrentVState = tovs_MDying;
			m_StepDone = true;
		}
	}
#ifdef __VOLEPESTICIDEON
	PesticideIngestion();
	ActOnPesticideDose();
#endif
}

//---------------------------------------------------------------------------

/**
\brief
Juvenile Male vole Step
*/
/**
The Step is one of the three timestep divisions. This is called repeatedly after BeginStep and before EndStep, until all voles report that they are done with Step. \n
\n Most of the behaviours are controlled by moving voles between behavioural states in Step (for other models this is also done in BeginStep and EndStep).
\n When a vole is done for the day it will signal this by setting m_StepDone==true. NB that a call to one behaviour may trigger a call to another behaviour on the next call to step inside the same timestep. In this way a daily cycle of activity can be undertaken.
*/
void Vole_JuvenileMale::Step() {
	if (m_StepDone || m_CurrentStateNo == -1) return;
	switch (CurrentVState)
	{
	case tovs_InitialState: //initial state
		CurrentVState = tovs_JuvenileExploration;
		break;
	case tovs_JuvenileExploration: // Juvenile Exploration
		st_JuvenileExplore();
		CurrentVState = tovs_MEvaluateExplore;
		m_StepDone = true;
		break;
	case tovs_MEvaluateExplore: // Eval&Explore
		CurrentVState = st_Eval_n_Explore();
		m_StepDone = true;
		break;
	case tovs_MMaturation:
		st_BecomeSubAdult();
		break;
	case tovs_MDying:
		FreeLocation();
		m_StepDone = true;
		break;
	default:
		m_OurLandscape->Warn("Vole_JuvenileMale::Step - unknown return error", nullptr);
		exit(1);
	}
}

//---------------------------------------------------------------------------

/**
\brief
Juvenile Male vole EndStep
*/
/**
The EndStep is one of the three timestep divisions. This is called once for each vole after BeginStep and Step. \n
The main function here is to remove voles that have died during step and otherwise to grow if not at max weight.
It also checks if the vole was killed due to human management and determines the potential territory size.
*/
void Vole_JuvenileMale::EndStep() {
	CheckManagement();
	if (CurrentVState == tovs_MDying) { st_Dying(); }
	else
	{
		m_Age++;
		m_Weight += growthperdayM;
		if (m_Age >= MinReproAgeM) CurrentVState = tovs_MMaturation;
		if (m_OurLandscape->SupplyDayInYear() == 1)
		{
			m_BornLastYear = true; // must be true if alive on 1st Jan.
		}
	}
}

//---------------------------------------------------------------------------

/**
\brief
JuvenileMale vole main territory assessment behaviour.
*/
/**
Evaluates the quality of his habitat and does some limited exploration in the surrounding area to see if
he can improve it by moving. \n
Returns 1 if it has died of extra dispersal mortality or starvation
*/
TTypeOfVoleState Vole_JuvenileMale::st_Eval_n_Explore(void) {
#ifdef __USINGTRAPLINES
	CheckTraps();
#endif
	const double Qual = CalculateCarryingCapacity(m_Location_x, m_Location_y, cfg_VoleDDepConst.value());
	const int NoFemales = m_OurPopulation->SupplyCountFemales(m_Location_x, m_Location_y, m_TerrRange);
	bool NoMates = false;
	if (NoFemales < 1) if (g_rand_uni_fnc() < g_NoFemalesMove) NoMates = true;
	//if (NoFemales<1) NoMates = true;
	const int OldMales = m_OurPopulation->SupplyInOlderTerr(m_Location_x, m_Location_y, m_Age, m_TerrRange); // returns -1 if a male has p_x,p_y in his territory and is older than p_Age else returns the number of females present
	if (Qual < m_MinJMVoleHabQual || OldMales < 1 || NoMates)
	//if ((Qual < m_MinJMVoleHabQual) || (OldMales < 1))
	{
		m_DispVector = g_random_fnc(8);
		// Not quality dependent dispersal, very directed
		return Dispersal(-1, g_random_fnc(MaleMovement[0]) + 1);
	}
	return Dispersal(Qual, 1 + g_random_fnc(MaleMovement[0]) + 1);
}

//---------------------------------------------------------------------------

/**
\brief
JuvenileMale vole dispersal behaviour.
*/
TTypeOfVoleState Vole_JuvenileMale::Dispersal(double p_OldQual, int p_Distance) {
	if (g_rand_uni_fnc() < g_extradispmort)
	{
		if (cfg_RecordVoleMort.value()) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_MDisp, 1);
		return tovs_MDying;
	}
	// Returns 1 for die, 0 for carry on
	if (m_DispVector == -1) m_DispVector = g_random_fnc(8); // Choose direction 0-7
	// Go that far in that direction (assuming it is possible to do that)
	const int oldx = m_Location_x;
	const int oldy = m_Location_y;
	MoveTo(m_DispVector, p_Distance, 10);
	//  Now we are there so what is the new quality
	// 1. Get the carrying capacity
	const double CC = CalculateCarryingCapacity(m_Location_x, m_Location_y, cfg_VoleDDepConst.value());
	if (CC <= p_OldQual)
	{
		// Don't want to move
		FreeLocation();
		m_Location_x = oldx;
		m_Location_y = oldy;
		m_OurPopulation->UpdateGuardMap(m_Location_x, m_Location_y, m_guard_cell_x, m_guard_cell_y);
		SetLocation();
	}
	return tovs_MEvaluateExplore;
}

//---------------------------------------------------------------------------
void Vole_JuvenileMale::st_BecomeSubAdult() {
	//cout << "\n I am in Vole_JuvenileMale::stBecomeSubAdult() and the locations are: ";
	//cout << m_Location_x << " , " << m_Location_y;
	const auto av = new struct_Vole_Adult;
	av->VPM = m_OurPopulation;
	av->L = m_OurLandscape;
	if (m_MyGenes.GetDirectFlag() == 0) av->m_dflag = false;
	else av->m_flag = true;
	if (m_MyGenes.GetGeneticFlag() == 0) av->m_gflag = false;
	else av->m_gflag = true;
	av->m_flag = true; // Used to signal pesticide effect to CreateObjects
	// Create the new voles (50% chance of male/female)
	av->PolyRefBorn = m_OurLandscape->SupplyPolyRef(m_Location_x, m_Location_y);
	av->ElemBorn = m_OurLandscape->SupplyElementType(m_Location_x, m_Location_y);
	av->VegBorn = m_OurLandscape->SupplyVegType(m_Location_x, m_Location_y);
	av->FatherId = m_FatherId;
	av->MotherId = m_MotherId;
	av->x = m_Location_x;
	av->y = m_Location_y;
	av->weight = m_Weight;
	// Do the genetics
	av->Genes = m_MyGenes;
	av->age = m_Age;
	av->new_Genes=new_Genes;
	av->new_Mates_Genes=new_Mates_Genes;
	//av->new_Genes.MitochondrialLine=new_Genes.MitochondrialLine;
	//av->new_Genes.YChromoLine=new_Genes.YChromoLine;
	av->xborn=xborn;
	av->yborn=yborn;
	av-> GenerationCount=GenerationCount;
	m_OurPopulation->CreateObjects(vob_Male, this, av, 1, filename_AlleleInput,true);
	// Remove the current object
	CurrentVState = tovs_MDying;
	FreeLocation();
	m_CurrentStateNo = -1;
	m_StepDone = true;
}

//---------------------------------------------------------------------------

/**
\brief
JuvenileMale vole death by external entity.
*/
/**
Response to external death event - most likely eaten by a explicitly modelled predator
*/
void Vole_JuvenileMale::OnKilled() {
	if (cfg_RecordVoleMort.value()) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_MPred, 1);
	CurrentVState = tovs_MDying;
}

//---------------------------------------------------------------------------

/**
\brief
Vole_JuvenileMale constructor.
*/
Vole_JuvenileMale::Vole_JuvenileMale(struct_Vole_Adult* p_aVoleStruct, const string& filename_AlleleInput, int i): Vole_Base(p_aVoleStruct, filename_AlleleInput, i) {
	//cout << "\n I am now in hte Vole_JuvenileMale constructor";
	//cout << "\n these are the genes after going to Juvenile male constructor";
	//p_aVoleStruct->new_Genes.PrintGenomes();
	//p_aVoleStruct->new_Mates_Genes.PrintGenomes();
	//cout << "\n I am now in hte jevenile Vole_Male  constructor and will print positions: ";
	//cout << p_aVoleStruct->x,p_aVoleStruct->y;
	m_Age = 14;
	m_Weight = WeanedWeight;
	m_Sex = true;
	m_TerrRange = m_MinMaleTerritorySize;
	m_MinTerrRange = m_MinMaleTerritorySize;
}

//---------------------------------------------------------------------------

void Vole_JuvenileMale::ReInit(struct_Vole_Adult* p_aVoleStruct, const string& filename_AlleleInput, int i) {
	Init(p_aVoleStruct, filename_AlleleInput, i);
	m_Age = 14;
	m_Weight = WeanedWeight;
	m_Sex = true;
	m_TerrRange = m_MinMaleTerritorySize;
	m_MinTerrRange = m_MinMaleTerritorySize;
}

//---------------------------------------------------------------------------

Vole_JuvenileMale::~Vole_JuvenileMale() {
	// Nothing to do
}

//---------------------------------------------------------------------------
#ifdef __VOLERODENTICIDEON
//#####################################################################################
//############################# RODENTICIDE HANDLING CODE ###############################
//#####################################################################################

void Vole_Base::RodenticideIngestion(void) {
	const double rodenticide = m_OurLandscape->SupplyRodenticide(m_Location_x, m_Location_y);
	if (rodenticide > 0)
	{
		// Todo effects
	}
	else return;
}
#endif

#ifdef __VOLEPESTICIDEON
//#####################################################################################
//############################# PESTICIDE HANDLING CODE ###############################
//#####################################################################################

void Vole_Base::PesticideIngestion( void )
{
	// If we are not running the pesticide module, then do nothing
	if ( !l_pest_enable_pesticide_engine.value()) return;
	BioDegradePesticide(); /** The extent of biodegredation is set by the config variable CfgFloat cfg_volepcidebiodegredrate - "VOLE_PCIDE_BIODEGREDATIONRATE" in the Config file. */
	/**
	Ingestion Multiplyer 1.39 \n
	Wt 0.025g					\n
	Ingested weight = 0.03475 kg \n
	...but we need this as dry weight, so multiply by 0.3 = 0.010425 \n
	Hence the vole consumes 0.010425 kg
	\n
	This value is used determine the inputs in mg, then we have have an option of assuming a fixed residue per unit
	vegetation, or we can alter the amount depending on the vegetation biomass currently present.
	*/
	double pest=0;
	pest =m_OurLandscape->SupplyPesticide( m_Location_x,m_Location_y, ppp_1); // ppp_1 is assumed. Units are mg/kg vegetation or mg/m2
	if (!cfg_pest_residue_or_rate.value())
	{
		/*
		// pest is is mg/m2 (corrected for wetweight)
		// veg is in g/m2
		// we need this in mg/kg body weight
		// mg/g veg = pest/veg
		// ingested mg (X) = pest/veg * 0.010425 * 1000
		// mg/kg = X/0.025  (vole weight in kg)
		pest*=g_SpeedyDivides[veg]; // assuming 70% is water
		*/
		double veg= 1 + m_OurLandscape->SupplyVegBiomass( m_Location_x,m_Location_y );  // Units are g dw/m2, 1 added to avoid divide by zero
		veg *= 3.33333; // assuming 70% water
		pest *= g_SpeedyDivides[(int)veg]; // mg/kg ww
	}
	//pest *= 0.010425; // mg ingested
	//pest *= 40; // = 1/0.025 // mg/kg bw ingested
	m_pesticideload += pest; // add this to the current body burden
}

void Vole_Base::ActOnPesticideDose()
{
	TTypesOfPesticide tp = m_OurLandscape->SupplyPesticideType();
	if ( tp == ttop_ModelinkPesticide21TWA ) ModelinkPesticide21TWA(m_pesticideload);
	else if (m_pesticideload > cfg_PesticideAccumulationThreshold.value()) /** If the body burden exceeds the trigger then call pesticide-specific actions by dose */
	{
		switch (tp)
		{
		case ttop_NoPesticide:
			break;
		case ttop_ReproductiveEffects:
			GeneralEndocrineDisruptor(m_pesticideload);
			break;
		case ttop_AcuteEffects:
			GeneralOrganophosphate(m_pesticideload); /** Calls the GeneralOrganophosphate action code */
			break;
		case ttop_Vinclozolin:
			Vinclozolin(m_pesticideload);
			break;
		case ttop_ModelinkPesticide:
			ModelinkPesticide();
			break;
		case ttop_GeneticDemo:
			GeneticDemoPesticide(m_pesticideload); /** Calls the GeneticDemoPesticide action code */
			break;
        default:
            exit(47);
		}
	}
}

//*************************************************************************************
//************************** PESTICIDE SPECIFIC EFFECT CODE ***************************
//*************************************************************************************

//-------------------------------------------------------------------------------------
//------------------------------ MODELINKPESTICIDE EFFECT CODE ------------------------
//-------------------------------------------------------------------------------------
void Vole_JuvenileMale::ModelinkPesticide()
{
	; // No impact
}

void Vole_JuvenileMale::ModelinkPesticide21TWA(double /* a_dose */)
{
	; // No impact
}

void Vole_JuvenileFemale::ModelinkPesticide()
{
	; // No impact
}

void Vole_JuvenileFemale::ModelinkPesticide21TWA(double /* a_dose */)
{
	; // No impact
}

void Vole_Male::ModelinkPesticide()
{
	; // No impact
}

void Vole_Male::ModelinkPesticide21TWA(double /* a_dose */)
{
	; // No impact
}

void Vole_Female::ModelinkPesticide()
{
	// Specify certain gestation days for the effects here
	if (((m_DaysUntilBirth>0) && (m_DaysUntilBirth<m_EndoCrineDisruptionGestationLength)) || (CurrentVState==tovs_Lactating) )
	{
		if (m_pesticideload > m_pesticideInfluenced3)
		{
			m_pesticideInfluenced3=m_pesticideload;
		}
	}
}

void Vole_Female::ModelinkPesticide21TWA(double /* a_dose */)
{
	m_pesticideloadarray[m_pesticideloadindex++]=m_pesticideload;
	if (m_pesticideloadindex==21) m_pesticideloadindex=0;
}

//-------------------------------------------------------------------------------------
//------------------------------ VINCLOZALIN EFFECT CODE ------------------------------
//-------------------------------------------------------------------------------------
void Vole_JuvenileMale::Vinclozolin(double /* a_dose */ )
{
	; // No impact
}

void Vole_JuvenileFemale::Vinclozolin(double /* a_dose */ )
{
	; // No impact
}

void Vole_Male::Vinclozolin(double /* a_dose */ )
{
	; // No impact
}

void Vole_Female::Vinclozolin(double /* a_dose */ )
{
	// May also wish to specify certain gestation days for the effects here
	if ((m_DaysUntilBirth>0) && (m_DaysUntilBirth<m_EndoCrineDisruptionGestationLength))
	{
		m_pesticideInfluenced1 = true;
	}
}

//-------------------------------------------------------------------------------------
//-------------------- GENERAL ENDOCRINE DISRUPTOR EFFECT CODE ------------------------
//-------------------------------------------------------------------------------------

void Vole_Female::GeneralEndocrineDisruptor(double /* a_dose */ )
{
	// May also wish to specify certain gestation days for the effects here
	if ((m_DaysUntilBirth>0) && (m_DaysUntilBirth< m_EndoCrineDisruptionGestationLength))
	{
		m_pesticideInfluenced1 = true;
	}
}

void Vole_Male::GeneralEndocrineDisruptor(double /* a_dose */ )
{
	; // No effect
}

void Vole_JuvenileMale::GeneralEndocrineDisruptor(double /* a_dose */ )
{
	; // No effect
}

void Vole_JuvenileFemale::GeneralEndocrineDisruptor(double /* a_dose */ )
{
	; // No effect
}



//-------------------------------------------------------------------------------------
//------------------------ GENERAL ORGANOPHOSPHATE EFFECT CODE ------------------------
//-------------------------------------------------------------------------------------
void Vole_JuvenileMale::GeneralOrganophosphate(double /* a_dose */)
{
	if ( g_rand_uni() > l_pest_daily_mort.value() ) return;
    CurrentVState=tovs_MDying;
    m_StepDone=true;
	if (cfg_RecordVoleMort.value() ) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_MPest,1);
	return;
}

void Vole_JuvenileFemale::GeneralOrganophosphate(double /* a_dose */)
{
	if ( g_rand_uni() > l_pest_daily_mort.value() ) return;
    CurrentVState=tovs_FDying;
    m_StepDone=true;
	if (cfg_RecordVoleMort.value() ) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_FPest,1);
	return;
}

void Vole_Male::GeneralOrganophosphate(double /* a_dose */)
{
	if ( g_rand_uni() > l_pest_daily_mort.value() ) return;
    CurrentVState=tovs_MDying;
    m_StepDone=true;
	if (cfg_RecordVoleMort.value() ) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_MPest,1);
	return;
}

void Vole_Female::GeneralOrganophosphate(double /* a_dose */)
{
	if ( g_rand_uni() > l_pest_daily_mort.value() ) return;
	CurrentVState=tovs_FDying;
    m_StepDone=true;
	if (cfg_RecordVoleMort.value() ) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_FPest,1);
	return;
}

void Vole_JuvenileMale::GeneticDemoPesticide(double /* a_dose */)
{
	if ( g_rand_uni() > l_pest_daily_mort.value() ) return;
	if ( cfg_ResistanceDominant.value()){ if ((SupplyMyAllele(3,0) == 1) || (SupplyMyAllele(3,1) == 1)) return; }
	  else

      {
          if ((SupplyMyAllele(3,0) == 1) && (SupplyMyAllele(3,1) == 1))
	{
			return;
	}
      }
	CurrentVState=tovs_MDying;
    m_StepDone=true;
	if (cfg_RecordVoleMort.value() ) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_MPest,1);
	return;
}

void Vole_JuvenileFemale::GeneticDemoPesticide(double /* a_dose */)
{
	if ( g_rand_uni() > l_pest_daily_mort.value() ) return;
	if ( cfg_ResistanceDominant.value()) { if ((SupplyMyAllele(3,0) == 1) || (SupplyMyAllele(3,1) == 1)) return; }
	  else if ((SupplyMyAllele(3,0) == 1) && (SupplyMyAllele(3,1) == 1))
	{
			return;
	}
	CurrentVState=tovs_FDying;
    m_StepDone=true;
	if (cfg_RecordVoleMort.value() ) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_FPest,1);
	return;
}

void Vole_Male::GeneticDemoPesticide(double /* a_dose */)
{
	if ( g_rand_uni() > l_pest_daily_mort.value() ) return; // No effect triggered due probability
	if ( cfg_ResistanceDominant.value()) { if ((SupplyMyAllele(3,0) == 1) || (SupplyMyAllele(3,1) == 1)) return; } // If we assume a dominant resistant allele then if either are set there is no effect
	  else if ((SupplyMyAllele(3,0) == 1) && (SupplyMyAllele(3,1) == 1))
	{
			return;
	}
	CurrentVState=tovs_MDying; // No save due to resistance so the vole dies
    m_StepDone=true;
	if (cfg_RecordVoleMort.value() ) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_MPest,1);
	return;
}

void Vole_Female::GeneticDemoPesticide(double /* a_dose */)
{
	if ( g_rand_uni() > l_pest_daily_mort.value() ) return; // No effect triggered due probability
	if ( cfg_ResistanceDominant.value()) { if ((SupplyMyAllele(3,0) == 1) || (SupplyMyAllele(3,1) == 1)) return; } // If we assume a dominant resistant allele then if either are set there is no effect
	  else if ((SupplyMyAllele(3,0) == 1) && (SupplyMyAllele(3,1) == 1))
	{
			return;
	}
	CurrentVState=tovs_FDying; // No save due to resistance so the vole dies
    m_StepDone=true;
	if (cfg_RecordVoleMort.value() ) m_OurPopulation->m_VoleRecordMort->ChangeData(tovmort_FPest,1);
	return;
}

//#####################################################################################
//########################### END PESTICIDE HANDLING CODE #############################
//#####################################################################################

#endif
