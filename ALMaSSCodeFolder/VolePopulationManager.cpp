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
<B>VolePopulationManager.cpp This file contains the source for the vole population manager class</B> \n
*/
/**
\file
 by Chris J. Topping \n
 Version of 28th Jan 2001 \n
 \n
 With additions as noted in: \n
 April 2006 \n
 November 2007 \n
 Doxygen formatted comments in May 2008 \n
*/
//---------------------------------------------------------------------------

#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <random>
#include "../Landscape/ls.h"
#include "../BatchALMaSS/PopulationManager.h"
#include "../BatchALMaSS/AOR_Probe.h"
#include "../BatchALMaSS/BinaryMapBase.h"
#include "../BatchALMaSS/MovementMap.h"
#include "../Vole/vole_all.h"
#include "../Vole/VolePopulationManager.h"
#include "../Vole/GeneticMaterial.h"
#include "../Vole/Vole_toletoc.h"

extern double MutationChance;

//---------------------------------------------------------------------------

/* External configuration variables related to outputs */
/** The first year for catastrophes */
extern CfgInt cfg_CatastropheEventStartYear;
/** The frequency of catastrophes */
extern CfgInt cfg_PmEventfrequency;
/** The size of catastrophes */
extern CfgInt cfg_PmEventsize;
/** Raw data for spatial output for calculating Ripley's spatial statistics */
extern CfgBool cfg_ReallyBigOutputMonthly_used;
extern CfgBool cfg_CfgRipleysOutputUsed;
/** Raw data for spatial output for calculating Ripley's spatial statistics with additional information about habitat and physiology */
extern CfgBool cfg_ReallyBigOutputUsed;
/** Extra mortality for dispersing individuals*/
extern CfgFloat cfg_extradispmort;
extern double g_extradispmort;
/** risk of daily mortality*/
extern double g_DailyMortChance;
extern double g_DailyMortChanceMaleTerr;
extern double g_NoFemalesMove;
extern CfgInt cfg_VoleDDepConst;
extern int g_sigAgeDiff;
extern int MinReproAgeM;
extern int MinReproAgeF;
extern CfgInt cfg_MinReproAgeM;
extern CfgInt cfg_MinReproAgeF;
extern unsigned char g_MaxAllele;

//---------------------------------------------------------------------------

// End of breeding season
static CfgInt cfg_vole_reprofinishdate("VOLE_REPROFINISHDATE", CFG_CUSTOM, 243); // September 1st
// VoleSampleFile for the vole
static CfgBool cfg_VoleSampleDataUsed("VOLE_SAMPLE_FILE_USED", CFG_CUSTOM, false);
//First sampling day for vole samplefile
static CfgInt cfg_VoleSampleFileDayOne("VOLE_SAMPLE_FILE_DAY_ONE", CFG_CUSTOM, 109);
//Second sampling day for vole samplefile
static CfgInt cfg_VoleSampleFileDayTwo("VOLE_SAMPLE_FILE_DAY_TWO", CFG_CUSTOM, 287);
//Whether a vole will move if no females present
static CfgFloat cfg_VoleFemaleMove("VOLE_FEMALEMOVE", CFG_CUSTOM, 0.0505);
// Radius whithin which the female picks a mate
static CfgInt cfg_MateRadius("FEMALE_MATE_RADIUS", CFG_CUSTOM, 250);
/** The unspecified background mortality rate */
static CfgFloat cfg_VoleBackgroundMort("VOLE_BCKMORT", CFG_CUSTOM, 0.0025); //changed from 0.0025 for testing purposes
/** The probability of moving to less favourable habitats when walking */
static CfgFloat cfg_MoveToLessFavourable("VOLEMOVETOLESSFAVOURALBLE", CFG_CUSTOM, 0.005);
/** The age difference before an old male can oust a younger one */
static CfgInt cfg_sigAgeDiff("VOLE_SIGAGEDIFF", CFG_CUSTOM, 30);
/** Toxicological variable for specific pesticide effects */
static CfgInt cfg_geneticproductfertilityeffect("VOLE_GENETICPRODUCTFERTILITYEFFECT", CFG_CUSTOM, 50);
/** Toxicological variable for specific pesticide effects */
static CfgFloat cfg_geneticsterilitychance("VOLE_GENETICSTERILITYCHANCE", CFG_CUSTOM, 0);
/** Toxicological variable for specific pesticide effects */
static CfgFloat cfg_f1sterilitychance("VOLE_FONESTERILITYCHANCE", CFG_CUSTOM, 0.0);
/** Record genetic data or not */
static CfgBool cfg_genetic_output("VOLE_GENETIC_OUTPUT", CFG_CUSTOM, true);
/** The number of each sex to sample for genetics */
static CfgInt cfg_genetetic_output_sample_size("VOLE_GENETIC_SAMPLE_SIZE", CFG_CUSTOM, 500);
// GeneticOutput variable for the first year the file should be produced
static CfgInt cfg_GeneticsResultOutputFirstYear("GENETICS_RESULT_OUTPUT_FIRST_YEAR", CFG_CUSTOM, 1); //default changed from 100000
// GeneticOutput variable for the start year in the first period
static CfgInt cfg_GeneticResultOutputSecondYear("GENETICS_RESULT_OUTPUT_SECOND_YEAR", CFG_CUSTOM, 100000);
// GeneticOutput variable for the start year in the second period
CfgBool cfg_GeneticsResultOutputUsed("GENETICS_RESULT_OUTPUT_USED", CFG_CUSTOM, true); //default changed from false
static CfgInt cfg_GeneticsResultOutputInterval_1("GENETICS_RESULT_OUTPUT_INTERVAL_ONE", CFG_CUSTOM, 1); //defualt changed from 100000
// GeneticOutput variable for the first interval the file should be produced
static CfgInt cfg_GeneticsResultOutputDay_1("GENETICS_RESULT_OUTPUT_DAY_ONE", CFG_CUSTOM, 30); //changed from 3 
// GeneticOutput variable for the second day of the year the file should be produced
static CfgInt cfg_GeneticsResultOutputDay_2("GENETICS_RESULT_OUTPUT_DAY_TWO", CFG_CUSTOM, 6);
// GeneticOutput variable for the when the first interval should end

/** Flag turning resistance output on or off */
static CfgBool cfg_VoleUseResistanceOuput("VOLE_RESISTANCEOUTPUT_ON_OFF", CFG_CUSTOM, false);
/** The day in the year to output the resitance gene frequency */
static CfgInt cfg_VoleResistanceDay("VOLE_RESISTANCEOUTPUT_DAY", CFG_CUSTOM, 365);
/** The starting frequency of '1' in chromosome 0 locus no 3 */
static CfgFloat cfg_ResistanceStartFrequency("VOLE_RESISTANCESTARTFREQ", CFG_CUSTOM, 0.01);

static CfgBool cfg_UseVoleTraplines("VOLE_TRAPLINES_ON_OFF", CFG_CUSTOM, false);
static CfgStr cfg_VoleTraplinesfile("VOLE_TRAPLINESFILE", CFG_CUSTOM, "VoleTraplines.txt");
static CfgInt cfg_VoleTrapResolution("VOLE_TRAPRESOLUTION", CFG_CUSTOM, 2);

CfgBool cfg_RecordVoleMort("VOLE_RECORDMORT_ON", CFG_CUSTOM, false);
static CfgStr cfg_VoleRecordMortFile("VOLE_RECORDMORTFILE", CFG_CUSTOM, "VoleMortalityCauses.txt");

CfgBool cfg_SexRatiosOutput_used("VOLE_SEXRATIOSOUTPUT_USED", CFG_CUSTOM, false);
static CfgStr cfg_SexRatiosOutput_filename("VOLE_SEXRATIOS_FILENAME", CFG_CUSTOM, "VoleSexRatios.txt");
static CfgInt cfg_SexRatiosOutput_interval("VOLE_SEXRATIOSOUTPUT_INTERVAL", CFG_CUSTOM, 1);
static CfgInt cfg_SexRatiosOutput_day("VOLE_SEXRATIOSOUTPUT_DAY", CFG_CUSTOM, -1); // -1 means daily
static CfgInt cfg_SexRatiosOutputFirstYear("VOLE_SEXRATIOSOUTPUT_FIRSTYEAR", CFG_CUSTOM, 1);
static CfgBool cfg_voleLandscapeGridOutputUsed("VOLE_LANDSCAPEGRIDOUTPUTUSED", CFG_CUSTOM, true);
static CfgInt cfg_voleLandscapeGridOutputDay("VOLE_LANDSCAPEGRIDOUTPUTDAY", CFG_CUSTOM, 1); //changed from 10000
static CfgBool cfg_volestartinoptimalonly("VOLE_STARTINOPTIMALONLY", CFG_CUSTOM, false);

/** This is specified as a global because it avoids costing time on look-ups. \n
 This is needed because of the extreme CPU intensive movement functions
*/
double MoveToLessFavourable;

/** How many male voles to start with */
static CfgInt cfg_vole_starting_numberM("VOLE_START_NO_M", CFG_CUSTOM, 10000);
/** How many female voles to start with */
static CfgInt cfg_vole_starting_numberF("VOLE_START_NO_F", CFG_CUSTOM, 10000);
/** The mean temperature over 14 days at which the grass is assummed to start to grow */
static CfgFloat cfg_GrassStartGrowth("VOLE_GRASSSTARTGROWTH", CFG_CUSTOM, 3.552);
/** The date before which repro is impossible */
static CfgFloat cfg_GrassStartGrowthDay("VOLE_GRASSSTARTGROWTHDAY", CFG_CUSTOM, 80);
/** Special output file toggle */
static CfgBool cfg_useagesexlocation("VOLE_USEAGESEXLOCATION", CFG_CUSTOM, false);

/** Litter sizes dependent upon month and year */
const float ReproTable_base[2][12] = {
	{0, 0, 0, 0, 4, 4, 4, 4, 3, 2, 3, 0}, // Born This year
	{0, 0, 4, 4, 5, 5, 5, 5, 4, 4, 0, 0} // Born Last year
}; // data from Andrea (1981)

/** Mutation probability */
static CfgFloat cfg_MutationChance("GENETICS_MUTATION_CHANCE", CFG_CUSTOM, 0.001);
/** The maximum allele value used */
static CfgInt cfg_MaxAllele("GENETICS_MAXALLELE", CFG_CUSTOM, 255);

extern int g_MaleReproductFinish;
extern CfgInt cfg_MinFemaleTerritorySize;
extern CfgInt cfg_MinMaleTerritorySize;
extern CfgBool cfg_ResistanceDominant;

static CfgStr params_for_beta0_cfg("PARAMS_FOR_BETA0", CFG_CUSTOM, "0.5 0.5");
static CfgStr params_for_beta1_cfg("PARAMS_FOR_BETA1", CFG_CUSTOM, "5.0 5.0");
static CfgStr landscape_info_cfg("LANDSCAPE_INFO", CFG_CUSTOM,"lallalalla");
static CfgStr QuadrantsDesired_cfg("QUADRANTS_DESIRED",CFG_CUSTOM,"four");
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//                      Vole_Population_Manager
//---------------------------------------------------------------------------

/**
Class for managing lists of all voles in the simulation and responsible for generating output.
Main functionality is the the base class Population_Manager
*/
Vole_Population_Manager::Vole_Population_Manager(Landscape* p_L) : Population_Manager(p_L, vob_foobar) {
	IDNumber = 0;
	thisYear = 0;
	YearsTotal = 0;
	m_VoleMap = new IDMap<TAnimal*>(p_L);
	AFreq = new AlleleFreq();
	m_GrowthStartDate = 9000; // delays this until we know the real one
	Init();
}
//---------------------------------------------------------------------------

/**
Vole_Population_Manager destructor \n
Only needed to close any potentially open output files
*/
Vole_Population_Manager::~Vole_Population_Manager() {
	// Save the impacted file
	//errno_t errNum;
	// Output files for the vole
	//fclose( m_GeneticsFile );
	//fclose( m_AlleleFreqsFile );
	//fclose( m_EasyPopRes );
	if (cfg_SexRatiosOutput_used.value()) CloseSexRatiosProbe();
	delete m_VoleMap;
	delete AFreq;
	// Vole traplines tidy-up
	if (cfg_UseVoleTraplines.value()) { delete m_Traplines; }
	if (cfg_RecordVoleMort.value())
	{
		m_VoleRecordMort->OPrint();
		delete m_VoleRecordMort;
	}
	if (cfg_useagesexlocation.value()) { (*m_VoleAgeSexLocationFile).close(); }
	if (cfg_VoleUseResistanceOuput.value()) CloseResistanceOutput();
}
//---------------------------------------------------------------------------

void Vole_Population_Manager::Init() {
	/** autom. called by constructor. \n
	Sets up output files and creates initial population
  
  */
	// Need to resize out vole habitat qualities vector to the correct size for the number of polygons
	// NB IF the landscape were to change the number of polygons due to morphing (not implemented currently), then
	// this code is broken.

	MakeAlleleInputFile(params_for_beta0_cfg.value(),params_for_beta1_cfg.value()); //make the file of allele inputs from the specified parameters
	const unsigned int p = m_TheLandscape->SupplyNumberOfPolygons();
	VoleHabitatBaseQualities.resize(p);
	m_is_paralleled = true;
	//// Initialise the vole grazing data
	//m_TheLandscape->ResetAllVoleGrazing();
	// Set up output files
	if (cfg_voleLandscapeGridOutputUsed.value())
	{
		ofstream ofile("VoleLandscapeGridData.txt", ios::out);
		ofile.close();
	}
	if (cfg_CfgRipleysOutputUsed.value()) { OpenTheRipleysOutputProbe(); }
	if (cfg_ReallyBigOutputUsed.value()) { OpenTheReallyBigProbe(); }
	else ReallyBigOutputPrb = nullptr;
	if (cfg_SexRatiosOutput_used.value()) OpenSexRatiosProbe();
	const FILE* GFile = fopen("VoleGeneticsOutput1.txt", "w");
	if (!GFile)
	{
		g_msg->Warn(WARN_FILE, "Population_Manager::Init(): ""Unable to open VoleGeneticsOutput1.txt", "");
		exit(1);
	}
	if (cfg_useagesexlocation.value()) { m_VoleAgeSexLocationFile = new ofstream("VoleAgeSexLocation.txt", ios::out); }
	// Pass some parameters to member variables (speed optimisation)
	g_MaleReproductFinish = cfg_vole_reprofinishdate.value();
	g_extradispmort = cfg_extradispmort.value();
	g_DailyMortChance = cfg_VoleBackgroundMort.value();
	g_NoFemalesMove = cfg_VoleFemaleMove.value();
	//m_BaseVoleDensity = 1/(cfg_MinFemaleTerritorySize.value() * 4.0) * cfg_VoleDDepConst.value();
	g_sigAgeDiff = cfg_sigAgeDiff.value();
	MinReproAgeM = cfg_MinReproAgeM.value();
	MinReproAgeF = cfg_MinReproAgeF.value();
#ifdef __VOLEPESTICIDEON
  /** Four memeber variables Used in pesticide simulations */
  m_impacted = 0;
  m_notimpacted = 0;
  m_geneticimpacted = 0;
  m_LittersLost = 0;
  m_LittersProduced = 0;
  FILE* impf;
  impf = fopen("Impacted.res","w");
  if (!impf) {
	  m_TheLandscape->Warn("Vole_Population_Manager Destructor","Could Not Open Impacted.Res File");
	  exit(0);
  }
  fclose(impf);
#endif
	/*
	 Create probabilities for ReproTable
  
	 ReproTable is two years each with 12 numbers representing the minimum litter size in the 0,1 position and in 2,3 positions then it is the number out of 100 that have one more than the first number\n
	 Currently this is set to no variation and litter sizes follow Andrea (1981).
	*/
	for (int i = 0; i < 12; i++)
	{
		ReproTable[0][i] = static_cast<int>(floor(ReproTable_base[0][i]));
		ReproTable[2][i] = static_cast<int>(floor(0.5 + ReproTable_base[0][i] * 100 - ReproTable[0][i] * 100));
		ReproTable[1][i] = static_cast<int>(floor(ReproTable_base[1][i]));
		ReproTable[3][i] = static_cast<int>(floor(0.5 + ReproTable_base[1][i] * 100 - ReproTable[1][i] * 100));
	}
	m_SimulationName = "Field Vole";
	// Create males and females
	const auto av = new struct_Vole_Adult;
	av->VPM = this;
	av->L = m_TheLandscape;
	MutationChance = cfg_MutationChance.value();
	g_MaxAllele = static_cast<unsigned char>(cfg_MaxAllele.value());

	GeneticMaterial AGene;
	for (int i = 0; i < cfg_vole_starting_numberM.value(); i++) //Initial creation of males
	{
		AGene.Initiation(AFreq);
		if (m_TheLandscape->SupplyPesticideType() != -1)
		// No big time penalty here, so we use this rather than the #ifdef
		{
			// If any pesticide is being used then we assume genetics is only used to control this - so make sure there is no
			// pesticide influenced genetetic code set (ie put a zero in 0 allele, of both chromosomes
			if (m_TheLandscape->SupplyPesticideType() == 5)
			{
				if (g_rand_uni_fnc() < cfg_ResistanceStartFrequency.value()) AGene.SetAllele(3, 1, 0);
				else AGene.SetAllele(3, 0, 0);
				if (g_rand_uni_fnc() < cfg_ResistanceStartFrequency.value()) AGene.SetAllele(3, 1, 1);
				else AGene.SetAllele(3, 0, 1);
			}
			AGene.SetAllele(0, 0, 0);
			AGene.SetAllele(0, 0, 1);
		}
		av->Genes = AGene;
		bool found = false;
		while (!found)
		{
			av->x = g_random_fnc(SimW);
			av->y = g_random_fnc(SimH);
			av->xborn= av->x;
			av->yborn= av->y;
			found = SuitableStartingLocation(av->x, av->y);
		}
		CreateObjectsInit(vob_Male, nullptr, av, 1, filename_AlleleInput, i); //i is passed on as because it's used as the unique paternal+maternal lineages for the new male.
	}

	for (int i = 0; i < cfg_vole_starting_numberF.value(); i++) //initial creation of females
	{
		AGene.Initiation(AFreq);
		if (m_TheLandscape->SupplyPesticideType() != -1)
		// No big time penalty here, so we use this rather than the #ifdef
		{
			// If any pesticide is being used then we assume genetics is only used to control this - so make sure there is no
			// pesticide influenced genetetic code set (ie put a zero in 0 allele, of both chromosomes
			if (m_TheLandscape->SupplyPesticideType() == 5)
			{
				if (g_rand_uni_fnc() < cfg_ResistanceStartFrequency.value()) AGene.SetAllele(3, 1, 0);
				else AGene.SetAllele(3, 0, 0);
				if (g_rand_uni_fnc() < cfg_ResistanceStartFrequency.value()) AGene.SetAllele(3, 1, 1);
				else AGene.SetAllele(3, 0, 1);
			}
			AGene.SetAllele(0, 0, 0);
			AGene.SetAllele(0, 0, 1);
		}
		av->Genes = AGene; // Their genes are set based on the Initiation (AlFreq) function
		bool found = false;
		while (!found)
		{
			av->x = g_random_fnc(SimW);
			av->y = g_random_fnc(SimH);
			av->xborn= av->x;
			av->yborn= av->y;
			found = SuitableStartingLocation(av->x, av->y);
		}
		CreateObjectsInit(vob_Female, nullptr, av, 1, filename_AlleleInput, (i+cfg_vole_starting_numberM.value())); //i+ nr of males is passed on because it's used as the unique paternal+maternal lineages for the new female.
	}
	delete av;

	// Load List of Animal Classes
	m_ListNames[0] = "Juvenile Male";
	m_ListNames[1] = "Juvenile Female";
	m_ListNames[2] = "Male";
	m_ListNames[3] = "Female";
	m_ListNameLength = vob_foobar;
	m_population_type = TOP_Vole;

	// Set up before step action sorts
	// sort w.r.t. x
	BeforeStepActions[0] = 0; // SortXIndex
	BeforeStepActions[1] = 0; // SortXIndex

	// Load State Names
	StateNames[tovs_InitialState] = "Initial State";
	//Males
	StateNames[tovs_JuvenileExploration] = "Juvenile Exploration";
	StateNames[tovs_MMaturation] = "Maturation";
	StateNames[tovs_MEvaluateExplore] = "Evaluate & Explore";
	StateNames[tovs_Infanticide] = "Infanticide";
	StateNames[tovs_MDying] = "Dying";
	//Females
	StateNames[tovs_FEvaluateExplore] = "Evaluate & Explore";
	StateNames[tovs_ReproBehaviour] = "Repro. Behaviour";
	StateNames[tovs_Lactating] = "Lactating";
	StateNames[tovs_GiveBirth] = "Give Birth";
	StateNames[tovs_FMaturation] = "Maturation";
	StateNames[tovs_Mating] = "Mating";
	StateNames[tovs_UpdateGestation] = "Update Gestation";
	StateNames[tovs_SpecialExplore] = "Special Explore";
	StateNames[tovs_FDying] = "Dying";

	MoveToLessFavourable = cfg_MoveToLessFavourable.value();
#ifdef __VOLEOUTPUTTESTFILES
	TestFile = fopen("TestFile.Txt", "w");
	if (!TestFile)
	{
		m_TheLandscape->Warn("Population_Manager::Population_Manager -", "Could Not Open TestFile.txt");
		exit(0);
	}
	TestFile2 = fopen("TestFile2.Txt", "w");
	m_StepSize = 1440; // Default - it is up to the individual animal models to
	if (!TestFile2)
	{
		m_TheLandscape->Warn("Population_Manager::Population_Manager -", "Could Not Open TestFile2.txt");
		exit(0);
	}
#endif
#ifdef __VOLEPESTICIDEON
  /** use for pesticide simulations only */
  m_geneticproductfertilityeffect = cfg_geneticproductfertilityeffect.value();
  /** use for pesticide simulations only */
  m_geneticsterilitychance = cfg_geneticsterilitychance.value();
  /** use for pesticide simulations only */
  m_f1sterilitychance = cfg_f1sterilitychance.value();
// Other output files
#endif
	char Nme[511];
	//strcpy( Nme, "GeneticsOutput.txt" );
	// m_GeneticsFile = fopen( Nme, "w" );
	//strcpy( Nme, "AlleleFreqs.txt" );
	//m_AlleleFreqsFile = fopen( Nme, "w" );
	//strcpy( Nme, "forEASYPOP.txt" );
	//m_EasyPopRes = fopen( Nme, "w" );

	// Vole traplines set-up
	if (cfg_UseVoleTraplines.value())
	{
		m_Traplines = new TrapLineMap(SimW, SimH, cfg_VoleTrapResolution.value(), cfg_VoleTraplinesfile.value());
	}
	// Vole record mortality set-up
	if (cfg_RecordVoleMort.value())
	{
		// The two ints in the constructor call refer to the number of int data and double data used respectively
		// The output will print ints first then doubles
		m_VoleRecordMort = new VoleSummaryOutput(cfg_VoleRecordMortFile.value(), m_TheLandscape, 14, 0);
		// In the string output allows us to provide titles:
		m_VoleRecordMort->OPrint(
			"Year\tDay\tMStarve\tFStarve\tMBck\tFBck\tMFarm\tFFarm\tMDisp\tFDisp\tMPred\tFPred\tMLife\tFLife\tMPest\tFPest\t");
		m_VoleRecordMort->OPrintEndl();
	}
	if (cfg_VoleUseResistanceOuput.value()) OpenResistanceOutput();
}
//---------------------------------------------------------------------------

bool Vole_Population_Manager::SuitableStartingLocation(const int x, const int y) const {
	/**
	* Two verisons of this exist, the first only lets voles start in good places, the second also allows fields and orchards.
	*/
	if (cfg_volestartinoptimalonly.value()) { return vole_tole_init_optimal(m_TheLandscape, x, y); }
	return vole_tole_init_friendly(m_TheLandscape, x, y);
}
//---------------------------------------------------------------------------

/** use for pesticide simulations only to determine the number of individual males currently impacted by pesticide
    Added April 2006
*/
void Vole_Population_Manager::ImpactedProbe() {
	FILE* MyFile = fopen("VoleImpactedProbe.txt", "a");
	if (MyFile == nullptr)
	{
		m_TheLandscape->Warn("Vole_Population_Manager::ImpactedProbe", "Could Not Open VoleImpactedProbe.txt File");
		exit(0);
	}
	// This just needs to trawl through the males and to count how many are affected by
	// a pesticide.  This function needs to be redefined for every new pesticide case
	// Current implementation is for Trine Dalkvists 'vinclozolin' work
	int dno = 0;
	int gno = 0;
	const unsigned size = GetLiveArraySize(vob_Male);
	for (unsigned j = 0; j < size; j++)
	{
		const auto vm = dynamic_cast<Vole_Male*>(SupplyAnimalPtr(vob_Male, j));
		if (vm->GetGeneticFlag() != 0) gno++;
		if (vm->GetDirectFlag() != 0) dno++;
	}
	const int tno = size;
	fprintf(MyFile, "%d\t%d\t%d\n", tno, dno, gno);
	fclose(MyFile);
}
//-------------------------------------------------------------------------------------------

void Vole_Population_Manager::DoFirst() {
	/**
	This method is called before the BeginStep of the main program loop by the base class. \n
	It controls when the grass starts to grow and therefore triggers start of reproduction but
	the primary role is as a safe location where outputs can be generated. At the time this is
	called there no chance that the voles are engaged in any kind of activity.
	*/
	/**
	* First job is to calculate a basic quality for each polygon. This is more efficient than getting the voles to do the
	* calculation each time they need to know.
	*/

	//for (int type = 0; type<vob_foobar; type++)
	//{
	//	int s = (int) TheArray[type].size();
	//	for (int i = 0; i<s; i++) m_TheLandscape->IncVoleGrazing(TheArray[type][i]->Supply_m_Location_x(), TheArray[type][i]->Supply_m_Location_y());
	//}
	//m_TheLandscape->CalcAllVoleGrazingDensity();
	unsigned sz = m_TheLandscape->SupplyNumberOfPolygons();
	for (unsigned p = 0; p < sz; p++)
	{
		double qual = AssessHabitat(p);
		//double volegrazing = m_TheLandscape->SupplyVoleGrazingDensityVector(p);
		//if (volegrazing<m_BaseVoleDensity) volegrazing = m_BaseVoleDensity;
		//VoleHabitatBaseQualities[p] = qual / (volegrazing/m_BaseVoleDensity);
		VoleHabitatBaseQualities[p] = qual;
	}

	int year = m_TheLandscape->SupplyYearNumber();
	if (year > -1)
	{
		g_DailyMortChance = cfg_VoleBackgroundMort.value();
		g_extradispmort = cfg_extradispmort.value();
	}
	else
	{
		g_DailyMortChance = 0;
		g_extradispmort = 0;
	}

	int today = m_TheLandscape->SupplyDayInYear();
	if (today == 0) { m_GrowthStartDate = 9000; }
	if (today > cfg_GrassStartGrowthDay.value() && m_GrowthStartDate == 9000)
	{
		if (m_TheLandscape->SupplyMeanTemp(m_TheLandscape->SupplyGlobalDate() - 7, 7) > cfg_GrassStartGrowth.value())
		{
			m_GrowthStartDate = today;
			printf("The start day is: %d\n", today);
		}
	}
	if (GetLiveArraySize(vob_Female) > 0)
	{
		if (today > SupplyGrowthStartDate() && today <= g_MaleReproductFinish)
		{
			dynamic_cast<Vole_Female*>(SupplyAnimalPtr(vob_Female, 0))->SetBreedingSeason(true);
			g_DailyMortChanceMaleTerr = g_DailyMortChance / 2.0;
		}
		else
		{
			dynamic_cast<Vole_Female*>(SupplyAnimalPtr(vob_Female, 0))->SetBreedingSeason(false);
			// Below could be used to alter mortality outside breeding season - currently unused (did not help POM fit)
			g_DailyMortChanceMaleTerr = g_DailyMortChance / 2.0;
		}
	}

	//** OUTPUT BELOW HERE **
	if (cfg_SexRatiosOutput_used.value())
	{
		if (year >= cfg_SexRatiosOutputFirstYear.value())
		{
			if (year % cfg_SexRatiosOutput_interval.value() == 0)
			{
				if (cfg_SexRatiosOutput_day.value() == today || cfg_SexRatiosOutput_day.value() == -1)
				{
					TheSexRatiosProbe();
				}
			}
		}
	}
	if (cfg_voleLandscapeGridOutputDay.value() == m_TheLandscape->SupplyDayInYear()) LandscapeQuadrantOutputProbe(
		m_TheLandscape->SupplyGlobalDate() - 365);
	if (cfg_voleLandscapeGridOutputDay.value() + 181 == m_TheLandscape->SupplyDayInYear()) LandscapeQuadrantOutputProbe(
		m_TheLandscape->SupplyGlobalDate() - 365);
	NoYoungKilledToday = 0;
	NoYoungKilledToday4 = 0;
	NoYoungKilledToday8 = 0;
	NoYoungKilledToday9 = 0;
	//AgeYoungKilledToday=0;
	YoungProducedToday = 0; // Zero the number of young
	JuvsProducedToday = 0; // Zero the number of juveniles
	if (today == 0)
	{
		// Annual Vole Mortality Record Operations
		if (cfg_RecordVoleMort.value()) m_VoleRecordMort->OPrint();
		if (cfg_RecordVoleMort.value()) m_VoleRecordMort->ResetData();
#ifdef __VOLEPESTICIDEON
	  // Save the impacted file
	  FILE* impf = fopen("Impacted.res","a");
	  if (!impf) {
		  m_TheLandscape->Warn("Vole_Population_Manager Destructor","Could Not Open Impacted.Res File");
		  exit(0);
	  }
	  fprintf(impf,"%d\t%d\t%d\t%d\t%d\n",m_notimpacted, m_impacted, m_geneticimpacted, m_LittersProduced, m_LittersLost);
	  fclose(impf);
	  m_impacted=0;
	  m_notimpacted=0;
	  m_geneticimpacted=0;
	  m_LittersLost=0;
	  m_LittersProduced = 0;
#endif
	}
	string QuadrantsDesired=QuadrantsDesired_cfg.value();

	if (RecordGeneticsToday(today, year, cfg_GeneticsResultOutputFirstYear.value(),
	                        cfg_GeneticsResultOutputInterval_1.value()))
	{
		if (QuadrantsDesired=="four"){
			Four_QuadrantBasedGeneticOutput(today, year, landscape_info_cfg.value());
		}else if (QuadrantsDesired=="nine"){
			Nine_QuadrantBasedGeneticOutput(today, year, landscape_info_cfg.value());
		}else if(QuadrantsDesired=="both"){
			Four_QuadrantBasedGeneticOutput(today, year, landscape_info_cfg.value());
			Nine_QuadrantBasedGeneticOutput(today, year, landscape_info_cfg.value());
		}else{
			cout<< "\n Please specify valid QuadrantsDesired (either 'four','nine' or 'both')!";
		}
		CompareGenomes(year, landscape_info_cfg.value());
	}
	if (today % 7==1){
		if (QuadrantsDesired=="four"){
			FourQuadrantsPopulationSizeProbe(today, year, landscape_info_cfg.value());
		}else if (QuadrantsDesired=="nine"){
			Nine_QuadrantBasedGeneticOutput(today, year, landscape_info_cfg.value());
		}else if(QuadrantsDesired=="both"){
			FourQuadrantsPopulationSizeProbe(today, year, landscape_info_cfg.value());
			NineQuadrantsPopulationSizeProbe(today, year, landscape_info_cfg.value());
		}else{
			cout << "\n Please specify valid QuadrantsDesired (either 'four', 'nine' or 'both')!";
		}
	}
	if(year%5==1 && today==1){
		GenerationCountOutput(year,today);
		LineagesOutput(year,today,1.0);
	}
	
	if (cfg_VoleUseResistanceOuput.value()) { if (cfg_VoleResistanceDay.value() == today) ResistanceOutput(); }

	bool terr = false;
	if (cfg_UseVoleTraplines.value())
	{
		int day = m_TheLandscape->SupplyGlobalDate();
		// Loop through all voles and see if any of them are in trap locations
		int s = static_cast<int>(GetLiveArraySize(vob_JuvenileMale));
		for (int i = 0; i < s; i++)
		{
			auto VJM = dynamic_cast<Vole_JuvenileMale*>(SupplyAnimalPtr(vob_JuvenileMale, i));
			if (VJM->SupplyInTrap())
			{
				InTrapPosition tp = VJM->SupplyTrapPosition();
				m_Traplines->Output(tp, day, vob_JuvenileMale, false, VJM->SupplyAge(), VJM->SupplyXBorn(),
				                    VJM->SupplyYBorn(), VJM->SupplyIDNo());
				VJM->SetFree();
			}
		}

		s = static_cast<int>(GetLiveArraySize(vob_JuvenileFemale));
		for (int i = 0; i < s; i++)
		{
			auto VJF = dynamic_cast<Vole_JuvenileFemale*>(SupplyAnimalPtr(vob_JuvenileFemale, i));
			if (VJF->SupplyInTrap())
			{
				InTrapPosition tp = VJF->SupplyTrapPosition();
				m_Traplines->Output(tp, day, vob_JuvenileFemale, false, VJF->SupplyAge(), VJF->SupplyXBorn(),
				                    VJF->SupplyYBorn(), VJF->SupplyIDNo());
				VJF->SetFree();
			}
		}

		s = static_cast<int>(GetLiveArraySize(vob_Female));
		for (int i = 0; i < s; i++)
		{
			auto VF = dynamic_cast<Vole_Female*>(SupplyAnimalPtr(vob_Female, i));
			if (VF->SupplyInTrap())
			{
				InTrapPosition tp = VF->SupplyTrapPosition();
				terr = VF->SupplyTerritorial();
				m_Traplines->Output(tp, day, vob_Female, terr, VF->SupplyAge(), VF->SupplyXBorn(), VF->SupplyYBorn(),
				                    VF->SupplyIDNo());
				VF->SetFree();
			}
		}

		s = static_cast<int>(GetLiveArraySize(vob_Male));
		for (int i = 0; i < s; i++)
		{
			auto VM = dynamic_cast<Vole_Male*>(SupplyAnimalPtr(vob_Male, i));
			if (VM->SupplyInTrap())
			{
				InTrapPosition tp = VM->SupplyTrapPosition();
				terr = VM->SupplyTerritorial();
				m_Traplines->Output(tp, day, vob_Male, terr, VM->SupplyAge(), VM->SupplyXBorn(), VM->SupplyYBorn(),
				                    VM->SupplyIDNo());
				VM->SetFree();
			}
		}
	}
	if (cfg_useagesexlocation.value()) { if (g_date->GetDayInMonth() == 1) TheAgeSexLocationProbe(); }
}
//---------------------------------------------------------------------------

/**
\brief Assess the quality of habitat at p_Polyref.
*/
double Vole_Population_Manager::AssessHabitat(const int p_Polyref) const {
	/**
	* Each unit of area (1m2) is assigned a score based on the polygon type. These scores fall currently into
	* a limited number of catagories, with 4, 3, 2.5, 2, 1, & 0 being the only scores currently used. 2.5 is only used by young forest
	* which is a hybrid between 2 & 3 represented by grassland and forest respectively. 4 denotes the best possible conditions.<br>
	* Once the score has been obtained it is modified by the digestability of the vegetation in terms of proportion of young green matter.
	* This modifies the value by multiplication with 0.7-1.0
	*/

	//double digestability = ( m_TheLandscape->SupplyVegDigestibilityVector(p_Polyref) + 0.2 );
	//double digestability = ( m_TheLandscape->SupplyVegDigestibilityVector(p_Polyref)+0.1) * (1.0/0.9);
	const double digestability = 1.0;
	const double score = vole_toletoc_asses_habitat_score(m_TheLandscape, p_Polyref);
	return score * digestability;
}
//---------------------------------------------------------------------------

bool Vole_Population_Manager::RecordGeneticsToday(const int p_today, int p_year, const int p_start_year,
                                                  const int p_interval) {
	if (p_today != 1) return false;
	if (p_year < p_start_year) return false;
	p_year -= p_start_year;
	if (p_year % p_interval != 0) return false;
	return true;
}

bool Vole_Population_Manager::SupplyOlderFemales(const unsigned p_x, const unsigned p_y, const unsigned p_Age,
                                                 const unsigned p_range) const {
	/** Returns false if there is an older female within the area p_x,p_y +/- range
	*/
	// Before checking the map remove ourselves so we don't count
	const auto caller = m_VoleMap->GetMapValue(p_x, p_y);
	m_VoleMap->ClearMapValue(p_x, p_y);
	// This is reset when the result is known
	int x = p_x - p_range;
	if (x < 0) x += SimW; // ensure we start in the simulation area!
	int y = p_y - p_range;
	if (y < 0) y += SimH;
	const int range_x = p_range * 2;
	const int range_y = p_range * 2;
	// create the extent variables
	const int xextent0 = x + range_x;
	const int yextent0 = y + range_y;
	const int xextent1 = x + range_x - SimW;
	const int yextent1 = y + range_y - SimH;
	// Create the looping variables needed
	// Create the looping variables needed
	int Dfinx;
	int Dfiny;
	int Afinx = 0; //unless the finx values for A are changed this stop
	int Afiny = 0; //the loop from executing
	const int Asty = 0;
	// int Astx,Dstx,Dtsy are always default so variables not used from them
	// NB Astx, Asty and Dstx are always 0, 0 & x respectively
	// Dsty is always y, Afiny is always yextent1 if it is used
	// Now create the loop values;
	if (xextent0 <= SimW) // No overlap with the eastern side
	{
		// Dstx, Dsty, Asty set by defaults
		//Astx & Afinx are not needed
		Dfinx = xextent0;
		// do we overlap the bottom?
		// Type B & D (overlap bottom, no overlap)
		if (yextent0 > SimH)
		{
			// Type B (overlap bottom only)
			Dfiny = SimH; // stop at the end
			Afiny = yextent1; // the overlap with the top
		}
		else Dfiny = yextent0;
	}
	else
	{
		// Type A & C overlap bottom & eastern edgdes
		if (yextent0 > SimH)
		{
			// relies on the default start for Asty, Astx, Dstx, Dsty
			Afinx = xextent1;
			Afiny = yextent1;
			Dfinx = SimW; // Stop at the end
			Dfiny = SimH;
		}
		else
		{
			// Type C overlap left edge only
			// Astx & Afiny are not needed here
			//Astx, Dstx, Dsty set by default
			Afinx = xextent1;
			Dfinx = SimW; // Stop at the end
			Dfiny = yextent0;
		}
	}
	Vole_Base* Ap;
	// A Loop
	for (int i = 0; i < Afinx; i++)
	{
		for (int j = Asty; j < Afiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (!Ap->SupplySex())
				{
					if (Ap->SupplyAge() >= p_Age)
					{
						m_VoleMap->SetMapValue(p_x, p_y, caller);
						return false;
					}
				}
			}
		}
		// C Loop
		for (int j = y; j < Dfiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (!Ap->SupplySex())
				{
					if (Ap->SupplyAge() >= p_Age)
					{
						m_VoleMap->SetMapValue(p_x, p_y, caller);
						return false;
					}
				}
			}
		}
	}
	// D Loop
	for (int i = x; i < Dfinx; i++)
	{
		for (int j = y; j < Dfiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (!Ap->SupplySex()) // female
				{
					if (Ap->SupplyAge() >= p_Age)
					{
						m_VoleMap->SetMapValue(p_x, p_y, caller);
						return false;
					}
				}
			}
		}
		// B Loop
		for (int j = 0; j < Afiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (!Ap->SupplySex()) // female
				{
					if (Ap->SupplyAge() >= p_Age)
					{
						m_VoleMap->SetMapValue(p_x, p_y, caller);
						return false;
					}
				}
			}
		}
	}
	// End of search algorithm
	m_VoleMap->SetMapValue(p_x, p_y, caller);
	return true;
}
//---------------------------------------------------------------------------

/**
  Counts vole postions on the map from a_x-a_range,p-y-a_range to a_x+p_size, a_y+a_range
*/
int Vole_Population_Manager::SupplyHowManyVoles(const unsigned a_x, const unsigned a_y, const unsigned a_range) const {
	int x = a_x - a_range;
	if (x < 0) x += SimW; // ensure we start in the simulation area!
	int y = a_y - a_range;
	if (y < 0) y += SimH;
	const int range_x = a_range * 2;
	const int range_y = a_range * 2;
	int Voles = 0;
	// create the extent variables
	const int xextent0 = x + range_x;
	const int yextent0 = y + range_y;
	const int xextent1 = x + range_x - SimW;
	const int yextent1 = y + range_y - SimH;
	// Create the looping variables needed
	// Create the looping variables needed
	int Dfinx;
	int Dfiny;
	int Afinx = 0; //unless the finx values for A are changed this stop
	int Afiny = 0; //the loop from executing
	// int Astx,Dstx,Dtsy are always default so variables not used from them
	// NB Astx, Asty and Dstx are always 0, 0 & x respectively
	// Dsty is always y, Afiny is always yextent1 if it is used
	// Now create the loop values;
	if (xextent0 <= SimW) // No overlap with the eastern side
	{
		// Dstx, Dsty, Asty set by defaults
		//Astx & Afinx are not needed
		Dfinx = xextent0;
		// do we overlap the bottom?
		// Type B & D (overlap bottom, no overlap)
		if (yextent0 > SimH)
		{
			// Type B (overlap bottom only)
			Dfiny = SimH; // stop at the end
			Afiny = yextent1; // the overlap with the top
		}
		else Dfiny = yextent0;
	}
	else
	{
		// Type A & C overlap bottom & eastern edgdes
		if (yextent0 > SimH)
		{
			// relies on the default start for Asty, Astx, Dstx, Dsty
			Afinx = xextent1;
			Afiny = yextent1;
			Dfinx = SimW; // Stop at the end
			Dfiny = SimH;
		}
		else
		{
			// Type C overlap left edge only
			// Astx & Afiny are not needed here
			//Astx, Dstx, Dsty set by default
			Afinx = xextent1;
			Dfinx = SimW; // Stop at the end
			Dfiny = yextent0;
		}
	}
	// A Loop
	for (int i = 0; i < Afinx; i++)
	{
		const int asty = 0;
		for (int j = asty; j < Afiny; j++) { if (m_VoleMap->GetMapValue(i, j) != nullptr) Voles++; }
		// C Loop
		for (int j = y; j < Dfiny; j++) { if (m_VoleMap->GetMapValue(i, j) != nullptr) Voles++; }
	}
	// D Loop
	for (int i = x; i < Dfinx; i++)
	{
		for (int j = y; j < Dfiny; j++) { if (m_VoleMap->GetMapValue(i, j) != nullptr) Voles++; }
		// B Loop
		for (int j = 0; j < Afiny; j++) { if (m_VoleMap->GetMapValue(i, j) != nullptr) Voles++; }
	}
	// End of search algorithm
	return Voles;
}
//---------------------------------------------------------------------------

void Vole_Population_Manager::ReproductionProbe() const {
	fprintf(YoungsFile, "%i   %i\n", YoungProducedToday, JuvsProducedToday);
}
//---------------------------------------------------------------------------

/** Prevents voles from mating with an individual on the other side of af barrier*/
bool Vole_Population_Manager::BarrierSearch(int F_x, int F_y, int M_x, int M_y) const {
	int poly2 = 200;

	if (F_x == M_x)
	{
		if (F_y == M_y) return true; // if they are on the same spot
		if (F_y < M_y)
		{
			// if they have the same x coordinate and different y coordinates and the female has the lowest y
			for (int j = F_y; j < M_y; j++)
			{
				int poly1 = m_TheLandscape->SupplyPolyRef(F_x, j);
				if (poly2 != poly1)
				{
					poly2 = poly1;
					if (vole_tole_assess_barrier(m_TheLandscape, poly1) == false) return false;
				}
			}
			return true;
		}
		// if (F_y > M_y) if they have the same x coordinate and different y coordinates and the male has the lowest y
		for (int j = M_y; j < F_y; j++)
		{
			int poly1 = m_TheLandscape->SupplyPolyRef(M_x, j);
			if (poly1 != poly2)
			{
				poly2 = poly1;
				if (vole_tole_assess_barrier(m_TheLandscape, poly1) == false) return false;
			}
		}
		return true;
	}
	if (F_y == M_y)
	{
		if (F_x < M_x)
		{
			// if they have the same y coordinate and different x coordinates and the female's x coordinate is the lowest
			for (int i = F_x; i < M_x; i++)
			{
				int poly1 = m_TheLandscape->SupplyPolyRef(i, F_y);
				if (poly1 != poly2)
				{
					poly2 = poly1;
					if (vole_tole_assess_barrier(m_TheLandscape, poly1) == false) return false;
				}
			}
			return true;
		}
		// (F_x > M_x) if they have the same y coordinate and different x coordinates and the male's x coordinate is the lowest
		for (int i = M_x; i < F_x; i++)
		{
			int poly1 = m_TheLandscape->SupplyPolyRef(i, M_y);
			if (poly1 != poly2)
			{
				poly2 = poly1;
				if (vole_tole_assess_barrier(m_TheLandscape, poly1) == false) return false;
			}
		}
		return true;
	}
	if (F_x < M_x)
	{
		// if the x coordinates differ and the females is the lowest
		if (F_y < M_y)
		{
			// if the y coordinates differ and the female has the lowest y coordinate
			unsigned diff_x = M_x - F_x;
			unsigned diff_y = M_y - F_y;
			if (diff_x <= diff_y)
			{
				// if the area between them is enlongated (y) or diagonal
				int j = F_y;

				for (int i = F_x; i < M_x; i++)
				{
					// the diagonal part
					int poly1 = m_TheLandscape->SupplyPolyRef(i, j);
					j++;
					if (poly1 != poly2)
					{
						poly2 = poly1;
						if (vole_tole_assess_barrier(m_TheLandscape, poly1) == false) return false;
					}
				}
				if (diff_x != diff_y)
				{
					// if not diagonal but enlongated and the x and y coordinates differ with the female having the lowest values
					int j_extent = F_y + diff_x;
					for (int j_y = j_extent; j_y < M_y; j_y++)
					{
						// the enlongated part
						int poly1 = m_TheLandscape->SupplyPolyRef(M_x, j_y);
						if (poly1 != poly2)
						{
							poly2 = poly1;
							if (vole_tole_assess_barrier(m_TheLandscape, poly1) == false) return false;
						}
					}
				}
				return true;
			}
			if (diff_x > diff_y)
			{
				// if not diagonal but widened and the x and y coordinates differ with the female having the lowest values
				int i = F_x;
				for (int j = F_y; j < M_y; j++)
				{
					// the diagonal part
					int poly1 = m_TheLandscape->SupplyPolyRef(i, j);
					i++;
					if (poly1 != poly2)
					{
						poly2 = poly1;
						if (vole_tole_assess_barrier(m_TheLandscape, poly1) == false) return false;
					}
				}
				int i_extent = diff_y + F_x;
				for (int i_x = i_extent; i_x < M_x; i_x++)
				{
					// the widened part
					int poly1 = m_TheLandscape->SupplyPolyRef(i_x, M_y);
					if (poly1 != poly2)
					{
						poly2 = poly1;
						if (vole_tole_assess_barrier(m_TheLandscape, poly1) == false) return false;
					}
				}
				return true;
			}
		}
		if (F_y > M_y)
		{
			// and if (F_x < M_x)
			int diff_x = M_x - F_x;
			int diff_y = F_y - M_y;
			if (diff_x <= diff_y)
			{
				// if the area between them is enlongated (y) or diagonal
				int j = F_y;
				for (int i = F_x; i < M_x; i++)
				{
					// the diagonal part
					int poly1 = m_TheLandscape->SupplyPolyRef(i, j);
					j--;
					if (poly1 != poly2)
					{
						poly2 = poly1;
						if (vole_tole_assess_barrier(m_TheLandscape, poly1) == false) return false;
					}
				}
				if (diff_x != diff_y)
				{
					// if not diagonal but enlongated and the x and y coordinates differ with the female having the lowest x value and highest y
					int j_extent = F_y - diff_x;
					for (int j_y = j_extent; j_y > M_y; j_y--)
					{
						// the enlongated part
						int poly1 = m_TheLandscape->SupplyPolyRef(M_x, j_y);
						if (poly1 != poly2)
						{
							poly2 = poly1;
							if (vole_tole_assess_barrier(m_TheLandscape, poly1) == false) return false;
						}
					}
				}
				return true;
			}
			// if (diff_x > diff_y) if not diagonal but widened and the x and y coordinates differ with the female having the lowest x and highest y
			int i = F_x;
			for (int j = F_y; j > M_y; j--)
			{
				// the diagonal part
				int poly1 = m_TheLandscape->SupplyPolyRef(i, j);
				i++;
				if (poly1 != poly2)
				{
					poly2 = poly1;
					if (vole_tole_assess_barrier(m_TheLandscape, poly1) == false) return false;
				}
			}
			int i_extent = diff_y + F_x;
			for (int i_x = i_extent; i_x < M_x; i_x++)
			{
				// the widened part
				int poly1 = m_TheLandscape->SupplyPolyRef(i_x, M_y);
				if (poly1 != poly2)
				{
					poly2 = poly1;
					if (vole_tole_assess_barrier(m_TheLandscape, poly1) == false) return false;
				}
			}
			return true;
		}
	}
	if (F_x > M_x)
	{
		if (F_y < M_y)
		{
			unsigned diff_x = F_x - M_x;
			unsigned diff_y = M_y - F_y;
			if (diff_x <= diff_y)
			{
				// if the area between them is enlongated (y) or diagonal
				int j = M_y;
				for (int i = M_x; i < F_x; i++)
				{
					// the diagonal part
					int poly1 = m_TheLandscape->SupplyPolyRef(i, j);
					j--;
					if (poly1 != poly2)
					{
						poly2 = poly1;
						if (vole_tole_assess_barrier(m_TheLandscape, poly1) == false) return false;
					}
				}

				if (diff_x != diff_y)
				{
					// if not diagonal but enlongated and the x and y coordinates differ with the female having the lowest values
					int j_extent = M_y - diff_x;
					for (int j_y = j_extent; j_y > F_y; j_y--)
					{
						// the enlongated part
						int poly1 = m_TheLandscape->SupplyPolyRef(F_x, j_y);
						if (poly1 != poly2)
						{
							poly2 = poly1;
							if (vole_tole_assess_barrier(m_TheLandscape, poly1) == false) return false;
						}
					}
				}
				return true;
			}
			// if (diff_x > diff_y) if not diagonal but widened and the x and y coordinates differ with the female having the lowest values
			int i = M_x;
			for (int j = M_y; j > F_y; j--)
			{
				int poly1 = m_TheLandscape->SupplyPolyRef(i, j);
				i++;
				if (poly1 != poly2)
				{
					poly2 = poly1;
					if (vole_tole_assess_barrier(m_TheLandscape, poly1) == false) return false;
				}
			}
			int i_extent = diff_y + M_x;
			for (int i_x = i_extent; i_x < F_x; i_x++)
			{
				// the widened part
				int poly1 = m_TheLandscape->SupplyPolyRef(i_x, F_y);
				if (poly1 != poly2)
				{
					poly2 = poly1;
					if (vole_tole_assess_barrier(m_TheLandscape, poly1) == false) return false;
				}
			}
			return true;
		}
		if (F_y > M_y)
		{
			// if (F_x > M_x)
			int diff_x = F_x - M_x;
			int diff_y = F_y - M_y;
			if (diff_x <= diff_y)
			{
				// if the area between them is enlongated (y) or diagonal
				int j = M_y;
				for (int i = M_x; i < F_x; i++)
				{
					// the diagonal part
					int poly1 = m_TheLandscape->SupplyPolyRef(i, j);
					j++;
					if (poly1 != poly2)
					{
						poly2 = poly1;
						if (vole_tole_assess_barrier(m_TheLandscape, poly1) == false) return false;
					}
				}
				if (diff_x != diff_y)
				{
					// if not diagonal but enlongated and the x and y coordinates differ with the female having the lowest x value and highest y
					int j_extent = M_y + diff_x;
					for (int j_y = j_extent; j_y < F_y; j_y++)
					{
						// the enlongated part
						int poly1 = m_TheLandscape->SupplyPolyRef(F_x, j_y);
						if (poly1 != poly2)
						{
							poly2 = poly1;
							if (vole_tole_assess_barrier(m_TheLandscape, poly1) == false) return false;
						}
					}
				}
				return true;
			}
			// if (diff_x > diff_y) if not diagonal but widened and the x and y coordinates differ with the female having the lowest x and highest y
			int i = M_x;
			for (int j = M_y; j < F_y; j++)
			{
				// the diagonal part
				int poly1 = m_TheLandscape->SupplyPolyRef(i, j);
				i++;
				if (poly1 != poly2)
				{
					poly2 = poly1;
					if (vole_tole_assess_barrier(m_TheLandscape, poly1) == false) return false;
				}
			}
			int i_extent = diff_y + M_x;
			for (int i_x = i_extent; i_x < F_x; i_x++)
			{
				// the widened part
				int poly1 = m_TheLandscape->SupplyPolyRef(i_x, F_y);
				if (poly1 != poly2)
				{
					poly2 = poly1;
					if (vole_tole_assess_barrier(m_TheLandscape, poly1) == false) return false;
				}
			}
			return true; //returns true if both the diagonal and widened part dont hit return false
		}
	}
	return true;
}

//-----------------------------------------------------------------------------

/**
Looks for males within p_Steps of p_x,p_y and returns the number of them. Pointers to these are saved in MList
*/
int Vole_Population_Manager::ListClosestMales(const int p_x, const int p_y, const int p_steps) {
	// First clear the MaleVoleList
	MList.clear();

	// as usual there are 4 possibilities of overlap
	// First convert centre co-rdinates to square co-ordinates

	int x = p_x - p_steps;
	if (x < 0) x += SimW; // ensure we start in the simulation area!
	int y = p_y - p_steps;
	if (y < 0) y += SimH;
	const int range_x = p_steps * 2;
	const int range_y = p_steps * 2;

	// create the extent variables
	const int xextent0 = x + range_x;
	const int yextent0 = y + range_y;
	const int xextent1 = x + range_x - SimW;
	const int yextent1 = y + range_y - SimH;
	// Create the looping variables needed
	// Create the looping variables needed
	int Dfinx;
	int Dfiny;
	int Afinx = 0; //unless the finx values for A are changed this stop
	int Afiny = 0; //the loop from executing
	const int Asty = 0;
	// int Astx,Dstx,Dtsy are always default so variables not used from them
	// NB Astx, Asty and Dstx are always 0, 0 & x respectively
	// Dsty is always y, Afiny is always yextent1 if it is used
	// Now create the loop values;
	if (xextent0 <= SimW) // No overlap with the eastern side
	{
		// Dstx, Dsty, Asty set by defaults
		//Astx & Afinx are not needed
		Dfinx = xextent0;
		// do we overlap the bottom?
		// Type B & D (overlap bottom, no overlap)
		if (yextent0 > SimH)
		{
			// Type B (overlap bottom only)
			Dfiny = SimH; // stop at the end
			Afiny = yextent1; // the overlap with the top
		}
		else Dfiny = yextent0;
	}
	else
	{
		// Type A & C overlap bottom & eastern edgdes
		if (yextent0 > SimH)
		{
			// relies on the default start for Asty, Astx, Dstx, Dsty
			Afinx = xextent1;
			Afiny = yextent1;
			Dfinx = SimW; // Stop at the end
			Dfiny = SimH;
		}
		else
		{
			// Type C overlap left edge only
			// Astx & Afiny are not needed here
			//Astx, Dstx, Dsty set by default
			Afinx = xextent1;
			Dfinx = SimW; // Stop at the end
			Dfiny = yextent0;
		}
	}
	Vole_Base* Ap;
	Vole_Male* AMale;
	int NoFound = 0;

	// A Loop
	for (int i = 0; i < Afinx; i++)
	{
		for (int j = Asty; j < Afiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (Ap->SupplySex()) // male
				{
					AMale = static_cast<Vole_Male*>(Ap);
					if (AMale->SupplyTerritorial())
					{
						MList.push_back(AMale);
						NoFound++;
					}
				}
			}
		}
		// C Loop
		for (int j = y; j < Dfiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (Ap->SupplySex()) // male
				{
					AMale = static_cast<Vole_Male*>(Ap);
					if (AMale->SupplyTerritorial())
					{
						MList.push_back(AMale);
						NoFound++;
					}
				}
			}
		}
	}
	// D Loop
	for (int i = x; i < Dfinx; i++)
	{
		for (int j = y; j < Dfiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (Ap->SupplySex()) // male
				{
					AMale = static_cast<Vole_Male*>(Ap);
					if (AMale->SupplyTerritorial())
					{
						MList.push_back(AMale);
						NoFound++;
					}
				}
			}
		}
		// B Loop
		for (int j = 0; j < Afiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (Ap->SupplySex()) // male
				{
					AMale = static_cast<Vole_Male*>(Ap);
					if (AMale->SupplyTerritorial())
					{
						MList.push_back(AMale);
						NoFound++;
					}
				}
			}
		}
	}
	return NoFound;
}
//---------------------------------------------------------------------------

/** Lists all females within p_steps of p_x & p_y and returns the number in the list */
int Vole_Population_Manager::ListClosestFemales(const int p_x, const int p_y, const int p_steps) {
	// First clear the FemaleVoleList
	FList.clear();
	// looks for a female within p_Steps of p_x,p_y
	// as usual there are 4 possibilities of overlap
	// First convert centre co-ordinates to square co-ordinates

	int x = p_x - p_steps;
	if (x < 0) x += SimW; // ensure we start in the simulation area!
	int y = p_y - p_steps;
	if (y < 0) y += SimH;
	const int range_x = p_steps * 2;
	const int range_y = p_steps * 2;

	// create the extent variables
	const int xextent0 = x + range_x;
	const int yextent0 = y + range_y;
	const int xextent1 = x + range_x - SimW;
	const int yextent1 = y + range_y - SimH;
	// Create the looping variables needed
	// Create the looping variables needed
	int Dfinx;
	int Dfiny;
	int Afinx = 0; //unless the finx values for A are changed this stop
	int Afiny = 0; //the loop from executing
	const int Asty = 0;
	// int Astx,Dstx,Dtsy are always default so variables not used from them
	// NB Astx, Asty and Dstx are always 0, 0 & x respectively
	// Dsty is always y, Afiny is always yextent1 if it is used
	// Now create the loop values;
	if (xextent0 <= SimW) // No overlap with the eastern side
	{
		// Dstx, Dsty, Asty set by defaults
		//Astx & Afinx are not needed
		Dfinx = xextent0;
		// do we overlap the bottom?
		// Type B & D (overlap bottom, no overlap)
		if (yextent0 > SimH)
		{
			// Type B (overlap bottom only)
			Dfiny = SimH; // stop at the end
			Afiny = yextent1; // the overlap with the top
		}
		else Dfiny = yextent0;
	}
	else
	{
		// Type A & C overlap bottom & eastern edgdes
		if (yextent0 > SimH)
		{
			// relies on the default start for Asty, Astx, Dstx, Dsty
			Afinx = xextent1;
			Afiny = yextent1;
			Dfinx = SimW; // Stop at the end
			Dfiny = SimH;
		}
		else
		{
			// Type C overlap left edge only
			// Astx & Afiny are not needed here
			//Astx, Dstx, Dsty set by default
			Afinx = xextent1;
			Dfinx = SimW; // Stop at the end
			Dfiny = yextent0;
		}
	}
	Vole_Female* AFemale;
	Vole_Base* Ap;
	int NoFound = 0;

	// A Loop
	for (int i = 0; i < Afinx; i++)
	{
		for (int j = Asty; j < Afiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (!Ap->SupplySex()) // female
				{
					AFemale = static_cast<Vole_Female*>(Ap);
					if (AFemale->SupplyTerritorial())
					{
						FList.push_back(AFemale);
						NoFound++;
					}
				}
			}
		}
		// C Loop
		for (int j = y; j < Dfiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (!Ap->SupplySex()) // female
				{
					AFemale = static_cast<Vole_Female*>(Ap);
					if (AFemale->SupplyTerritorial())
					{
						FList.push_back(AFemale);
						NoFound++;
					}
				}
			}
		}
	}
	// D Loop
	for (int i = x; i < Dfinx; i++)
	{
		for (int j = y; j < Dfiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (!Ap->SupplySex()) // female
				{
					AFemale = static_cast<Vole_Female*>(Ap);
					if (AFemale->SupplyTerritorial())
					{
						FList.push_back(AFemale);
						NoFound++;
					}
				}
			}
		}
		// B Loop
		for (int j = 0; j < Afiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (!Ap->SupplySex()) // female
				{
					AFemale = static_cast<Vole_Female*>(Ap);
					if (AFemale->SupplyTerritorial())
					{
						FList.push_back(AFemale);
						NoFound++;
					}
				}
			}
		}
	}
	return NoFound;
}
//---------------------------------------------------------------------------

Vole_Female* Vole_Population_Manager::FindClosestFemale(const int p_x, const int p_y, const int p_steps) const {
	/** looks for the closest female within p_Steps of p_x,p_y */
	// as usual there are 4 possibilities of overlap
	// First convert centre co-rdinates to square co-ordinates

	int x = p_x - p_steps;
	if (x < 0) x += SimW; // ensure we start in the simulation area!
	int y = p_y - p_steps;
	if (y < 0) y += SimH;
	const int range_x = p_steps * 2;
	const int range_y = p_steps * 2;

	// create the extent variables
	const int xextent0 = x + range_x;
	const int yextent0 = y + range_y;
	const int xextent1 = x + range_x - SimW;
	const int yextent1 = y + range_y - SimH;
	// Create the looping variables needed
	// Create the looping variables needed
	int Dfinx;
	int Dfiny;
	int Afinx = 0; //unless the finx values for A are changed this stop
	int Afiny = 0; //the loop from executing
	const int Asty = 0;
	// int Astx,Dstx,Dtsy are always default so variables not used from them
	// NB Astx, Asty and Dstx are always 0, 0 & x respectively
	// Dsty is always y, Afiny is always yextent1 if it is used
	// Now create the loop values;
	if (xextent0 <= SimW) // No overlap with the eastern side
	{
		// Dstx, Dsty, Asty set by defaults
		//Astx & Afinx are not needed
		Dfinx = xextent0;
		// do we overlap the bottom?
		// Type B & D (overlap bottom, no overlap)
		if (yextent0 > SimH)
		{
			// Type B (overlap bottom only)
			Dfiny = SimH; // stop at the end
			Afiny = yextent1; // the overlap with the top
		}
		else Dfiny = yextent0;
	}
	else
	{
		// Type A & C overlap bottom & eastern edgdes
		if (yextent0 > SimH)
		{
			// relies on the default start for Asty, Astx, Dstx, Dsty
			Afinx = xextent1;
			Afiny = yextent1;
			Dfinx = SimW; // Stop at the end
			Dfiny = SimH;
		}
		else
		{
			// Type C overlap left edge only
			// Astx & Afiny are not needed here
			//Astx, Dstx, Dsty set by default
			Afinx = xextent1;
			Dfinx = SimW; // Stop at the end
			Dfiny = yextent0;
		}
	}
	Vole_Female* AFemale;
	Vole_Female* Found = nullptr;
	int dist, disty;
	int FoundDist = SimW; // too big
	Vole_Base* Ap;

	// A Loop
	for (int i = 0; i < Afinx; i++)
	{
		for (int j = Asty; j < Afiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (!Ap->SupplySex()) // female
				{
					AFemale = static_cast<Vole_Female*>(Ap);
					if (AFemale->SupplyTerritorial())
					{
						dist = abs(p_x - i); // remove the signed bit
						disty = abs(p_y - j);
						dist += disty;
						if (dist < FoundDist)
						{
							Found = AFemale;
							FoundDist = dist;
						}
					}
				}
			}
		}
		// C Loop
		for (int j = y; j < Dfiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (!Ap->SupplySex()) // female
				{
					AFemale = static_cast<Vole_Female*>(Ap);
					if (AFemale->SupplyTerritorial())
					{
						dist = abs(p_x - i); // remove the signed bit
						disty = abs(p_y - j);
						dist += disty;
						if (dist < FoundDist)
						{
							Found = AFemale;
							FoundDist = dist;
						}
					}
				}
			}
		}
	}
	// D Loop
	for (int i = x; i < Dfinx; i++)
	{
		for (int j = y; j < Dfiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (!Ap->SupplySex()) // female
				{
					AFemale = static_cast<Vole_Female*>(Ap);
					if (AFemale->SupplyTerritorial())
					{
						dist = abs(p_x - i); // remove the signed bit
						disty = abs(p_y - j);
						dist += disty;
						if (dist < FoundDist)
						{
							Found = AFemale;
							FoundDist = dist;
						}
					}
				}
			}
		}
		// B Loop
		for (int j = 0; j < Afiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (!Ap->SupplySex()) // female
				{
					AFemale = static_cast<Vole_Female*>(Ap);
					if (AFemale->SupplyTerritorial())
					{
						dist = abs(p_x - i); // remove the signed bit
						disty = abs(p_y - j);
						dist += disty;
						if (dist < FoundDist)
						{
							Found = AFemale;
							FoundDist = dist;
						}
					}
				}
			}
		}
	}
	return Found;
}
//---------------------------------------------------------------------------
Vole_Male* Vole_Population_Manager::FindOutsideRadiusMale(const int p_x, const int p_y) //**TD**
{
	/** Returns a pointer to a random male - he does not even have to have a territory */
	int x = p_x - cfg_MateRadius.value(); // 'left most' x
	if (x < 0) x += SimW; // ensure we start within the landscape!
	int y = p_y - cfg_MateRadius.value(); // 'top' y
	if (y < 0) y += SimH;
	const int range_x = cfg_MateRadius.value() * 2; // width of the search area
	const int range_y = cfg_MateRadius.value() * 2; // hight of the search area

	// create the extent variables - if out of bounds
	const int xextent0 = x + range_x; // 'right most' x
	const int yextent0 = y + range_y; // 'lowest'(bottom) y
	const int xextent1 = x + range_x - SimW; // incase we are out of bounds
	const int yextent1 = y + range_y - SimH;

	// Create the looping variables needed
	int Dfinx;
	int Dfiny;
	int Afinx = 0; //unless the finx values for A are changed this stop
	int Afiny = 0; //the loop from executing
	//int Asty=0;

	// Dfinx, Dfiny are the stop x and y values for 'the mate picking area'
	// Afinx, Afiny are the 'overlapping area' if 'the mate picking area' goes out off bounds
	// Set default at 0 and only changes value if out off bounds

	// Now create the loop values;
	// Type A (no overlap with the eastern bound)
	if (xextent0 <= SimW)
	{
		//Afinx are not needed
		Dfinx = xextent0; // because no overlap with the right most bound
		// do we overlap the bottom (lowest y)?
		// Type B & D (overlap bottom, no overlap)
		if (yextent0 > SimH)
		{
			// Type B (overlap bottom only)
			Dfiny = SimH; // stop at the end
			Afiny = yextent1; // the overlap with the top because off bounds
		}
		else Dfiny = yextent0; // if no overlap
	}
	else
	{
		// Type A & C eastern bound & overlap bottom
		if (yextent0 > SimH)
		{
			// relies on the default start for Afiny, Afinx, Dfinx, Dfiny
			Afinx = xextent1; // Afinx gets a value other than 0 = overlapping value
			Afiny = yextent1; //
			Dfinx = SimW; // Stop at the end
			Dfiny = SimH;
		}
		else
		{
			// Type C overlap west bound only
			//Afiny, Dfinx, Dfiny set by default
			Afinx = xextent1;
			Dfinx = SimW; // Stop at the end
			Dfiny = yextent0;
		}
	}

	Vole_Male* Found = nullptr;
	const int size = static_cast<int>(GetLiveArraySize(vob_Male));

	do
	{
		const int test = g_random_fnc(size);
		TAnimal* Ap = SupplyAnimalPtr(vob_Male, test); // Males only
		Found = static_cast<Vole_Male*>(Ap); // Male pointer to TAnimal[0][test]
		const int territory = Found->SupplyTerritorial();
		if (Found->GetCurrentStateNo() == -1) Found = nullptr; // If end of timestep reached skip check
		{
			const int i = Found->Supply_m_Location_x();
			const int j = Found->Supply_m_Location_y();

			if (Afinx == 0) // no right overlap
			{
				if (Afiny == 0) // no bottom overlap
				{
					if (x <= i && i <= Dfinx && (y <= j && j <= Dfiny)) { Found = nullptr; }
				}

				if (Afiny != 0) // overlap bottom
				{
					if (x <= i && i <= Dfinx && ((0 <= j && j <= Afiny) || (y <= j && j <= Dfiny))) { Found = nullptr; }
				}
			}
			if (Afinx != 0) // overlap right
			{
				if (Afiny == 0)
				{
					if (y <= j && j <= Dfiny && ((0 <= i && i <= Afinx) || (x <= i && i <= Dfinx))) { Found = nullptr; }
				}

				if (Afiny != 0) //overlap bottom
				{
					if (0 <= i && i <= Afinx && ((0 <= j && j <= Afiny) || (y <= j && j <= Dfiny))) { Found = nullptr; }
					if (x <= i && i <= Dfinx && ((0 <= j && j <= Afiny) || (y <= j && j <= Dfiny))) { Found = nullptr; }
				}
			}
		}
		if (territory == false) Found = nullptr;
	}
	while (Found == nullptr && size > 0);
	return Found;
}

//--------------------------------------------------------------------------------------

Vole_Male* Vole_Population_Manager::FindWithinRadiusMale(const int p_x, const int p_y) const
// ***TD***
{
	/** looks for males within cfg_MateRadius of p_x,p_y and returns a list of all males*/

	const auto vbl = new vector<Vole_Male*>;

	// as usual there are 4 possibilities of overlap
	// First convert centre co-rdinates to square co-ordinates
	// To create start x and y values for 'the mate picking area'
	int x = p_x - cfg_MateRadius.value(); // 'left most' x
	if (x < 0) x += SimW; // ensure we start within the landscape!
	int y = p_y - cfg_MateRadius.value(); // 'top' y
	if (y < 0) y += SimH;
	const int range_x = cfg_MateRadius.value() * 2; // width of the search area
	const int range_y = cfg_MateRadius.value() * 2; // hight of the search area

	// create the extent variables - if out of bounds
	const int xextent0 = x + range_x; // 'right most' x
	const int yextent0 = y + range_y; // 'lowest'(bottom) y
	const int xextent1 = x + range_x - SimW; // incase we are out of bounds (positive if we are)
	const int yextent1 = y + range_y - SimH;

	// Create the looping variables needed
	int Dfinx;
	int Dfiny;
	int Afinx = 0; //unless the finx values for A are changed this stop
	int Afiny = 0; //the loop from executing
	//int Asty=0;

	// Dfinx, Dfiny are the stop x and y values for 'the mate picking area'
	// Afinx, Afiny are the 'overlapping area' if 'the mate picking area' goes out off bounds
	// Set default at 0 and only changes value if out off bounds

	// Now create the loop values;
	// Type A (no overlap with the eastern bound)
	if (xextent0 <= SimW)
	{
		//Afinx are not needed
		Dfinx = xextent0; // because no overlap with the right most bound
		// do we overlap the bottom (lowest y)?
		// Type B & D (overlap bottom, no overlap)
		if (yextent0 > SimH)
		{
			// Type B (overlap bottom only)
			Dfiny = SimH; // stop at the end
			Afiny = yextent1; // the overlap with the top because off bounds
		}
		else Dfiny = yextent0; // if no overlap bottom
	}
	else
	{
		// Type A & C eastern bound & overlap bottom
		if (yextent0 > SimH)
		{
			// relies on the default start for Afiny, Afinx, Dfinx, Dfiny
			Afinx = xextent1; // Afinx gets a value other than 0 = overlapping value
			Afiny = yextent1; //
			Dfinx = SimW; // Stop at the end
			Dfiny = SimH;
		}
		else
		{
			// Type C overlap west bound only
			//Afiny, Dfinx, Dfiny set by default
			Afinx = xextent1;
			Dfinx = SimW; // Stop at the end
			Dfiny = yextent0;
		}
	}
	Vole_Male* AMale;
	Vole_Male* Found = nullptr;
	Vole_Base* Ap;

	// A Loop
	for (int i = 0; i < Afinx; i++) // Afinx and y will be 0 if no overlap
	{
		for (int j = 0; j < Afiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j)); // Tjeck to see if any voles are at the location
			if (Ap)
			{
				if (Ap->SupplySex()) // male
				{
					AMale = static_cast<Vole_Male*>(Ap);
					if (AMale->SupplyTerritorial()) { vbl->push_back(AMale); }
				}
			}
		}
		// C Loop
		for (int j = y; j < Dfiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j)); // Tjeck to see if any voles are at the location
			if (Ap)
			{
				if (Ap->SupplySex()) // male
				{
					AMale = static_cast<Vole_Male*>(Ap);
					if (AMale->SupplyTerritorial()) { vbl->push_back(AMale); }
				}
			}
		}
	}
	// D Loop
	for (int i = x; i < Dfinx; i++)
	{
		for (int j = y; j < Dfiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j)); // Tjeck to see if any voles are at the location
			if (Ap)
			{
				if (Ap->SupplySex()) // male
				{
					AMale = static_cast<Vole_Male*>(Ap);
					if (AMale->SupplyTerritorial()) { vbl->push_back(AMale); }
				}
			}
		}
		// B Loop
		for (int j = 0; j < Afiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j)); // Tjeck to see if any voles are at the location
			if (Ap)
			{
				if (Ap->SupplySex()) // male
				{
					AMale = static_cast<Vole_Male*>(Ap);
					if (AMale->SupplyTerritorial()) { vbl->push_back(AMale); }
				}
			}
		}
	}

	const int l = vbl->size();
	if (l == 0) { return Found = nullptr; }
	const int i = g_random_fnc(l);
	Found = vbl->at(i);
	return Found;
}
//--------------------------------------------------------------------------------------

Vole_Male* Vole_Population_Manager::FindClosestMale(const int p_x, const int p_y, const int p_steps) const {
	/** looks for the closest male within p_Steps of p_x,p_y */
	// as usual there are 4 possibilities of overlap
	// First convert centre co-rdinates to square co-ordinates
	int x = p_x - p_steps;
	if (x < 0) x += SimW; // ensure we start in the simulation area!
	int y = p_y - p_steps;
	if (y < 0) y += SimH;
	const int range_x = p_steps * 2;
	const int range_y = p_steps * 2;

	// create the extent variables
	const int xextent0 = x + range_x;
	const int yextent0 = y + range_y;
	const int xextent1 = x + range_x - SimW;
	const int yextent1 = y + range_y - SimH;
	// Create the looping variables needed
	// Create the looping variables needed
	int Dfinx;
	int Dfiny;
	int Afinx = 0; //unless the finx values for A are changed this stop
	int Afiny = 0; //the loop from executing
	const int Asty = 0;
	// int Astx,Dstx,Dtsy are always default so variables not used from them
	// NB Astx, Asty and Dstx are always 0, 0 & x respectively
	// Dsty is always y, Afiny is always yextent1 if it is used
	// Now create the loop values;
	if (xextent0 <= SimW) // No overlap with the eastern side
	{
		// Dstx, Dsty, Asty set by defaults
		//Astx & Afinx are not needed
		Dfinx = xextent0;
		// do we overlap the bottom?
		// Type B & D (overlap bottom, no overlap)
		if (yextent0 > SimH)
		{
			// Type B (overlap bottom only)
			Dfiny = SimH; // stop at the end
			Afiny = yextent1; // the overlap with the top
		}
		else Dfiny = yextent0;
	}
	else
	{
		// Type A & C overlap bottom & eastern edgdes
		if (yextent0 > SimH)
		{
			// relies on the default start for Asty, Astx, Dstx, Dsty
			Afinx = xextent1;
			Afiny = yextent1;
			Dfinx = SimW; // Stop at the end
			Dfiny = SimH;
		}
		else
		{
			// Type C overlap left edge only
			// Astx & Afiny are not needed here
			//Astx, Dstx, Dsty set by default
			Afinx = xextent1;
			Dfinx = SimW; // Stop at the end
			Dfiny = yextent0;
		}
	}
	Vole_Male* AMale;
	Vole_Male* Found = nullptr;
	int dist, disty;
	int FoundDist = SimW; // too big
	Vole_Base* Ap;
	// A Loop
	for (int i = 0; i < Afinx; i++)
	{
		for (int j = Asty; j < Afiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (Ap->SupplySex()) // male
				{
					{
						AMale = static_cast<Vole_Male*>(Ap);
						if (AMale->SupplyTerritorial())
						{
							dist = abs(p_x - i); // remove the signed bit
							disty = abs(p_y - j);
							dist += disty;

							if (dist < FoundDist)
							{
								const bool Barrier = BarrierSearch(p_x, p_y, i, j);
								if (Barrier == true)
								{
									Found = AMale;
									FoundDist = dist;
								}
							}
						}
					}
				}
			}
		}
		// C Loop
		for (int j = y; j < Dfiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (Ap->SupplySex()) // male
				{
					{
						AMale = static_cast<Vole_Male*>(Ap);
						if (AMale->SupplyTerritorial())
						{
							dist = abs(p_x - i); // remove the signed bit
							disty = abs(p_y - j);
							dist += disty;

							if (dist < FoundDist)
							{
								const bool Barrier = BarrierSearch(p_x, p_y, i, j);
								if (Barrier == true)
								{
									Found = AMale;
									FoundDist = dist;
								}
							}
						}
					}
				}
			}
		}
	}
	// D Loop
	for (int i = x; i < Dfinx; i++)
	{
		for (int j = y; j < Dfiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (Ap->SupplySex()) // male
				{
					{
						AMale = static_cast<Vole_Male*>(Ap);
						if (AMale->SupplyTerritorial())
						{
							dist = abs(p_x - i); // remove the signed bit
							disty = abs(p_y - j);
							dist += disty;
							if (dist < FoundDist)
							{
								const bool Barrier = BarrierSearch(p_x, p_y, i, j);
								if (Barrier == true)
								{
									Found = AMale;
									FoundDist = dist;
								}
							}
						}
					}
				}
			}
		}
		// B Loop
		for (int j = 0; j < Afiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (Ap->SupplySex()) // male
				{
					AMale = static_cast<Vole_Male*>(Ap);
					{
						if (AMale->SupplyTerritorial())
						{
							dist = abs(p_x - i); // remove the signed bit
							disty = abs(p_y - j);
							dist += disty;
							if (dist < FoundDist)
							{
								const bool Barrier = BarrierSearch(p_x, p_y, i, j);
								if (Barrier == true)
								{
									Found = AMale;
									FoundDist = dist;
								}
							}
						}
					}
				}
			}
		}
	}
	return Found;
}
//---------------------------------------------------------------------------

Vole_Male* Vole_Population_Manager::FindRandomMale() {
	/** Returns a pointer to a random male - he does not even have to have a territory */
	Vole_Male* Found = nullptr;
	const int size = static_cast<int>(GetLiveArraySize(vob_Male));
	if (size < 1) return nullptr;
	do
	{
		const int test = g_random_fnc(size);
		TAnimal* Ap = SupplyAnimalPtr(vob_Male, test);
		Found = static_cast<Vole_Male*>(Ap);
		if (Found->GetCurrentStateNo() == -1) Found = nullptr; //if not alive
	}
	while (Found == nullptr && size > 0);
	return Found;
}
//---------------------------------------------------------------------------
int Vole_Population_Manager::SupplyCountFemales(const unsigned p_x, const unsigned p_y,
                                                const unsigned p_TerrRange) const {
	/** returns -1 if a male has p_x,p_y in his territory and is older than p_Age else returns the number of females present
	*/

	// Before checking the map remove ourselves so we don't count
	//PointerInt c=VoleMap->GetMapValue(p_x,p_y);
	//TAnimal* caller=(TAnimal*) c;
	TAnimal* caller = m_VoleMap->GetMapValue(p_x, p_y);
	m_VoleMap->ClearMapValue(p_x, p_y);
	int x = p_x - p_TerrRange;
	if (x < 0) x += SimW; // ensure we start in the simulation area!
	int y = p_y - p_TerrRange;
	if (y < 0) y += SimH;
	const int range_x = p_TerrRange * 2;
	const int range_y = p_TerrRange * 2;
	int Females = 0;

	// create the extent variables
	const int xextent0 = x + range_x;
	const int yextent0 = y + range_y;
	const int xextent1 = x + range_x - SimW;
	const int yextent1 = y + range_y - SimH;
	// Create the looping variables needed
	// Create the looping variables needed
	int Dfinx;
	int Dfiny;
	int Afinx = 0; //unless the finx values for A are changed this stop
	int Afiny = 0; //the loop from executing
	const int Asty = 0;
	// int Astx,Dstx,Dtsy are always default so variables not used from them
	// NB Astx, Asty and Dstx are always 0, 0 & x respectively
	// Dsty is always y, Afiny is always yextent1 if it is used
	// Now create the loop values;
	if (xextent0 <= SimW) // No overlap with the eastern side
	{
		// Dstx, Dsty, Asty set by defaults
		//Astx & Afinx are not needed
		Dfinx = xextent0;
		// do we overlap the bottom?
		// Type B & D (overlap bottom, no overlap)
		if (yextent0 > SimH)
		{
			// Type B (overlap bottom only)
			Dfiny = SimH; // stop at the end
			Afiny = yextent1; // the overlap with the top
		}
		else Dfiny = yextent0;
	}
	else
	{
		// Type A & C overlap bottom & eastern edgdes
		if (yextent0 > SimH)
		{
			// relies on the default start for Asty, Astx, Dstx, Dsty
			Afinx = xextent1;
			Afiny = yextent1;
			Dfinx = SimW; // Stop at the end
			Dfiny = SimH;
		}
		else
		{
			// Type C overlap left edge only
			// Astx & Afiny are not needed here
			//Astx, Dstx, Dsty set by default
			Afinx = xextent1;
			Dfinx = SimW; // Stop at the end
			Dfiny = yextent0;
		}
	}
	Vole_Base* Ap;

	// A Loop
	for (int i = 0; i < Afinx; i++)
	{
		for (int j = Asty; j < Afiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (!Ap->SupplySex()) // female
				{
					if (Ap->SupplyTerritorial()) { Females++; }
				}
			}
		}
		// C Loop
		for (int j = y; j < Dfiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (!Ap->SupplySex()) // female
				{
					if (Ap->SupplyTerritorial()) { Females++; }
				}
			}
		}
	}
	// D Loop
	for (int i = x; i < Dfinx; i++)
	{
		for (int j = y; j < Dfiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (!Ap->SupplySex()) // female
				{
					if (Ap->SupplyTerritorial()) { Females++; }
				}
			}
		}
		// B Loop
		for (int j = 0; j < Afiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (!Ap->SupplySex()) // female
				{
					if (Ap->SupplyTerritorial()) { Females++; }
				}
			}
		}
	}
	// End of search algorithm
	m_VoleMap->SetMapValue(p_x, p_y, caller);
	return Females; // returns the number of adult females
}
//---------------------------------------------------------------------------

int Vole_Population_Manager::SupplyInOlderTerr(const unsigned p_x, const unsigned p_y, const unsigned p_Age,
                                               const unsigned p_Range) const {
	/** returns -1 if a male has p_x,p_y in his territory and is older than p_Age
	*/

	// Before checking the map remove ourselves so we don't count
	//PointerInt c=VoleMap->GetMapValue(p_x,p_y);
	//TAnimal* caller=(TAnimal*) c;
	int state = 0;
	TAnimal* caller = m_VoleMap->GetMapValue(p_x, p_y);
	m_VoleMap->ClearMapValue(p_x, p_y);
	int x = p_x - p_Range;
	if (x < 0) x += SimW; // ensure we start in the simulation area!
	int y = p_y - p_Range;
	if (y < 0) y += SimH;
	const int range_x = p_Range * 2;
	const int range_y = p_Range * 2;

	// create the extent variables
	const int xextent0 = x + range_x;
	const int yextent0 = y + range_y;
	const int xextent1 = x + range_x - SimW;
	const int yextent1 = y + range_y - SimH;
	// Create the looping variables needed
	// Create the looping variables needed
	int Dfinx;
	int Dfiny;
	int Afinx = 0; //unless the finx values for A are changed this stop
	int Afiny = 0; //the loop from executing
	const int Asty = 0;
	// int Astx,Dstx,Dtsy are always default so variables not used from them
	// NB Astx, Asty and Dstx are always 0, 0 & x respectively
	// Dsty is always y, Afiny is always yextent1 if it is used
	// Now create the loop values;
	if (xextent0 <= SimW) // No overlap with the eastern side
	{
		// Dstx, Dsty, Asty set by defaults
		//Astx & Afinx are not needed
		Dfinx = xextent0;
		// do we overlap the bottom?
		// Type B & D (overlap bottom, no overlap)
		if (yextent0 > SimH)
		{
			// Type B (overlap bottom only)
			Dfiny = SimH; // stop at the end
			Afiny = yextent1; // the overlap with the top
		}
		else Dfiny = yextent0;
	}
	else
	{
		// Type A & C overlap bottom & eastern edgdes
		if (yextent0 > SimH)
		{
			// relies on the default start for Asty, Astx, Dstx, Dsty
			Afinx = xextent1;
			Afiny = yextent1;
			Dfinx = SimW; // Stop at the end
			Dfiny = SimH;
		}
		else
		{
			// Type C overlap left edge only
			// Astx & Afiny are not needed here
			//Astx, Dstx, Dsty set by default
			Afinx = xextent1;
			Dfinx = SimW; // Stop at the end
			Dfiny = yextent0;
		}
	}
	Vole_Male* AMale;
	Vole_Base* Ap;

	// A Loop
	for (int i = 0; i < Afinx; i++)
	{
		for (int j = Asty; j < Afiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (Ap->SupplySex()) // male
				{
					AMale = static_cast<Vole_Male*>(Ap);
					if (AMale->SupplyTerritorial())
					{
						if (AMale->SupplyAge() >= p_Age)
						{
							if (AMale->SupplyAge() - g_sigAgeDiff > p_Age || g_rand_uni_fnc() < 0.5)
							{
								m_VoleMap->SetMapValue(p_x, p_y, caller);
								state = -1;
								return state; // No Good
							}
						}
					}
				}
				else state++;
			}
		}
		// C Loop
		for (int j = y; j < Dfiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (Ap->SupplySex()) // male
				{
					AMale = static_cast<Vole_Male*>(Ap);
					if (AMale->SupplyTerritorial())
					{
						if (AMale->SupplyAge() >= p_Age)
						{
							if (AMale->SupplyAge() - g_sigAgeDiff > p_Age || g_rand_uni_fnc() < 0.5)
							{
								m_VoleMap->SetMapValue(p_x, p_y, caller);
								state = -1;
								return state; // No Good
							}
						}
					}
				}
				else state++;
			}
		}
	}
	// D Loop
	for (int i = x; i < Dfinx; i++)
	{
		for (int j = y; j < Dfiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (Ap->SupplySex()) // male
				{
					AMale = static_cast<Vole_Male*>(Ap);
					if (AMale->SupplyTerritorial())
					{
						if (AMale->SupplyAge() >= p_Age)
						{
							if (AMale->SupplyAge() - g_sigAgeDiff > p_Age || g_rand_uni_fnc() < 0.5)
							{
								m_VoleMap->SetMapValue(p_x, p_y, caller);
								state = -1;
								return state; // No Good
							}
						}
					}
				}
				else state++;
			}
		}
		// B Loop
		for (int j = 0; j < Afiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap)
			{
				if (Ap->SupplySex()) // male
				{
					AMale = static_cast<Vole_Male*>(Ap);
					if (AMale->SupplyTerritorial())
					{
						if (AMale->SupplyAge() >= p_Age)
						{
							if (AMale->SupplyAge() - g_sigAgeDiff > p_Age || g_rand_uni_fnc() < 0.5)
							{
								m_VoleMap->SetMapValue(p_x, p_y, caller);
								state = -1;
								return state; // No Good
							}
						}
					}
				}
				else state++;
			}
		}
	}
	// End of search algorithm
	m_VoleMap->SetMapValue(p_x, p_y, caller);
	return state; // No Males it return the number of adult females
}
//---------------------------------------------------------------------------

bool Vole_Population_Manager::InSquare(const int p_x, const int p_y, const int p_sqx, const int p_sqy,
                                       const int p_range) const {
	/** Determines if p_x,p_y is the in the square denoted by p_sqx,p_sqy + p_range */

	const int x_extent = p_sqx + p_range;
	const int y_extent = p_sqy + p_range;
	if (x_extent >= SimW)
	{
		if (y_extent >= SimH) // overlaps TR corner of sim area
		// Must test four rectangles
		/*
		  1) p_sqx to <m_SimW, p_sq_y to m_SimW (TR)
		  2) 0 to x_extent-m_SimW p_sq_y to m_SimH (TL)
		  3) 0 to x_extent-m_SimW, 0 to y_extent-m_SimH (BL)
		  4) p_sq_x to <m_SimW, 0 to y_extent-m_SimH (BR)
		*/
		{
			// 1 Top right square (limited by SimAreaHeight & SimAreaWidth
			if (p_x >= p_sqx && p_y >= p_sqy) return true;
			// 2 Top Left Square (limited by 0,SimAreaHeight)
			if (p_x < x_extent - SimW && p_y > p_sqy) return true;
			// 3 Bottom Left square (limited by 0,0)
			if (p_x < x_extent - SimW && p_y < y_extent - SimH) return true;
			// Bottom Right square (limited by SimAreaWidth,0)
			if (p_x >= p_sqx && p_y < y_extent - SimH) return true;
		}
		else // Overlaps the west edge of the sim area
		{
			if (p_y >= p_sqy && p_y < y_extent)
			{
				// y is in square
				if (p_x >= p_sqx) return true;
				if (p_x < x_extent - SimW) return true;
			}
		}
	}
	else
	{
		if (y_extent >= SimH) // overlaps top of simulation area
		{
			if (p_x >= p_sqx && p_x < x_extent)
			{
				// x is OK
				if (p_y >= p_sqy) return true;
				if (p_y < y_extent - SimH) return true;
			}
		}
		else // square does not overlap end of simulation area
		{
			if (p_x >= p_sqx && p_x < x_extent && p_y >= p_sqy && p_y < y_extent) return true;
		}
	}
	return false; // not in square
}
//---------------------------------------------------------------------------

vector<Vole_Base*>* Vole_Population_Manager::SupplyVoleList(const unsigned p_x, const unsigned p_y,
                                                            const unsigned p_Range) const {
	/** returns a list of all voles in p_x,p_y, p_range square */

	const auto vbl = new vector<Vole_Base*>;
	// This is reset when the result is known
	int x = p_x;
	if (x < 0) x += SimW; // ensure we start in the simulation area!
	int y = p_y;
	if (y < 0) y += SimH;
	const int range_x = p_Range;
	const int range_y = p_Range;
	// create the extent variables
	const int xextent0 = x + range_x;
	const int yextent0 = y + range_y;
	const int xextent1 = x + range_x - SimW;
	const int yextent1 = y + range_y - SimH;
	// Create the looping variables needed
	// Create the looping variables needed
	int Dfinx;
	int Dfiny;
	int Afinx = 0; //unless the finx values for A are changed this stop
	int Afiny = 0; //the loop from executing
	const int Asty = 0;
	// int Astx,Dstx,Dtsy are always default so variables not used from them
	// NB Astx, Asty and Dstx are always 0, 0 & x respectively
	// Dsty is always y, Afiny is always yextent1 if it is used
	// Now create the loop values;
	if (xextent0 <= SimW) // No overlap with the eastern side
	{
		// Dstx, Dsty, Asty set by defaults
		//Astx & Afinx are not needed
		Dfinx = xextent0;
		// do we overlap the bottom?
		// Type B & D (overlap bottom, no overlap)
		if (yextent0 > SimH)
		{
			// Type B (overlap bottom only)
			Dfiny = SimH; // stop at the end
			Afiny = yextent1; // the overlap with the top
		}
		else Dfiny = yextent0;
	}
	else
	{
		// Type A & C overlap bottom & eastern edgdes
		if (yextent0 > SimH)
		{
			// relies on the default start for Asty, Astx, Dstx, Dsty
			Afinx = xextent1;
			Afiny = yextent1;
			Dfinx = SimW; // Stop at the end
			Dfiny = SimH;
		}
		else
		{
			// Type C overlap left edge only
			// Astx & Afiny are not needed here
			//Astx, Dstx, Dsty set by default
			Afinx = xextent1;
			Dfinx = SimW; // Stop at the end
			Dfiny = yextent0;
		}
	}
	Vole_Base* Ap;

	// A Loop
	for (int i = 0; i < Afinx; i++)
	{
		for (int j = Asty; j < Afiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap) { vbl->push_back(Ap); }
		}
		// C Loop
		for (int j = y; j < Dfiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap) { vbl->push_back(Ap); }
		}
	}
	// D Loop
	for (int i = x; i < Dfinx; i++)
	{
		for (int j = y; j < Dfiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap) { vbl->push_back(Ap); }
		}
		// B Loop
		for (int j = 0; j < Afiny; j++)
		{
			Ap = static_cast<Vole_Base*>(m_VoleMap->GetMapValue(i, j));
			if (Ap) { vbl->push_back(Ap); }
		}
	}
	return vbl;
}
//---------------------------------------------------------------------------

void Vole_Population_Manager::SendMessage(const TTypeOfVoleMessage p_message_type, const unsigned p_x,
                                          const unsigned p_y, const unsigned p_range, const bool p_sex,
                                          double /*p_Weight*/ /*, unsigned  p_IDNo */) {
	/** Passes a message to recipients. In this case the only one we use is infanticide sent to all females within an area */
	if (p_message_type == tovm_Infanticide)
	{
		if (p_sex == false) // female
		{
			for (unsigned i = 0; i < GetLiveArraySize(vob_Female); i++)
			{
				const auto AFemale = static_cast<Vole_Female*>(SupplyAnimalPtr(vob_Female, i));
				{
					// is it in the square defined by p_x,p_y & p_range
					const unsigned x = AFemale->SupplyX();
					const unsigned y = AFemale->SupplyY();
					// Need to know if female is within p_range of p_x, p_y

					unsigned dx = abs(static_cast<int>(p_x - x));
					unsigned dy = abs(static_cast<int>(p_y - y));
					if (dx > SimWH) dx = SimW - dx;
					if (dy > SimHH) dy = SimH - dy;
					if (dx <= p_range && dy <= p_range)
					{
						//				 if (p_Weight>AFemale->SupplyWeight())
						//				 {
						AFemale->OnInfanticideAttempt();
						//InfanticideOutput(p_ActualAge,p_IDNo);
						//				 }
					}
				}
			}
		}
		else
		{
			m_TheLandscape->Warn("Vole_Population_Manager::SendMessage Error", "Wrong sex specified for infanticide");
		}
	}
	else { m_TheLandscape->Warn("Vole_Population_Manager::SendMessage Error", "Unknown message"); }
}

//---------------------------------------------------------------------------

/**
Creates 'number' of vole objects of the type ob_type using 'as' for the base data
*/
void Vole_Population_Manager::CreateObjects(const VoleObject a_obType, TAnimal* /* a_pvo */, struct_Vole_Adult* a_as,
	const int a_number, const string& filename_AlleleInput, bool isOffspring) {
	for (int i = 0; i < a_number; i++)
	{
		if (a_obType == vob_JuvenileMale)
		{
			if(isOffspring==true){ //just to make sure we don't override its inherited genes, but maybe unnecessary as inits are done in CreateObjectsInit anyway...
				auto new_JMale = new Vole_JuvenileMale(a_as, ""); //pass an empty string as allele freq file to not override inherited genes
				PushIndividual(a_obType, new_JMale);
				m_VoleMap->SetMapValue(a_as->x, a_as->y, new_JMale);
				IncLiveArraySize(a_obType);
			}else{
				auto new_JMale = new Vole_JuvenileMale(a_as, filename_AlleleInput); //If its not an offspring, it gets to make its own genes, so gets passed a file with allele freqs
				PushIndividual(a_obType, new_JMale);
				m_VoleMap->SetMapValue(a_as->x, a_as->y, new_JMale);
				IncLiveArraySize(a_obType);
			}

#ifdef __VOLEPESTICIDEON
			if (as->m_dflag)
			{  // Chromo 1. direct effect
				new_JMale->SetDirectFlag();
				new_JMale->SetGeneticFlag();
				if (g_rand_uni() >= m_f1sterilitychance) new_JMale->SetFertile(true);
				else new_JMale->SetFertile(false); //
				new_JMale->SetPesticideInfluenced2(true);
			}
			else
			{
				if (as->m_gflag)
				{  // Chromo 0, genetic effect
					new_JMale->UnsetDirectFlag();
					new_JMale->SetGeneticFlag();
					if (g_rand_uni() >= m_geneticsterilitychance) new_JMale->SetFertile(true); //100 means all are sterile, 0 means none are
					else new_JMale->SetFertile(false);
					new_JMale->SetPesticideInfluenced2(true);
				}
				else
				{
					new_JMale->UnsetGeneticFlag();
					new_JMale->UnsetDirectFlag();
					new_JMale->SetFertile(true);
				}
			}
#endif
		}
		
		if (a_obType == vob_JuvenileFemale)
		{
			if (isOffspring==true){//just to make sure we don't override its inherited genes, but maybe unnecessary as inits are done in CreateObjectsInit anyway...
				auto new_JFemale = new Vole_JuvenileFemale(a_as, "");  //pass an empty string as allele freq file to not override inherited genes
				PushIndividual(a_obType, new_JFemale);
				m_VoleMap->SetMapValue(a_as->x, a_as->y, new_JFemale);
				IncLiveArraySize(a_obType);

			}else{
				auto new_JFemale = new Vole_JuvenileFemale(a_as, filename_AlleleInput); //If its not an offspring, it gets to make its own genes, so gets passed a file with allele freqs
				PushIndividual(a_obType, new_JFemale);
				m_VoleMap->SetMapValue(a_as->x, a_as->y, new_JFemale);
				IncLiveArraySize(a_obType);
			}
			


#ifdef __VOLEPESTICIDEON
			new_JFemale->SetMaturityDelay(as->misc_use);
#endif
		}
		if (a_obType == vob_Male)
		{
			if (isOffspring==true){ //just to make sure we don't override its inherited genes, but maybe unnecessary as inits are done in CreateObjectsInit anyway...
				auto new_Male = new Vole_Male(a_as, ""); //pass an empty string as allele freq file to not override inherited genes
				PushIndividual(a_obType, new_Male);
				m_VoleMap->SetMapValue(a_as->x, a_as->y, new_Male);
				IncLiveArraySize(a_obType);
			}else{
				auto new_Male = new Vole_Male(a_as, filename_AlleleInput); //If its not an offspring, it gets to make its own genes, so gets passed a file with allele freqs
				PushIndividual(a_obType, new_Male);
				m_VoleMap->SetMapValue(a_as->x, a_as->y, new_Male);
				IncLiveArraySize(a_obType);
			}
			

#ifdef __VOLEPESTICIDEON
			if (as->m_dflag)
			{  // Chromo 1. direct effect
				new_Male->SetDirectFlag();
				new_Male->SetGeneticFlag();
				if (g_rand_uni() >= m_f1sterilitychance) new_Male->SetFertile(true);
				else new_Male->SetFertile(false); //
				new_Male->SetPesticideInfluenced2(true);
			}
			else
			{
				if (as->m_gflag)
				{  // Chromo 0, genetic effect
					new_Male->UnsetDirectFlag();
					new_Male->SetGeneticFlag();
					if (g_rand_uni() >= m_geneticsterilitychance) new_Male->SetFertile(true); //100 means all are sterile, 0 means none are
					else new_Male->SetFertile(false);
					new_Male->SetPesticideInfluenced2(true);
				}
				else
				{
					new_Male->UnsetGeneticFlag();
					new_Male->UnsetDirectFlag();
					new_Male->SetFertile(true);
				}
			}
#endif
		}
		if (a_obType == vob_Female)
		{
			if (isOffspring==true){ //just to make sure we don't override its inherited genes, but maybe unnecessary as inits are done in CreateObjectsInit anyway...
				const auto new_Female = new Vole_Female(a_as, ""); //pass an empty string as allele freq file to not override inherited genes
				PushIndividual(a_obType, new_Female);
				m_VoleMap->SetMapValue(a_as->x, a_as->y, new_Female);
				IncLiveArraySize(a_obType);
			}else{
				const auto new_Female = new Vole_Female(a_as, filename_AlleleInput); //If its not an offspring, it gets to make its own genes, so gets passed a file with allele freqs
				PushIndividual(a_obType, new_Female);
				m_VoleMap->SetMapValue(a_as->x, a_as->y, new_Female);
				IncLiveArraySize(a_obType);
			}
		}
	}
}
//---------------------------------------------------------------------------

/**
Creates 'number' of vole objects of the type ob_type using 'as' for the base data for use at the beginning of a simulation. \n
A number of the attributes are set at defaults or randomised
*/
void Vole_Population_Manager::CreateObjectsInit(const VoleObject a_obType, TAnimal* a_pvo,
                                                struct_Vole_Adult* a_voleStruct, const int a_number, const string& filename_AlleleInput, int id) {
	struct_Vole_Adult* as = nullptr;
	const int Year = m_TheLandscape->SupplyYearNumber();

	for (int i = 0; i < a_number; i++)
	{
		if (a_obType == vob_Male)
		{
			as = a_voleStruct;
			as->FatherId = 10000000;
			as->MotherId = 10000000;
			as->FatherStateAtBirth = 2;
			as->new_Genes.EraseGenomes(); // probably unnecessary, but make sure a new vole has no prior genes
			as->new_Mates_Genes.EraseGenomes(); // probably unnecessary, but make sure a new vole has no prior mate's genes
			as->GenerationCount=0; //an initialized vole is generation 0
			//as->new_Genes.InitializeTestingGenomesMale(); //will give all males a genome only with 0 and one only with 1 if activated, for testing purposes
			const auto new_Male = new Vole_Male(as, filename_AlleleInput, id); //i is passed on as id because it's used as the unique paternal+maternal lineages for the new male.

			new_Male->SetWeight(40.0);
			new_Male->Setm_Mature();
			new_Male->Set_Age(g_random_fnc(500));
			new_Male->Set_BirthYear(Year);
			new_Male->SetFertile(true);

			PushIndividual(a_obType, new_Male);
		}
		if (a_obType == vob_Female)
		{
			as = a_voleStruct;
			as->FatherId = 10000000;
			as->MotherId = 10000000;
			as->FatherStateAtBirth = 2;
			as->new_Genes.EraseGenomes(); // probably unnecessary, but make sure a new vole has no prior genes
			as->new_Mates_Genes.EraseGenomes();// probably unnecessary, but make sure a new vole has no prior mate's genes
			as->GenerationCount=0; //an initialized vole is generation 0
			//as->new_Genes.InitializeTestingGenomesFemale(); //will give all females a genome only with 2 and one only with 3 if activated for testing purposes
			const auto new_Female = new Vole_Female(as, filename_AlleleInput, id); //i is passed on as id because it's used as the unique paternal+maternal lineages for the new female
			new_Female->SetWeight(40.0);
			new_Female->Setm_Mature();
			new_Female->Set_Age(1); // (random(500))
			new_Female->Set_BirthYear(Year);
			//cout << "\n in create objects init after making vole, here are genes:";
			//new_Female->new_Genes.PrintGenomes();
			PushIndividual(a_obType, new_Female);
		}
		IncLiveArraySize(a_obType);
	}
}
//---------------------------------------------------------------------------

void Vole_Population_Manager::Catastrophe(void) {
	/** This version simply alters populations after 1st January - it is very dangerous to
	add individuals in many of the models so beware!!!! */

	// First do we have the right day?
	const int today = m_TheLandscape->SupplyDayInYear();
	if (today != 1) return;
	// First do we have the right year?
	const int year = m_TheLandscape->SupplyYearNumber() - cfg_CatastropheEventStartYear.value();
	if (year % cfg_PmEventfrequency.value() != 0) return;

	// assignment to get rid of warning
	// Now if the % decrease is higher or lower than 100 we need to do different things
	int esize = cfg_PmEventsize.value();
	if (esize < 100)
	{
		unsigned size2 = GetLiveArraySize(vob_Male);
		for (unsigned j = 0; j < size2; j++)
		{
			if (g_random_fnc(100) > esize)
			{
				const auto mv = dynamic_cast<Vole_Male*>(SupplyAnimalPtr(vob_Male, j));
				mv->CurrentVState = tovs_MDying; // Kill it
			}
		}
		size2 = GetLiveArraySize(vob_Female);
		for (unsigned j = 0; j < size2; j++)
		{
			if (g_random_fnc(100) > esize)
			{
				const auto fv = dynamic_cast<Vole_Female*>(SupplyAnimalPtr(vob_Female, j));
				fv->CurrentVState = tovs_FDying; // Kill it
			}
		}
	}
	else if (esize > 100)
	{
		Vole_Base* VB = nullptr;
		// This is a tricky thing to do because we need to duplicate voles, but dare not mess
		// mate pointers etc up.
		// This also requires a copy method in the target vole
		// esize also needs translating  120 = 20%, 200 = 100%
		if (esize < 200)
		{
			esize -= 100;
			for (int i = vob_JuvenileMale; i <= static_cast<int>(vob_Female); i++)
			{
				const unsigned size2 = GetLiveArraySize(i);
				for (unsigned j = 0; j < size2; j++)
				{
					if (g_random_fnc(100) < esize)
					{
						VB = dynamic_cast<Vole_Base*>(SupplyAnimalPtr(i, j));
						VB->CopyMyself(static_cast<VoleObject>(i)); // Duplicate it
					}
				}
			}
		}
		else
		{
			esize -= 100;
			esize /= 100; // this will throw away fractional parts so will get 1, 2, 3  from 200, 350 400
			for (int i = vob_JuvenileMale; i <= static_cast<int>(vob_Female); i++)
			{
				const unsigned size2 = GetLiveArraySize(i);
				for (unsigned j = 0; j < size2; j++)
				{
					for (int e = 0; e < esize; e++)
					{
						VB = dynamic_cast<Vole_Base*>(SupplyAnimalPtr(i, j));
						VB->CopyMyself(static_cast<VoleObject>(i)); // Duplicate it
					}
				}
			}
		}
	}
	else return; // No change so do nothing
}
//-----------------------------------------------------------------------------

void Vole_Population_Manager::LandscapeQuadrantOutputProbe(const int a_day) {
	/** Output file facility added in January 2013 */
	TAnimal* VB;
	int gridcount[256];
	/**
	* This output splits the landscape up into 16x16 cells which gives landscapes of 625x625m grids for 10x10 km. <br>
	* All voles are counted in the landspape and assigned to a grid. The sum of all individials in the grid is then output to the output file separated by tabs and prefixed by the global simulation day. <br>
	* NB this only works for square landscapes.
	*/
	const int width = m_TheLandscape->SupplySimAreaWidth();
	for (int i = 0; i < 256; i++) gridcount[i] = 0;
	const double sqwidth = width / 16.0;
	const unsigned int totalYM = GetLiveArraySize(vob_JuvenileMale);
	for (unsigned j = 0; j < totalYM; j++) //juvenile males
	{
		int x{0}, y{0};
		SupplyLocXY(vob_JuvenileMale, j, x, y);
		const int gx = static_cast<int>(floor(x / sqwidth));
		const int gy = static_cast<int>(floor(y / sqwidth));
		gridcount[gx + gy * 16]++;
	}
	const unsigned int totalYF = GetLiveArraySize(vob_JuvenileFemale);
	for (unsigned j = 0; j < totalYF; j++) //juvenile females
	{
		int x{0}, y{0};
		SupplyLocXY(vob_JuvenileFemale, j, x, y);

		const int gx = static_cast<int>(floor(x / sqwidth));
		const int gy = static_cast<int>(floor(y / sqwidth));
		gridcount[gx + gy * 16]++;
	}
	const unsigned int totalAM = GetLiveArraySize(vob_Male);
	for (unsigned j = 0; j < totalAM; j++) //adult males
	{
		int x{0}, y{0};
		SupplyLocXY(vob_Male, j, x, y);
		const int gx = static_cast<int>(floor(x / sqwidth));
		const int gy = static_cast<int>(floor(y / sqwidth));
		gridcount[gx + gy * 16]++;
	}
	const unsigned int totalAF = GetLiveArraySize(vob_Female);
	for (unsigned j = 0; j < totalAF; j++) //adult females
	{
		int x{0}, y{0};
		SupplyLocXY(vob_Female, j, x, y);
		const int gx = static_cast<int>(floor(x / sqwidth));
		const int gy = static_cast<int>(floor(y / sqwidth));
		gridcount[gx + gy * 16]++;
	}
	// Do the output
	/* Open the output file and append */
	ofstream ofile("VoleLandscapeGridData.txt", ios::app);
	ofile << a_day << '\t' << totalYM << '\t' << totalYF << '\t' << totalAM << '\t' << totalAF << '\t';
	for (int i = 0; i < 255; i++) { ofile << gridcount[i] << '\t'; }
	ofile << gridcount[255] << endl;
	ofile.close();
}

void Vole_Population_Manager::TheAOROutputProbe() { m_AOR_Probe->DoProbe(vob_Female); }
//-----------------------------------------------------------------------------

void Vole_Population_Manager::TheRipleysOutputProbe(ofstream* a_prb) {
	Vole_Female* FS;
	const unsigned int totalF = GetLiveArraySize(vob_Female);
	int x, y;
	const int w = m_TheLandscape->SupplySimAreaWidth();
	const int h = m_TheLandscape->SupplySimAreaWidth();
	*a_prb << 0 << '\t' << w << '\t' << 0 << '\t' << h << '\t' << totalF << endl;
	for (unsigned j = 0; j < totalF; j++) //adult females
	{
		SupplyLocXY(vob_Female, j, x, y);

		*a_prb << x << '\t' << y << endl;
	}
	a_prb->flush();
}
//-----------------------------------------------------------------------------
void Vole_Population_Manager::TheAgeSexLocationProbe() {
	Vole_Base* VB;
	const int year = m_TheLandscape->SupplyYearNumber();
	const int day = m_TheLandscape->SupplyDayInYear();
	int x, y;
	const char sex[4] = {'m', 'f', 'M', 'F'};

	/**
	The code is separated into two loops and calls duplicated. This is because there may be sex-specific calls placed
	here which would be imcompatible with a single data input/output method.\n
	This is called the Really Big Probe because 10yrs of daily output on 10x10km landscape can easily produce a >500MB text file.
	Use carefully. \n
	*/
	// Check all males
	for (int vob = 0; vob <= static_cast<int>(vob_Female); vob++)
	{
		const unsigned int total = GetLiveArraySize(vob);
		for (unsigned j = 0; j < total; j++)
		{
			Population_Manager::SupplyLocXY(vob, j, x, y);

			const int Age = VB->SupplyAge();
			*m_VoleAgeSexLocationFile << year << '\t' << day << '\t' << sex[vob] << '\t' << Age << '\t' << x << '\t' <<
				y << endl;
		}
	}
}
//-----------------------------------------------------------------------------

/**
An output facility. This method can be re-written to provide any data necessary on a daily or annual basis
e.g. for pattern oriented modelling purposes. Outputs can be altered and added to the print statement as necessary.
*/
void Vole_Population_Manager::TheReallyBigOutputProbe() {
	const int year = m_TheLandscape->SupplyYearNumber();
	const int day = m_TheLandscape->SupplyDayInYear();
	int x, y;
	/**
	The code is separated into two loops and calls duplicated. This is because there may be sex-specific calls placed
	here which would be imcompatible with a single data input/output method.\n
	This is called the Really Big Probe because 10yrs of daily output on 10x10km landscape can easily produce a >500MB text file.
	Use carefully. \n
	*/
	// Check all males
	unsigned int totalM = GetLiveArraySize(vob_Male);
	for (unsigned j = 0; j < totalM; j++)
	{
		const auto MV = dynamic_cast<Vole_Male*>(SupplyAnimalPtr(vob_Male, j));
		SupplyLocXY(vob_Male, j, x, y);

		const int poly = m_TheLandscape->SupplyPolyRef(x, y);
		const int ele = m_TheLandscape->SupplyElementType(poly);
		const int vegt = m_TheLandscape->SupplyVegType(poly);
		const int Age = MV->SupplyAge();
		const int Ter = MV->SupplyTerritorial();
		*ReallyBigOutputPrb << year << "\t" << day << "\t" << x << "\t" << y << "\t adM\t" << "0\t" << poly << "\t" <<
			ele << "\t" << vegt << "\t" << Age << "\t" << Ter << "\n";
	}
	totalM = GetLiveArraySize(vob_JuvenileMale);
	for (unsigned j = 0; j < totalM; j++)
	{
		const auto JMV = dynamic_cast<Vole_JuvenileMale*>(SupplyAnimalPtr(vob_JuvenileMale, j));
		SupplyLocXY(vob_JuvenileMale, j, x, y);

		const int poly = m_TheLandscape->SupplyPolyRef(x, y);
		const int ele = m_TheLandscape->SupplyElementType(poly);
		const int vegt = m_TheLandscape->SupplyVegType(poly);
		const int Age = JMV->SupplyAge();
		const int Ter = JMV->SupplyTerritorial();
		*ReallyBigOutputPrb << year << "\t" << day << "\t" << x << "\t" << y << "\t juvM\t" << "0\t" << poly << "\t" <<
			ele << "\t" << vegt << "\t" << Age << "\t" << Ter << "\n";
	}
	// Do the same for females
	unsigned int totalF = GetLiveArraySize(vob_Female);
	for (unsigned j = 0; j < totalF; j++) //adult females
	{
		const auto FV = dynamic_cast<Vole_Female*>(SupplyAnimalPtr(vob_Female, j));
		SupplyLocXY(vob_Female, j, x, y);
		const int poly = m_TheLandscape->SupplyPolyRef(x, y);
		const int ele = m_TheLandscape->SupplyElementType(poly);
		const int vegt = m_TheLandscape->SupplyVegType(poly);
		const int Age = FV->SupplyAge();
		const int Ter = FV->SupplyTerritorial();
		*ReallyBigOutputPrb << year << "\t" << day << "\t" << x << "\t" << y << "\t juvM\t" << "0\t" << poly << "\t" <<
			ele << "\t" << vegt << "\t" << Age << "\t" << Ter << "\n";
	}
	totalF = GetLiveArraySize(vob_JuvenileFemale);
	for (unsigned j = 0; j < totalF; j++) //adult females
	{
		const auto JFV = dynamic_cast<Vole_JuvenileFemale*>(SupplyAnimalPtr(vob_JuvenileFemale, j));
		SupplyLocXY(vob_JuvenileFemale, j, x, y);
		x = JFV->Supply_m_Location_x();
		y = JFV->Supply_m_Location_y();
		const int poly = m_TheLandscape->SupplyPolyRef(x, y);
		const int ele = m_TheLandscape->SupplyElementType(poly);
		const int vegt = m_TheLandscape->SupplyVegType(poly);
		const int Age = JFV->SupplyAge();
		const int Ter = JFV->SupplyTerritorial();
		*ReallyBigOutputPrb << year << "\t" << day << "\t" << x << "\t" << y << "\t juvM\t" << "0\t" << poly << "\t" <<
			ele << "\t" << vegt << "\t" << Age << "\t" << Ter << "\n";
	}
	ReallyBigOutputPrb->flush();
}
//-----------------------------------------------------------------------------

/**
Open the sex ratio probe
*/
bool Vole_Population_Manager::OpenSexRatiosProbe() {
	SexRatiosPrb = fopen(cfg_SexRatiosOutput_filename.value(), "w");
	if (!SexRatiosPrb)
	{
		g_msg->Warn(WARN_FILE, "Population_Manager::OpenSexRatiosProbe(): ""Unable to open probe file",
		            cfg_SexRatiosOutput_filename.value());
		exit(1);
	}
	fprintf(SexRatiosPrb,
	        "Year\tDay\tSubMales\tAdMalesThisYear\tAdMalesLatYear\tSubFemales\tAdFemalesThisYear\tAdFemalesLatYear\tJuvMales\tJuvFemales\tTotalMales\tTotalFemales\n");
	return true;
}
//-----------------------------------------------------------------------------

/**
Close the sex ratio probe
*/
void Vole_Population_Manager::CloseSexRatiosProbe() {
	if (SexRatiosPrb) fclose(SexRatiosPrb);
	SexRatiosPrb = nullptr;
}
//-----------------------------------------------------------------------------

/**
This does the same as the ReallyBigProbe, but digests the data to produce a daily output instead of
individual vole data. This only works for the whole population not a subset.
*/
void Vole_Population_Manager::TheSexRatiosProbe() {
	const int year = m_TheLandscape->SupplyYearNumber();
	const int day = m_TheLandscape->SupplyDayInYear();
	fprintf(SexRatiosPrb, "%d\t%d\t", year, day);
	// Create our counters
	int Ad1 = {0};
	int Ad2 = {0};
	int Sub = {0};
	for (int v = vob_Male; v <= static_cast<int>(vob_Female); v++)
	{
		const unsigned int total = GetLiveArraySize(v);
		for (unsigned j = 0; j < total; j++)
		{
			const auto VB = dynamic_cast<Vole_Base*>(SupplyAnimalPtr(v, j));
			//TTypesOfLandscapeElement tole = VB->SupplyElemType();
			//if ((tole == tole_Orchard) || (tole == tole_NaturalGrassDry) )
			{
				const bool mature = static_cast<int>(VB->SupplyMature());
				const bool born_lastyear = VB->SupplyBornLastYear();
				if (!mature) Sub++;
				else if (born_lastyear) Ad2++;
				else Ad1++;
			}
		}
		fprintf(SexRatiosPrb, "%d\t%d\t%d\t", Sub, Ad1, Ad2);
		Sub = 0;
		Ad1 = 0;
		Ad2 = 0;
	}
	fprintf(SexRatiosPrb, "%d\t%d\t%d\t%d\n", static_cast<int>(GetLiveArraySize(vob_JuvenileMale)),
	        static_cast<int>(GetLiveArraySize(vob_JuvenileFemale)), static_cast<int>(GetLiveArraySize(vob_Male)),
	        GetLiveArraySize(vob_Female));
}
//-----------------------------------------------------------------------------
TrapLineMap::TrapLineMap(const unsigned int a_width, const unsigned int a_height, const unsigned int a_resolution,
                         const char* a_file) : BinaryMapBase(a_width, a_height, a_resolution, 2) { Init(a_file); }

TrapLineMap::~TrapLineMap() { fclose(m_ofile); }

/**
* Clears the map then reads the list of co-ordinates from the text file, saves these in m_TrapList.
* At the same time this creates the trap map by putting 1s at each co-ordinate.
*/
void TrapLineMap::Init(const char* a_inifile) {
	ClearMap();
	FILE* inpfile = fopen(a_inifile, "r");
	if (!inpfile)
	{
		g_msg->Warn(WARN_FILE, " TrapLineMap::Init Unable to open file ", a_inifile);
		exit(1);
	}
	APoint pt;
	unsigned input, x, y;
	fscanf(inpfile, "%d\n", &input);
	m_noTraps = input;
	for (unsigned i = 0; i < m_noTraps; i++)
	{
		fscanf(inpfile, "%d\t%d\n", &x, &y);
		pt.m_x = x;
		pt.m_y = y;
		m_TrapCoords.push_back(pt);
		SetValue(x, y, 1);
	}
	fclose(inpfile);
	// Open the output file
	m_ofile = fopen("VoleTrapCounts.txt", "w");
	if (!m_ofile)
	{
		g_msg->Warn(WARN_FILE, " TrapLineMap::Init Unable to open output file ", "VoleTrapCounts.txt");
		exit(1);
	}
	fprintf(m_ofile, "Day\tx\ty\tElement\tVegType\tSex\tTerritorial\tAge\tBornX\tBornY\tID no\n");
}

void TrapLineMap::Output(const InTrapPosition& a_tp, const int a_day, const int a_sex, const bool a_terr,
                         const int a_age, const int a_bx, const int a_by, const int a_ID) const {
	char terr;
	if (a_terr) terr = 'Y';
	else terr = 'N';
	fprintf(m_ofile, "%d\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\t%d\t%d\n", a_day, a_tp.m_x, a_tp.m_y, a_tp.m_EleType,
	        a_tp.m_VegType, a_sex, terr, a_age, a_bx, a_by, a_ID);
	fflush(m_ofile);
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

VoleSummaryOutput::VoleSummaryOutput(const char* a_filename, Landscape* a_land, const int a_numdataINT,
                                     const int a_numdataDOUBLE) {
	m_landscape = a_land;
	OpenOutput(a_filename);
	m_ndInt = a_numdataINT;
	m_ndDouble = a_numdataDOUBLE;
	ResetData();
}
VoleSummaryOutput::~VoleSummaryOutput() { CloseOutput(); }
void VoleSummaryOutput::OPrint() {
	*m_File << m_landscape->SupplyYearNumber() << '\t' << m_landscape->SupplyDayInYear();
	for (int i = 0; i < m_ndInt; i++) *m_File << '\t' << m_dataI[i];
	for (int i = 0; i < m_ndDouble; i++) *m_File << '\t' << m_dataD[i];
	*m_File << endl;
}
void VoleSummaryOutput::OPrint(const int a_value) { *m_File << a_value << '\t'; }
void VoleSummaryOutput::OPrint(const double a_value) { *m_File << a_value << '\t'; }
void VoleSummaryOutput::OPrint(const char* a_value) { *m_File << a_value << '\t'; }
void VoleSummaryOutput::OPrintEndl() { *m_File << endl; }
void VoleSummaryOutput::OpenOutput(const char* a_filename) { m_File = new ofstream(a_filename, ios::out); }
void VoleSummaryOutput::CloseOutput() const {
	m_File->close();
	delete m_File;
}
void VoleSummaryOutput::ResetData() {
	for (int i = 0; i < m_ndInt; i++) m_dataI[i] = 0;
	for (int i = 0; i < m_ndDouble; i++) m_dataD[i] = 0.0;
}
void VoleSummaryOutput::ChangeData(const int a_data, const int a_value) { m_dataI[a_data] += a_value; }
void VoleSummaryOutput::ChangeData(const int a_data, const double a_value) { m_dataD[a_data] += a_value; }
//-----------------------------------------------------------------------------

/** Used to output vole genetics at user-defined dates*/
void Vole_Population_Manager::GeneticsOutputFile(const unsigned listindex) {
	FILE* vfile = fopen("GeneticsData.txt", "a");
	if (vfile == nullptr)
	{
		m_TheLandscape->Warn("Vole_Population_Manager::GeneticsOutputFile", "Could Not Open GeneticsData.txt File");
		exit(0);
	}
	const int Y = m_TheLandscape->SupplyYearNumber();
	const int M = m_TheLandscape->SupplyMonth();
	const int D = m_TheLandscape->SupplyDayInYear();
	const unsigned size = GetLiveArraySize(listindex);
	const unsigned TotalSize = GetLiveArraySize(vob_Male) + GetLiveArraySize(vob_Female);
	const int w = m_TheLandscape->SupplySimAreaWidth();
	const int h = m_TheLandscape->SupplySimAreaWidth();
	fprintf(vfile, "\n");
	fprintf(vfile, "%d\t %d\t %d\t %d\t %d\t %d\n", 0, w, 0, h, size, TotalSize);
	fprintf(vfile, "\n");
	fprintf(vfile,
	        "%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\n",
	        "Year", "Mo", "Day", "ID", "Sex", "Age", "Birth Yr", "Terr", "Ma", "TotY", "x", "y", "Poly", "Ele", "Vege",
	        "GA", "DA", "Bx", "By", "Bpol", "Bele", "Bveg");
	for (unsigned j = 0; j < size; j++)
	{
		const auto MV = dynamic_cast<Vole_Base*>(SupplyAnimalPtr(vob_Male, j));
		int x{0}, y{0};
		SupplyLocXY(vob_Male, j, x, y);

		const int poly = m_TheLandscape->SupplyPolyRef(x, y);
		const int ele = m_TheLandscape->SupplyElementType(poly);
		const int vegt = m_TheLandscape->SupplyVegType(poly);
		const int GA = MV->GetGeneticFlag();
		const int DA = MV->GetDirectFlag();

		const unsigned Age = MV->SupplyAge();
		const int BY = MV->SupplyBirthYear();
		const int Mature = MV->SupplyMature();
		const int Terr = MV->SupplyTerritorial();
		const int TotYoung = 0;

		const int BornX = MV->SupplyXBorn();
		const int BornY = MV->SupplyYBorn();
		const int BornPoly = MV->SupplyPolyRefBorn();
		const int BornEle = MV->SupplyElemBorn();
		const int BornVegt = MV->SupplyVegBorn();
		const int sex = MV->SupplySex();
		int ID = MV->SupplyIDNo();

		fprintf(vfile, "%d\t", Y); //Year
		fprintf(vfile, "%d\t", M); // month letter
		fprintf(vfile, "%d\t", D);

		fprintf(vfile, "%d\t", ID);
		fprintf(vfile, "%d\t", sex); // sex
		fprintf(vfile, "%u\t", Age);
		fprintf(vfile, "%d\t", BY); //BirthYear
		fprintf(vfile, "%d\t", Terr);
		fprintf(vfile, "%d\t", Mature);
		fprintf(vfile, "%d\t", TotYoung);
		fprintf(vfile, "%d\t", x); // coordinate
		fprintf(vfile, "%d\t", y); // coordinate
		fprintf(vfile, "%d\t", poly); // polygon
		fprintf(vfile, "%d\t", ele); // element type
		fprintf(vfile, "%d\t", vegt); // vegetation type
		fprintf(vfile, "%d\t", GA); // Genetic affected
		fprintf(vfile, "%d\t", DA); // Direct affected
		fprintf(vfile, "%d\t", BornX); //X-coordinat for birth location
		fprintf(vfile, "%d\t", BornY); //Y-coordinat for birth location
		fprintf(vfile, "%d\t", BornPoly);
		fprintf(vfile, "%d\t", BornEle);
		fprintf(vfile, "%d\t", BornVegt);

		for (int i = 0; i < 16; i++)
		{
			for (int jj = 0; jj < 2; jj++)
			{
				const uint32 allele = MV->SupplyMyAllele(i, jj);
				fprintf(vfile, "%u\t", allele);
			}
		}
		fprintf(vfile, "\n");
		ID++;
	}
	fclose(vfile);
}
//-----------------------------------------------------------------------------

/**
Output recording the genetic information for each individual in the population of listindex type
*/
void Vole_Population_Manager::GeneticsResultsOutput(FILE* ofile, const unsigned listindex) {
	const char month[13] = {'0', 'J', 'F', 'M', 'A', 'm', 'j', 'u', 'a', 'S', 'O', 'N', 'D'};
	int ID = 0;
	const int Y = m_TheLandscape->SupplyYear();
	const int M = m_TheLandscape->SupplyMonth();
	const unsigned size = GetLiveArraySize(listindex);
	if (size > 0)
	{
		if (listindex == 0)
		{
			for (unsigned j = 0; j < size; j++)
			{
				const auto MV = dynamic_cast<Vole_Male*>(SupplyAnimalPtr(vob_Male, j));
				int x{0}, y{0};
				SupplyLocXY(vob_Male, j, x, y);

				// Outputs genetic results
				// Year \t Month \t sex \t individual code \t XX\tXX\tXX\tXX\tXX\tXX\tXX\tXXn .... 32 loci x 2 alleles
				fprintf(ofile, "%d\t", ID);
				fprintf(ofile, "M");
				fprintf(ofile, "%c", month[M]);
				fprintf(ofile, "%d\t", Y);
				fprintf(ofile, "%d\t", x);
				fprintf(ofile, "%d\t", y);
				for (int g = 0; g < 32; g++)
				{
					int allele = 1 + MV->SupplyAllele(g, 0);
					fprintf(ofile, "%d\t", allele);
					allele = 1 + MV->SupplyAllele(g, 1);
					fprintf(ofile, "%d\t", allele);
				}
				fprintf(ofile, "\n");
				ID++;
			}
		}
		else
		{
			for (unsigned j = 0; j < size; j++)
			{
				const auto FV = dynamic_cast<Vole_Female*>(SupplyAnimalPtr(vob_Female, j));
				int x{0}, y{0};
				SupplyLocXY(vob_Female, j, x, y);

				// Outputs genetic results
				// Year \t Month \t sex \t individual code \t XX\tXX\tXX\tXX\tXX\tXX\tXX\tXXn .... 32 loci x 2 alleles
				fprintf(ofile, "%d\t", ID);
				fprintf(ofile, "F");
				fprintf(ofile, "%c", month[M]);
				fprintf(ofile, "%d\t", Y);
				fprintf(ofile, "%d\t", x);
				fprintf(ofile, "%d\t", y);
				for (int g = 0; g < 32; g++)
				{
					int allele = 1 + FV->SupplyAllele(g, 0);
					fprintf(ofile, "%d\t", allele);
					allele = 1 + FV->SupplyAllele(g, 1);
					fprintf(ofile, "%d\t", allele);
				}
				fprintf(ofile, "\n");
				ID++;
			}
		}
	}
}
//-----------------------------------------------------------------------------

void Vole_Population_Manager::OpenResistanceOutput() {
	m_VoleResistanceOutputFile = new ofstream("VoleResistanceOutput.txt", ios::out);
	if (!m_VoleResistanceOutputFile->is_open())
	{
		m_TheLandscape->Warn("Vole_Population_Manager::GeneticsOutputFile",
		                     "Could Not Open VoleResistanceOutput.txt File");
		exit(0);
	}
	*m_VoleResistanceOutputFile << "year" << '\t' << "day" << '\t' << "Resistgene" << '\t' << "Neutralgene" << '\t' <<
		"Frequency" << '\t' << "Population Size" << endl;
}
//-----------------------------------------------------------------------------

void Vole_Population_Manager::CloseResistanceOutput() const {
	if (!m_VoleResistanceOutputFile->is_open()) m_VoleResistanceOutputFile->close();
	delete m_VoleResistanceOutputFile;
}
//-----------------------------------------------------------------------------

void Vole_Population_Manager::ResistanceOutput() {
	int Rgene = 0;
	int Ngene = 0;
	const int Y = m_TheLandscape->SupplyYearNumber();
	const int D = m_TheLandscape->SupplyDayInYear();
	for (int listindex = 0; listindex < vob_foobar; listindex++)
	{
		const unsigned size = GetLiveArraySize(listindex);
		for (unsigned j = 0; j < size; j++)
		{
			const auto new_VB = dynamic_cast<Vole_Base*>(SupplyAnimalPtr(listindex, j));
			const uint32 allele1 = new_VB->SupplyMyAllele(3, 0);
			const uint32 allele2 = new_VB->SupplyMyAllele(3, 1);
			if (cfg_ResistanceDominant.value())
			{
				if (allele1 == 1 || allele2 == 1) Rgene++;
				else Ngene++;
			}
			else
			{
				if (allele1 == 1 && allele2 == 1) Rgene++;
				else Ngene++;
			}
		}
	}
	const double freq = static_cast<double>(Rgene) / static_cast<double>(Rgene + Ngene);
	*m_VoleResistanceOutputFile << Y << '\t' << D << '\t' << Rgene << '\t' << Ngene << '\t' << freq << '\t' << Rgene +
		Ngene << endl;
}
//---------------------------------------------------------------------------
// I start here (anastasia) 
//This only works for 4 quadrants (a 2*2 Quadrant grid) :)) But the function outputs allele frequencies, heterozyg.,FIS,FST, and positions of voles categorized by the 4 quadrants
void Vole_Population_Manager::Four_QuadrantBasedGeneticOutput(const int a_day, const int year, const string landscape_info) {
    TAnimal* VB;

	string DivisionDirection;
	DivisionDirection=landscape_info.substr(33, 2); //adjust depending on the naming of the landscape

	int GridCellsNr=9; // Makes 9*9 grid cells, but that is reduced to 4 quadrants later, only put odd numbers!
	int DivisionPoint=GridCellsNr/2; // middle point (will be uses as index for buffer zone)
	int ExcludedVolesCount=0; // To count voles in buffer zones
	int VoleSum=0; 

    const int width = m_TheLandscape->SupplySimAreaWidth(); // get landscape width
	const double sqwidth = width/GridCellsNr; // length of a gridcell

    //making a new file for GeneticsAndPosition for every year, so these titles need to be sent to the new file every year
	string genetics_out_file = "GeneticsAndPosition" + to_string(year) + ".txt";
	ofstream genetics_out(genetics_out_file);
	genetics_out << "year,day,type,nowx,nowy,nowquadrant,division_direction,"; //making the title for the MANY columns
	for (int z = 0; z < ChromosomeCounting * LocusCounting; ++z) {
        // Append "locusZ" to the output stream, the title for the column
        genetics_out << "locus" << z << ",";
    }

    //these files are for appending every year, which is why we only write the "header" in the first year, not every year
    ofstream position_out("BirthPositionAndPositionNow.txt",ios::app);
    ofstream genetics_count_quadrants("GeneticsInQuadrants.txt", ios::app);
    ofstream FSTout("FST_output_matrix_style.txt", ios::app);
	ofstream FSToutForAnalysis("FST_output_dataset_style.txt", ios::app);
	ofstream FISoutForAnalysis("FIS_output_dataset_style.txt", ios::app);
	ofstream HetOut("Heterozyg_output.txt",ios::app);

	if (year==1){ //writing headers/titles in the first year
		FSToutForAnalysis << "year,poprow,popcol,FST,division_direction \n";
		FISoutForAnalysis << "year,pop,FIS,division_direction \n";
        position_out << "year,day,type,birthx,birthy,nowx,nowy,birthquadrant,nowquadrant,division_direction\n";
        genetics_count_quadrants <<"year,day,quadrant,chromosome,locus,base,division_direction\n";
		HetOut << "\nyear,pop,avr_het_exp,avr_het";
	}

    //Structures to hold Allele counts and heterozygosity counts
	vector<vector<vector<vector<double>>>> AlleleCountInQuadrants((4), vector<vector<vector<double>>>(ChromosomeCounting, vector<vector<double>>(LocusCounting, vector<double>(DiffBaseCounting,0))));
	vector<vector<vector<double>>> HetCount((4), vector<vector<double>>(ChromosomeCounting, vector<double>(LocusCounting,0.0)));
	vector<vector<double>> HetAndExpHetPopwise((4), vector<double>(2,-2.0));
	
    //Going through the juvenile males to get their postions, alleles and heterozygosity
	const unsigned int totalYM = GetLiveArraySize(vob_JuvenileMale);
	VoleSum+=totalYM;
	for (unsigned j = 0; j < totalYM; j++){	
		auto JM = dynamic_cast<Vole_JuvenileMale*>(SupplyAnimalPtr(vob_JuvenileMale,j));
		int nowx= JM -> Supply_m_Location_x(); //getting his position
		int nowy= JM -> Supply_m_Location_y();

		int gx = static_cast<int>(floor(nowx / sqwidth)); //getting the small quadrant he is in
		int gy = static_cast<int>(floor(nowy / sqwidth));

		int GridCellIndex; //let's put him in a big quadrant now, and exclude any voles in the buffer zones.
		GridCellIndex=Four_AssignQuadrant(gx,gy,DivisionPoint);
		if (GridCellIndex==100){
			ExcludedVolesCount++;
			continue;
		}
        int nowquadrant=GridCellIndex;

        int birthx= JM->xborn; //now let's find out where he was born
		int birthy= JM->yborn;

		int birthquadrant;

        int birthx_q = static_cast<int>(floor(birthx / sqwidth)); //put his birthplace in the small quadrant
		int birthy_q = static_cast<int>(floor(birthy / sqwidth));
		birthquadrant=Four_AssignQuadrant(birthx_q,birthy_q,DivisionPoint);
		if (birthquadrant==100){
			ExcludedVolesCount++;
			continue;
		}

        //writing to position file where he is and where he was born and what type he is for the file about birthplace and current postion
        position_out << year << "," << a_day << "," << "JuvenileMale" <<"," << birthx << "," << birthy << "," << nowx << "," << nowy << "," << birthquadrant << "," << nowquadrant << "," << DivisionDirection << "\n";

		//now let's look at his genes
		vector<vector<int>> Genome0= JM->new_Genes.Genome0; //get both of his genomes (because he's a diploid boy ;) )
		vector<vector<int>> Genome1= JM->new_Genes.Genome1;

        //the first part of the gene printout (the rest will be "filled out in the following loop")
		genetics_out <<"\n" << year <<"," << a_day << "," << "JuvenileMale" << "," << nowx << ","<< nowy << "," << nowquadrant << "," <<DivisionDirection << ","; 

		for (int a = 0; a < ChromosomeCounting; a++){
			for (int b = 0; b < LocusCounting; b++){
				//"collecting" bases for writing genome to genetics_out file.
				int base0=Genome0[a][b];
				int base1=Genome1[a][b];
				//counting heterozygosity in the subpop in the quadrant
				if (base0!=base1){
					HetCount[GridCellIndex][a][b]+=1;
				}
                //the following is for counting allele frequencies in the quadrant
				AlleleCountInQuadrants[GridCellIndex][a][b][base0]+=1;
				AlleleCountInQuadrants[GridCellIndex][a][b][base1]+=1;
				//the last part of the genetics_out printout that started before the loop!
				genetics_out << base0 << base1 << ",";
			}
		} 

	}

	//now we'll go through all the juvenile females and record their positions and allele frequencies
    const unsigned int totalYF = GetLiveArraySize(vob_JuvenileFemale);
	VoleSum+=totalYF;
	for (unsigned j = 0; j < totalYF; j++){	
		auto JF = dynamic_cast<Vole_JuvenileFemale*>(SupplyAnimalPtr(vob_JuvenileFemale,j));
		int nowx= JF -> Supply_m_Location_x(); //getting her position
		int nowy= JF -> Supply_m_Location_y();

		int gx = static_cast<int>(floor(nowx / sqwidth)); //getting the small quadrant she is in
		int gy = static_cast<int>(floor(nowy / sqwidth));

		int GridCellIndex; //let's put her in a big quadrant now, and exclude any voles in the buffer zones.
		GridCellIndex=Four_AssignQuadrant(gx,gy,DivisionPoint);
		if (GridCellIndex==100){
			ExcludedVolesCount++;
			continue;
		}
        int nowquadrant=GridCellIndex;

        int birthx= JF->xborn; //now let's find out where she was born
		int birthy= JF->yborn;

		int birthquadrant;

        int birthx_q = static_cast<int>(floor(birthx / sqwidth)); //put her birthplace in the small quadrant
		int birthy_q = static_cast<int>(floor(birthy / sqwidth));
		birthquadrant=Four_AssignQuadrant(birthx_q,birthy_q,DivisionPoint);
		if (birthquadrant==100){
			ExcludedVolesCount++;
			continue;
		}

        //writing to position file where she is and where she was born and what type she is..
        position_out << year << "," << a_day << "," << "JuvenileFemale" <<"," << birthx << "," << birthy << "," << nowx << "," << nowy << "," << birthquadrant << "," << nowquadrant << "," << DivisionDirection << "\n";

		//now let's look at her genes
		vector<vector<int>> Genome0= JF->new_Genes.Genome0; //getting both of her genomes (because she's a diploid girl :) )
		vector<vector<int>> Genome1= JF->new_Genes.Genome1;

        //the first part of the gene printout (the rest will be "filled out in the following loop")
		genetics_out <<"\n" << year <<"," << a_day << "," << "JuvenileFemale" << "," << nowx << ","<< nowy << "," << nowquadrant << "," << DivisionDirection <<","; 

		for (int a = 0; a < ChromosomeCounting; a++){
			for (int b = 0; b < LocusCounting; b++){
				//"collecting" bases for writing genome to genetics_out file.
				int base0=Genome0[a][b];
				int base1=Genome1[a][b];
				//counting heterozygosity in the subpop in the quadrant
				if (base0!=base1){
					HetCount[GridCellIndex][a][b]+=1;
				}
                //the following is for counting allele frequencies in the quadrant
				AlleleCountInQuadrants[GridCellIndex][a][b][base0]+=1;
				AlleleCountInQuadrants[GridCellIndex][a][b][base1]+=1;
				//the last part of the genetics_out printout that started before the loop!
				genetics_out << base0 << base1 << ",";
			}
		} 

	}

	//now we go through the adult females and record their positions and genetics
    const unsigned int totalAF = GetLiveArraySize(vob_Female);
	VoleSum+=totalAF;
	for (unsigned j = 0; j < totalAF; j++){	
		auto VF = dynamic_cast<Vole_Female*>(SupplyAnimalPtr(vob_Female,j));
		int nowx= VF -> Supply_m_Location_x(); //getting her position
		int nowy= VF -> Supply_m_Location_y();

		int gx = static_cast<int>(floor(nowx / sqwidth)); //getting the small quadrant she is in
		int gy = static_cast<int>(floor(nowy / sqwidth));

		int GridCellIndex; //let's put her in a big quadrant now, and exclude any voles in the buffer zones.
		GridCellIndex=Four_AssignQuadrant(gx,gy,DivisionPoint);
		if (GridCellIndex==100){
			ExcludedVolesCount++;
			continue;
		}
        int nowquadrant=GridCellIndex;

        int birthx= VF->xborn; //now let's find out where she was born
		int birthy= VF->yborn;

		int birthquadrant;

        int birthx_q = static_cast<int>(floor(birthx / sqwidth)); //put her birthplace in the small quadrant
		int birthy_q = static_cast<int>(floor(birthy / sqwidth));
		birthquadrant=Four_AssignQuadrant(birthx_q,birthy_q,DivisionPoint);
		if (birthquadrant==100){
			ExcludedVolesCount++;
			continue;
		}

        //writing to position file where she is and where she was born and what type she is..
        position_out << year << "," << a_day << "," << "AdultFemale" <<"," << birthx << "," << birthy << "," << nowx << "," << nowy << "," << birthquadrant << "," << nowquadrant << "," << DivisionDirection << "\n";

		//now let's look at her genes
		vector<vector<int>> Genome0= VF->new_Genes.Genome0; //get both genomes because she's a diploid lady :) 
		vector<vector<int>> Genome1= VF->new_Genes.Genome1;

        //the first part of the gene printout (the rest will be "filled out in the following loop")
		genetics_out <<"\n" << year <<"," << a_day << "," << "AdultFemale" << "," << nowx << ","<< nowy << "," << nowquadrant << "," << DivisionDirection << ","; 

		for (int a = 0; a < ChromosomeCounting; a++){
			for (int b = 0; b < LocusCounting; b++){
				//"collecting" bases for writing genome to genetics_out file.
				int base0=Genome0[a][b];
				int base1=Genome1[a][b];
				//counting heterozygosity in the subpop in the quadrant
				if (base0!=base1){
					HetCount[GridCellIndex][a][b]+=1;
				}
                //the following is for counting allele frequencies in the quadrant
				AlleleCountInQuadrants[GridCellIndex][a][b][base0]+=1;
				AlleleCountInQuadrants[GridCellIndex][a][b][base1]+=1;
				//the last part of the genetics_out printout that started before the loop!
				genetics_out << base0 << base1 << ",";
			}
		} 

	}
	
	//Now for looking through the final group, adult males:
    const unsigned int totalAM = GetLiveArraySize(vob_Male);
	VoleSum+=totalAM;
	for (unsigned j = 0; j < totalAM; j++){	
		auto VM = dynamic_cast<Vole_Male*>(SupplyAnimalPtr(vob_Male,j));
		int nowx= VM -> Supply_m_Location_x(); //getting his position
		int nowy= VM -> Supply_m_Location_y();

		int gx = static_cast<int>(floor(nowx / sqwidth)); //getting the small quadrant he is in
		int gy = static_cast<int>(floor(nowy / sqwidth));

		int GridCellIndex; //let's put him in a big quadrant now, and exclude any voles in the buffer zones.
		GridCellIndex=Four_AssignQuadrant(gx,gy,DivisionPoint);
		if (GridCellIndex==100){
			ExcludedVolesCount++;
			continue;
		}
        int nowquadrant=GridCellIndex;

        int birthx= VM->xborn; //now let's find out where he was born
		int birthy= VM->yborn;

		int birthquadrant;

        int birthx_q = static_cast<int>(floor(birthx / sqwidth)); //put his birthplace in the small quadrant
		int birthy_q = static_cast<int>(floor(birthy / sqwidth));
		GridCellIndex=Four_AssignQuadrant(birthx_q,birthy_q,DivisionPoint);
		if (GridCellIndex==100){
			ExcludedVolesCount++;
			continue;
		}

        //writing to position file where he is and where he was born and what type he is..
        position_out << year << "," << a_day << "," << "AdultMale" <<"," << birthx << "," << birthy << "," << nowx << "," << nowy << "," << birthquadrant << "," << nowquadrant <<"," << DivisionDirection << "\n";

		//now let's look at his genes
		vector<vector<int>> Genome0= VM->new_Genes.Genome0; //get both genomes because he's a diploid guy ;) 
		vector<vector<int>> Genome1= VM->new_Genes.Genome1;

        //the first part of the gene printout (the rest will be "filled out in the following loop")
		genetics_out <<"\n" << year <<"," << a_day << "," << "AdultMale" << "," << nowx << ","<< nowy << "," << nowquadrant << "," << DivisionDirection << ","; 

		for (int a = 0; a < ChromosomeCounting; a++){
			for (int b = 0; b < LocusCounting; b++){
				//"collecting" bases for writing genome to genetics_out file.
				int base0=Genome0[a][b];
				int base1=Genome1[a][b];
				//counting heterozygosity in the subpop in the quadrant
				if (base0!=base1){
					HetCount[GridCellIndex][a][b]+=1;
				}
                //the following is for counting allele frequencies in the quadrant
				AlleleCountInQuadrants[GridCellIndex][a][b][base0]+=1;
				AlleleCountInQuadrants[GridCellIndex][a][b][base1]+=1;
				//the last part of the genetics_out printout that started before the loop!
				genetics_out << base0 << base1 <<",";
			}
		} 

	}
    // Output AlleleCountInQuadrants data
    for (int i = 0; i < 4; ++i) {
        for (int a = 0; a < ChromosomeCounting; ++a) {
            for (int b = 0; b < LocusCounting; ++b) {
                for (int c = 0; c < DiffBaseCounting; ++c) {
					genetics_count_quadrants << "\n" << year << "," << a_day << "," << i << "," << a << "," << b << "," << c << ",";
                    genetics_count_quadrants << AlleleCountInQuadrants[i][a][b][c] << "," << DivisionDirection;
                }
            }
        }
    }
    genetics_count_quadrants.close(); // now allele counts per quadrant have been output

	//Now that we have the allele frequencies pr quadrant we can "easily" calculate F_ST and F_IS between the quadrants :)
	vector<vector<double>> FST_Matrix(4, vector<double>(4,0.0)); //To hold the F_ST values 
	vector<double> FISVector(4,0); //To hold the F_IS VALUES
	
	for (int i = 0; i < 4; ++i) { //looping through the 4 big quadrants and their populations to calc. F_ST (between two pops at a time)
		vector<vector<vector<double>>> HetFreq((4), vector<vector<double>>(ChromosomeCounting, vector<double>(LocusCounting, 0.0))); //each pop needs it's own heterozyg freq. holder
		for (int j = 0; j < 4; ++j){
			if (i < j){
				double FIS_Sum0=0.0; //FIS of the ith population
				double FIS_Sum1=0.0; //FIS of the jth population
				double FIS_Average0=0.0; // to be calculates from sum
				double FIS_Average1=0.0; // to be calculates from sum

				double FST_Sum=0; //FST of the ith and jth population compared/together
				double ExpHetPop0Sum=0; // for accumulating all the expected hets for the genome for pop0 
				double ExpHetPop1Sum=0; // for accumulating all the expected hets for the genome for pop1
				double ExpHetPooledSum=0; // for accumulating all the expected hets for the genome for pop0 +pop1 poolrf
				double ExpHetAvr0And1Sum=0; // for the weighted average of ExpHetPop0Sum and ExpHetPop1Sum (after they themselved have been averaged, of course)
				double HetFreqPop0Sum=0; // for accumulating all the obs hets for the genome for pop0 
				double HetFreqPop1Sum=0; // for accumulating all the obs hets for the genome for pop1
				double FST_Average; //FST for the pair to be held by this variable

                vector<vector<vector<double>>> AlleleCountInQuadrantsP(ChromosomeCounting, vector<vector<double>>(LocusCounting, vector<double>(DiffBaseCounting, 0))); //holder vector for allele count in ith+ jth population pooled together
               
                //allele counts need to be transformed into frequencies to calculate FST and FIS
				vector<vector<vector<double>>> AlleleFrequencyInQuadrants0(ChromosomeCounting, vector<vector<double>>(LocusCounting, vector<double>(DiffBaseCounting, 0))); //holder vector for allele freqs for ith population
				vector<vector<vector<double>>> AlleleFrequencyInQuadrants1(ChromosomeCounting, vector<vector<double>>(LocusCounting, vector<double>(DiffBaseCounting, 0))); //holder vector for allele freqs for jth population
				vector<vector<vector<double>>> AlleleFrequencyInQuadrantsP(ChromosomeCounting, vector<vector<double>>(LocusCounting, vector<double>(DiffBaseCounting, 0))); //holder vector for alle freqs for ith + jth population pooled
				int TotalAllelesAtLocus0=0; //for total alleles for a locus in pop i, needs init here to be available at right level
				int TotalAllelesAtLocus1=0; //for total alleles for a locus in pop j
				
				for (int a = 0; a < ChromosomeCounting; ++a) { //to get frequncies we need to get total number (so add together counts and assign allele sum to pooled pop)
					for (int b = 0; b < LocusCounting; ++b) {
						TotalAllelesAtLocus0=0; //for total alleles in pop i
						TotalAllelesAtLocus1=0; //for total alleles in pop j 
						for (int c = 0; c < DiffBaseCounting; ++c) { // loop through the bases
							int BaseCount0= AlleleCountInQuadrants[i][a][b][c];
							TotalAllelesAtLocus0=TotalAllelesAtLocus0+BaseCount0; //count how many alleles were found all in all at this locus (2* animal count..)

							int BaseCount1= AlleleCountInQuadrants[j][a][b][c];
							TotalAllelesAtLocus1=TotalAllelesAtLocus1+BaseCount1;

							AlleleCountInQuadrantsP[a][b][c]=AlleleCountInQuadrants[j][a][b][c]+AlleleCountInQuadrants[i][a][b][c]; // alleles for the pooled pop.
						}
						int IndividualsPop0=TotalAllelesAtLocus0/2; //To have nr of individuals in each pop directly
						int IndividualsPop1=TotalAllelesAtLocus1/2;

						for (int c = 0; c < DiffBaseCounting; ++c){ //Calculate those allele frequencies! For the ith, jth and pooled populations
							if (TotalAllelesAtLocus0==0.0 ){
								AlleleFrequencyInQuadrants0[a][b][c]=0.0; // if there is no population, there is no genetic diff; FST has to be zero
							}else{
								AlleleFrequencyInQuadrants0[a][b][c]=AlleleCountInQuadrants[i][a][b][c]/TotalAllelesAtLocus0;
							}
							if(TotalAllelesAtLocus1==0.0){
								AlleleFrequencyInQuadrants1[a][b][c]=0.0;
							}else{
								AlleleFrequencyInQuadrants1[a][b][c]=AlleleCountInQuadrants[j][a][b][c]/TotalAllelesAtLocus1;
							}
							if (TotalAllelesAtLocus0==0.0 && TotalAllelesAtLocus1==0.0){
								AlleleFrequencyInQuadrantsP[a][b][c]=0.0;
							}else{
								AlleleFrequencyInQuadrantsP[a][b][c]=AlleleCountInQuadrantsP[a][b][c]/(TotalAllelesAtLocus0+TotalAllelesAtLocus1);
							}
						}

                        //Calculating the heterozygosity frequencies for the ith and jth population (for FIS calculation)
                        for (int a = 0; a < ChromosomeCounting; ++a) {
							for (int b = 0; b < LocusCounting; ++b) {
								HetFreq[i][a][b]=HetCount[i][a][b]/IndividualsPop0; 
								HetFreq[j][a][b]=HetCount[j][a][b]/IndividualsPop1;
							}
						}

						//now that we have the Allele Frequencies in the quadrants we can calculate the expected heterozygosity at the locus (needed for FST and FIS)
						double F_ST_ThisLocus=0.0;
                        //this allows for 4 different bases, but these days I only use 2!
						double ExpHetThisLocusPop0 = 1 - ((AlleleFrequencyInQuadrants0[a][b][0] * AlleleFrequencyInQuadrants0[a][b][0]) + 
														(AlleleFrequencyInQuadrants0[a][b][1] * AlleleFrequencyInQuadrants0[a][b][1]) + 
														(AlleleFrequencyInQuadrants0[a][b][2] * AlleleFrequencyInQuadrants0[a][b][2]) + 
														(AlleleFrequencyInQuadrants0[a][b][3] * AlleleFrequencyInQuadrants0[a][b][3]));

						double ExpHetThisLocusPop1 = 1 - ((AlleleFrequencyInQuadrants1[a][b][0] * AlleleFrequencyInQuadrants1[a][b][0]) + 
														(AlleleFrequencyInQuadrants1[a][b][1] * AlleleFrequencyInQuadrants1[a][b][1]) + 
														(AlleleFrequencyInQuadrants1[a][b][2] * AlleleFrequencyInQuadrants1[a][b][2]) + 
														(AlleleFrequencyInQuadrants1[a][b][3] * AlleleFrequencyInQuadrants1[a][b][3]));

						double ExpHetThisLocusPopP = 1 - ((AlleleFrequencyInQuadrantsP[a][b][0] * AlleleFrequencyInQuadrantsP[a][b][0]) + 
														(AlleleFrequencyInQuadrantsP[a][b][1] * AlleleFrequencyInQuadrantsP[a][b][1]) + 
														(AlleleFrequencyInQuadrantsP[a][b][2] * AlleleFrequencyInQuadrantsP[a][b][2]) + 
														(AlleleFrequencyInQuadrantsP[a][b][3] * AlleleFrequencyInQuadrantsP[a][b][3]));

					
						//add to FST and FIS sums for later averaging
						ExpHetPop0Sum+=ExpHetThisLocusPop0; //for FST calc and FIS calc
						ExpHetPop1Sum+=ExpHetThisLocusPop1;
						
						HetFreqPop0Sum+=HetFreq[i][a][b]; // for FIS calc
						HetFreqPop1Sum+=HetFreq[j][a][b];
						
						ExpHetPooledSum+=ExpHetThisLocusPopP;

						ExpHetAvr0And1Sum+=(ExpHetThisLocusPop0*IndividualsPop0+ExpHetThisLocusPop1*IndividualsPop1)/(IndividualsPop0+IndividualsPop1); //weighted by pop sizes
					}
				}
				//Now for calculating a genome wide FST average for ith and jth pop together, and FIS average for ith and jth pop seperately
				if (ExpHetPooledSum!=0.0){ 
					FST_Average = 1 - (((ExpHetAvr0And1Sum / (static_cast<double>(ChromosomeCounting * LocusCounting))) / (ExpHetPooledSum / (static_cast<double>(ChromosomeCounting * LocusCounting)))));
				}else{
					FST_Average=0.0; // if no hets are expected in the pooled pop, there is no diff; FST must be zero.
				}
				// Now for calculating FIS
				HetAndExpHetPopwise[i][0]=(ExpHetPop0Sum / (static_cast<double>(ChromosomeCounting * LocusCounting))); //expected het average pop i
				HetAndExpHetPopwise[j][0]=(ExpHetPop1Sum / (static_cast<double>(ChromosomeCounting * LocusCounting))); //expected het average pop j
				HetAndExpHetPopwise[i][1]=(HetFreqPop0Sum / (static_cast<double>(ChromosomeCounting * LocusCounting))); //expected het average pop i
				HetAndExpHetPopwise[j][1]=(HetFreqPop1Sum / (static_cast<double>(ChromosomeCounting * LocusCounting))); //expected het average pop j
				//calculating FIS for both populations (ith and jth)
				if (ExpHetPop0Sum!=0.0){
					FIS_Average0 = 1 - ((HetFreqPop0Sum / (static_cast<double>(ChromosomeCounting * LocusCounting))) / (ExpHetPop0Sum / (static_cast<double>(ChromosomeCounting * LocusCounting))));
				}else{
					FIS_Average0=0.0; // if no hets are expected in the pop, there is no diff; FIS must be zero.
				}
				if (ExpHetPop1Sum!=0.0){
					FIS_Average1 = 1 - ((HetFreqPop1Sum / (static_cast<double>(ChromosomeCounting * LocusCounting))) / (ExpHetPop1Sum / (static_cast<double>(ChromosomeCounting * LocusCounting))));
				}else{
					FIS_Average1= 0.0;// if no hets are expected in the pop, there is no diff; FIS must be zero.
				}
				
				if (TotalAllelesAtLocus0!=0 && TotalAllelesAtLocus1!=0){
					FST_Matrix[i][j]=FST_Average; // put FST average in FST matrix
					FST_Matrix[j][i]=FST_Average; // put same FST average at other side of diagonal :)
				}else{
					FST_Matrix[i][j]=-2.0; // if there is no individual in one quadrant, set to invalid number to signal that
					FST_Matrix[j][i]=-2.0;
				}

				if (TotalAllelesAtLocus0!=0){
					FISVector[i]=FIS_Average0; //FIS average of ith pop into ith postition in FIS vector
				}else{
					FISVector[i]=-2.0; //if no pop, there's no FIS
				}

				if (TotalAllelesAtLocus1!=0){
					FISVector[j]=FIS_Average1; //FIS average of jth pop into jth postition in FIS vector
				}else{
	 				FISVector[j]=-2.0;//if no pop, there's no FIS
				}
			}
		}
    }
    //write to FST and FIS output files from FST matrix and FIS vector.
	FSTout << year << "," << DivisionDirection << "\n";
	for (size_t i = 0; i < FST_Matrix.size(); ++i) {
		for (size_t j = 0; j < FST_Matrix[i].size(); ++j) {
			FSTout << std::setprecision(5) << FST_Matrix[i][j] << " "; // write to the matrix style FST output
			FSToutForAnalysis << year <<","<< i << "," << j << "," << FST_Matrix[i][j] <<"," << DivisionDirection <<endl; //write to the dateset style FST output
		}
		FSTout << endl;
	}
	for (size_t i = 0; i < FISVector.size(); ++i) {
            FISoutForAnalysis << year << "," << i << "," << FISVector[i] <<"," << DivisionDirection << endl; // Write each FIS element to the file, followed by a newline
        }
	//also write heterozygosity (obs and exp) to a file
    for (int i = 0; i < HetAndExpHetPopwise.size(); ++i) {
        HetOut <<"\n" << year << "," << i << ",";
        for (int k = 0; k < HetAndExpHetPopwise[i].size(); ++k) {
            HetOut << HetAndExpHetPopwise[i][k];
            // Check if this is the last element in the row
            if (k != HetAndExpHetPopwise[i].size() - 1) {
                // Add comma if it's not the last element to make a nice csv format!
                HetOut << ",";
            }
        }
	}
};

// This version is for getting output of genetics of voles in a 3*3 quadrant grid!! 
// More in depth comments for the same concept can be found in Four_QuadrantsBasesGeneticOutput
// But the function outputs allele frequencies, heterozyg.,FIS,FST, and positions of voles categorized by the 4 quadrants
void Vole_Population_Manager::Nine_QuadrantBasedGeneticOutput(const int a_day, const int year, const string landscape_info) {
    TAnimal* VB;

	string DivisionDirection;
	DivisionDirection=landscape_info.substr(33, 2); // adjust this based on the naming of your landscapes


	int GridCellsNr=9;
	int DivisionPoint=GridCellsNr/2; // so 3 smaller quadrants will be in each big quadrant
	int ExcludedVolesCount=0;
	int VoleSum=0;

    const int width = m_TheLandscape->SupplySimAreaWidth(); // get width of area
	const double sqwidth = width/GridCellsNr;

    //making a new file for GeneticsAndPosition for every year. so needs new titles/headers each year :) 
	string genetics_out_file = "GeneticsAndPosition" + to_string(year) + ".txt";
	ofstream genetics_out(genetics_out_file);
	genetics_out << "year,day,type,nowx,nowy,nowquadrant,division_direction,"; //making the title for the MANY columns
	for (int z = 0; z < ChromosomeCounting * LocusCounting; ++z) {
        // Append "locusZ" to the output stream, the title for the column
        genetics_out << "locus" << z << ",";
    }

    //these files are for appending every year, which is why we only write the "header" in the first year
    ofstream position_out("BirthPositionAndPositionNow.txt",ios::app);
    ofstream genetics_count_quadrants("GeneticsInQuadrants.txt", ios::app);
    ofstream FSTout("FST_output_matrix_style.txt", ios::app);
	ofstream FSToutForAnalysis("FST_output_dataset_style.txt", ios::app);
	ofstream FISoutForAnalysis("FIS_output_dataset_style.txt", ios::app);
	ofstream HetOut("Heterozyg_output.txt",ios::app);

	if (year==1){ //only write headers/titles in the first year for most files
		FSToutForAnalysis << "year,poprow,popcol,FST,division_direction \n";
		FISoutForAnalysis << "year,pop,FIS,division_direction \n";
        position_out << "year,day,type,birthx,birthy,nowx,nowy,birthquadrant,nowquadrant,division_direction\n";
        genetics_count_quadrants <<"year,day,quadrant,chromosome,locus,base,division_direction\n";
		HetOut << "\nyear,pop,avr_het_exp,avr_het";
	}

    //Structures to hold Allele counts and heterozygosity counts
	vector<vector<vector<vector<double>>>> AlleleCountInQuadrants((9), vector<vector<vector<double>>>(ChromosomeCounting, vector<vector<double>>(LocusCounting, vector<double>(DiffBaseCounting,0))));
	vector<vector<vector<double>>> HetCount((9), vector<vector<double>>(ChromosomeCounting, vector<double>(LocusCounting,0.0)));
	vector<vector<double>> HetAndExpHetPopwise((9), vector<double>(2,-2.0));
	
    //Going through the juvenile males to get their postions, alleles and heterozygosity
	const unsigned int totalYM = GetLiveArraySize(vob_JuvenileMale);
	VoleSum+=totalYM;
	for (unsigned j = 0; j < totalYM; j++){	
		auto JM = dynamic_cast<Vole_JuvenileMale*>(SupplyAnimalPtr(vob_JuvenileMale,j));
		int nowx= JM -> Supply_m_Location_x(); //getting his position
		int nowy= JM -> Supply_m_Location_y();

		int gx = static_cast<int>(floor(nowx / sqwidth)); //getting the small quadrant he is in
		int gy = static_cast<int>(floor(nowy / sqwidth));

		int GridCellIndex; //let's put him in a big quadrant now, and exclude any voles in the buffer zones.
		GridCellIndex=Nine_AssignQuadrant(gx,gy,DivisionPoint);
        int nowquadrant=GridCellIndex;

        int birthx= JM->xborn; //now let's find out where he was born
		int birthy= JM->yborn;

		int birthquadrant;

        int birthx_q = static_cast<int>(floor(birthx / sqwidth)); //put his birthplace in the small quadrant
		int birthy_q = static_cast<int>(floor(birthy / sqwidth));
		birthquadrant=Nine_AssignQuadrant(birthx_q,birthy_q,DivisionPoint);
        //writing to position file where he is and where he was born and what type he is..
        position_out << year << "," << a_day << "," << "JuvenileMale" <<"," << birthx << "," << birthy << "," << nowx << "," << nowy << "," << birthquadrant << "," << nowquadrant << "," << DivisionDirection << "\n";

		vector<vector<int>> Genome0= JM->new_Genes.Genome0; //now let's look at his genes
		vector<vector<int>> Genome1= JM->new_Genes.Genome1;

        //the first part of the gene printout (the rest will be "filled out in the following loop")
		genetics_out <<"\n" << year <<"," << a_day << "," << "JuvenileMale" << "," << nowx << ","<< nowy << "," << nowquadrant << "," <<DivisionDirection << ","; 

		for (int a = 0; a < ChromosomeCounting; a++){
			for (int b = 0; b < LocusCounting; b++){
				//"collecting" bases for writing genome to genetics_out file.
				int base0=Genome0[a][b];
				int base1=Genome1[a][b];
				//counting heterozygosity in the subpop in the quadrant
				if (base0!=base1){
					HetCount[GridCellIndex][a][b]+=1;
				}
                //the following is for counting allele frequencies in the quadrant
				AlleleCountInQuadrants[GridCellIndex][a][b][base0]+=1;
				AlleleCountInQuadrants[GridCellIndex][a][b][base1]+=1;
				//the last part of the genetics_out printout that started before the loop!
				genetics_out << base0 << base1 << ",";
			}
		} 

	}

	//Going through the juvenile females to get their postions, alleles and heterozygosity
    const unsigned int totalYF = GetLiveArraySize(vob_JuvenileFemale);
	VoleSum+=totalYF;
	for (unsigned j = 0; j < totalYF; j++){	
		auto JF = dynamic_cast<Vole_JuvenileFemale*>(SupplyAnimalPtr(vob_JuvenileFemale,j));
		int nowx= JF -> Supply_m_Location_x(); //getting her position
		int nowy= JF -> Supply_m_Location_y();

		int gx = static_cast<int>(floor(nowx / sqwidth)); //getting the small quadrant she is in
		int gy = static_cast<int>(floor(nowy / sqwidth));

		int GridCellIndex; //let's put her in a big quadrant now, and exclude any voles in the buffer zones.
		GridCellIndex=Nine_AssignQuadrant(gx,gy,DivisionPoint);
        int nowquadrant=GridCellIndex;

        int birthx= JF->xborn; //now let's find out where she was born
		int birthy= JF->yborn;

		int birthquadrant;

        int birthx_q = static_cast<int>(floor(birthx / sqwidth)); //put her birthplace in the small quadrant
		int birthy_q = static_cast<int>(floor(birthy / sqwidth));
		birthquadrant=Nine_AssignQuadrant(birthx_q,birthy_q,DivisionPoint);

        //writing to position file where she is and where she was born and what type she is..
        position_out << year << "," << a_day << "," << "JuvenileFemale" <<"," << birthx << "," << birthy << "," << nowx << "," << nowy << "," << birthquadrant << "," << nowquadrant << "," << DivisionDirection << "\n";

		vector<vector<int>> Genome0= JF->new_Genes.Genome0; //now let's look at her genes
		vector<vector<int>> Genome1= JF->new_Genes.Genome1;

        //the first part of the gene printout (the rest will be "filled out in the following loop")
		genetics_out <<"\n" << year <<"," << a_day << "," << "JuvenileFemale" << "," << nowx << ","<< nowy << "," << nowquadrant << "," << DivisionDirection <<","; 

		for (int a = 0; a < ChromosomeCounting; a++){
			for (int b = 0; b < LocusCounting; b++){
				//"collecting" bases for writing genome to genetics_out file.
				int base0=Genome0[a][b];
				int base1=Genome1[a][b];
				//counting heterozygosity in the subpop in the quadrant
				if (base0!=base1){
					HetCount[GridCellIndex][a][b]+=1;
				}
                //the following is for counting allele frequencies in the quadrant
				AlleleCountInQuadrants[GridCellIndex][a][b][base0]+=1;
				AlleleCountInQuadrants[GridCellIndex][a][b][base1]+=1;
				//the last part of the genetics_out printout that started before the loop!
				genetics_out << base0 << base1 << ",";
			}
		} 

	}
	//Going through the adult females to get their postions, alleles and heterozygosity
    const unsigned int totalAF = GetLiveArraySize(vob_Female);
	VoleSum+=totalAF;
	for (unsigned j = 0; j < totalAF; j++){	
		auto VF = dynamic_cast<Vole_Female*>(SupplyAnimalPtr(vob_Female,j));
		int nowx= VF -> Supply_m_Location_x(); //getting her position
		int nowy= VF -> Supply_m_Location_y();

		int gx = static_cast<int>(floor(nowx / sqwidth)); //getting the small quadrant she is in
		int gy = static_cast<int>(floor(nowy / sqwidth));

		int GridCellIndex; //let's put her in a big quadrant now, and exclude any voles in the buffer zones.
		GridCellIndex=Nine_AssignQuadrant(gx,gy,DivisionPoint);
        int nowquadrant=GridCellIndex;

        int birthx= VF->xborn; //now let's find out where she was born
		int birthy= VF->yborn;

		int birthquadrant;

        int birthx_q = static_cast<int>(floor(birthx / sqwidth)); //put her birthplace in the small quadrant
		int birthy_q = static_cast<int>(floor(birthy / sqwidth));
		birthquadrant=Nine_AssignQuadrant(birthx_q,birthy_q,DivisionPoint);

        //writing to position file where she is and where she was born and what type she is..
        position_out << year << "," << a_day << "," << "AdultFemale" <<"," << birthx << "," << birthy << "," << nowx << "," << nowy << "," << birthquadrant << "," << nowquadrant << "," << DivisionDirection << "\n";

		vector<vector<int>> Genome0= VF->new_Genes.Genome0; //now let's look at her genes
		vector<vector<int>> Genome1= VF->new_Genes.Genome1;

        //the first part of the gene printout (the rest will be "filled out in the following loop")
		genetics_out <<"\n" << year <<"," << a_day << "," << "AdultFemale" << "," << nowx << ","<< nowy << "," << nowquadrant << "," << DivisionDirection << ","; 

		for (int a = 0; a < ChromosomeCounting; a++){
			for (int b = 0; b < LocusCounting; b++){
				//"collecting" bases for writing genome to genetics_out file.
				int base0=Genome0[a][b];
				int base1=Genome1[a][b];
				//counting heterozygosity in the subpop in the quadrant
				if (base0!=base1){
					HetCount[GridCellIndex][a][b]+=1;
				}
                //the following is for counting allele frequencies in the quadrant
				AlleleCountInQuadrants[GridCellIndex][a][b][base0]+=1;
				AlleleCountInQuadrants[GridCellIndex][a][b][base1]+=1;
				//the last part of the genetics_out printout that started before the loop!
				genetics_out << base0 << base1 << ",";
			}
		} 

	}
	//Going through the adult males to get their postions, alleles and heterozygosity
    const unsigned int totalAM = GetLiveArraySize(vob_Male);
	VoleSum+=totalAM;
	for (unsigned j = 0; j < totalAM; j++){	
		auto VM = dynamic_cast<Vole_Male*>(SupplyAnimalPtr(vob_Male,j));
		int nowx= VM -> Supply_m_Location_x(); //getting his position
		int nowy= VM -> Supply_m_Location_y();

		int gx = static_cast<int>(floor(nowx / sqwidth)); //getting the small quadrant he is in
		int gy = static_cast<int>(floor(nowy / sqwidth));

		int GridCellIndex; //let's put him in a big quadrant now, and exclude any voles in the buffer zones.
		GridCellIndex=Nine_AssignQuadrant(gx,gy,DivisionPoint);
        int nowquadrant=GridCellIndex;

        int birthx= VM->xborn; //now let's find out where he was born
		int birthy= VM->yborn;

		int birthquadrant;

        int birthx_q = static_cast<int>(floor(birthx / sqwidth)); //put his birthplace in the small quadrant
		int birthy_q = static_cast<int>(floor(birthy / sqwidth));
		birthquadrant=Nine_AssignQuadrant(birthx_q,birthy_q,DivisionPoint);

        //writing to position file where he is and where he was born and what type he is..
        position_out << year << "," << a_day << "," << "AdultMale" <<"," << birthx << "," << birthy << "," << nowx << "," << nowy << "," << birthquadrant << "," << nowquadrant <<"," << DivisionDirection << "\n";

		vector<vector<int>> Genome0= VM->new_Genes.Genome0; //now let's look at his genes
		vector<vector<int>> Genome1= VM->new_Genes.Genome1;

        //the first part of the gene printout (the rest will be "filled out in the following loop")
		genetics_out <<"\n" << year <<"," << a_day << "," << "AdultMale" << "," << nowx << ","<< nowy << "," << nowquadrant << "," << DivisionDirection << ","; 

		for (int a = 0; a < ChromosomeCounting; a++){
			for (int b = 0; b < LocusCounting; b++){
				//"collecting" bases for writing genome to genetics_out file.
				int base0=Genome0[a][b];
				int base1=Genome1[a][b];
				//counting heterozygosity in the subpop in the quadrant
				if (base0!=base1){
					HetCount[GridCellIndex][a][b]+=1;
				}
                //the following is for counting allele frequencies in the quadrant
				AlleleCountInQuadrants[GridCellIndex][a][b][base0]+=1;
				AlleleCountInQuadrants[GridCellIndex][a][b][base1]+=1;
				//the last part of the genetics_out printout that started before the loop!
				genetics_out << base0 << base1 <<",";
			}
		} 

	}
    // Output AlleleCountInQuadrants data
    for (int i = 0; i < 9; ++i) {
        for (int a = 0; a < ChromosomeCounting; ++a) {
            for (int b = 0; b < LocusCounting; ++b) {
                for (int c = 0; c < DiffBaseCounting; ++c) {
					genetics_count_quadrants << "\n" << year << "," << a_day << "," << i << "," << a << "," << b << "," << c << ",";
                    genetics_count_quadrants << AlleleCountInQuadrants[i][a][b][c] << "," << DivisionDirection;
                }
            }
        }
    }
    genetics_count_quadrants.close();

	//Now that we have the allele frequencies pr quadrant we can "easily" calculate F_ST and F_IS between the quadrants :)
	vector<vector<double>> FST_Matrix(9, vector<double>(9,0.0)); //To hold the F_ST values //INCONSISTENT:0.0 or -1??
	vector<double> FISVector(9,-1); //To hold the F_IS VALUES
	
	for (int i = 0; i < 9; ++i) { //looping through the 9 big quadrants and their populations to calc. F_ST (between two pops at a time)
		vector<vector<vector<double>>> HetFreq((9), vector<vector<double>>(ChromosomeCounting, vector<double>(LocusCounting, 0.0))); //each pop needs it's own heterozyg freq. holder
		for (int j = 0; j < 9; ++j){
			if (i < j){
				double FIS_Sum0=0.0; //FIS of the ith population
				double FIS_Sum1=0.0; //FIS of the jth population
				double FIS_Average0=0.0; 
				double FIS_Average1=0.0;

				double FST_Sum=0; //FST of the ith and jth population compared/together
				double ExpHetPop0Sum=0;
				double ExpHetPop1Sum=0;
				double ExpHetPooledSum=0;
				double ExpHetAvr0And1Sum=0;
				double HetFreqPop0Sum=0;
				double HetFreqPop1Sum=0;
				double FST_Average;

                vector<vector<vector<double>>> AlleleCountInQuadrantsP(ChromosomeCounting, vector<vector<double>>(LocusCounting, vector<double>(DiffBaseCounting, 0))); //holder vector for allele count in ith+ jth population pooled together
               
                //allele counts need to be transformed into frequencies to calculate FST and FIS
				vector<vector<vector<double>>> AlleleFrequencyInQuadrants0(ChromosomeCounting, vector<vector<double>>(LocusCounting, vector<double>(DiffBaseCounting, 0))); //holder vector for allele freqs for ith population
				vector<vector<vector<double>>> AlleleFrequencyInQuadrants1(ChromosomeCounting, vector<vector<double>>(LocusCounting, vector<double>(DiffBaseCounting, 0))); //holder vector for allele freqs for jth population
				vector<vector<vector<double>>> AlleleFrequencyInQuadrantsP(ChromosomeCounting, vector<vector<double>>(LocusCounting, vector<double>(DiffBaseCounting, 0))); //holder vector for alle freqs for ith + jth population pooled
				int TotalAllelesAtLocus0=0; //for total alleles for a locus in pop i, needs init here to be available at right level
				int TotalAllelesAtLocus1=0; //for total alleles for a locus in pop j
				
				for (int a = 0; a < ChromosomeCounting; ++a) { //to get frequncies we need to get total number (so add together counts and assign allele sum to pooled pop)
					for (int b = 0; b < LocusCounting; ++b) {
						TotalAllelesAtLocus0=0; //for total alleles in pop i
						TotalAllelesAtLocus1=0; //for total alleles in pop j 
						for (int c = 0; c < DiffBaseCounting; ++c) {
							//cout << "c is " << c; 
							int BaseCount0= AlleleCountInQuadrants[i][a][b][c];
							TotalAllelesAtLocus0=TotalAllelesAtLocus0+BaseCount0;

							int BaseCount1= AlleleCountInQuadrants[j][a][b][c];
							TotalAllelesAtLocus1=TotalAllelesAtLocus1+BaseCount1;

							AlleleCountInQuadrantsP[a][b][c]=AlleleCountInQuadrants[j][a][b][c]+AlleleCountInQuadrants[i][a][b][c];
						}
						int IndividualsPop0=TotalAllelesAtLocus0/2; //To have nr of individuals in each pop directly
						int IndividualsPop1=TotalAllelesAtLocus1/2;

						for (int c = 0; c < DiffBaseCounting; ++c){ //Calculate those allele frequencies! For the ith, jth and pooled populations
							if (TotalAllelesAtLocus0==0.0 ){
								AlleleFrequencyInQuadrants0[a][b][c]=0.0;
							}else{
								AlleleFrequencyInQuadrants0[a][b][c]=AlleleCountInQuadrants[i][a][b][c]/TotalAllelesAtLocus0;
							}
							if(TotalAllelesAtLocus1==0.0){
								AlleleFrequencyInQuadrants1[a][b][c]=0.0;
							}else{
								AlleleFrequencyInQuadrants1[a][b][c]=AlleleCountInQuadrants[j][a][b][c]/TotalAllelesAtLocus1;
							}
							if (TotalAllelesAtLocus0==0.0 && TotalAllelesAtLocus1==0.0){
								AlleleFrequencyInQuadrantsP[a][b][c]=0.0;
							}else{
								AlleleFrequencyInQuadrantsP[a][b][c]=AlleleCountInQuadrantsP[a][b][c]/(TotalAllelesAtLocus0+TotalAllelesAtLocus1);
							}
						}

                        //Calculating the heterozygosity frequencies for the ith and jth population (for FIS calculation)
                        for (int a = 0; a < ChromosomeCounting; ++a) {
							for (int b = 0; b < LocusCounting; ++b) {
								HetFreq[i][a][b]=HetCount[i][a][b]/IndividualsPop0; 
								HetFreq[j][a][b]=HetCount[j][a][b]/IndividualsPop1;
							}
						}

						//now that we have the Allele Frequencies in the quadrants we can calculate the expected heterozygosity at the locus (needed for FST and FIS)
						double F_ST_ThisLocus=0.0;
                        //this allows for 4 different bases, but these days I only use 2.
						double ExpHetThisLocusPop0 = 1 - ((AlleleFrequencyInQuadrants0[a][b][0] * AlleleFrequencyInQuadrants0[a][b][0]) + 
														(AlleleFrequencyInQuadrants0[a][b][1] * AlleleFrequencyInQuadrants0[a][b][1]) + 
														(AlleleFrequencyInQuadrants0[a][b][2] * AlleleFrequencyInQuadrants0[a][b][2]) + 
														(AlleleFrequencyInQuadrants0[a][b][3] * AlleleFrequencyInQuadrants0[a][b][3]));

						double ExpHetThisLocusPop1 = 1 - ((AlleleFrequencyInQuadrants1[a][b][0] * AlleleFrequencyInQuadrants1[a][b][0]) + 
														(AlleleFrequencyInQuadrants1[a][b][1] * AlleleFrequencyInQuadrants1[a][b][1]) + 
														(AlleleFrequencyInQuadrants1[a][b][2] * AlleleFrequencyInQuadrants1[a][b][2]) + 
														(AlleleFrequencyInQuadrants1[a][b][3] * AlleleFrequencyInQuadrants1[a][b][3]));

						double ExpHetThisLocusPopP = 1 - ((AlleleFrequencyInQuadrantsP[a][b][0] * AlleleFrequencyInQuadrantsP[a][b][0]) + 
														(AlleleFrequencyInQuadrantsP[a][b][1] * AlleleFrequencyInQuadrantsP[a][b][1]) + 
														(AlleleFrequencyInQuadrantsP[a][b][2] * AlleleFrequencyInQuadrantsP[a][b][2]) + 
														(AlleleFrequencyInQuadrantsP[a][b][3] * AlleleFrequencyInQuadrantsP[a][b][3]));

					
						//add to FST and FIS sums for later averaging
						ExpHetPop0Sum+=ExpHetThisLocusPop0;
						ExpHetPop1Sum+=ExpHetThisLocusPop1;
						
						HetFreqPop0Sum+=HetFreq[i][a][b];
						HetFreqPop1Sum+=HetFreq[j][a][b];
						
						ExpHetPooledSum+=ExpHetThisLocusPopP;
						ExpHetAvr0And1Sum+=(ExpHetThisLocusPop0*IndividualsPop0+ExpHetThisLocusPop1*IndividualsPop1)/(IndividualsPop0+IndividualsPop1);
					}
				}
				//Now for calculating a genome wide FST average for ith and jth pop together, and FIS average for ith and jth pop seperately
				if (ExpHetPooledSum!=0.0){
					FST_Average = 1 - (((ExpHetAvr0And1Sum / (static_cast<double>(ChromosomeCounting * LocusCounting))) / (ExpHetPooledSum / (static_cast<double>(ChromosomeCounting * LocusCounting)))));
				}else{
					FST_Average=0.0;
				}

				HetAndExpHetPopwise[i][0]=(ExpHetPop0Sum / (static_cast<double>(ChromosomeCounting * LocusCounting))); //expected het average pop i
				HetAndExpHetPopwise[j][0]=(ExpHetPop1Sum / (static_cast<double>(ChromosomeCounting * LocusCounting))); //expected het average pop j
				HetAndExpHetPopwise[i][1]=(HetFreqPop0Sum / (static_cast<double>(ChromosomeCounting * LocusCounting))); //expected het average pop i
				HetAndExpHetPopwise[j][1]=(HetFreqPop1Sum / (static_cast<double>(ChromosomeCounting * LocusCounting))); //expected het average pop j
				//calculating FIS for both populations (ith and jth)
				if (ExpHetPop0Sum!=0.0){
					FIS_Average0 = 1 - ((HetFreqPop0Sum / (static_cast<double>(ChromosomeCounting * LocusCounting))) / (ExpHetPop0Sum / (static_cast<double>(ChromosomeCounting * LocusCounting))));
				}else{
					FIS_Average0=0.0;
				}
				if (ExpHetPop1Sum!=0.0){
					FIS_Average1 = 1 - ((HetFreqPop1Sum / (static_cast<double>(ChromosomeCounting * LocusCounting))) / (ExpHetPop1Sum / (static_cast<double>(ChromosomeCounting * LocusCounting))));
				}else{
					FIS_Average1= 0.0;
				}
				
				if (TotalAllelesAtLocus0!=0 && TotalAllelesAtLocus1!=0){
					FST_Matrix[i][j]=FST_Average; // put FST average in FST matrix
					FST_Matrix[j][i]=FST_Average; // put same FST average at other side of diagonal :)
				}else{
					FST_Matrix[i][j]=-2.0; // if there is no individual in one quadrant, set to invalid number to signal that
					FST_Matrix[j][i]=-2.0;
				}

				if (TotalAllelesAtLocus0!=0){
					FISVector[i]=FIS_Average0; //FIS average of ith pop into ith postition in FIS vector
				}else{
					FISVector[i]=-2.0; //if no pop, there's no FIS
				}

				if (TotalAllelesAtLocus1!=0){
					FISVector[j]=FIS_Average1; //FIS average of jth pop into jth postition in FIS vector
				}else{
	 				FISVector[j]=-2.0;//if no pop, there's no FIS
				}
			}
		}
    }
    //write to FST and FIS output files from FST matrix and FIS vector.
	FSTout << year << "," << DivisionDirection << "\n";
	for (size_t i = 0; i < FST_Matrix.size(); ++i) {
		for (size_t j = 0; j < FST_Matrix[i].size(); ++j) {
			FSTout << std::setprecision(5) << FST_Matrix[i][j] << " ";
			FSToutForAnalysis << year <<","<< i << "," << j << "," << FST_Matrix[i][j] <<"," << DivisionDirection <<endl;
		}
		FSTout << endl;
	}
	for (size_t i = 0; i < FISVector.size(); ++i) {
            FISoutForAnalysis << year << "," << i << "," << FISVector[i] <<"," << DivisionDirection << endl; // Write each element to the file, followed by a newline
        }
    for (int i = 0; i < HetAndExpHetPopwise.size(); ++i) {
        HetOut <<"\n" << year << "," << i << ",";
        for (int k = 0; k < HetAndExpHetPopwise[i].size(); ++k) {
            HetOut << HetAndExpHetPopwise[i][k];
            // Check if this is the last element in the row
            if (k != HetAndExpHetPopwise[i].size() - 1) {
                // Add comma if it's not the last element
                HetOut << ",";
            }
        }
	}
};

//this function makes the allele frequency input file for the vole genome generation. To be called when pop manager is being initialized.
//right now only allows the bases 0 and 1. It allows the use of 2 different beta distributions for the first and second half of chromosomes.
//The parameters should be specified in the config file (or elswhere)
//the input strings should look like this for example "0.5 0.5" and "5.0 5.0".
void Vole_Population_Manager::MakeAlleleInputFile(const string params_for_beta0, const string params_for_beta1){
	ofstream AlleleOut("AlleleFreqs.txt"); //we will print to this file, from which we will read when initializing the voles
	
    for (int i = 0; i < ChromosomeCounting/2; ++i) { //make frequencies according to first beta dist parameters
        for (int j = 0; j < LocusCounting; ++j) {
			//Now allowing only 0 or 1, the vector structure is to allow random placement of the alle freq into 0 or 1, more useful if more alleles allowed
			vector<double> prob_holder(2, 0.0);
			auto seed = std::chrono::steady_clock::now().time_since_epoch().count();
			probability_distribution p2("BETA", params_for_beta0); //make a beta dist to draw allele freqs from
			double result0=p2.Get(); //Draw an allele frequency
          	double result1=(1-result0); //Two alleles allowed so the other frequency is just 1-freq.
			
			probability_distribution p3("DISCRETE", "0.5 0.5"); // Distribution to draw position (base) that get the frequency assinged. more useful if you have more than 2 alleles allowed
			int result0Position;
			int result1Position;

			while (true) { //make sure they get a different position
				result0Position = p3.Geti(); //drawing either 0 or 1 at equal chance.
				result1Position = p3.Geti();
				if (result0Position != result1Position) {
					break;
				}
			}

			prob_holder[result0Position]=result0; // assign frequencies to their new position
			prob_holder[result1Position]=result1;

			for (int k = 0; k < 2; ++k){ //writing the found frequencies to the file
				AlleleOut << prob_holder[k];
				if (k!=1){
					AlleleOut << ","; // only putting comma if not the last freq
				}
				if (k==1){
					AlleleOut << "\n";
				}
			}
		}
	}
	//make frequencies according to last beta distribution parameters, otherwise exactly same concept..
	for (int i = ChromosomeCounting/2; i < ChromosomeCounting; ++i) { 
        for (int j = 0; j < LocusCounting; ++j) {
			std::vector<double> prob_holder(2, 0.0);
			auto seed = std::chrono::steady_clock::now().time_since_epoch().count();
			probability_distribution p2("BETA", params_for_beta0);
			double result0=p2.Get(); 
          	double result1=(1-result0);

			probability_distribution p3("DISCRETE", "0.5 0.5");
			int result0Position;
			int result1Position;

			while (true) {
				result0Position = p3.Geti();
				result1Position = p3.Geti();
				if (result0Position != result1Position) {
					break;
				}
			}
			prob_holder[result0Position]=result0;
			prob_holder[result1Position]=result1;
			for (int k = 0; k < 2; ++k){
				AlleleOut << prob_holder[k];
				if (k!=1){
					AlleleOut << ",";
				}
				if (k==1){
					AlleleOut << "\n";
				}
			}
		}
	}
};

//this function compares genomes by analysing the resulting output from the functs Four_QuadratBasedGeneticOutput and Nine_QuadratBasedGeneticOutput
void Vole_Population_Manager::CompareGenomes(const int year, const string landscape_info){
	
	string DivisionDirection;
	DivisionDirection=landscape_info.substr(33, 2); //update this depending on the naming of your landscapes!

	string genetics_file = "GeneticsAndPosition" + to_string(year) + ".txt"; // A file for genes was made each year, now get the one from the current year
	ofstream DistDiffFile;

	if (year == 1) {
		DistDiffFile.open("DistanceAndGeneticDifference.txt", ios::out); //new file if it's year one, else append
	} else {
		DistDiffFile.open("DistanceAndGeneticDifference.txt", ios::app);
	}

	ifstream genetics_file_in(genetics_file); //we are reading in the most current genetics file

	if (year==1){  // put headers/titles in distance vs similarity file if it's the first year.
		DistDiffFile << "year, distance, quadrant_of_0, quadrant_of_1, similarity,division_direction,x0,y0,x1,y1\n"; 
	}

    vector<string> lines; // for holding the lines from the genetics. Probably bad memory use.. But next person reading this: you can probably fix it ;) 

    // Read all lines from the file
    string line;
    bool isFirstLine = true;
	while (getline(genetics_file_in, line)) {
		if (isFirstLine) { 
			isFirstLine = false;
			continue;  // Skip the first line as it is the titles
		}
		lines.push_back(line); // fill each line into "lines" holder.
	}
	if (lines.size() > 500) { // I only want at most 500*500 comparisons! 
        unsigned seed = static_cast<unsigned>(std::time(nullptr)); // Seed using current time
        std::shuffle(lines.begin(), lines.end(), std::default_random_engine(seed)); // shuffle order of lines
        lines.resize(500); // resize: get random lines because of the shuffeling before ;) 
    }

    // Compare each line with every other line
    for (size_t i = 0; i < lines.size(); ++i) {
        for (size_t j = i + 1; j < lines.size(); ++j) {
			double SimilarityCount=0.0;
			stringstream ss0(lines[i]); // Create stringstream for the current two lines 
			stringstream ss1(lines[j]); 
        	string token0; // holders for tokens extracted from the lines
			string token1;
			double AverageSimilarity=-1.0;
			double DistanceBetweenIndividuals=-1;
			int ElementCounter=0;

			//positions of the animals compared
			int nowx0; 
			int nowx1;
			int nowy0;
			int nowy1;
			int quadrant0;
			int quadrant1;
			


        	while (getline(ss0, token0, ',') && (getline(ss1, token1, ',' ))) {
				if (ElementCounter == 3) { //get x value position of animals
					nowx0 = stoi(token0);
					nowx1 = stoi(token1);
				}
				if (ElementCounter == 4) { //get y value position of animals
					nowy0 = stoi(token0);
					nowy1 = stoi(token1);
				}

				if (ElementCounter == 5) { //get the quadrants the animals were found in
					quadrant0 = stoi(token0);
					quadrant1 = stoi(token1);
				}
				
				
				if (ElementCounter>6){ // now it's the alleles from the loci being extracted!
					std::vector<int> digits_token0; //holder for the digits.
					std::vector<int> digits_token1;
			
					for (char digit : token0) {
						digits_token0.push_back(digit);
					}
				
					for (char digit : token1) {
						digits_token1.push_back(digit );
					}
					
					//counting how many alleles the two animals have in common at the locus
					if (digits_token0[0]==digits_token0[1] && digits_token1[0]==digits_token1[1] && digits_token0[0]!=digits_token1[0] ){
						//Do nothing if they have nothing in common
					}else if (digits_token0[0]==digits_token1[0] && digits_token0[1]==digits_token1[1] ||
					digits_token0[1]==digits_token1[0] && digits_token0[0]==digits_token1[1] ){// if they have both in common, add 2
						SimilarityCount+=2;
					}else{
						SimilarityCount++; // else (that means they have 1 in common in this case: add 1)
					}
				}
				ElementCounter++; //increase to keep iterating
			}
			AverageSimilarity=SimilarityCount/(2*ChromosomeCounting*LocusCounting); // get average similariy
			DistanceBetweenIndividuals = sqrt(pow((nowx0 - nowx1), 2) + pow((nowy0 - nowy1), 2)); // calc euclidean distance
			//now output for the two individuals compared
			DistDiffFile << year << "," << DistanceBetweenIndividuals << "," << quadrant0 << "," << quadrant1 << "," << AverageSimilarity <<"," << DivisionDirection<<","<<nowx0 << ","<<nowy0<<","<<nowx1 << ","<<nowy1<< "\n"; 

		}
	}
};
//This function gives you the population sizes pr quadrant (here four in a 2*2 grid), divided into male and females
void Vole_Population_Manager::FourQuadrantsPopulationSizeProbe(const int a_day, const int year, const string landscape_info) {
    TAnimal* VB;

	string DivisionDirection;

	DivisionDirection=landscape_info.substr(33, 2); //adjust it based on the naming of your landscapes

	int GridCellsNr=9;
	int DivisionPoint=GridCellsNr/2; // 4 will be the buffer zone
	vector<double> ExcludedVolesCount(2,0.0); //excluded males (0) and females (1)
	int VoleSum=0;

    const int width = m_TheLandscape->SupplySimAreaWidth(); // get size of each landscape
	const double sqwidth = width/GridCellsNr;// set size of each small quadrant

	ofstream pop_in_quadrants_out("PopulationInQuadrants.txt", ios::app); // first year, first day, write titles/headers
	if (a_day==1 && year==1){
		pop_in_quadrants_out << "day,quadrant0pop_males,quadrant0pop_females," 
                      << "quadrant1pop_males,quadrant1pop_females," 
                      << "quadrant2pop_males,quadrant2pop_females," 
                      << "quadrant3pop_males,quadrant3pop_females," 
                      << "excluded_males,excluded_females," 
                      << "sum_without_excluded_males,sum_without_excluded_females," 
                      << "sum_with_excluded_males,sum_with_excluded_females," 
                      << "sum_with_excluded_all,sum_without_excluded_all," 
                      << "division_direction" 
                      << "\n";
	}
	//4 quadrants, and a place for excluded, sum without excluded, sum with excluded, each having a slot for male and female.
	vector<vector<double>> PopInQuadrants(7, vector<double>(2, 0.0));
    //Going through the juvenile males to get their numbers
	const unsigned int totalYM = GetLiveArraySize(vob_JuvenileMale);
	VoleSum+=totalYM;
	for (unsigned j = 0; j < totalYM; j++){	
		auto JM = dynamic_cast<Vole_JuvenileMale*>(SupplyAnimalPtr(vob_JuvenileMale,j));
		int nowx= JM -> Supply_m_Location_x(); //getting his position
		int nowy= JM -> Supply_m_Location_y();

		int gx = static_cast<int>(floor(nowx / sqwidth)); //getting the small quadrant he is in
		int gy = static_cast<int>(floor(nowy / sqwidth));

		int GridCellIndex; //let's put him in a big quadrant now, and exclude any voles in the buffer zones.
		GridCellIndex=Four_AssignQuadrant(gx, gy, DivisionPoint); //place the vole in the big quadrant
		if (GridCellIndex!=100){ //if the voles was not exclued, count him.
			int nowquadrant=GridCellIndex;
			PopInQuadrants[nowquadrant][0]+=1; //add one male to pop of quadrant
		}else{
			ExcludedVolesCount[0]++;
		}
        
	}

	//going through juvenile females to get their positions, quadrants, numbers
    const unsigned int totalYF = GetLiveArraySize(vob_JuvenileFemale);
	VoleSum+=totalYF;
	for (unsigned j = 0; j < totalYF; j++){	
		auto JF = dynamic_cast<Vole_JuvenileFemale*>(SupplyAnimalPtr(vob_JuvenileFemale,j));
		int nowx= JF -> Supply_m_Location_x(); //getting her position
		int nowy= JF -> Supply_m_Location_y();

		int gx = static_cast<int>(floor(nowx / sqwidth)); //getting the small quadrant she is in
		int gy = static_cast<int>(floor(nowy / sqwidth));

		int GridCellIndex; //let's put her in a big quadrant now, and exclude any voles in the buffer zones.
		GridCellIndex=Four_AssignQuadrant(gx, gy, DivisionPoint); //place the vole in the big quadrant
		if (GridCellIndex!=100){ //if the voles was not exclued, count her.
			int nowquadrant=GridCellIndex;
			PopInQuadrants[nowquadrant][1]+=1; //add one female to pop of quadrant
		}else{
			ExcludedVolesCount[1]++;
		}
	}
	//going through adult females to get their positions, quadrants, numbers
    const unsigned int totalAF = GetLiveArraySize(vob_Female);
	VoleSum+=totalAF;
	for (unsigned j = 0; j < totalAF; j++){	
		auto VF = dynamic_cast<Vole_Female*>(SupplyAnimalPtr(vob_Female,j));
		int nowx= VF -> Supply_m_Location_x(); //getting her position
		int nowy= VF -> Supply_m_Location_y();

		int gx = static_cast<int>(floor(nowx / sqwidth)); //getting the small quadrant she is in
		int gy = static_cast<int>(floor(nowy / sqwidth));

		int GridCellIndex; //let's put her in a big quadrant now, and exclude any voles in the buffer zones.
		GridCellIndex=Four_AssignQuadrant(gx, gy, DivisionPoint); //place the vole in the big quadrant
		if (GridCellIndex!=100){ //if the voles was not exclued, count her.
			int nowquadrant=GridCellIndex;
			PopInQuadrants[nowquadrant][1]+=1; //add one female to pop of quadrant
		}else{
			ExcludedVolesCount[1]++;
		}
	}
	//going through adjult males to get their positions, quadrants, numbers
    const unsigned int totalAM = GetLiveArraySize(vob_Male);
	VoleSum+=totalAM;
	for (unsigned j = 0; j < totalAM; j++){	
		auto VM = dynamic_cast<Vole_Male*>(SupplyAnimalPtr(vob_Male,j));
		int nowx= VM -> Supply_m_Location_x(); //getting his position
		int nowy= VM -> Supply_m_Location_y();

		int gx = static_cast<int>(floor(nowx / sqwidth)); //getting the small quadrant he is in
		int gy = static_cast<int>(floor(nowy / sqwidth));

		int GridCellIndex; //let's put him in a big quadrant now, and exclude any voles in the buffer zones.
		GridCellIndex=Four_AssignQuadrant(gx, gy, DivisionPoint); //place the vole in the big quadrant
		if (GridCellIndex!=100){ //if the voles was not exclued, count him.
			int nowquadrant=GridCellIndex;
			PopInQuadrants[nowquadrant][0]+=1; //add one male to pop of quadrant
		}else{
			ExcludedVolesCount[0]++;
		}
	}
	PopInQuadrants[4][0]=ExcludedVolesCount[0]; //excluded males
	PopInQuadrants[4][1]=ExcludedVolesCount[1]; //excluded females
	PopInQuadrants[5][0]=PopInQuadrants[0][0]+PopInQuadrants[1][0]+PopInQuadrants[2][0]+PopInQuadrants[3][0]; //total males in quadrants
	PopInQuadrants[5][1]=PopInQuadrants[0][1]+PopInQuadrants[1][1]+PopInQuadrants[2][1]+PopInQuadrants[3][1]; //total females in quadrants
	PopInQuadrants[6][0]=PopInQuadrants[5][0]+PopInQuadrants[4][0]; //sum males with excluded voles
	PopInQuadrants[6][1]=PopInQuadrants[5][1]+PopInQuadrants[4][1]; //sum females with excluded voles
	double SumAllWithExcluded=PopInQuadrants[6][0]+ PopInQuadrants[6][1]; // sum of all quadrants + excluded
	double SumAllWithoutExcluded=PopInQuadrants[5][0]+PopInQuadrants[5][1]; //sum of all quadrants and not the excluded (from buffer zones)

	int day= ((year-1)*365)+a_day;
	pop_in_quadrants_out << day << ",";
	for (int i = 0; i < PopInQuadrants.size(); ++i) { //for each quadrant and the excluded, output males and females
        pop_in_quadrants_out<< PopInQuadrants[i][0] << "," << PopInQuadrants[i][1] <<",";
    }
	// after that, output all the summarizing elements 
	pop_in_quadrants_out << SumAllWithExcluded << "," << SumAllWithoutExcluded << "," << DivisionDirection << "\n";

};
//counts animals in 9 quadrants (3*3 grid) with no buffer zones. Divided into male + female.
void Vole_Population_Manager::NineQuadrantsPopulationSizeProbe(const int a_day, const int year, const string landscape_info) {
    TAnimal* VB;

	string DivisionDirection;

	DivisionDirection=landscape_info.substr(33, 2); //change this depending on the naming of your landscapes


	int GridCellsNr=9;
	int DivisionPoint=GridCellsNr/2;
	vector<double> ExcludedVolesCount(2,0.0); //excluded males (0) and females (1)
	int VoleSum=0;

    const int width = m_TheLandscape->SupplySimAreaWidth(); // get width of landscape
	const double sqwidth = width/GridCellsNr; //width of a small quadrant..

	ofstream pop_in_quadrants_out("PopulationInQuadrants.txt", ios::app);
	if (a_day==1 && year==1){ //first day, first year: output the MANY titles to the output file.
		pop_in_quadrants_out << "day,";
		for (int i = 0; i < 9; ++i) {
			pop_in_quadrants_out << "quadrant" << i << "pop_males,quadrant" << i << "pop_females,"; // for each quadrant
		}
		pop_in_quadrants_out << "excluded_males,excluded_females,"
							<< "sum_without_excluded_males,sum_without_excluded_females,"
							<< "sum_with_excluded_males,sum_with_excluded_females,"
							<< "sum_with_excluded_all,sum_without_excluded_all,"
							<< "division_direction"
							<< "\n";
	}
	//9 quadrants, and a place for excluded, sum without excluded, sum with excluded, each having a slot for male and female.
	vector<vector<double>> PopInQuadrants(12, vector<double>(2, 0.0));
    //Going through the juvenile males to get their postions, alleles and heterozygosity
	const unsigned int totalYM = GetLiveArraySize(vob_JuvenileMale);
	VoleSum+=totalYM;
	//going through the juvenile males and getting their position.
	for (unsigned j = 0; j < totalYM; j++){	
		auto JM = dynamic_cast<Vole_JuvenileMale*>(SupplyAnimalPtr(vob_JuvenileMale,j));
		int nowx= JM -> Supply_m_Location_x(); //getting his position
		int nowy= JM -> Supply_m_Location_y();

		int gx = static_cast<int>(floor(nowx / sqwidth)); //getting the small quadrant he is in
		int gy = static_cast<int>(floor(nowy / sqwidth));

		int GridCellIndex; //let's put him in a big quadrant now	
		GridCellIndex=Nine_AssignQuadrant(gx,gy,DivisionPoint); //assign him to a quadrant
        int nowquadrant=GridCellIndex;
		PopInQuadrants[nowquadrant][0]+=1; //add one male to pop of quadrant
	}

	//now going through the juvenile females
    const unsigned int totalYF = GetLiveArraySize(vob_JuvenileFemale);
	VoleSum+=totalYF;
	for (unsigned j = 0; j < totalYF; j++){	
		auto JF = dynamic_cast<Vole_JuvenileFemale*>(SupplyAnimalPtr(vob_JuvenileFemale,j));
		int nowx= JF -> Supply_m_Location_x(); //getting her position
		int nowy= JF -> Supply_m_Location_y();

		int gx = static_cast<int>(floor(nowx / sqwidth)); //getting the small quadrant she is in
		int gy = static_cast<int>(floor(nowy / sqwidth));

		int GridCellIndex; //let's put her in a big quadrant now
		GridCellIndex=Nine_AssignQuadrant(gx,gy,DivisionPoint); //assign her to a quadrant
        int nowquadrant=GridCellIndex;
		PopInQuadrants[nowquadrant][1]+=1; //add a female

	}
	//now going through the adult females
    const unsigned int totalAF = GetLiveArraySize(vob_Female);
	VoleSum+=totalAF;
	for (unsigned j = 0; j < totalAF; j++){	
		auto VF = dynamic_cast<Vole_Female*>(SupplyAnimalPtr(vob_Female,j));
		int nowx= VF -> Supply_m_Location_x(); //getting her position
		int nowy= VF -> Supply_m_Location_y();

		int gx = static_cast<int>(floor(nowx / sqwidth)); //getting the small quadrant she is in
		int gy = static_cast<int>(floor(nowy / sqwidth));

		int GridCellIndex; //let's put her in a big quadrant now
		GridCellIndex=Nine_AssignQuadrant(gx,gy,DivisionPoint);
        int nowquadrant=GridCellIndex;
		PopInQuadrants[nowquadrant][1]+=1; //add a female

	}
	//now we'll go through the adult males
    const unsigned int totalAM = GetLiveArraySize(vob_Male);
	VoleSum+=totalAM;
	for (unsigned j = 0; j < totalAM; j++){	
		auto VM = dynamic_cast<Vole_Male*>(SupplyAnimalPtr(vob_Male,j));
		int nowx= VM -> Supply_m_Location_x(); //getting his position
		int nowy= VM -> Supply_m_Location_y();

		int gx = static_cast<int>(floor(nowx / sqwidth)); //getting the small quadrant he is in
		int gy = static_cast<int>(floor(nowy / sqwidth));

		int GridCellIndex; //let's put him in a big quadrant now, and exclude any voles in the buffer zones.
		GridCellIndex=Nine_AssignQuadrant(gx,gy,DivisionPoint);
        int nowquadrant=GridCellIndex;
		PopInQuadrants[nowquadrant][0]+=1; //add a male
	}
	PopInQuadrants[9][0]=ExcludedVolesCount[0]; //excluded males (will be zero in this setup of no buffer zones)
	PopInQuadrants[9][1]=ExcludedVolesCount[1]; //excluded females (will be zero in this setup of no buffer zones)

	for (int i = 0; i <= 8; ++i) { // Iterate from 0 to 8
		PopInQuadrants[10][0] += PopInQuadrants[i][0]; // Accumulate the sum for total males in quadrants
		PopInQuadrants[10][1]+= PopInQuadrants[i][1]; // Accumulate the sum for total females in quadrants
	}
	PopInQuadrants[11][0]=PopInQuadrants[10][0]+PopInQuadrants[9][0]; //sum males with excluded voles
	PopInQuadrants[11][1]=PopInQuadrants[10][1]+PopInQuadrants[9][1]; //sum females with excluded voles
	double SumAllWithExcluded=PopInQuadrants[11][0]+ PopInQuadrants[11][1]; //sum of males and females in all quadrants+ the nonexistant buffer zones
	double SumAllWithoutExcluded=PopInQuadrants[10][0]+PopInQuadrants[10][1]; //sum of males and females in all quadrants
	
	//now: let's output resulting counts of males and females in each quadrant. 
	int day= ((year-1)*365)+a_day;
	pop_in_quadrants_out << day << ",";
	for (int i = 0; i < PopInQuadrants.size(); ++i) {
        pop_in_quadrants_out<< PopInQuadrants[i][0] << "," << PopInQuadrants[i][1] <<",";
    }
	//and output the sums
	pop_in_quadrants_out << SumAllWithExcluded << "," << SumAllWithoutExcluded << "," << DivisionDirection << "\n";
};
//This function outputs the amount of (maternal) generations have passed for each vole
void Vole_Population_Manager::GenerationCountOutput(const int year,const int a_day) {
    TAnimal* VB;
	ofstream GenerationOut("GenerationCount.txt", ios::app); //make the output file
	if (year==1 && a_day==1){ //if it's the first year and first day, put headers/titles into output file
		GenerationOut << "year,generation\n";
	}
	// go through juvenile males and "ask them" what their generation count is
	const unsigned int totalYM = GetLiveArraySize(vob_JuvenileMale);
	for (unsigned j = 0; j < totalYM; j++){
		float SampleOrNot=g_rand_uni_fnc();
		if (SampleOrNot>0.9){
			auto JM = dynamic_cast<Vole_JuvenileMale*>(SupplyAnimalPtr(vob_JuvenileMale,j));
			int Generation= JM -> GenerationCount;
			GenerationOut <<year<<","<<Generation<<"\n"; 
		}
	}

	// go through juvenile females and "ask them" what their generation count is
    const unsigned int totalYF = GetLiveArraySize(vob_JuvenileFemale);
	for (unsigned j = 0; j < totalYF; j++){
		float SampleOrNot=g_rand_uni_fnc();
		if (SampleOrNot>0.9){
			auto JF = dynamic_cast<Vole_JuvenileFemale*>(SupplyAnimalPtr(vob_JuvenileFemale,j));
			int Generation= JF -> GenerationCount;
			GenerationOut <<year<<","<<Generation<<"\n"; 
		}
	}

	// go through adult females and "ask them" what their generation count is
    const unsigned int totalAF = GetLiveArraySize(vob_Female);
	for (unsigned j = 0; j < totalAF; j++){	
		float SampleOrNot=g_rand_uni_fnc();
		if (SampleOrNot>0.9){
			auto VF = dynamic_cast<Vole_Female*>(SupplyAnimalPtr(vob_Female,j));
			int Generation= VF -> GenerationCount;
			GenerationOut <<year<<","<<Generation<<"\n";
		}
	}
	// go through adult males and "ask them" what their generation count is
    const unsigned int totalAM = GetLiveArraySize(vob_Male);
	for (unsigned j = 0; j < totalAM; j++){	
		float SampleOrNot=g_rand_uni_fnc();
		if (SampleOrNot>0.9){
			auto VM = dynamic_cast<Vole_Male*>(SupplyAnimalPtr(vob_Male,j));
			int Generation= VM -> GenerationCount;
			GenerationOut <<year<<","<<Generation<<"\n";
		}
	}

};
// This function outputs the maternal and paternal lineage ID for each vole at the beginning of the year at a specified sample rate..
void Vole_Population_Manager::LineagesOutput(const int year,const int a_day, const float sample_rate) {
    TAnimal* VB;
	ofstream GenerationOut("LineagesOut.txt", ios:: app); //make output file
	if (year==1 && a_day==1){ // if the year and day are 1, put headers/titles into the output file
		GenerationOut << "year,mito_line,ychrom_line\n";
	}
	//go through the juvenile males and "ask" them about their maternal and paternal lineage
	const unsigned int totalYM = GetLiveArraySize(vob_JuvenileMale);
	for (unsigned j = 0; j < totalYM; j++){
		float SampleOrNot=g_rand_uni_fnc();//check whether to sample or not
		if (SampleOrNot > (1 - sample_rate)){
			auto JM = dynamic_cast<Vole_JuvenileMale*>(SupplyAnimalPtr(vob_JuvenileMale,j));
			int mito= JM -> new_Genes.MitochondrialLine;
			int ychr = JM -> new_Genes. YChromoLine;
			GenerationOut <<year<<","<<mito<<","<<ychr<<"\n"; 
		}
	}

	//go through the juvenile females and "ask" them about their maternal and paternal lineage
    const unsigned int totalYF = GetLiveArraySize(vob_JuvenileFemale);
	for (unsigned j = 0; j < totalYF; j++){
		float SampleOrNot=g_rand_uni_fnc();//check whether to sample or not
		if (SampleOrNot > (1 - sample_rate)){
			auto JF = dynamic_cast<Vole_JuvenileFemale*>(SupplyAnimalPtr(vob_JuvenileFemale,j));
			int mito= JF -> new_Genes.MitochondrialLine;
			int ychr = JF -> new_Genes. YChromoLine;
			GenerationOut <<year<<","<<mito<<","<<ychr<<"\n"; 
		}
	}

	//go through the adult females and "ask" them about their maternal and paternal lineage
    const unsigned int totalAF = GetLiveArraySize(vob_Female);
	for (unsigned j = 0; j < totalAF; j++){	
		float SampleOrNot=g_rand_uni_fnc();//check whether to sample or not
		if (SampleOrNot > (1 - sample_rate)){
			auto VF = dynamic_cast<Vole_Female*>(SupplyAnimalPtr(vob_Female,j));
			int mito= VF -> new_Genes.MitochondrialLine;
			int ychr = VF -> new_Genes. YChromoLine;
			GenerationOut <<year<<","<<mito<<","<<ychr<<"\n"; ;
		}
	}
	
	//go through the adult males and "ask" them about their maternal and paternal lineage
    const unsigned int totalAM = GetLiveArraySize(vob_Male);
	for (unsigned j = 0; j < totalAM; j++){	
		float SampleOrNot=g_rand_uni_fnc();//check whether to sample or not
		if (SampleOrNot > (1 - sample_rate)){
			auto VM = dynamic_cast<Vole_Male*>(SupplyAnimalPtr(vob_Male,j));
			int mito= VM -> new_Genes.MitochondrialLine;
			int ychr = VM -> new_Genes. YChromoLine;
			GenerationOut <<year<<","<<mito<<","<<ychr<<"\n";
		}
	}
};

int Vole_Population_Manager::Four_AssignQuadrant(int gx, int gy, int DivisionPoint){
	int GridCellIndex;
	if (gx<DivisionPoint){
		if (gy<DivisionPoint){
			GridCellIndex=0;
		}else if(gy>DivisionPoint){
			GridCellIndex=2;
		}else{
			GridCellIndex=100; //indicator number that a vole was excluded
		}
	}else if(gx>DivisionPoint){
		if(gy<DivisionPoint){
			GridCellIndex=1;
		}else if(gy>DivisionPoint){
			GridCellIndex=3;
		}else{
			GridCellIndex=100; //indicator number that a vole was excluded
		}
	}else{
		GridCellIndex=100; //indicator number that a vole was excluded
	}
	return(GridCellIndex);
}

int Vole_Population_Manager::Nine_AssignQuadrant(int gx, int gy, int DivisionPoint){
	int GridCellIndex;
	if (gx<(DivisionPoint-1)){ //assignment made for a 3*3 grid :)
		if (gy<(DivisionPoint-1)){
			GridCellIndex=0;
		}else if (gy>(DivisionPoint+1)){
			GridCellIndex=6;
		}else{
			GridCellIndex=3;
		}
	}else if(gx>(DivisionPoint+1)){
		if (gy<(DivisionPoint-1)){
			GridCellIndex=2;
		}else if (gy>(DivisionPoint+1)){
			GridCellIndex=8;
		}else{
			GridCellIndex=5;
		}
	}else{
		if (gy<(DivisionPoint-1)){
			GridCellIndex=1;
		}else if (gy>(DivisionPoint+1)){
			GridCellIndex=7;
		}else{
			GridCellIndex=4;
		}
	}
	return(GridCellIndex);
}

