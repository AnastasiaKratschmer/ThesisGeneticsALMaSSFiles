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
extern std::string filename_AlleleInput;
extern int ChromosomeCounting;
extern int LocusCounting;
extern int DiffBaseCounting;

#ifndef VolePopulationManagerH
#define VolePopulationManagerH

//---------------------------------------------------------------------------

class struct_Vole_Adult;
class TAnimal;
class Vole_Male;
class Vole_Female;
class AlleleFreq1616;
class AlleleFreq;
class APoint;

//------------------------------------------------------------------------------

typedef vector<Vole_Male*> Vole_MaleList;
typedef vector<Vole_Female*> Vole_FemaleList;

/**
\brief
Types of vole message
*/
typedef enum {
	tovm_Infanticide = 0
}TTypeOfVoleMessage;

/**
\brief
Types of vole mortality
*/
typedef enum {
	tovmort_MStarve = 0,
	tovmort_FStarve,
	tovmort_MBck,
	tovmort_FBck,
	tovmort_MFarm,
	tovmort_FFarm,
	tovmort_MDisp,
	tovmort_FDisp,
	tovmort_MPred,
	tovmort_FPred,
	tovmort_MLife,
	tovmort_FLife,
	tovmort_MPest,
	tovmort_FPest
} TTypeOfVoleMortality;

/**
\brief
A class for simulation virtual traplines
*/
class TrapLineMap : public BinaryMapBase {
public:
	/** \brief TrapLineMap constructor */
	TrapLineMap(unsigned int a_width, unsigned int a_height, unsigned int a_resolution, const char* a_file);
	/** \brief TrapLineMap destructor */
	~TrapLineMap();
	/** \brief Is there a trap at this x,y? */
	bool IsTrap(const unsigned a_x, const unsigned a_y) {
		if (GetValue(a_x, a_y) == 1) return true;
		return false;
	}
	void Output(const InTrapPosition& a_tp, int a_day, int a_sex, bool a_terr, int a_age, int a_bx, int a_by,
	            int a_ID) const;

protected:
	/** \brief Reads in the trap coords, creates the map and sets up the output files */
	void Init(const char* a_inifile);
	/** \brief The total number of traps */
	unsigned m_noTraps;
	/** \brief List of trap co-ordinates */
	vector<APoint> m_TrapCoords;
	FILE* m_ofile;
};

/**
\brief
A base class for summary outputs
*/
class VoleSummaryOutput {
public:
	VoleSummaryOutput(const char* a_filename, Landscape* a_land, int a_numdataINT, int a_numdataDOUBLE);
	virtual ~VoleSummaryOutput();
	virtual void OPrint();
	virtual void OPrint(int a_value);
	virtual void OPrint(double a_value);
	virtual void OPrint(const char* a_value);
	virtual void OPrintEndl();
	void ResetData();
	void ChangeData(int a_data, int a_value);
	void ChangeData(int a_data, double a_value);

protected:
	ofstream* m_File;
	Landscape* m_landscape;
	int m_ndInt;
	int m_ndDouble;
	int m_dataI[100]; // This is way more than we should ever need
	double m_dataD[100]; // This is way more than we should ever need
	void OpenOutput(const char* a_filename);
	void CloseOutput() const;
};

/**
\brief
The class to handle all vole population related matters
*/
class Vole_Population_Manager : public Population_Manager {
protected:
	int thisYear, YearsTotal;
	void DoFirst() override;
	void Catastrophe() override;
	int m_GrowthStartDate;
	int m_male1;
	int m_female1;
	int m_total1;

	int m_male2;
	int m_female2;
	int m_total2;

	int m_male3;
	int m_female3;
	int m_total3;
	//double m_BaseVoleDensity;
	TrapLineMap* m_Traplines;
	vector<double> VoleHabitatBaseQualities;
	bool OpenSexRatiosProbe();
	void CloseSexRatiosProbe();
	bool SuitableStartingLocation(int a_x, int a_y) const;
	double AssessHabitat(int p_Polyref) const;

#ifdef __VOLEPESTICIDEON
public:
	void AddToImpacted() {m_impacted++;}
    void AddToNotImpacted() {m_notimpacted++;}
    void AddToGeneticImpacted() {m_geneticimpacted++;}
	void incLittersLost() { m_LittersLost++; }
	void incLittersProduced() { m_LittersProduced++; }
protected:
	int m_impacted; // counter for voles hit and affected by pesticide
    int m_notimpacted; // counter for voles not hit and affected by pesticide
    int m_geneticimpacted; // counter for voles not hit and affected by pesticide
	double m_geneticsterilitychance;
	double m_f1sterilitychance;
	int m_LittersLost;
	int m_LittersProduced;
#endif
	ofstream* m_VoleAgeSexLocationFile;
	ofstream* m_VoleResistanceOutputFile;

public:
	VoleSummaryOutput* m_VoleRecordMort;
	FILE* SexRatiosPrb;
	/** \brief used to get a specific vole allele from outside the population manager */
	uint32 GetVoleAllele(const int a_list, const int a_vole, const int a_locus, const int a_chromosome) {
		const auto new_VB = dynamic_cast<Vole_Base*>(SupplyAnimalPtr(a_list, a_vole));
		return new_VB->SupplyMyAllele(a_locus, a_chromosome);
	}
	void ImpactedProbe() override;
	bool IsTrap(int p_x, int p_y) const {
		m_TheLandscape->CorrectCoords(p_x, p_y);
		return m_Traplines->IsTrap(p_x, p_y);
	}
	//virtual void Catastrophe(int a_mort);
	void TheAgeSexLocationProbe();
	void TheAOROutputProbe() override;
	void TheRipleysOutputProbe(ofstream* a_prb) override;
	void TheReallyBigOutputProbe() override;
	virtual void TheSexRatiosProbe();
	virtual void LandscapeQuadrantOutputProbe(int a_day);
	void AddToFrag1Male() { m_male1++; } // ***TD***
	void AddToFrag1Female() { m_female1++; }
	void AddToFrag2Male() { m_male2++; }
	void AddToFrag2Female() { m_female2++; }
	void AddToFrag3Male() { m_male3++; }
	void AddToFrag3Female() { m_female3++; }
	int ReproTable[4][12]; // Filled in by init function
	unsigned IDNumber;
	/** Used to store lists of females for further use */
	Vole_FemaleList FList;
	/** Used to store lists of males for further use */
	Vole_MaleList MList;
	Vole_Population_Manager(Landscape* p_L);
	~Vole_Population_Manager() override;
	void CreateObjects(VoleObject a_obType, TAnimal* a_pvo, struct_Vole_Adult* a_as, int a_number, const string& filename_AlleleInput, bool isOffspring=false);
	void CreateObjectsInit(VoleObject a_obType, TAnimal* a_pvo, struct_Vole_Adult* a_voleStruct, int a_number, const string& filename_AlleleInput, int id=-1);
	virtual void Init();
	static bool RecordGeneticsToday(int p_today, int p_year, int p_start_year, int p_interval);
	double GetHabitatQuality(const int a_index) const { return VoleHabitatBaseQualities[a_index]; }
	bool SupplyOlderFemales(unsigned p_x, unsigned p_y, unsigned p_Age, unsigned p_range) const;
	int SupplyHowManyVoles(unsigned a_x, unsigned a_y, unsigned a_range) const;
	int SupplyGrowthStartDate() const { return m_GrowthStartDate; }
	int SupplyInOlderTerr(unsigned p_x, unsigned p_y, unsigned p_Age, unsigned p_Range) const;
	int SupplyCountFemales(unsigned p_x, unsigned p_y, unsigned p_TerrRange) const;
	void SendMessage(TTypeOfVoleMessage p_message, unsigned p_x, unsigned p_y, unsigned p_range, bool p_sex,
	                 double p_Weight /*, unsigned p_IDNo*/);
	vector<Vole_Base*>* SupplyVoleList(unsigned x, unsigned y, unsigned range) const;
	bool InSquare(int p_x, int p_y, int p_sqx, int p_sqy, int p_range) const;
	int ListClosestFemales(int p_x, int p_y, int p_steps);
	int ListClosestMales(int p_x, int p_y, int p_steps);
	Vole_Female* FindClosestFemale(int p_x, int p_y, int p_steps) const;
	Vole_Male* FindClosestMale(int p_x, int p_y, int p_steps) const;
	Vole_Male* FindRandomMale();
	void AddToYoung(const int young) { YoungProducedToday += young; }
	Vole_Male* FindWithinRadiusMale(int p_x, int p_y) const;
	Vole_Male* FindOutsideRadiusMale(int p_x, int p_y); // ***TD***
	bool BarrierSearch(int F_x, int F_y, int M_x, int M_y) const; // ***TD***
	void AddToYoungYr(const int young) { YoungProducedToday += young; }
	void AddToJuvs(const int juvs) { JuvsProducedToday += juvs; }
	void AddToNoYoungInfanticideCount(const int m_NoOfYoung, const int m_YoungAge) {
		// ***TD***
		NoYoungKilledToday += m_NoOfYoung;
		if (m_YoungAge < 5) { NoYoungKilledToday4 += m_NoOfYoung; }
		else if (m_YoungAge < 9) { NoYoungKilledToday8 += m_NoOfYoung; }
		else { NoYoungKilledToday9 += m_NoOfYoung; }
	}
	void ReproductionProbe() const;
	/* THIS FUNCTION IS OUTDATED; BUT MAYBE A REPLACEMENT IS NEEDED!!
	void TheGeneticProbe(unsigned listindex,int ProbeNo,unsigned &SubSize);
	*/
	void GeneticsOutputFile(unsigned listindex);
	void GeneticsResultsOutput(FILE* ofile, unsigned listindex) override;
	void GeneticsResultsOutputNEW(FILE* ofile, unsigned listindex);
	/** \brief Opens the output file ready for resistance results */
	void OpenResistanceOutput();
	/** \brief Closes the output file ready for resistance results */
	void CloseResistanceOutput() const;
	/** \brief Resistance results file output  */
	void ResistanceOutput();
	void Four_QuadrantBasedGeneticOutput(const int a_day, const int year, const string landscape_info);
	void Nine_QuadrantBasedGeneticOutput(const int a_day, const int year, const string landscape_info);
	void FourQuadrantsPopulationSizeProbe(const int today, const int year, const string landscape_info);
	void NineQuadrantsPopulationSizeProbe(const int a_day, const int year, const string landscape_info);
	void MakeAlleleInputFile(const string params_for_beta0, const string params_for_beta1);
	void CompareGenomes(const int year,const string landscape_info);
	void GenerationCountOutput(const int year,const int a_day);
	void LineagesOutput(const int year,const int a_day, const float sample_rate);
	int Four_AssignQuadrant(int gx, int gy, int DivisionPoint);
	int Nine_AssignQuadrant(int gx, int gy, int Divisionpoint);

	// Attributes
	AlleleFreq* AFreq;
	int YoungProducedToday, JuvsProducedToday;
	int NoYoungKilledToday, NoYoungKilledToday4, NoYoungKilledToday8, NoYoungKilledToday9;
	FILE* YoungsFile;
	IDMap<TAnimal*>* m_VoleMap;
	int m_geneticproductfertilityeffect;
};

//------------------------------------------------------------------------------
#endif
