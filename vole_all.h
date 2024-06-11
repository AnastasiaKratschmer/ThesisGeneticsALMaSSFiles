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
<B>Vole_all.h This header file contains the code for all vole lifestage classes</B>\n
*/
/**
\file 
 by Chris J. Topping \n
 Version of 01st Nov 2000 \n
 With additions as noted in: \n
 April 2006 \n
 November 2007 \n
 Doxygen formatted comments in May 2008 \n
*/

//---------------------------------------------------------------------------
#ifndef Vole_allH
#define Vole_allH

#include "../Vole/GeneticMaterial.h"

class MovementMap16;
extern std::string filename_AlleleInput;
extern int ChromosomeCounting;
extern int LocusCounting;
extern int DiffBaseCounting;



//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

/**
\brief
Vole behavioural states
*/
typedef enum {
	tovs_InitialState = 0,
	//Males
	tovs_JuvenileExploration,
	tovs_MMaturation,
	tovs_MEvaluateExplore,
	tovs_Infanticide,
	tovs_MDying,
	//Females
	tovs_FEvaluateExplore,
	tovs_ReproBehaviour,
	tovs_Lactating,
	tovs_GiveBirth,
	tovs_FMaturation,
	tovs_Mating,
	tovs_UpdateGestation,
	tovs_SpecialExplore,
	tovs_FDying,
}TTypeOfVoleState;


typedef enum {
	vdisp_CarryOn = 0,
	vdisp_Die,
	vdisp_Mature,
	vdisp_Infanticide
}VoleDispersalReturns;
//------------------------------------------------------------------------------

// Forward Declarations

class Vole_Population_Manager;
class Landscape;
class AlleleFreq1616;
class GeneticMaterial1616;
class GeneticMaterial256_16;
class AlleleFreq;
class GeneticMaterial;

//------------------------------------------------------------------------------

typedef enum {
	vob_JuvenileMale = 0,
	vob_JuvenileFemale,
	vob_Male,
	vob_Female,
	vob_foobar
} VoleObject;

/**
\brief
A class for storing the position of the trap the vole is in
*/
class InTrapPosition : public AnimalPosition {
public:
	bool m_inAtrap;
};

/**
\brief
A struct for passing data to create a new vole
*/
class struct_Vole_Adult {
public:
	int x;
	int y;
	int age;
	double weight;
	int xborn;
	int yborn;
	int GenerationCount;
	int PolyRefBorn;
	int ElemBorn;
	int VegBorn;
	int BirthYear;
	int FatherId;
	int MotherId;
	int FatherStateAtBirth;
	Landscape* L;
	Vole_Population_Manager* VPM;
	//GeneticMaterial256_16 Genes256_16; 
	GeneticMaterial Genes;
	bool m_flag;
	bool m_gflag;
	bool m_dflag;
	double misc_use; // used to pass a parameter e.g. for pesticide use
	//GeneticMaterialNEW new_Genes=;
	//GeneticMaterialNEW new_Mates_Genes=;
	GeneticMaterialNEW new_Genes = GeneticMaterialNEW();
	GeneticMaterialNEW new_Mates_Genes = GeneticMaterialNEW();
	//GeneticMaterialNEW new_Genes;
	//GeneticMaterialNEW new_Mates_Genes;
};

//------------------------------------------------------------------------------
/**
\brief
Base class for voles - all vole objects are descended from this class. 
*/
class Vole_Base : public TAnimal, public struct_Vole_Adult {
	// Attributes
//public:
//	GeneticMaterialNEW new_Genes;
//	GeneticMaterialNEW new_Mates_Genes;
protected:
	/** The year of birth*/
	int m_BirthYear;
	/**
	\brief
	A flag set if the female was born the year before.
	*/
	bool m_BornLastYear;
	/** The minimum territory range*/
	unsigned int m_MinTerrRange;
	/** The id number of the mother*/
	unsigned m_MotherId;
	/** The id number of the fater*/
	unsigned m_FatherId;
	/** Death course*/
	int m_Death;
	/** The size of their territory (radius of a square) */
	int m_TerrRange;
	/** Their sex Male==true Female==false */
	bool m_Sex;
	/** Whether they are mature or not */
	bool m_Mature;
	/** Their age in days */
	int m_Age;
	/** Their x location at birth*/
	int m_XBorn; //*TD*
	/** Their y location at birth*/
	int m_YBorn; //*TD*
	/** Their polygon ref at birth location*/
	int m_PolyRefBorn;
	/** The element type at at birth location*/
	int m_ElemBorn; //*TD*
	/** The vegetation type at birth location*/
	int m_VegBorn; //*TD*
	/** Their lifespan remaining (unless killed by external events) */
	int m_LifeSpan;
	/** Their weight in grams */
	double m_Weight;
	/**	\brief	Flag indicating the fertility state (true means fertile)	*/
	bool m_fertile;
#ifdef __VoleStarvationDays
	/** How many days they have been starving */
    int m_StarvationDays;   
#endif
	// No of young produced;
	int m_NoOfYoungTotal;
	/** The current dispersal direction */
	int m_DispVector;
	/** Do they have a terriory? */
	bool m_Have_Territory;
	/** Their reserves - in days that they can survive without food */
	int m_Reserves;
	/** Their individual ID number */
	unsigned IDNo;
	/** The size of simulation landscape */
	int m_SimH, m_SimW;
	/** Their genes */
	//GeneticMaterial256_16 MyGenes;
	GeneticMaterial m_MyGenes;
	InTrapPosition m_intrappos;
	/** Maximum territory size male*/
	static unsigned int m_MaxMaleTerritorySize;
	/** Maximum territory size female*/
	static unsigned int m_MaxFemaleTerritorySize;
	/** Minimum territory size male*/
	static unsigned int m_MinMaleTerritorySize;
	/** Minimum territory size female*/
	static unsigned int m_MinFemaleTerritorySize;
	/** Minimum acceptable habitat quality - assumes that minimum hab qual for survival is 2 * min female terr size */
	static double m_MinFVoleHabQual;
	/** Minimum acceptable habitat juvenile male quality - assumes that minimum hab qual for survival is 2 * min male terr size */
	static double m_MinJMVoleHabQual;
	/** Minimum acceptable habitat quality - assumes that minimum hab qual for survival is 2 * min female terr size */
	static double m_MinMVoleHabQual;
	/** only used to create an increasing territory size with age */
	static double m_MaleTerritoryRangeSlope;
	/** only used to create an increasing territory size with age */
	static double m_FemaleTerritoryRangeSlope;
	// The values below were fitted by trial and error but await a detailed sensitivity analysis when spatial data is available
	/** Habitat is v.good above this */
	static double m_FHabQualThreshold3;
	/** Habitat is OK above this */
	static double m_FHabQualThreshold2;
	/** Habitat is v.bad below this */
	static double m_FHabQualThreshold1;
	/** Habitat is v.good above this */
	static double m_MHabQualThreshold3;
	/** Habitat is OK above this */
	static double m_MHabQualThreshold2;
	/** Habitat is v.bad below this */
	static double m_MHabQualThreshold1;
	/** Local storage for whether it is breeding season */
	static bool m_BreedingSeason;
#ifdef __VoleStarvationDays
	/** The number of days without suitable food the vole will live - this primarily entails risk of predation whilst moving rather than actual startvation.\n
	The value 21 was fitted based on expert knowledge applied to evaluating individual vole dispersal movements. It has a large influence on colonisation ability!*/
	static int m_MaxStarvationDays;  
#endif

#ifdef __VOLEPESTICIDEON
public:
	double SupplyPesticideLoad() { return m_pesticideload; }
	double SupplyBioDegradeRate() { return m_pesticideBioDegradeRate; }
	void SetPesticideInfluenced1(bool a_bool) { m_pesticideInfluenced1 = a_bool; }
	void SetPesticideInfluenced2(bool a_bool) { m_pesticideInfluenced2 = a_bool; }
	bool GetPesticideInfluenced1() { return m_pesticideInfluenced1; }
	bool GetPesticideInfluenced2() { return m_pesticideInfluenced2; }
protected:
	bool m_pesticideInfluenced1;
	bool m_pesticideInfluenced2;
	double m_pesticideInfluenced3;
	double m_pesticideload;
	double m_pesticideBioDegradeRate;
	void AddToPesticidLoad(double a_pestamount) { m_pesticideload+=a_pestamount; }
	void BioDegradePesticide() { m_pesticideload *= m_pesticideBioDegradeRate; }
	void ClearPesticidLoad() { m_pesticideload = 0.0; }
	void SetBioDegradeRate(double a_rate) { m_pesticideBioDegradeRate = a_rate; }
	virtual void PesticideIngestion( void );
	virtual void ActOnPesticideDose( void );
	virtual void ModelinkPesticide() { ; }
	virtual void ModelinkPesticide21TWA(double /* a_dose */ ) { ; }
	virtual void Vinclozolin(double /* a_dose */ ) { ; }
	virtual void GeneralOrganophosphate(double /* a_dose */) { ; }
	virtual void GeneralEndocrineDisruptor(double /* a_dose */) { ; }	
	virtual void GeneticDemoPesticide(double /* a_dose */){ ; }
#endif
#ifdef __VOLERODENTICIDEON
	virtual void RodenticideIngestion(void);
#endif

public:
	// Functions
	Vole_Base(struct_Vole_Adult* a_AVoleStruct_ptr,const string& filename_AlleleInput, int i=-1);
	~Vole_Base() override;
	virtual void Init(struct_Vole_Adult* a_AVoleStruct_ptr, const string& filename_AlleleInput, int i=-1);
	virtual void ReInit(struct_Vole_Adult* a_AVoleStruct_ptr, const string& filename_AlleleInput, int i=-1);

	void BeginStep() override {
	};

	void Step() override {
	};
	void EndStep() override;
	void st_Dying();

	/** \brief Set Breeding Season flag */
	void SetBreedingSeason(bool a_flag) { m_BreedingSeason = a_flag; }
	/** Set our weight */
	void SetWeight(double W) { m_Weight = W; }
	/** Become adult */
	void Setm_Mature() { m_Mature = true; }
	/** Set year of birth */
	void Set_BirthYear(int BirthYear) { m_BirthYear = BirthYear; }

	/** Set mother ID */
	void Set_MotherId(unsigned MotherIdNo) { m_MotherId = MotherIdNo; }; // ***TD***
	/** Set father ID */
	void Set_FatherId(unsigned FatherIdNo) { m_FatherId = FatherIdNo; }; // ***TD***
	/** Sets the number of young produced by the female */
	void Set_NoYoungTot(int a_NoOfYoung) { m_NoOfYoungTotal += a_NoOfYoung; }; // ***TD***
	/** Set x-coordinate of birth location */
	void Set_XBorn(int a_Location_x) { m_XBorn = a_Location_x; }; // ***TD***
	/** Set y-coordinate of birth location */
	void Set_YBorn(int a_Location_y) { m_YBorn = a_Location_y; }; // ***TD***
	/** Set element type of birth location */
	void Set_ElemBorn(int a_Location_x, int a_Location_y) // ***TD***
	{
		TTypesOfLandscapeElement const m_EBorn = m_OurLandscape->SupplyElementType(a_Location_x, a_Location_y);
		m_ElemBorn = m_OurLandscape->BackTranslateEleTypes(m_EBorn);
	}

	/** Set vegetation type of birth location */
	void Set_VegBorn(int a_Location_x, int a_Location_y) // ***TD***
	{
		TTypesOfVegetation const m_VBorn = (m_OurLandscape->SupplyVegType(a_Location_x, a_Location_y));
		m_VegBorn = (m_OurLandscape->BackTranslateVegTypes(m_VBorn));
	}

	/** Set polygonref of birth location */
	void Set_PolyRefBorn(int a_Location_x, int a_Location_y) { m_PolyRefBorn = (m_OurLandscape->SupplyPolyRef(a_Location_x, a_Location_y)); }; // ***TD***
	/** Set our age */
	void Set_Age(int Age) { m_Age = Age; }

	/** Get our current vole state */
	int WhatState() override { return CurrentVState; }
	/** \brief Were we born this year? */
	bool SupplyBornLastYear() { return m_BornLastYear; }
	/** Tell whether we have a territory */
	bool SupplyTerritorial() { return m_Have_Territory; } // ***TD***
	/** Tell father ID */
	int SupplyFatherId() { return m_FatherId; }; // ***TD***
	/** Tell mother ID*/
	int SupplyMotherId() { return m_MotherId; }; // ***TD***
	/** Tell our sex */
	bool SupplySex() { return m_Sex; };
	/** Tell our birth year */
	int SupplyBirthYear() { return m_BirthYear; }; // ***TD***
	/** Tell our number of young produced */
	int SupplyTotNoYoung() { return m_NoOfYoungTotal; }; // ***TD***
	/** Tell our x-coordinate at birth */
	int SupplyXBorn() { return m_XBorn; }; // ***TD***
	/** Tell our x-coordinate at birth */
	int SupplyYBorn() { return m_YBorn; }; // ***TD***
	/** Tell our polygon ref at birth location */
	int SupplyPolyRefBorn() { return m_PolyRefBorn; }; // ***TD***
	/** Tell our elementype at birth location */
	int SupplyElemBorn() { return m_ElemBorn; }; // ***TD***
	/** Provide our current elementype */
	TTypesOfLandscapeElement SupplyElemType() { return m_OurLandscape->SupplyElementType(m_Location_x, m_Location_y); };
	/** Tell our vegetation type at birth location */
	int SupplyVegBorn() { return m_VegBorn; }; // ***TD***
	/** Tell our territory range */
	int SupplyTerrRange() { return m_TerrRange; }; // ***TD***
	/** Tell our weight */
	double SupplyWeight() { return m_Weight; };
	/** Tell our ID number */
	int SupplyIDNo() { return IDNo; }; // ***TD***
	/** Tell if mature */
	bool SupplyMature() { return m_Mature; };
	/** Tell the cause of death */
	int SupplyDeathCause() { return m_Death; }; // ***TD***
	/** Tell our age */
	unsigned SupplyAge() { return m_Age; };
	/** Tell our x coordinate */
	unsigned SupplyX() { return m_Location_x; };
	/** Tell our y coordinate */
	unsigned SupplyY() { return m_Location_y; };
	/** Are we in a trap? */
	bool SupplyInTrap() { return m_intrappos.m_inAtrap; }
	/** Get the trap location */
	InTrapPosition SupplyTrapPosition() { return m_intrappos; }
	/** Release vole from trap */
	void SetFree() { m_intrappos.m_inAtrap = false; }
	/** Genetic functionality */
	int SupplyHomoZyg() { return m_MyGenes.HomozygosityCount(); }
	/** Genetic functionality */
	int SupplyHeteroZyg() { return m_MyGenes.HeterozygosityCount(); }
	/** Genetic functionality */
	int SupplyAllele(int locus, int allele) { return m_MyGenes.GetAllele(locus, allele); }
	uint32 SupplyMyAllele(int i, int j) { return m_MyGenes.GetAllele(i, j); } //*TD
	/** Genetic functionality */
	int GetGeneticFlag() { return m_MyGenes.GetGeneticFlag(); }
	/** Genetic functionality */
	int GetDirectFlag() { return m_MyGenes.GetDirectFlag(); }
	/** Genetic functionality */
	void SetGeneticFlag() { m_MyGenes.SetGeneticFlag(); }
	/** Genetic functionality */
	void SetDirectFlag() { m_MyGenes.SetDirectFlag(); }
	/** Genetic functionality */
	void UnsetGeneticFlag() { m_MyGenes.UnsetGeneticFlag(); }
	/** Genetic functionality */
	void UnsetDirectFlag() { m_MyGenes.UnsetDirectFlag(); }
	/** Genetic functionality */
	GeneticMaterial SupplyGenes() { return m_MyGenes; }

	virtual void OnKilled() {
	};
	virtual bool MortalityTest();
	/** Our current behavioural state */
	TTypeOfVoleState CurrentVState;
	Vole_Population_Manager* m_OurPopulation;
	void CopyMyself(VoleObject a_vole);
	/**
	\brief
	Set the male vole fertility 
	*/
	void SetFertile(bool f) {
		/**
		Primarily used in ecotoxicological simulations where toxic effects may render the vole sterile.
		*/
		m_fertile = f;
	}

	/**
	\brief
	Get the male vole fertility 
	*/
	bool GetFertile() {
		/**
		Primarily used in ecotoxicological simulations where toxic effects may render the vole sterile.
		*/
		return m_fertile;
	}

protected:
	double CalculateCarryingCapacity(int x, int y, int a_ddep) const;
	//double CalculateCarryingCapacity(int x,int y,int &p_stand_x,int &p_stand_y);
	int MoveQuality(int p_x, int p_y) const;
	void MoveTo(int p_Vector, int p_Distance, int iterations);
	void DoWalking(int p_Distance, int& p_Vector, int& vx, int& vy) const;
	void DoWalkingCorrect(int p_Distance, int& p_Vector, int& vx, int& vy) const;
	void Escape(int p_Vector, int p_Distance);
	void CheckTraps();

	virtual void SetLocation() {
	};

	virtual void FreeLocation() {
	};
	virtual bool GetLocation(int /*unused*/, int /*unused*/) { return false; };
};

//------------------------------------------------------------------------------

/**
\brief
The class for juvenile male voles
*/
/**
Contains all the behaviour specific to the male vole. Only st_Infanticide and st_JuvenileExplore are specific to the male, other behaviours differ only in details from the female.
*/
class Vole_JuvenileMale : public Vole_Base {
public:
	Vole_JuvenileMale(struct_Vole_Adult* p_aVoleStruct,const string& filename_AlleleInput, int i =-1);
	~Vole_JuvenileMale() override;
	void ReInit(struct_Vole_Adult* p_aVoleStruct, const string& filename_AlleleInput, int i =-1) override;
	bool OnFarmEvent(FarmToDo event) override;
	void BeginStep() override;
	void Step() override;
	void EndStep() override;
	void OnKilled() override;

protected:
	void SetLocation() override;
	void FreeLocation() override;
	bool GetLocation(int px, int py) override;
	TTypeOfVoleState Dispersal(double p_OldQual, int p_Distance);
	void DetermineTerritorySize();
	inline bool CanFeed();
	void st_JuvenileExplore(void);
	void st_BecomeSubAdult(void);
	TTypeOfVoleState st_Eval_n_Explore(void);
#ifdef __VOLEPESTICIDEON
	virtual void ModelinkPesticide();
	virtual void ModelinkPesticide21TWA(double a_dose);
	virtual void Vinclozolin(double a_dose);
	virtual void GeneralOrganophosphate(double a_dose);
	virtual void GeneralEndocrineDisruptor(double /* a_dose */);
	virtual void GeneticDemoPesticide(double /* a_dose */);

#endif
};

//------------------------------------------------------------------------------

/**
\brief
The class for male voles
*/
/**
Contains all the behaviour specific to the male vole. Only st_Infanticide and st_JuvenileExplore are specific to the male, other behaviours differ only in details from the female.
*/
class Vole_Male : public Vole_JuvenileMale {
public:
	Vole_Male(struct_Vole_Adult* p_aVoleStruct,const string& filename_AlleleInput, int i =-1);
	void ReInit(struct_Vole_Adult* p_aVoleStruct,const string& filename_AlleleInput, int i =-1) override;
	~Vole_Male() override;
	//virtual bool OnFarmEvent(FarmToDo event);
	//virtual void BeginStep();
	void Step() override;
	void EndStep() override;

protected:
	VoleDispersalReturns Dispersal(double p_OldQual, int p_Distance);
	void DetermineTerritorySize();
	inline bool CanFeed() const;
	int st_Maturation(void) const;
	void st_Infanticide(void) const;
	bool MortalityTest() override;
	TTypeOfVoleState st_Eval_n_Explore(void);
#ifdef __VOLEPESTICIDEON
	virtual void ModelinkPesticide();
	virtual void ModelinkPesticide21TWA(double a_dose);
	virtual void Vinclozolin(double a_dose);
	virtual void GeneralOrganophosphate(double a_dose);
	virtual void GeneralEndocrineDisruptor(double /* a_dose */);
	virtual void GeneticDemoPesticide(double /* a_dose */);
#endif
};

//------------------------------------------------------------------------------

/**
\brief
The class for female voles
*/
/**
Contains all the behaviour specific to the female vole. The differences between the male and female are primarily in female reproductive behaviour, but there are small differences in other behaviours requiring re-implementation of many of the behaviours (e.g. dispersal).
*/
class Vole_JuvenileFemale : public Vole_Base {
public:
	Vole_JuvenileFemale(struct_Vole_Adult* p_aVoleStruct,const string& filename_AlleleInput, int i =-1);
	void ReInit(struct_Vole_Adult* p_aVoleStruct,const string& filename_AlleleInput, int i =-1) override;
	~Vole_JuvenileFemale() override;
	bool OnFarmEvent(FarmToDo event) override;
	void BeginStep() override;
	void Step() override;
	void EndStep() override;
	void OnKilled() override;

protected:
	// Methods
	int Dispersal(double p_OldQual, int p_Distance);
	int st_Evaluate_n_Explore();
	void st_BecomeSubAdult(void);
	//unsigned HowManyOldFemales(int p_x, int p_y);
	void SetLocation() override;
	void FreeLocation() override;
	bool GetLocation(int px, int py) override;
#ifdef __VOLEPESTICIDEON
	double m_maturitydelay;
	virtual void ModelinkPesticide();
	virtual void ModelinkPesticide21TWA(double a_dose);
	virtual void Vinclozolin(double a_dose);
	virtual void GeneralOrganophosphate(double a_dose);
	virtual void GeneralEndocrineDisruptor(double /* a_dose */);
	virtual void GeneticDemoPesticide(double /* a_dose */);
public:
	void SetMaturityDelay(double a_delay) { m_maturitydelay = a_delay; }
#endif
};

//---------------------------------------------------------------------------

/**
\brief
The class for female voles
*/
/**
Contains all the behaviour specific to the female vole. The differences between the male and female are primarily in female reproductive behaviour, but there are small differences in other behaviours requiring re-implementation of many of the behaviours (e.g. dispersal).
*/
class Vole_Female : public Vole_JuvenileFemale {
public:
	Vole_Female(struct_Vole_Adult* p_aVoleStruct,const string& filename_AlleleInput, int i =-1);
	void ReInit(struct_Vole_Adult* p_aVoleStruct,const string& filename_AlleleInput, int i =-1) override;
	~Vole_Female() override;
	void Step() override;
	int SupplyNoOfYoung() { return m_NoOfYoung; };

	int SupplyMateId() { return m_MatesIdNo; };
	int SupplyMateSB() { return m_MateLive; };

	unsigned SupplyYoungAge() { return m_YoungAge; };
	bool SupplyPregnant() { return m_Pregnant; };

	void Set_MateState(int MateState) { m_MateLive = MateState; };
	void OnInfanticideAttempt();

protected:
	//Attributes

	int m_MateLive;
	int m_MatesIdNo;
	/**
	\brief
	The number of young in the current litter (if one).
	*/
	int m_NoOfYoung;
	/**
	\brief
	A flag indicating whether pregnant or not.
	*/
	bool m_Pregnant;
	/**
	\brief
	A counter counting down gestation days.
	*/
	int m_DaysUntilBirth;
	/**
	\brief
	The age of current litter in days.
	*/
	unsigned m_YoungAge;
	/**
	\brief
	The DNA passed from the male on mating.
	*/
	GeneticMaterial m_MatesGenes;

	// Methods
	TTypeOfVoleState st_ReproBehaviour() const;
	TTypeOfVoleState st_UpdateGestation();
	TTypeOfVoleState st_BecomeReproductive();
	TTypeOfVoleState st_GiveBirth();
	TTypeOfVoleState st_Lactating();
	int st_Evaluate_n_Explore();
	int st_Special_Explore();
	TTypeOfVoleState st_Mating();
#ifdef __VOLEPESTICIDEON
	unsigned int m_pesticideloadindex;
	double m_pesticideloadarray[21];
	virtual void ModelinkPesticide();
	virtual void ModelinkPesticide21TWA(double a_dose);
	virtual void Vinclozolin(double a_dose);
	virtual void GeneralOrganophosphate(double a_dose);
	virtual void GeneralEndocrineDisruptor(double /* a_dose */);
	static int m_EndoCrineDisruptionGestationLength;
	virtual void GeneticDemoPesticide(double /* a_dose */);
#endif
};

//---------------------------------------------------------------------------
#endif
