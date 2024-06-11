#include "../Vole/Vole_toletoc.h"



int vole_tole_move_quality(Landscape* m_TheLandscape, int x, int y)
{
	int poly_index = m_TheLandscape->SupplyPolyRefIndex(x, y);
	int score = -9999;
	TTypesOfLandscapeElement tole = m_TheLandscape->SupplyElementTypeFromVector(poly_index);
	switch (tole)
	{
	case tole_NaturalGrassDry: // 110
	case tole_NaturalGrassWet: // 110
	case tole_UnknownGrass:
	case tole_GreenFallow:
	case tole_BeetleBank:
		score = 4;
		break;
	case tole_PermPasture: // 35
	case tole_OPermPasture:
	case tole_PermPastureLowYield: // 35
	case tole_OPermPastureLowYield:
	case tole_NaturalFarmGrass:
	case tole_PermPasturePigs:
	case tole_OPermPasturePigs:
	case tole_OtherPermCrop:
		if (m_TheLandscape->SupplyGrazingPressureVector(poly_index) > 0) score = 1; else score = 4;
		break;
	case tole_PermPastureTussocky:
	case tole_PermPastureTussockyWet:
		if (m_TheLandscape->SupplyGrazingPressureVector(poly_index) > 0) score = 2; else score = 4;
		break;
	case tole_Orchard: // 56
	case tole_OOrchard:
	case tole_BushFruit:
	case tole_OBushFruit:
	case tole_Vineyard:
	case tole_OliveGrove:
	case tole_WalnutPlantation:
	case tole_AlmondPlantation:
	case tole_StonePineForest:
		// Quality depends on when it was mown
		if (m_TheLandscape->SupplyJustMownVector(poly_index)) score = 0; else score = 4;
		break;
	case tole_MownGrassStrip: // 58
		// Quality depends on when it was mown
		if (m_TheLandscape->SupplyJustMownVector(poly_index)) score = 1; else score = 4;
		break;
	case tole_OrchardBand: // 57
		score = 4 - (int)(0.6 * m_TheLandscape->SupplyJustSprayedVector(poly_index));
		break;
	case tole_RoadsideSlope:
	case tole_Railway: // 118
	case tole_Pipeline:
	case tole_FieldBoundary: // 160
	case tole_WoodlandMargin: 
	case tole_FlowerStrip:
	case tole_RoadsideVerge: // 13
	case tole_HedgeBank:
	case tole_PermanentSetaside:
	case tole_YoungForest: // 60
	case tole_OFarmYoungForest:
	case tole_FarmYoungForest:
	case tole_Vildtager:
	case tole_Wasteland:
	case tole_ForestAisle:
	case tole_SolarPanel:
	case tole_FlowersPerm:
	case tole_WaterBufferZone:
	case tole_FarmBufferZone:
		score = 4;
		break;
	case tole_Saltmarsh: // 95
	case tole_Marsh: // 95
	case tole_Scrub: // 70
	case tole_Heath:
	case tole_Hedges: // 130 (internal ALMaSS representation for Hedges)
	case tole_RiversidePlants: // 98
	case tole_RiceField:
	case tole_AsparagusPerm:
	case tole_OAsparagusPerm: 
	case tole_PitDisused: // 75
	case tole_MixedForest: // 60
	case tole_DeciduousForest: // 40
	case tole_MontadoCorkOak:
	case tole_MontadoHolmOak:
	case tole_MontadoMixed:
	case tole_AgroForestrySystem:
	case tole_FarmForest:
	case tole_OFarmForest:
	case tole_CorkOakForest:
	case tole_HolmOakForest:
	case tole_OtherOakForest:
	case tole_ChestnutForest:
	case tole_EucalyptusForest:
	case tole_InvasiveForest:
	case tole_SwampForest:
	case tole_Copse:
		score = 3;
		break;
	case tole_RiversideTrees: // 97
	case tole_ConiferousForest: // 50
	case tole_MaritimePineForest:
	case tole_ChristmasTrees:
	case tole_OChristmasTrees:
	case tole_BuiltUpWithParkland:
	case tole_Parkland:
	case tole_AmenityGrass:
	case tole_WoodyEnergyCrop:
	case tole_OEnergyCrop:
	case tole_EnergyCrop:
	case tole_IndividualTree:
		score = 2;
		break;
	case tole_Garden: //11
	case tole_Track: // 123
	case tole_SmallRoad: // 122
	//case tole_LargeRoad: // 121
	case tole_BareRock:
	case tole_Saltpans:
	case tole_UrbanNoVeg:
	case tole_UrbanVeg:
	case tole_PlantNursery:
	case tole_UrbanPark:
	case tole_SandDune:
	case tole_Churchyard:	//204
	case tole_Airport:
	case tole_Portarea:
	case tole_Stream:
		score = 1;
		break;
	case tole_LargeRoad:
		double CanItCrossTheRoad;
		CanItCrossTheRoad= g_rand_uni_fnc();
		if ( CanItCrossTheRoad > 0.9){
			score=1;
		}else{
			score=-1;
		}
		break;
	case tole_MetalledPath:
	case tole_Carpark:
	case tole_HeritageSite:
	case tole_Building: // 5
	case tole_ActivePit: // 115
	case tole_RefuseSite:
	case tole_Freshwater: // 90
	case tole_FishFarm:
	case tole_Pond:
	case tole_River: // 96
	case tole_RiverBed:
	case tole_Saltwater: // 80
	case tole_Coast: // 100
	case tole_StoneWall: // 15
	case tole_Fence: //225
	case tole_WindTurbine:
	case tole_Pylon:
	case tole_DrainageDitch:
	case tole_Canal:
	case tole_FarmFeedingGround: 
	case tole_MushroomPerm:
	case tole_Chameleon:
	case tole_Missing:

		score = -1;
		break;
	case tole_Field: // 20 & 30
	case tole_UnsprayedFieldMargin:
	{
		TTypesOfCrops CType = m_TheLandscape->SupplyCropType(poly_index);
		switch (CType) {
		case toc_OPermanentGrassLowYield:
		case toc_PermanentGrassLowYield:
		case toc_PermanentGrassTussocky:
		case toc_OSetAside:
		case toc_SetAside:
		case toc_PermanentSetAside:
		case toc_YoungForestCrop:
		case toc_OYoungForestCrop: /** todo Should the orchards be here?*/
			score = 4;
			break;
		case toc_OSeedGrass1:
		case toc_SeedGrass1:
		case toc_OSeedGrass2:
		case toc_SeedGrass2:
		case toc_OCloverGrassSilage1:
		case toc_CloverGrassGrazed1:
		case toc_CloverGrassGrazed2:
		case toc_OCloverGrassGrazed1:
		case toc_OCloverGrassGrazed2:
		case toc_OPermanentGrassGrazed:
		case toc_PermanentGrassGrazed:
		case toc_GrassGrazed1:
		case toc_GrassGrazed2:
		case toc_GrassGrazedExtensive:
		case toc_GrassGrazedLast:
		case toc_CloverGrassGrazed3:
		case toc_OCloverGrassGrazed3:
			score = 3 - m_TheLandscape->SupplyGrazingPressureVector(poly_index);
			break;
		case toc_AsparagusEstablishedPlantation:
		case toc_Beans:
		case toc_Beans_Whole:
		case toc_Beet:
		case toc_BushFruit:
		case toc_Cabbage:
		case toc_CabbageSpring:
		case toc_Carrots:
		case toc_CarrotsSpring:
		case toc_CatchCropPea:
		case toc_CorkOak:
		case toc_DummyCropPestTesting:
		case toc_FarmForest:
		case toc_FieldPeas:
		case toc_FieldPeasSilage:
		case toc_FieldPeasStrigling:
		case toc_FodderBeet:
		case toc_FodderGrass:
		case toc_FodderLucerne1:
		case toc_FodderLucerne2:
		case toc_GenericCatchCrop:
		case toc_GrazingPigs:
		case toc_Horticulture:
		case toc_Maize:
		case toc_MaizeSilage:
		case toc_MaizeSpring:
		case toc_MaizeStrigling:
		case toc_MixedVeg:
		case toc_OAsparagusEstablishedPlantation:
		case toc_Oats:
		case toc_OBarleyPeaCloverGrass:
		case toc_OBeans:
		case toc_OBeans_Whole:
		case toc_OBushFruit:
		case toc_OCabbage:
		case toc_OCarrots:
		case toc_OFarmForest:
		case toc_OFieldPeas:
		case toc_OFieldPeasSilage:
		case toc_OFirstYearDanger:
		case toc_OFodderBeet:
		case toc_OFodderGrass:
		case toc_OGrazingPigs:
		case toc_OLentils:
		case toc_OliveGrove:
		case toc_OLupines:
		case toc_OMaize:
		case toc_OMaizeSilage:
		case toc_OMixedVeg:
		case toc_OOats:
		case toc_OOrchApple:
		case toc_OOrchardCrop:
		case toc_OOrchCherry:
		case toc_OOrchOther:
		case toc_OOrchPear:
		case toc_OPotatoes:
		case toc_OPotatoesIndustry:
		case toc_OPotatoesSeed:
		case toc_OrchApple:
		case toc_OrchardCrop:
		case toc_OrchCherry:
		case toc_OrchOther:
		case toc_OrchPear:
		case toc_ORyeGrass:
		case toc_OSBarleySilage:
		case toc_OSetAside_Flower:
		case toc_OSpringBarley:
		case toc_OSpringBarleyCloverGrass:
		case toc_OSpringBarleyExtensive:
		case toc_OSpringBarleyPeaCloverGrass:
		case toc_OSpringBarleyPigs:
		case toc_OSpringBarleySilage:
		case toc_OSpringRape:
		case toc_OSpringRye:
		case toc_OSpringWheat:
		case toc_OStarchPotato:
		case toc_OSugarBeet:
		case toc_OTriticale:
		case toc_OVegSeeds:
		case toc_OWinterBarley:
		case toc_OWinterBarleyExtensive:
		case toc_OWinterRape:
		case toc_OWinterRye:
		case toc_OWinterWheat:
		case toc_OWinterWheatUndersown:
		case toc_OWinterWheatUndersownExtensive:
		case toc_PlantNursery:
		case toc_Potatoes:
		case toc_PotatoesIndustry:
		case toc_PotatoesSeed:
		case toc_PotatoesSpring:
		case toc_Ryegrass:
		case toc_ORyegrass:
		case toc_Sorghum:
		case toc_SpringBarley:
		case toc_SpringBarleyCloverGrass:
		case toc_SpringBarleyPeaCloverGrass:
		case toc_SpringBarleySeed:
		case toc_SpringBarleySilage:
		case toc_SpringRape:
		case toc_SpringRye:
		case toc_SpringWheat:
		case toc_StarchPotato:
		case toc_SugarBeet:
		case toc_Sunflower:
		case toc_Triticale:
		case toc_Tulips:
		case toc_Turnip:
		case toc_VegSeeds:
		case toc_Vineyards:
		case toc_WinterBarley:
		case toc_WinterRape:
		case toc_WinterRye:
		case toc_WinterTriticale:
		case toc_WinterWheat:
		case toc_YellowLupin:
		case toc_Fallow:
		case toc_Unmanaged:
		{
			double cover = m_TheLandscape->SupplyVegCoverVector(poly_index);
			double height = m_TheLandscape->SupplyVegHeightVector(poly_index);
			if ((cover > 0.80) && (height > 40)) score = 3;
			else if ((cover < 0.50) || (height < 10)) score = 0;
			else score = 2;
		}
			break;
		
		case toc_Foobar:
		case toc_None:
		default: 
			//   Unknown crop type
			g_msg->Warn("vole_tole_move_quality - Unknown toc type: ", static_cast<int>(CType));
			exit(1);
			break;
		}
	}
	break;
	case tole_Foobar: // 999 !! type unknown - should not happen
	default:
		static char errornum[20];
		sprintf(errornum, "%d", m_TheLandscape->SupplyElementTypeFromVector(poly_index));
		m_TheLandscape->Warn("Vole_Base:MoveQuality: Unknown tole_type",
			errornum);
		exit(1);

	}
	return score; 
}

int vole_tole_init_optimal(Landscape* m_TheLandscape, int x, int y)
{
	TTypesOfLandscapeElement tole = m_TheLandscape->SupplyElementType(x, y);
	switch (tole)
	{
	case tole_FieldBoundary:
	case tole_PermPastureLowYield:
	case tole_OPermPastureLowYield:
	case tole_PermPastureTussocky:
	case tole_PermPastureTussockyWet:
	case tole_PermanentSetaside:
	case tole_NaturalGrassDry:
	case tole_PermPasture:
	case tole_OPermPasture:
	case tole_YoungForest:
	case tole_FarmYoungForest:
	case tole_OFarmYoungForest:
	case tole_OrchardBand:
	case tole_FlowerStrip:
	case tole_FlowersPerm:
	case tole_NaturalFarmGrass:
	case tole_NaturalGrassWet:
	//case tole_RoadsideVerge: /*\todo some were already commented out */
	//case tole_Railway:
	//case tole_Scrub:
	//case tole_PitDisused:
	//case tole_OFarmForest://Jordans
	//case tole_FarmForest://Jordans
	//case tole_Garden:
	//case tole_HedgeBank:
	case tole_BeetleBank:
		return true;
	
	default:
		return false;
	}
}

int vole_tole_init_friendly(Landscape* m_TheLandscape, int x, int y)
{
	
	TTypesOfLandscapeElement tole = m_TheLandscape->SupplyElementType(x, y);
	
	switch (tole)
	{
	case tole_Orchard:
	case tole_Vineyard:
	case tole_OliveGrove:
	case tole_OrchardBand:
	case tole_MownGrassStrip:
	case tole_RoadsideVerge:
	case tole_RoadsideSlope:
	case tole_Railway:
	case tole_FlowerStrip: //commented out because they are also in init_optimal - can/should they be in both places? - Both
	case tole_FieldBoundary:
	case tole_Scrub:
	case tole_Field:
	case tole_PermPastureLowYield:
	case tole_OPermPastureLowYield:
	case tole_PermPastureTussocky:
	case tole_PermanentSetaside:
	case tole_PermPasture:
	case tole_NaturalGrassDry:
	case tole_NaturalGrassWet:
	case tole_PitDisused:
	case tole_YoungForest:
	case tole_FarmYoungForest:
	case tole_Garden:
	case tole_HedgeBank:
	case tole_BeetleBank:
	case tole_SolarPanel:
	case tole_BushFruit: //Jordan
	case tole_Churchyard: //Jordan
	case tole_Copse: //Jordan
	case tole_ForestAisle: //Jordan
	case tole_GreenFallow: //Jordan
	case tole_Heath: //Jordan
	case tole_Hedges: //Jordan
	case tole_IndividualTree: //Jordan
	case tole_OBushFruit: //Jordan
	case tole_OChristmasTrees: //Jordan
	case tole_ChristmasTrees:
	case tole_OEnergyCrop: //Jordan
	case tole_OOrchard: //Jordan
	case tole_OPermPasture: //Jordan
	case tole_OtherPermCrop: //Jordan
	case tole_Parkland: //Jordan
	case tole_PermPastureTussockyWet: //Jordan
	case tole_RiversidePlants: //Jordan
	case tole_RiversideTrees: //Jordan
	case tole_UnknownGrass: //Jordan
	case tole_UnsprayedFieldMargin: //Jordan
	case tole_Wasteland: //Jordan
	case tole_WaterBufferZone: //Jordan
	case tole_WoodlandMargin: //Jordan
	case tole_Saltmarsh: // 95
	case tole_Marsh: // 95
	case tole_Vildtager:
	case tole_FarmBufferZone:
		return true;
	default:
		return false;
	}

}


double vole_toletoc_asses_habitat_score(Landscape* m_TheLandscape, int p_Polyref)
{
	TTypesOfLandscapeElement ElementType = m_TheLandscape->SupplyElementTypeFromVector(p_Polyref);
	double score = -9999;
	double Cover;

	switch (ElementType)
	{
	case tole_Railway: // 118
	case tole_FlowerStrip:
	case tole_FlowersPerm: //Jordan
	case tole_FieldBoundary: // 160
	case tole_RoadsideVerge: // 13
	case tole_RoadsideSlope:
	case tole_HedgeBank:
	case tole_Hedges: // 130 (internal ALMaSS representation for Hedges)
	case tole_BeetleBank:
	case tole_NaturalGrassWet:
	case tole_NaturalFarmGrass: //Jordan
	case tole_SolarPanel: //Jordan
	case tole_WaterBufferZone: //Jordan
	case tole_NaturalGrassDry: // 110
	case tole_FarmBufferZone:
		if (m_TheLandscape->SupplyVegHeightVector(p_Polyref) <= 5.0) score = 2.25;
		else score = 3.0;
		break;
	case tole_OtherPermCrop:
	case tole_PermPasture: // 35
	case tole_OPermPasture:
	case tole_PermPasturePigs:
	case tole_OPermPasturePigs:
		if (m_TheLandscape->SupplyGrazingPressureVector(p_Polyref) > 0) score = 2.0; else score = 2.7;
		break;
	case tole_PermPastureLowYield:
	case tole_OPermPastureLowYield:
		if (m_TheLandscape->SupplyGrazingPressureVector(p_Polyref) > 0) score = 2.1; else score = 2.8;
		break;
	case tole_PermPastureTussocky:
	case tole_PermPastureTussockyWet: //Jordan
		if (m_TheLandscape->SupplyGrazingPressureVector(p_Polyref) > 0) score = 2.25; else score = 3.0;
		break;
	case tole_Orchard: // 56
	case tole_OOrchard: //Jordan
	case tole_BushFruit:
	case tole_OBushFruit:
	case tole_Vineyard:
	case tole_OliveGrove:
	case tole_AlmondPlantation:
	case tole_WalnutPlantation:

		// Quality depends on when it was mown
	{
		if (m_TheLandscape->SupplyVegCoverVector(p_Polyref) < 0.8) score = 2.5;
		else if (m_TheLandscape->SupplyVegHeightVector(p_Polyref) <= 40) score = 2.5;
		else score = 3.0;
	}
		if (m_TheLandscape->SupplyJustMownVector(p_Polyref)) score = 2.0;
		break;
	case tole_MownGrassStrip: // 58
		// Quality depends on when it was mown
		if (m_TheLandscape->SupplyJustMownVector(p_Polyref)) score = 1.0; else score = 1.5;
		break;
	case tole_OrchardBand: // 57
		// If spraying herbicide then quality depends on time since spraying. SupplyJustSprayed returns a counter that counts down to zero after spraying.
		// score = 3.0 - (0.5*m_TheLandscape->SupplyJustSprayedVector(p_Polyref)); // For herbicide
		score = 3.0;
		break;
	case tole_PermanentSetaside:
	case tole_YoungForest: // 60
	case tole_FarmYoungForest://Jordans
	case tole_OFarmYoungForest:
	case tole_ForestAisle:
		score = 2.2;
		break;
	case tole_Marsh: // 95
	case tole_Scrub: // 70
	case tole_Vildtager:
	case tole_Heath:
	case tole_RiversidePlants: // 98
	case tole_RiceField:
	case tole_PitDisused: // 75
	case tole_UnknownGrass:
	case tole_GreenFallow:
	case tole_Wasteland:
	case tole_AsparagusPerm:
	case tole_OAsparagusPerm:
		score = 1.5;
		break;
	case tole_MixedForest: // 60
	case tole_DeciduousForest: // 40
	case tole_MontadoCorkOak:
	case tole_MontadoHolmOak:
	case tole_MontadoMixed:
	case tole_AgroForestrySystem:
	case tole_OFarmForest:
	case tole_FarmForest:
	case tole_CorkOakForest:
	case tole_HolmOakForest:
	case tole_OtherOakForest:
	case tole_ChestnutForest:
	case tole_EucalyptusForest:
	case tole_InvasiveForest:
	case tole_SwampForest:
	case tole_RiversideTrees: // 97
	case tole_ConiferousForest: // 50
	case tole_ChristmasTrees:
	case tole_OChristmasTrees:
	case tole_MaritimePineForest:
	case tole_StonePineForest:
	case tole_BuiltUpWithParkland:
	case tole_Parkland:
	case tole_Copse:
	case tole_AmenityGrass:
	case tole_MetalledPath:  //202
	case tole_WoodyEnergyCrop:// 59
	case tole_EnergyCrop:
	case tole_OEnergyCrop:
	case tole_WoodlandMargin:
	case tole_IndividualTree:
		score = 1;
		break;
	case tole_Garden: //11
	case tole_Track: // 123
	case tole_SmallRoad: // 122
	case tole_LargeRoad: // 121
	case tole_BareRock:
	case tole_UrbanNoVeg:
	case tole_UrbanPark:
	case tole_UrbanVeg:
	case tole_SandDune:
	case tole_Churchyard:
	case tole_HeritageSite:
	case tole_Saltmarsh:
	case tole_PlantNursery:
	case tole_WindTurbine:
	case tole_Pylon:
	case tole_Saltpans:
	case tole_Pipeline:
		score = 0;
		break;
	case tole_Building: // 5
	case tole_Carpark:
	case tole_ActivePit: // 115
	case tole_Freshwater: // 90
	case tole_Pond:
	case tole_Stream:
	case tole_Saltwater: // 80
	case tole_Coast: // 100
	case tole_StoneWall: // 15
	case tole_Fence: //225
	case tole_DrainageDitch:
	case tole_RefuseSite:
	case tole_Canal:
	case tole_River:
	case tole_Airport:
	case tole_Portarea:
	case tole_MushroomPerm:
	case tole_Missing:
	case tole_FishFarm:
	case tole_RiverBed:
	case tole_Chameleon:
	case tole_FarmFeedingGround:
		score = -1;
		break;
	case tole_Field: // 20 & 30
	case tole_UnsprayedFieldMargin:
	{
		TTypesOfCrops CType = m_TheLandscape->SupplyCropType(p_Polyref);
		Cover = m_TheLandscape->SupplyVegCoverVector(p_Polyref);
		if (Cover < 0.20) score = 0;
		else switch (CType)
		{
		case toc_SeedGrass1:
		case toc_SeedGrass2:
		case toc_OSeedGrass1:
		case toc_OSeedGrass2:
		case toc_ORyeGrass:
		case toc_ORyegrass:
			score = 2.25;
			break;
		case toc_Beans_Whole:
		case toc_OBeans:
		case toc_Beans:
		case toc_OBeans_Whole:
		case toc_OAsparagusEstablishedPlantation:
		case toc_AsparagusEstablishedPlantation:
		case toc_CatchCropPea:
		case toc_FieldPeas:
		case toc_FieldPeasStrigling:
		case toc_OLentils:
		case toc_FieldPeasSilage:
		case toc_WinterBarley:
		case toc_OWinterBarley:
		case toc_WinterRye:
		case toc_OWinterRye:
		case toc_OSpringRye:
		case toc_Ryegrass:
		case toc_SpringRye:
		case toc_OSpringBarleyCloverGrass:
		case toc_OSpringBarleyExtensive:
		case toc_OSpringBarleyPeaCloverGrass:
		case toc_OSpringBarleySilage:
		case toc_OSpringRape:
		case toc_Sorghum:
		case toc_OSpringWheat:
		case toc_OStarchPotato:
		case toc_OSugarBeet:
		case toc_OVegSeeds:
		case toc_OWinterBarleyExtensive:
		case toc_OWinterWheatUndersownExtensive:
		case toc_Oats:
		case toc_WinterWheat:
		case toc_Triticale:
		case toc_WinterTriticale:
		case toc_OTriticale:
		case toc_WinterRape:
		case toc_SpringRape:
		case toc_OSpringBarley:
		case toc_SpringBarley:
		case toc_SpringBarleySeed:
		case toc_SpringBarleyCloverGrass:
		case toc_SpringBarleyPeaCloverGrass:
		case toc_SpringBarleySilage:
		case toc_SpringWheat:
		case toc_OBarleyPeaCloverGrass:
		case toc_OSBarleySilage:
		case toc_OWinterWheatUndersown:
		case toc_OWinterWheat:
		case toc_OFieldPeas:
		case toc_OFieldPeasSilage:
		case toc_OSpringBarleyPigs:
		case toc_OWinterRape:
		case toc_Maize:
		case toc_OMaize:
		case toc_MaizeSpring:
		case toc_MaizeStrigling:
		case toc_MaizeSilage:
		case toc_OMaizeSilage:
		case toc_OOats:
		case toc_GenericCatchCrop:
			score = 1;
			break;
		case toc_OGrazingPigs:
		case toc_Beet:
		case toc_FodderBeet:
		case toc_SugarBeet:
		case toc_OFodderBeet:
		case toc_Carrots:
		case toc_OCarrots:
		case toc_CarrotsSpring:
		case toc_Potatoes:
		case toc_OPotatoes:
		case toc_PotatoesIndustry:
		case toc_OPotatoesIndustry:
		case toc_OPotatoesSeed:
		case toc_PotatoesSeed:
		case toc_VegSeeds:
		case toc_PotatoesSpring:
		case toc_StarchPotato:
		case toc_OMixedVeg:
		case toc_MixedVeg:
		case toc_Turnip:
		case toc_OCabbage:
		case toc_Cabbage:
		case toc_CabbageSpring:
		case toc_Sunflower:
		case toc_Tulips:
		case toc_OLupines:
		case toc_YellowLupin:
		case toc_PlantNursery:
		case toc_Horticulture:
		case toc_DummyCropPestTesting:
		case toc_OFirstYearDanger:
		case toc_Fallow:
			score = 0;
			break;
			//  Needs some more thought
		case toc_OCloverGrassSilage1:
		case toc_FodderGrass:
		case toc_CloverGrassGrazed1:
		case toc_CloverGrassGrazed2:
		case toc_OCloverGrassGrazed1:
		case toc_OCloverGrassGrazed2:
		case toc_OCloverGrassGrazed3:
		case toc_CloverGrassGrazed3:
		case toc_OFodderGrass:
		case toc_FodderLucerne1:
		case toc_FodderLucerne2:
		case toc_GrassGrazed1:
		case toc_GrassGrazed2:
		case toc_GrassGrazedExtensive:
		case toc_GrassGrazedLast:
		case toc_GrazingPigs:
			if (m_TheLandscape->SupplyGrazingPressureVector(p_Polyref) > 0) score = 1;
			else score = 2;
			break;
		case toc_OPermanentGrassGrazed:
		case toc_PermanentGrassGrazed:
		case toc_PermanentGrassLowYield:
		case toc_OPermanentGrassLowYield:
			if (m_TheLandscape->SupplyGrazingPressureVector(p_Polyref) > 0) score = 2.2;
			else score = 2.6;
			break;
		case toc_PermanentGrassTussocky:
			if (m_TheLandscape->SupplyGrazingPressureVector(p_Polyref) > 0) score = 2.25;
			else score = 3.0;
			break;
		case toc_PermanentSetAside:
		case toc_OSetAside_Flower:
		case toc_OSetAside:
		case toc_SetAside:
		case toc_Unmanaged:
			score = 3.0;
			break;
		case toc_OBushFruit:
		case toc_BushFruit:
		case toc_OrchardCrop:
		case toc_CorkOak:
		case toc_OFarmForest:
		case toc_FarmForest:
		case toc_OliveGrove:
		case toc_OOrchApple:
		case toc_OrchApple:
		case toc_OOrchardCrop:
		case toc_OOrchCherry:
		case toc_OrchCherry:
		case toc_OOrchOther:
		case toc_OrchOther:
		case toc_OOrchPear:
		case toc_OrchPear:
		case toc_Vineyards:
		case toc_YoungForestCrop:
		case toc_OYoungForestCrop:
			// Quality depends on when it was mown
			if (m_TheLandscape->SupplyVegCoverVector(p_Polyref) < 0.8) score = 2.5;
			else if (m_TheLandscape->SupplyVegHeightVector(p_Polyref) <= 40) score = 2.5;
			else score = 3.0;
			break;
		case toc_Foobar:
		case toc_None:
		default: 
			//   Unknown crop type
			g_msg->Warn("vole_toletoc_asses_habitat_score - Unknown toc type: ", static_cast<int>(CType) );
			exit(1);
			break;
		}
	}
		break;
	
	case tole_Foobar: // 999 !! type unknown - should not happen
	default:  // NOLINT(clang-diagnostic-covered-switch-default)
		static char errornum[20];
		sprintf(errornum, "%d", ElementType);
		m_TheLandscape->Warn("vole_toletoc_asses_habitat_score - Unknown tole_type: ", errornum);
		exit(1);
	}

	return score;

}

bool vole_tole_assess_barrier(Landscape* m_TheLandscape, int p_Polyref)
{
	TTypesOfLandscapeElement Elem;
	Elem = m_TheLandscape->SupplyElementType(p_Polyref);

	switch (Elem) {
// Roads, stream, garden, churchyards etc. and all types of forest
	case tole_Churchyard:
	case tole_Saltmarsh:
	case tole_Stream:
	case tole_HeritageSite:
	//case tole_LargeRoad: //moved to true
	case tole_Airport:
	case tole_Portarea:
	case tole_Coast:
	case tole_AgroForestrySystem:
	case tole_OFarmForest:
	case tole_FarmForest:
	case tole_DeciduousForest:
	case tole_ConiferousForest:
	case tole_OChristmasTrees:
	case tole_MaritimePineForest:	
	case tole_ChristmasTrees: //Jordan
	case tole_MontadoCorkOak:
	case tole_MontadoHolmOak:
	case tole_MontadoMixed:
	case tole_CorkOakForest:
	case tole_HolmOakForest:
	case tole_OtherOakForest:
	case tole_ChestnutForest:	
	case tole_MixedForest: //Jordan
	case tole_InvasiveForest:
	case tole_SwampForest:	
	case tole_Railway: //Jordan
	//case tole_SmallRoad: //Moved to true
	case tole_UrbanNoVeg: //Jordan
	case tole_UrbanVeg:
	case tole_UrbanPark:
	case tole_Garden: //11
	case tole_Track: // 123	
	case tole_AmenityGrass:
	case tole_Building:
	case tole_Canal:
	case tole_Copse:
	case tole_EucalyptusForest:
	case tole_FishFarm:
	case tole_Freshwater:
	case tole_Pond:
	case tole_River:
	case tole_Saltwater:
	case tole_StonePineForest:
	case tole_WindTurbine:
	case tole_RefuseSite:
		return false;

	case tole_ActivePit:
	case tole_AlmondPlantation:
	case tole_AsparagusPerm:
	case tole_BareRock:
	case tole_BeetleBank:
	case tole_BuiltUpWithParkland:
	case tole_BushFruit:
	case tole_Carpark:
	case tole_Chameleon:
	case tole_DrainageDitch:
	case tole_EnergyCrop:
	case tole_FarmBufferZone:
	case tole_FarmFeedingGround:
	case tole_FarmYoungForest:
	case tole_Fence:
	case tole_Field:
	case tole_FieldBoundary:
	case tole_FlowersPerm:
	case tole_FlowerStrip:
	case tole_Foobar:
	case tole_ForestAisle:
	case tole_GreenFallow:
	case tole_Heath:
	case tole_HedgeBank:
	case tole_Hedges:
	case tole_IndividualTree:
	case tole_Marsh:
	case tole_MetalledPath:
	case tole_Missing:
	case tole_MownGrassStrip:
	case tole_MushroomPerm:
	case tole_NaturalFarmGrass:
	case tole_NaturalGrassDry:
	case tole_NaturalGrassWet:
	case tole_OAsparagusPerm:
	case tole_OBushFruit:
	case tole_OEnergyCrop:
	case tole_OFarmYoungForest:
	case tole_OliveGrove:
	case tole_OOrchard:
	case tole_OPermPasture:
	case tole_OPermPastureLowYield:
	case tole_OPermPasturePigs:
	case tole_Orchard:
	case tole_OrchardBand:
	case tole_OtherPermCrop:
	case tole_Parkland:
	case tole_PermanentSetaside:
	case tole_PermPasture:
	case tole_PermPastureLowYield:
	case tole_PermPasturePigs:
	case tole_PermPastureTussocky:
	case tole_PermPastureTussockyWet:
	case tole_Pipeline:
	case tole_PitDisused:
	case tole_PlantNursery:
	case tole_Pylon:
	case tole_RiceField:
	case tole_RiverBed:
	case tole_RiversidePlants:
	case tole_RiversideTrees:
	case tole_RoadsideSlope:
	case tole_RoadsideVerge:
	case tole_Saltpans:
	case tole_SandDune:
	case tole_Scrub:
	case tole_SolarPanel:
	case tole_StoneWall:
	case tole_UnknownGrass:
	case tole_UnsprayedFieldMargin:
	case tole_Vildtager:
	case tole_Vineyard:
	case tole_WalnutPlantation:
	case tole_Wasteland:
	case tole_WaterBufferZone:
	case tole_WoodlandMargin:
	case tole_WoodyEnergyCrop:
	case tole_YoungForest:
	case tole_SmallRoad:
	case tole_LargeRoad:

		return true;

	default:
		//   Unknown crop type
		g_msg->Warn("bool vole_tole_assess_barrier - Unknown tole type: ", static_cast<int>(Elem));
		exit(1);


	}

	

}	