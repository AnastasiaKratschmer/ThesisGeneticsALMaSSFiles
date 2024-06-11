#ifndef Vole_toletoc_H
#define Vole_toletoc_H

#include "../Landscape/ls.h"

int vole_tole_move_quality(Landscape* m_TheLandscape, int x, int y);
int vole_tole_init_optimal(Landscape* m_TheLandscape, int x, int y);
int vole_tole_init_friendly(Landscape* m_TheLandscape, int x, int y);
double vole_toletoc_asses_habitat_score(Landscape* m_TheLandscape, int p_Polyref);
bool vole_tole_assess_barrier(Landscape* m_TheLandscape, int p_Polyref);
#endif