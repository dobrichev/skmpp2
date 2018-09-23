
#define _CRT_SECURE_NO_DEPRECATE
//#define ZHOU_OLD
#include <sys/timeb.h>
#include "main.h"  // tab0  

uint64_t p_cptg[40], p_cpt1g[20], p_cpt2g[20];


extern ZHOU    zhou[50], zhou_i, zhou_solve;
extern ZH_GLOBAL zh_g;
extern SGO sgo;

ofstream  fout1,fout2,fout3,fout4,
fout_diam, fout_pearl, fout_l45, fout_l65, fout_solved, fout_unsolved;

#include "solver_step.h"
FINPUT finput;
#include "go_0_cpp.h"
// updated file for test upload
