
extern PM_GO pm_go;
void Go_c110(){// template serate mode
    if (!sgo.finput_name) return;
	long long cptg[20];
	memset(cptg, 0, sizeof cptg);
	cout << "Go_110 entry " << sgo.finput_name << " input" << endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()){
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	char ze[200]; ze[81] = 0;
	uint32_t npuz = 0;
	while (finput.GetPuzzle(ze)){
		npuz++;
		if (npuz < sgo.vx[0]) continue;
		if(pm_go.opprint2)cout << finput.ze << "to process npuz=" << npuz << endl;
		zh_g.npuz = npuz;
		if (zh_g.Go_InitSolve(ze)) {
			cout << finput.ze << "invalid or multiple solutions npuz=" << npuz << endl;
			continue;
		}
		pm_go.SolveSerate110();
		if (sgo.vx[1] && npuz >= sgo.vx[1]) break;
	}
}

//========================================
const char * zh_g_cpt[10] = { "npuz", "guess", "close_d ", "upd1 ", "upd2 ", "fupd ", "hpair ", "htripl ", " ", " " };
#include "go_0xx_cpp.h"
//#include "go_1xx_cpp.h"
void Go_sol_1xx(){
	cout << "command 1xx command=" << sgo.command << endl;
	pm_go.opprint =sgo.bfx[9];
	pm_go.opprint2 = sgo.bfx[8];
	switch (sgo.command){
	case 110: Go_c110(); break;// template solve serate mode
	//case 111: Go_c111(); break;// solve/split quick serate mode
	//case 199: Go_c199(); break;// current test
	}
	cerr << "back go-sol_1xx" << endl;

}

void Go_0( ){
	// open  outputs files 1;2;3 output +_filex.txt
	if (sgo.foutput_name){
		char zn[200];
		strcpy(zn, sgo.foutput_name);
		int ll = (int)strlen(zn);
		strcpy(&zn[ll], "_file1.txt");
		fout1.open(zn);  
		zn[ll + 5] = '2'; fout2.open(zn);
		zn[ll + 5] = '3'; fout3.open(zn);
	}
	switch (sgo.command / 100){
	case 1: Go_sol_1xx(); break;
	}
	cerr << "go_0 return" << endl;
}
