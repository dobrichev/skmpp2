
#define _CRT_SECURE_NO_DEPRECATE
#include "main.h"
#include "Zhn_cpp.h"

// catching time as seconds+millis  (seconds since year 1970)
int32_t GetTimeMillis() {
	struct _timeb tbuf;
	_ftime64_s(&tbuf); 
	return ((int32_t)(1000 * tbuf.time) + tbuf.millitm);
}
// builing an appropriate message depending on the elapsed time te-ts
void PrintTime(int32_t ts,int32_t te){
	uint32_t dt=te-ts,dtmil=dt%1000,dts=dt/1000,dth=dts/3600;   dth=dth%1000;
	cerr << endl<<"total elapsed time "; 
    uint32_t dtm=dts/60; dts=dts%60 ,   dth=dtm/60, dtm=dtm%60;
    if(dth) cerr <<dth<<"h "; 
	if(dth || dtm) cerr <<dtm<<"m "; 
	cerr	<<dts <<"s ";
	if(dtmil<10) cerr << "00"; else  if(dtmil<100) cerr << '0';
	cerr <<dtmil<<"ms "<<endl;   return;
}

void PrintTimeCout(int32_t ts, int32_t te){
	uint32_t dt = te - ts, dtmil = dt % 1000, dts = dt / 1000, dth = dts / 3600;   dth = dth % 1000;
	cout << endl << "total elapsed time ";
	uint32_t dtm = dts / 60; dts = dts % 60, dth = dtm / 60, dtm = dtm % 60;
	if (dth) cout << dth << "h ";
	if (dth || dtm) cout << dtm << "m ";
	cout << dts << "s ";
	if (dtmil<10) cout << "00"; else  if (dtmil<100) cout << '0';
	cout << dtmil << "ms " << endl;   return;
}


int Search_ccd(const char * ww)
{	// List of 2 char command, 9 commands
	const char * ccd[]={"-i" ,    // input name including extension
				  "-o" ,	//	output or second filename
				  "-c",  // main command option + check boxes
				  "-v" ,  // value   0 to 9 default 0
				  "-b" ,  // bit field 0 to 9 default is 0
				  "-s",  // strings 0 to 9
	};   // puzzle processing specificities
	char wt[4]; 
	strncpy_s(wt,4,ww,2);
	wt[2]=0;
	for(int i=0;i<6;i++)
		if(!strcmp(wt,ccd[i])) 
			return i;
	return -1;
}

SGO sgo;
extern void Go_0();
extern FINPUT finput;
int main(int narg, char *argv[]) {
	cerr << "mainstart" << endl;
	int32_t tdeb=GetTimeMillis();
	char * finput_name=0,*foutput_name=0,* ww;
	char * s_strings[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };//optionnal 10 strings

	uint32_t command = 0, 
		vx[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, //integers 0 to 9
		bfx[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };// bit fields 

	for(int i=0;i<narg;i++)	{
		ww = argv[i];
		int ir=Search_ccd(ww);
		if(ir<0) continue;
		if (ir == 3){// -vn-xxxx
			if (ww[3] - '-') continue; //must be -vn-  with n 0_9
			int ind=ww[2] - '0';
			if (ind < 0 || ind>9)continue;
			vx[ind] = atoi(&ww[4]);
			continue;
		}	
		else if (ir == 4){//  -bn- followed by a bit field bits rigth to left max 8 bits
			if (ww[3] - '-') continue; 
			int ind = ww[2] - '0';
			if (ind < 0 || ind>9)continue;
			int length = (int)strlen(&ww[4]);
			if (length > 8)continue; 
			for (int i = 0; i < length; i++) if (ww[4 + i] == '1') bfx[ind] |= (1 << i);
		}
		else if (ir == 5){//  -sn- followed by a string
			if (ww[3] - '-') continue; //must be -vn-  with n 0_9
			int ind = ww[2] - '0';
			if (ind < 0 || ind>9)continue;
			s_strings[ind] = &ww[4];
		}
		else{
			switch (ir)	{
			case 0: finput_name = &ww[2]; 	break;  // -i
			case 1: foutput_name = &ww[2]; 	break;  // -o
			case 2: command = atoi(&ww[2]);	break; //-c
			}// end command  
		}
	}// end loop on options
	if(finput_name) cerr <<" file1 (input) " << finput_name<<endl;
	if(foutput_name) cerr <<" file2 (output) " << foutput_name<<endl;
	cerr << "command " << command<<endl;
	// set zh_g general commands
	zh_g.modeguess =(int) bfx[9];
	if (!(zh_g.modeguess & 1))zh_g.modeguess = 0;
	zh_g.maxindex = (int)vx[7];
	// store command line parameters 
	sgo.command = command;
	sgo.bfx = bfx;
	sgo.finput_name = finput_name;
	sgo.foutput_name = foutput_name;
	sgo.s_strings = s_strings;
	sgo.vx = vx;
	
	Go_0();
	cerr << " print cout time "  << endl;

	int32_t tfin=GetTimeMillis();
    PrintTimeCout(tdeb,tfin);
	PrintTime(tdeb, tfin);
	return 0;
}

uint64_t p_cptg[40], p_cpt1g[20], p_cpt2g[20];


extern ZHOU    zhou[50];
extern ZH_GLOBAL zh_g;
extern SGO sgo;

ofstream  fout1,fout2;

#include "solver_step.h"

extern PM_GO pm_go;
void Go_c110(){// template serate mode
    if (!sgo.finput_name) return;
	long long cptg[20];
	memset(cptg, 0, sizeof cptg);
	cout << "Go_110 entry " << sgo.finput_name << " input" << endl;
	FINPUT finput;
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
	}
	switch (sgo.command / 100){
	case 1: Go_sol_1xx(); break;
	}
	cerr << "go_0 return" << endl;
}


