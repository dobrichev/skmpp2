

//========================================
const char * zh_g_cpt[10] = { "npuz", "guess", "close_d ", "upd1 ", "upd2 ",
"fupd ", "hpair ", "htripl ", " ", " " };
#include "go_0xx_cpp.h"
#include "go_1xx_cpp.h"
#include "go_2xx_cpp.h"
#include "go_4xx_cpp.h"
 void Go_0xx(){
		Go_c0();
}
void Go_sol_1xx(){
	cout << "command 1xx command=" << sgo.command << endl;
	pm_go.opprint =sgo.bfx[9];
	pm_go.opprint2 = sgo.bfx[8];
	switch (sgo.command){
	case 110: Go_c110(); break;// template solve serate mode
	case 111: Go_c111(); break;// solve/split quick serate mode
	case 199: Go_c199(); break;// current test 
	}
	cerr << "back go-sol_1xx" << endl;

}
void Go_gen_2xx(){
	cout << "command 2xx command=" << sgo.command << endl;
	switch (sgo.command){
	case 200: Go_c200(); break;// split the entry file in files 1;2;3
	case 201: Go_c201(); break;// change n clues or 1 to n clues
	case 202: Go_c202(); break;// gen symmetry of given
	case 210:  Go_c210(); break;// create a seed on a pattern
	case 211:  Go_c211(); break;// create a seed gen interim file
	case 212:  Go_c212(); break;// create a seed gen interim file
	case 217:  Go_c217(); break;// restore a data base (v1=n v2=guesses)
	case 218:  Go_c218(); break;// extract played seeds
	case 219:  Go_c219(); break;// restore a data base
	}
}
void Go_can_3xx(){

}
/* subtask v0 for task 400
0 add sequence
1 add string0
2 add nclues
9 add stcd puz to entry
10 '.' for empty cell
11 erase '"' in entry
12 cut entry to v1
15 mantext in output
16 mintext in output
21 extract 81 character starting in v1
22 first v1 puzzles in output 1 others output2
23 sampling start v1 one every v2
40 count digits
41 count given per band
42 count digits per band
43 count given per unit
44count digits per unit
*/

void Go_misc_4xx(){
	cout << "command 4xx command=" << sgo.command << endl;
	switch (sgo.command){
	case 400: Go_c400(); break;// small tasks  see subtask v0
	case 401: Go_c401(); break;// .dat to .txt  
	case 402: Go_c402(); break;// morph rows cols diag  s1 s2 v1
	case 440: Go_c440(); break;// parse game submissions  
	case 445: Go_c445(); break;// split the entry file on int param 

	case 480: Go_c480(); break;// add compressed clues to entry (game data base)
	case 481: Go_c480(); break;// check/update game data base
	}

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
	case 0: Go_0xx(); break;
	case 1: Go_sol_1xx(); break;
	case 2: Go_gen_2xx(); break;
	case 3: Go_can_3xx(); break;
	case 4: Go_misc_4xx(); break;
	}
	cerr << "go_0 return" << endl;
}
