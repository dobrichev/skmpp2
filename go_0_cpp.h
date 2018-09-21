

//========================================
char * zh_g_cpt[10] = { "npuz", "guess", "close_d ", "upd1 ", "upd2 ",
"fupd ", "hpair ", "htripl ", " ", " " };
#include "go_0xx.cpp"
#include "go_1xx.cpp"
#include "go_2xx.cpp"
#include "go_4xx.cpp"
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
void Go_misc_4xx(){
	cout << "command 4xx command=" << sgo.command << endl;
	switch (sgo.command){
		/*
  	case 200: Go_c200(); break;// split the entry file in files 1;2;3
	case 201: Go_c201(); break;// change n clues or 1 to n clues
	case 203: Go_c201_old(); break;// change n clues or 1 to n clues
	case 202: Go_c202(); break;// gen symmetry of given
	case 210:  Go_c210(); break;// create a seed on a pattern
	case 211:  Go_c211(); break;// create a seed gen interim file
	case 212:  Go_c212(); break;// create a seed gen interim file
	case 217:  Go_c217(); break;// restore part of a data base
	case 218:  Go_c218(); break;// extract played seeds
	case 219:  Go_c219(); break;// restore a data base
      */
	case 440: Go_c440(); break;// parse game submissions  
	case 445: Go_c445(); break;// split the entry file on int param 

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
	/*
	// open  outputs files 1;2;3 output +_filex.txt
	zn[ll+5]='2';
	if(output2.OpenFO(zn))
	return 2; // 1 if error open	*/
	switch (sgo.command / 100){
	case 0: Go_0xx(); break;
	case 1: Go_sol_1xx(); break;
	case 2: Go_gen_2xx(); break;
	case 3: Go_can_3xx(); break;
	case 4: Go_misc_4xx(); break;
	}
	cerr << "go_0 return" << endl;
}
