#include "sk_t.h"
#include "solver_step.h"
#include "Zhn.h"
//#include "Zhtables.cpp" // also describes the pattern

extern PM_GO pm_go;
extern ZHOU zhou_i;

/*
LastCell=10,				///< last cell in row column box
SingleBox=12,               ///< single in box
Single_R_C=15,              ///< single in row or column
Single_after_Locked=17,		///< locked in box  clearing row/col ?? giving a fix??
HiddenPair_single=20,		///< hidden pair, hidden fix
NakedSingle=23,				///< cell one candidate
HiddenTriplet_single=25,    ///< Hidden triplet, fix
Locked_box=26,				///< locked in box, no fix
Locked_RC=28,				///< locked in row/col  no fix
NakedPair=30,               ///< 2 cells containing 2 digits
XWing=32,                   ///< XWing
HiddenPair=34,              ///< 2 digits locked in 2 cell
Naked_triplet=36,           ///< 3 cells containing 3 digits
swordfish=38,               ///< swordfish
HiddenTriplet=40,           ///< 3 digits in 3 cells
XYWing=42,					///< XYWing
XYZWing=44,					///< XYZWing
*/

void ZH_GLOBAL::Debug(){
	cout << "zh_g debug digits maps" << endl;
	for (int i = 0; i < 9; i++) cout <<zhou_i.FD[i][1].bf.u32[3]+1;
	cout << "\t";
	for (int i = 0; i < 9; i++) cout << x3_dmap_inv[i] + 1;
	cout << endl << "cells maps" << endl;
	for (int i = 0; i < 81; i++)cout << x3_cmap[i] + 1;
	cout << endl << "list of given" << endl;
	for (int i = 0; i < ngiven; i++)
		cout << (int)(tgiven[i].u8[1] + 1) << "_" << (int)(tgiven[i].u8[0] + 1) << "  ";
	cout << endl;
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

void ZHOU::Debug(int all){
	//	if(0) return;
	cout << "DEBUG index=" << (int)index << " nbsol=" << zh_g.nsol << " unsolved=" << cells_unsolved.Count()
		<< " ndigits=" << ndigits << " unsolved digits 0" 
        <<oct<< unsolved_digits <<dec<< endl;
	char zi[82];  SetKnown(zi);
	cout <<zi<<" known rows digits"<<endl;
	if(!all) return;

	cout<<"map per digit"<<endl;
	for(int ib=0;ib<3;ib++) {
		for(int ir=0;ir<3;ir++){
			for(int idig=0;idig<9;idig++){
				uint32_t vf = FD[idig][0].bf.u32[ib],
					wr = (vf >> (9 * ir)) & 0x1ff, digit = FD[idig][1].bf.u32[3] + 1;
				for(int k=0;k<9;k++){
					if(wr & (1<<k))		cout <<digit;
					else 		cout<<".";
					if(k==2 || k==5) 	cout << " ";
				}
				cout<< "  ";
			}
			cout << endl; //end of row
		}
		cout << endl; // end of block
	}
	cout << endl; // end of map per digit

}
int ZHOU::GetAllDigits(int cell){
	int ir = 0, xcell = C_To128[cell];;
	for (int i = 0; i < 9; i++) if (FD[i][0].On(xcell))ir |= (1 << FD[i][1].bf.u32[3]);
	return ir;
}
void ZHOU::ImageCandidats() {// only active digits ??
	int dig_cells[81];for(int i=0;i<81;i++) dig_cells[i]=GetAllDigits(i);
	int i, j, l, lcol[9], tcol = 0,ncand=0;
	cout <<"PM map "<<endl<<endl;
	for(i = 0; i < 9; i++) {  // attention ici i indice colonne
		lcol[i] = 2;    // 2  mini tous chiffres impos�s
		for(j = 0; j < 9; j++) {
			l = _popcnt32(dig_cells[9 * j + i]);
			if(l > lcol[i])       lcol[i] = l;
		}
		tcol += lcol[i];
	}
	for(i = 0; i < 9; i++) {
		if((i == 3) || (i == 6))cout <<"|";
		cout <<(char)('A' + i)<<Blancs(lcol[i], 1);
	}
	cout << endl;
	for(i = 0; i < 9; i++) { // maintenant indice ligne
		if((i == 3) || (i == 6)) {
			for(int ix = 0; ix < (tcol + 10); ix++)       cout <<(char)'-';
			cout << endl;
		}
		for(j = 0; j < 9; j++) {
			if ((j == 3) || (j == 6))cout << "|";
			int cell = 9 * i + j, digs = dig_cells[cell], ndigs = _popcnt32(digs);
			ncand += ndigs;
			for (int id = 0; id < 9; id++)if (digs & (1 << id))
				cout << id + 1;
			cout <<Blancs(lcol[j] + 1 - ndigs, 1);
		} // end for j
		cout << endl;
	} // end for i
	cout << endl;

}

//void ZHOU::ImageCandidats_b3() {// only active digits ??
//	int dig_cells[81];for(int i=54;i<81;i++) dig_cells[i]=GetAllDigits(i);
//	int i, j, l, lcol[9], tcol = 0,ncand=0;
//	cout <<"PM map bande 3"<<endl<<endl;
//	for(i = 0; i < 9; i++) {  // attention ici i indice colonne
//		lcol[i] = 2;    // 2  mini tous chiffres impos�s
//		for(j = 6; j < 9; j++) {
//			l = _popcnt32(dig_cells[9 * j + i]);
//			if(l > lcol[i])       lcol[i] = l;
//		}
//		tcol += lcol[i];
//	}
//	for (i = 0; i < 9; i++) {
//		if ((i == 3) || (i == 6))cout << "|";
//		cout << (char)('A' + i) << Blancs(lcol[i], 1);
//	}
//	cout << endl;
//	for (i = 6; i < 9; i++) { // maintenant indice ligne
//		for (j = 0; j < 9; j++) {
//			if ((j == 3) || (j == 6))cout << "|";
//			int cell = 9 * i + j, digs = dig_cells[cell], ndigs = _popcnt32(digs);
//			ncand += ndigs;
//			for (int id = 0; id < 9; id++)if (digs & (1 << id))
//				cout << id + 1;
//			cout << Blancs(lcol[j] + 1 - ndigs, 1);
//		} // end for j
//		cout << endl;
//	} // end for i
//	cout << endl;
//
//}






void ZHOU::DebugDigit(int digit){
	cout << "DEBUG index=" << index << " digit=" << digit+1  << endl;
	for (int ib = 0; ib<3; ib++) {
		for (int ir = 0; ir<3; ir++){
			unsigned vf = FD[digit][0].bf.u32[ib],
				wr = (vf >> (9 * ir)) & 0x1ff;
			for (int k = 0; k<9; k++){
				if (wr & (1 << k))		cout << digit + 1;
				else 		cout << ".";
				if (k == 2 || k == 5) 	cout << " ";
			}
			cout << endl; //end of row
		}
		cout << endl; // end of block
	}
	cout << endl; // end of map per digit

}

