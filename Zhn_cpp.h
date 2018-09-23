/*
Based on code posted to <http://forum.enjoysudoku.com/3-77us-solver-2-8g-cpu-testcase-17sodoku-t30470.html>
by user zhouyundong_2012.
The copyright is not specified.
*/
//====================================
//#include "sk_t.h"  // tab0  
//#include "Zhn.h"
#include "Zhtables_cpp.h" // also describes the pattern
//#define DIAG

ZH_GLOBAL zh_g;

ZHOU zhou[50]; // must host main brute force plus minimality analysis and recursive generation
ZHOU zhou_i,zhou_ip,//zhou_i================== initial for a GAME
    zhou_solve;// basis to solve a puzzle using elimination logic
//============================= ZH_GLOBAL code
ZH_GLOBAL::ZH_GLOBAL(){
	diag = 0;
	modeguess = 1;
	goadduas = ngiven=0;
	zsol = pat = puzfinal = 0; // no solution unless required buy the user
	// init Tblgame
	for (int i = 0; i < 3; i++)init_3x.bf.u32[i] = BIT_SET_27;
	init_3x.bf.u32[3] = 0;
	init_digit = init_3x;
	init_digit.bf.u32[3] = 0777;//all rows unsolved
	for (int i = 0; i < 9; i++)	{
		zhou_i.FD[i][0] = init_digit;
		zhou_i.FD[i][1] = init_3x;
	}
	zhou_i.cells_unsolved = init_3x;
	zhou_i.ndigits = 9;
	zhou_i.index = 0;
	zhou_i.box_hidden_pair = 0;
	zhou_i.unsolved_digits = 0777;
	NoMorph();
	if (0){
		cout << "gl0bal debug" << endl;
		Debug();		zhou_i.Debug(1);
	}
}
void ZH_GLOBAL::NoMorph(){
	for (int i = 0; i < 9; i++){
		x3_dmap_inv[i] = i;
		zhou_i.FD[i][1].bf.u32[3] = i;
	}
	for (int i = 0; i < 81; i++)x3_cmap[i] = i;
}
void ZH_GLOBAL::MorphPat(char * ze){// sort entry to optimize brute force
	char zdiag[81], *source;
	int count[6], *sourcei, zei[81], *zdiagi = C_transpose_d;
	__stosd((unsigned long*)count, 0, 6);
	for (int i = 0; i < 81; i++) {
		if (ze[i] - '.'){
			int band = i / 27, stack = C_stack[i] + 3;
			++count[band];			++count[stack];
		}
		zdiag[C_transpose_d[i]] = ze[i];
		zei[i] = i;
	}
	int imax = 0, nmax = 0;// find biggest count band stack
	for (int ib = 0; ib < 6; ib++)if (count[ib]>nmax){	nmax = count[ib];	imax = ib;	}
	// morph the puzzle to the high band in band
	if (imax > 2){
		source = zdiag;		sourcei = zdiagi;
		__movsd((unsigned long*)count, (unsigned long *) &count[3], 3);
	}
	else { source = ze; sourcei = zei; }
	int tsort[3], w;// sort bands increasing order of count
	tsort[0] = count[0] << 8;
	tsort[1] = 1 | (count[1] << 8);
	tsort[2] = 2 | (count[2] << 8);
	for (int i = 0; i < 2; i++) for (int j = i + 1; j < 3; j++)
		if (tsort[i]>tsort[j]){ w = tsort[i]; tsort[i] = tsort[j]; tsort[j] = w; }
	int ib1 = tsort[0] & 3, ib2 = tsort[1] & 3, ib3 = tsort[2] & 3;
	puz[81]=0,
	memcpy(puz, &source[27 * ib1], 27);
	memcpy(&puz[27], &source[27 * ib2], 27);
	memcpy(&puz[54], &source[27 * ib3], 27);
	__movsd((unsigned long*)x3_cmap, (unsigned long *)&sourcei[27 * ib1], 27);
	__movsd((unsigned long*)&x3_cmap[27], (unsigned long *)&sourcei[27 * ib2], 27);
	__movsd((unsigned long*)&x3_cmap[54], (unsigned long *)&sourcei[27 * ib3], 27);

}
//td is 8bits cell + 8 bits digits
void ZH_GLOBAL::Morph_digits(int morph){// using the given entry
	ngiven = 0;
	unsigned long count[9]; __stosd(count, 0, 9);
	for (int i = 0; i < 81; i++){
		register int c = puz[i];
		if (c<'1' || c>'9') continue;
		c -= '1';
		count[c]++;
		tgiven[ngiven++].u16 = i | (c << 8);
	}
	if (morph){// if morph asked, map digits on increasing count.
		int tsort[9], w;// sort bands increasing order of count
		for (int i = 0; i < 9; i++)tsort[i] = (count[i] << 8) | i;
		for (int i = 0; i < 8; i++) for (int j = i + 1; j < 9; j++)
			if (tsort[i]>tsort[j]){ w = tsort[i]; tsort[i] = tsort[j]; tsort[j] = w; }
		for (int i = 0; i < 9; i++){
			register int ii = tsort[i] & 15;
			zhou_i.FD[i][1].bf.u32[3] = ii;
			x3_dmap_inv[ii] = i;
		}
		// map the given to new order
		for (int ic = 0; ic < ngiven; ic++)
			tgiven[ic].u8[1] = x3_dmap_inv[tgiven[ic].u8[1]];
	}
}
void ZH_GLOBAL::Map_Morph_digits(GINT16 * td, int nc){//applying x3_cmap to the cells

}
void ZH_GLOBAL::Morph_digits(GINT16 * td, int nc){// must be in line with the morphed pattern
	ngiven = 0;
	unsigned long count[9]; __stosd(count, 0, 9);
	for (int it = 0; it < nc; it++){
		register int cell = td[it].u8[0], dig = td[it].u8[1];
		count[dig]++;
	}
	int tsort[9], w;// sort bands increasing order of count
	for (int i = 0; i < 9; i++)tsort[i] = (count[i] << 8) | i;
	for (int i = 0; i < 8; i++) for (int j = i + 1; j < 9; j++)
		if (tsort[i]>tsort[j]){ w = tsort[i]; tsort[i] = tsort[j]; tsort[j] = w; }
	for (int i = 0; i < 9; i++){
		register int ii = tsort[i] & 15;
		zhou_i.FD[i][1].bf.u32[3] = ii;
		x3_dmap_inv[ii] = i;
	}
	// map the given to new order
	for (int ic = 0; ic < ngiven; ic++)
		tgiven[ic].u8[1] = x3_dmap_inv[tgiven[ic].u8[1]];
}
int ZH_GLOBAL::InitSudoku(){ return zhou[0].InitSudoku(tgiven,ngiven); }
int ZH_GLOBAL::Go_InitSudoku(char * ze){
	MorphPat(ze);
	if (diag)cout << puz << " morphed puzzle" << endl;
	Morph_digits(1);
	//if (diag)Debug();
	return InitSudoku();
}
int ZH_GLOBAL::Go_InitSudoku_NoMorph(char * ze){
	ze[81] = 0;
	strcpy(puz, ze);
	NoMorph();
	ngiven = 0;
	int digs = 0;
    for (int i = 0; i < 81; i++){
		register int c = puz[i];
		if (c<'1' || c>'9') continue;
		c -= '1';
		digs |= 1 << c;
		tgiven[ngiven++].u16 = i | (c << 8);
	}
	if (_popcnt32(digs) < 8) return 1; // don't accept less than 8 digits given
	return InitSudoku();
}
int ZH_GLOBAL::Go_InitSolve(char * ze){
	zsol = stdfirstsol;
	InitCount(1);
	ze[81] = 0;
	strcpy(puz, ze);
	NoMorph();
	ngiven = 0;
	int digs = 0;
	for (int i = 0; i < 81; i++){
		register int c = puz[i];
		if (c<'1' || c>'9') continue;
		c -= '1';
		digs |= 1 << c;
		tgiven[ngiven++].u16 = i | (c << 8);
	}
	if (_popcnt32(digs) < 8) return 1; // don't accept less than 8 digits given
	if( InitSudoku()) return 1;
	zhou_solve = zhou[0];
	zhou[0].ComputeNext();
	if (nsol != 1) return 1;
	for (int i = 0; i < 81; i++)zerobased_sol[i] = stdfirstsol[i] - '1';
	__stosq((unsigned long long *)locked_nacked_brc_done[0].bf.u64, 0, 6);
	__stosd((unsigned long *)row_col_x2[0], 0, 18);

	return 0;
}

int ZH_GLOBAL::Go_InitSolve(GINT16 * td, int nc){
	zsol = stdfirstsol;
	InitCount(1);
	strcpy(puz, empty_puzzle);
	NoMorph();
	ngiven = nc;
	for (int i = 0; i < nc; i++){
		puz[td[i].u8[0]] = td[i].u8[1] + '1';
		tgiven[i] = td[i];
	}
	if (InitSudoku()) return 1;
	zhou_solve = zhou[0];
	zhou[0].ComputeNext();
	if (nsol != 1) return 1;
	for (int i = 0; i < 81; i++)zerobased_sol[i] = stdfirstsol[i] - '1';
	__stosq((unsigned long long *)locked_nacked_brc_done[0].bf.u64, 0, 6);
	__stosd((unsigned long *)row_col_x2[0], 0, 18);

	return 0;
}
void ZH_GLOBAL::ValidPuzzle(ZHOU * z){
	if (zsol && (!nsol)){// store the first solution
		z->SetKnown(zsol);
		for (int i = 0; i < 9; i++)digit_sol[i] = z->FD[i][0];
	}
	nsol++;
}



//============================= zhou code

inline void ZHOU::Copy(ZHOU & o){
	*this = o;
	index++;
}
inline void ZHOU::Assign(int digit, int cell, int xcell){
	FD[digit][0] &= AssignMask_Digit[cell];
	cells_unsolved.Clear(xcell);
}
void ZHOU::Setcell(int cell){ 
	SetaCom(zh_g.zerobased_sol[cell], cell, C_To128[cell]);
}

void ZHOU::CleanCellForDigits(int cell, int digits){
	unsigned long digit;
	while (_BitScanForward(&digit, digits)){
		digits ^= 1 << digit;
		ClearCandidate_c(digit, cell);
	}
}
int ZHOU::CleanCellsForDigits(BF128 &  cells, int digits){
	int iret = 0;
	unsigned long digit;
	while (_BitScanForward(&digit, digits)){
		digits ^= 1 << digit;
		BF128 clean = FD[digit][0] & cells;
		if (clean.isNotEmpty()){
			iret = 1;
			FD[digit][0] -= clean;
		}
	}
	return iret;
}
int ZHOU::Clean(PM3X & elims){// return 1 if something has been done
	int iret = 0;
	for (int idig = 0; idig < 9; idig++){
		BF128 w = FD[idig][0] & elims.pmdig[idig];
		if (w.isNotEmpty()){
			iret = 1;
			FD[idig][0] -= w;
		}
	}
	return iret;
}

#define UPD_AND_CL(W,X)	last_assigned = W+1;cur_assigned = 0;\
B = ~(A & TblRowMask[S]);\
cells_unsolved.bf.u32[X] &= B;\
FD[0][0].bf.u32[X] &= B; FD[1][0].bf.u32[X] &= B; FD[2][0].bf.u32[X] &= B; \
FD[3][0].bf.u32[X] &= B; FD[4][0].bf.u32[X] &= B; FD[5][0].bf.u32[X] &= B; \
FD[6][0].bf.u32[X] &= B; FD[7][0].bf.u32[X] &= B; FD[8][0].bf.u32[X] &= B; }\
FD[W][0].bf.u32[X] = FD[W][1].bf.u32[X] = A;

#define UPD_012(W,X,Y,Z)A = FD[W][0].bf.u32[X];\
if (A != FD[W][1].bf.u32[X]){\
Shrink = (TblShrinkMask[A & 0x1FF]\
| TblShrinkMask[(A >> 9) & 0x1FF] << 3\
| TblShrinkMask[(A >> 18)] << 6);\
if ((A &= TblComplexMask[Shrink]) == 0)  return 0;\
S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);\
B = FD[W][0].bf.u32[Y];\
FD[W][0].bf.u32[Z] &= TblMaskSingle[S] & TblMaskDouble[S | ((B | (B >> 9) | (B >> 18)) & 0x1FF)];\
B = FD[W][0].bf.u32[Z];\
FD[W][0].bf.u32[Y] &= TblMaskSingle[S] & TblMaskDouble[S | ((B | (B >> 9) | (B >> 18)) & 0x1FF)];\
S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];

#define UPD_ONE_DIGIT(W,V) if (cur_assigned > W+1)goto exit_digits;\
	if (FD[W][0].bf.u64[0] != FD[W][1].bf.u64[0]\
|| FD[W][0].bf.u32[2] != FD[W][1].bf.u32[2]){\
r_free = FD[W][0].bf.u32[3];\
UPD_012(W, 0, 1, 2)	if ((r_free & 7) != S){\
r_free &= 0770 | S;	UPD_AND_CL(W, 0)}\
UPD_012(W, 1, 0, 2)	if (((r_free >> 3) & 7) != S){\
r_free &= 0707 | (S << 3);	UPD_AND_CL(W, 1)}\
UPD_012(W, 2, 0, 1)	if (((r_free >> 6) & 7) != S){\
r_free &= 077 | (S << 6);	UPD_AND_CL(W, 2)}\
FD[W][0].bf.u32[3] = r_free;if (!r_free){\
unsolved_digits &= V;ndigits--;}}


int ZHOU::Update(){
	if (zh_g.diag)cout << "Update index=" << index << endl;
	register uint32_t Shrink = 1, r_free, B, A, S, last_assigned = 0,cur_assigned;
loop_upd:
	zh_g.cpt[3]++;
	cur_assigned = last_assigned; last_assigned = 0;
	UPD_ONE_DIGIT(8, 0377) UPD_ONE_DIGIT(7, 0577) UPD_ONE_DIGIT(6, 0677)
		//if (zh_g.diag){ cout << "apr�s 987" << endl; Debug(); }
		//if (cur_assigned > 5)goto exit_digits;
	UPD_ONE_DIGIT(5, 0737) UPD_ONE_DIGIT(4, 0757) UPD_ONE_DIGIT(3, 0767)
		//if (zh_g.diag){ cout << "apr�s 654" << endl; Debug(0); }
	//if (cur_assigned > 3)goto exit_digits;
	UPD_ONE_DIGIT(2, 0773) UPD_ONE_DIGIT(1, 0775) UPD_ONE_DIGIT(0, 0776)
	exit_digits:
	if (zh_g.diag>1){ cout << "avant test loop last_assigned=" << last_assigned << endl;
	Debug(1); ImageCandidats();
	}
	if (last_assigned) goto loop_upd;// nothing to do in the next cycle
	return 1;
}
int ZHOU::Upd1(int digit){
	//if (1) return 1;
	if (zh_g.diag)cout << "Upd index=" << index << " update digit=" <<digit << endl;
	register uint32_t Shrink = 1, r_free , B, A, S; 
	int cur_assigned = 0;
digitloop:
	int last_assigned = 0;
	switch (digit){
	case 0:UPD_ONE_DIGIT(0, 0776) break;
	case 1:UPD_ONE_DIGIT(1, 0775) break;
	case 2:UPD_ONE_DIGIT(2, 0773)  break;
	case 3:UPD_ONE_DIGIT(3, 0767)  break;
	case 4:UPD_ONE_DIGIT(4, 0757)  break;
	case 5:UPD_ONE_DIGIT(5, 0737)  break;
	case 6:UPD_ONE_DIGIT(6, 0677)  break;
	case 7:UPD_ONE_DIGIT(7, 0577)  break;
	case 8:UPD_ONE_DIGIT(8, 0377) break;
	}
	exit_digits:
		if (last_assigned) goto digitloop;// nothing to do in the next cycle
	return 1;
}
int ZHOU::InitSudoku(GINT16 * t, int n){// if morph, done before
	BF128 Digit_cell_Assigned[9];
	__stosq((uint64_t*)Digit_cell_Assigned, 0, 18);
	*this = zhou_i;
	for (int ic = 0; ic < n; ic++)   {
		int digit = t[ic].u8[1], cell = t[ic].u8[0], xcell = C_To128[cell];
		if (FD[digit][0].Off(xcell))  return 1;// check not valid entry
		Assign(digit, cell, xcell);
		Digit_cell_Assigned[digit].Set(xcell);
	}
	BF128 w = cells_unsolved; w.bf.u32[3] = ~0;
	for (int i = 0; i<9; i++)  FD[i][0] &= w | Digit_cell_Assigned[i];
	return 0;
}
int ZHOU::InitSudoku(char * zpuz, int morph){
	if (!morph){
		zh_g.NoMorph();
		__movsb((unsigned char *)zh_g.puz, (unsigned char *)zpuz, 82);
	}
	else zh_g.MorphPat(zpuz);
	zh_g.Morph_digits(morph);
	return InitSudoku(zh_g.tgiven, zh_g.ngiven);
}
char * ZHOU::SetKnown(char * zs){
	strcpy(zs, empty_puzzle);
	BF128  fd;
	for (int idigit = 0; idigit < 9; idigit++){
		fd = FD[idigit][0]; // unsolved or partially solved, still in zhou
		int digit = FD[idigit][1].bf.u32[3];
		if (fd.bf.u32[3] == 0777) continue;
		for (int ib = 0; ib < 3; ib++){// 3 blocs per digit
			int arows = (fd.bf.u32[3] >> (3 * ib)) & 7;
			if (arows == 7) continue; // not assigned
			unsigned int band = fd.bf.u32[ib];
			for (int j = 0; j<3; j++) if (!(arows & (1 << j))) {
				int row = (band >> TblMult9[j]) & 0x1ff;
				unsigned long  irow;
				_BitScanForward(&irow, row);
				int	cell = Tblstartblock[ib] + TblMult9[j] + irow;
				zs[cell] = digit + '1';
			}
		}
	}
	return zs;
}

void ZHOU::SetaCom(int digit, int cell, int xcell){ // single in cell
//cout<<"index"<<(int)index << "setacom d;c;xc\t" << digit + 1 << "\t" << cell << "\t" << xcell << endl;
if (zh_g.diag){
cout << zh_g.cpt[0] << "setacom digit=" << digit + 1 << " cell=" << cell << endl;
ImageCandidats();
}
BF128 *Fd = FD[digit];
*Fd &= AssignMask_Digit[cell];
cells_unsolved.Clear(xcell);
BF128 * RF = FD[8];
for (; RF >= FD[0]; RF -= 2)RF->Clear(xcell);
Fd->setBit(xcell); // restore bit for digit assigned
}

int ZHOU::Isvalid(){ // usually after init 2 steps
	zh_g.nsol = 0; zh_g.lim = 1; ComputeNext();  return zh_g.nsol;
}
int ZHOU::CheckValidityQuick(char *puzzle){
	zh_g.nsol = 0; zh_g.lim = 1;
	if (zh_g.Go_InitSudoku_NoMorph(puzzle)) return 0;
	//if (zh_g.Go_InitSudoku(puzzle)) return 0;;
	//if (ApplySingleOrEmptyCells())	return 0; // pas bon
	if (0){
		Debug(1);
		ImageCandidats();
		cout << "call computenext() unsolved digits= 0"<<oct
			<< unsolved_digits<<dec << endl;
		//return 0;
	}

	ComputeNext();
	return zh_g.nsol;
}
int ZHOU::IsMinimale(GINT16 * to, int no){// assumed checked valid before
	zh_g.lim = 0;
	GINT16 td[81], cx;
	int nd;
	for (int i = 0; i < no; i++){
		nd = 0;
		for (int i2 = 0; i2 < no; i2++){
			if (i2 - i)td[nd++] = to[i2];
			else cx = to[i2];
		}
		zh_g.nsol = 0;
		InitSudoku(td, nd); // always valid
		ClearCandidate_c(cx.u8[1], cx.u8[0]);// clean the valid known
		//cout << "check minimale step" << endl;
		//ImageCandidats();
		ComputeNext();
		if (!zh_g.nsol) return 0; //could check also not valid
	}
	return 1;
}

#include "Zhn_doc_debug_cpp.h"
/*Zhn Control of flow
-b9- use of hidden pairs triplets zhg_modeguess
 abc  a use pairs  b use triplets c use pairs not hidden pair
 -v7- max index for hidden pairs ...
*/

int ZHOU::FullUpdate(){
	if (zh_g.nsol > zh_g.lim) return 0;
	zh_g.cpt[5]++;
	while (1){
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (zh_g.diag>1)	{
			cout << "look for singles in cells zh_g.cpt[5]=" << zh_g.cpt[5] << endl;
			//ImageCandidats();
		}

		if (ApplySingleOrEmptyCells())	return 0; // locked empty cell or conflict singles in cells
		if (zh_g.diag>1 && zh_g.single_applied)	{
			cout << "full update boucle  update zh_g.cpt[5]=" << zh_g.cpt[5] << endl;
			Debug(1);
			ImageCandidats();
		}
		if (!zh_g.single_applied)	break;
	}
	if (zh_g.diag)		cout<<"index="<<index << " exit fullupdate  updatecompte full=" 
		<< zh_g.cpt[5]	<< " compte upd=" << zh_g.cpt[3] << endl;
	return 1;
}
void ZHOU::Guess(){// low chance to be solved here
	if (zh_g.nsol > zh_g.lim) return;
	if (cells_unsolved.isEmpty()){
		if (zh_g.diag)cout <<endl<< "valid guess1\n" << endl;
		zh_g.ValidPuzzle(this);
		return;
	}
	if (index<=zh_g.maxindex && zh_g.modeguess){
		//cout << "first entry guess see more " << endl;
		int ir = SolveHiddenPair_Box_Row();
		//Debug(1);
		if (zh_g.diag)ImageCandidats();
		if (ir){ComputeNext();		return; }
		if (zh_g.modeguess & 2){
			ir = SolveHiddenTriplet_Box_Row();
			if (ir){ ComputeNext();		return; }
					}
		if (zh_g.modeguess & 4){
			if (zh_g.pairs_naked.isNotEmpty()){
				zh_g.cpt[1]++;
				GuessBivalueInCell(zh_g.pairs_naked);
				return; 
			}
		}
	}
	zh_g.cpt[1]++;
	if (zh_g.diag)		cout << "entry guess index=" << (int)index << " ndigits=" << ndigits << endl;
	if (zh_g.pairs.isNotEmpty()){ GuessBivalueInCell(zh_g.pairs);	return; }
	if (zh_g.diag)		cout << " guess pas biv cell " << endl;
	if (GuessHiddenBivalue()) return;
	if (zh_g.diag)		cout << " guess pas hidden biv " << endl;
	if (GuessHiddenTriplet()) return;
	int xcell = cells_unsolved.getFirst96(),cell=From_128_To_81[xcell];
	for (int idig = 0; idig < 9; idig++){
		if (FD[idig][0].On(xcell)){// one valid digit		
			ZHOU * mynext = this + 1; // start next guess
			mynext->Copy(*this);
			mynext->SetaCom(idig, cell, xcell);
			mynext->ComputeNext();
		}
	}
}

#define NAKED64(X) 	Map=FD[X][0].bf.u64[0];R3|=R2&Map;R2|=R1&Map;R1|=Map;
#define NAKED32(X) 	Map=FD[X][0].bf.u32[2];R3|=R2&Map;R2|=R1&Map;R1|=Map;


int ZHOU::ApplySingleOrEmptyCells_Band3(){
	register uint32_t  Map = FD[0][0].bf.u32[2], R3 = FD[1][0].bf.u32[2],
		R2 = Map & R3, R1 = Map | R3;// digits 12
	Map = FD[2][0].bf.u32[2]; R3 = R2&Map; R2 |= R1&Map; R1 |= Map;	
	NAKED32(3) NAKED32(4) NAKED32(5) NAKED32(6) NAKED32(7) NAKED32(8)

	if (cells_unsolved.bf.u32[2] & (~R1)) return 1; // locked
	R1 &= ~R2;	R1 &= cells_unsolved.bf.u32[2]; // forget solved seen as singles
	if (R1) zh_g.single_applied = 1;
	else{	zh_g.pairs.bf.u32[2] = R2& (~R3);	return 0;	}
	unsigned long res;
	while (_BitScanForward(&res, R1)){// usually a very small number of cells to assign
		R2 = 1<< res; // switch to the bit value
		R1 &= ~R2;  // clear the bit
		if (R2 & FD[0][0].bf.u32[2]){ Assign(0, res + 54, res + 64); goto nextr1; }
		if (R2 & FD[1][0].bf.u32[2]){ Assign(1, res + 54, res + 64); goto nextr1; }
		if (R2 & FD[2][0].bf.u32[2]){ Assign(2, res + 54, res + 64); goto nextr1; }
		if (R2 & FD[3][0].bf.u32[2]){ Assign(3, res + 54, res + 64); goto nextr1; }
		if (R2 & FD[4][0].bf.u32[2]){ Assign(4, res + 54, res + 64); goto nextr1; }
		if (R2 & FD[5][0].bf.u32[2]){ Assign(5, res + 54, res + 64); goto nextr1; }
		if (R2 & FD[6][0].bf.u32[2]){ Assign(6, res + 54, res + 64); goto nextr1; }
		if (R2 & FD[7][0].bf.u32[2]){ Assign(7, res + 54, res + 64); goto nextr1; }
		if (R2 & FD[8][0].bf.u32[2]){ Assign(8, res + 54, res + 64); goto nextr1; }
		return 1; //conflict with a previous cell assign
	nextr1:;
	}
	return 0;// not locked 
}
int ZHOU::ApplySingleOrEmptyCells_B12(){
	register uint64_t  Map = FD[0][0].bf.u64[0], R3 = FD[1][0].bf.u64[0],
		R2 = Map & R3, R1 = Map | R3; // digits 12
		Map = FD[2][0].bf.u64[0]; R3 = R2&Map; R2 |= R1&Map; R1 |= Map;
	NAKED64(3) NAKED64(4) NAKED64(5) NAKED64(6) NAKED64(7) NAKED64(8)

		if (cells_unsolved.bf.u64[0] & (~R1)) return 1; // locked
	R1 &= ~R2;	R1 &= cells_unsolved.bf.u64[0]; // forget solved seen as singles
	if (R1) zh_g.single_applied = 1;
	else{		zh_g.pairs.bf.u64[0] = R2& (~R3);		return 0;	}
	unsigned long res;
	while (_BitScanForward64(&res, R1)){// usually a very small number of cells to assign
		R2= 1; R2 <<= res; // switch to the bit value
		R1 &= ~R2;  // clear the bit
		R3 = From_128_To_81[res];
		if (R2 & FD[0][0].bf.u64[0]){ Assign(0, (int)R3, res); goto nextr1; }
		if (R2 & FD[1][0].bf.u64[0]){ Assign(1, (int)R3, res); goto nextr1; }
		if (R2 & FD[2][0].bf.u64[0]){ Assign(2, (int)R3, res); goto nextr1; }
		if (R2 & FD[3][0].bf.u64[0]){ Assign(3, (int)R3, res); goto nextr1; }
		if (R2 & FD[4][0].bf.u64[0]){ Assign(4, (int)R3, res); goto nextr1; }
		if (R2 & FD[5][0].bf.u64[0]){ Assign(5, (int)R3, res); goto nextr1; }
		if (R2 & FD[6][0].bf.u64[0]){ Assign(6, (int)R3, res); goto nextr1; }
		if (R2 & FD[7][0].bf.u64[0]){ Assign(7, (int)R3, res); goto nextr1; }
		if (R2 & FD[8][0].bf.u64[0]){ Assign(8, (int)R3, res); goto nextr1; }

		return 1; //conflict with a previous cell assugn
	nextr1:;
	}
	return 0;// not locked 
}
int ZHOU::ApplySingleOrEmptyCells(){
	zh_g.single_applied = 0;
	if (ApplySingleOrEmptyCells_Band3()) return 1;
	return ApplySingleOrEmptyCells_B12();
}

inline void ZHOU::GuessBivalueInCell(BF128 & wc){// priority to highest digits
	unsigned long res;
	register uint64_t R3 = wc.bf.u64[0];
	if (_BitScanForward64(&res, R3)){// first pair is ok to go
		goto gopair;
	}
	R3 = wc.bf.u64[1];
	if (_BitScanForward64(&res, R3)){
		res += 64;
		goto gopair;
	}
	return; // no bi value in cell  can not be
gopair: // do digit update for the 2 digits
	int cell = From_128_To_81[res],tdig[2],ndig=0;
	for (int idig = 0; idig < 9; idig++) 
		if (FD[idig][0].On(res))tdig[ndig++]=idig;
	ZHOU * mynext = this + 1; // start next guess
	mynext->Copy(*this);
	mynext->SetaCom(tdig[0], cell, res);
	mynext->Upd1( tdig[0]); // try to save one update cycle
	mynext->ComputeNext();
	SetaCom(tdig[1], cell, res);
	Upd1( tdig[1]);// try to save one update cycle
	ComputeNext();
}
int ZHOU::GuessHiddenBivalue(){// look for a single or a hidden pair in row or box
	if (zh_g.diag)		cout << "guess hidden biv" << endl;
	uint32_t hidden;
	int idig;
	int dcell, dxcell;
	for (idig = 0; idig <9; idig++){// priority to high digits last done
		BF128 & fd = FD[idig][0];
		register int Rows;
		if (!(Rows = fd.bf.u32[3])) continue;
		if (Rows & 7){//try row band1
			dcell = dxcell = 0;
			register int  band = fd.bf.u32[0];
			hidden = band & 0777;			if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 0777000;		if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 0777000000;		if (_popcnt32(hidden) == 2)goto exitok;
		}
		if (Rows & 070)	{// try row band2
			dcell = 27; dxcell = 32;
			register int  band = fd.bf.u32[1];
			hidden = band & 0777;			if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 0777000;		if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 0777000000;		if (_popcnt32(hidden) == 2)goto exitok;
		}
		if (Rows & 0700)	{//try row band3
			dcell = 54; dxcell = 64;
			register int  band = fd.bf.u32[2];
			hidden = band & 0777;			if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 0777000;		if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 0777000000;		if (_popcnt32(hidden) == 2)goto exitok;
		}
	}
	// no bi value row, try bi value box
	for (idig = 0; idig <9; idig++){// priority to high digits last done
		BF128 & fd = FD[idig][0];
		register int Rows;
		if (!(Rows = fd.bf.u32[3])) continue;
		if (Rows & 7){//try bow band1
			dcell = dxcell = 0;
			register int  band = fd.bf.u32[0];
			hidden = band & 07007007;		if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 070070070;		if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 0700700700;		if (_popcnt32(hidden) == 2)goto exitok;
		}
		if (Rows & 070)	{// try box band2
			dcell = 27; dxcell = 32;
			register int  band = fd.bf.u32[1];
			hidden = band & 07007007;		if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 070070070;		if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 0700700700;		if (_popcnt32(hidden) == 2)goto exitok;
		}
		if (Rows & 0700)	{//try bow band3
			dcell = 54; dxcell = 64;
			register int  band = fd.bf.u32[2];
			hidden = band & 07007007;		if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 070070070;		if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 0700700700;		if (_popcnt32(hidden) == 2)goto exitok;
		}
	}
	return 0;
exitok:
	unsigned long res;
	_BitScanForward(&res, hidden);
	hidden ^= 1 << res; //clear bit
	ZHOU * mynext = this + 1; // start next guess
	mynext->Copy(*this);
	mynext->SetaCom(idig, res + dcell, res + dxcell);
	mynext->Upd1(idig);
	mynext->ComputeNext();
	if (zh_g.diag)cout << "second digit biv last" << endl;
	_BitScanForward(&res, hidden);
	SetaCom(idig, res + dcell, res + dxcell);
	Upd1(idig);
	ComputeNext();
	return 1;
}
int ZHOU::GuessHiddenTriplet(){// look for a triplet in row or box
	if (zh_g.diag)		cout << "guess hidden triplet" << endl;
	register uint32_t hidden;
	int idig;
	for (idig= 0; idig <9; idig++){// priority to high digits
		BF128 & fd = FD[idig][0];
		register int rows = fd.bf.u32[3] >> 6;
		if (!rows)continue;// solved in band3
		register int  band = fd.bf.u32[2], maskrow = 0777, maskbox = 07007007;
		for (int i = 0; i < 3; i++, rows >>= 1, maskrow <<= 9, maskbox <<= 3){// band or box
			if (rows & 1){
				hidden = band &maskrow;		if (_popcnt32(hidden) ==3)goto exitok;
			}
			hidden = band &maskbox;		if (_popcnt32(hidden) == 3)goto exitok;
		}
	}
	return 0;
exitok:
	int ndig = 3;
	uint32_t w = hidden;
	unsigned long res;
	while (_BitScanForward(&res, w)){
		w ^= 1 << res; //clear bit
		if (--ndig){//not last cell
			ZHOU * mynext = this + 1; // start next guess
			mynext->Copy(*this);
			mynext->SetaCom(idig, res + 54, res + 64);
			mynext->ComputeNext();
		}
		else{// this is the last do it in the same place
			SetaCom(idig, res + 54, res + 64);
			ComputeNext();
			return 1;
		}
	}
	return 1;// should never be reached
}


int ZHOU::PartialInitSudoku(GINT16 * t, int n){// if morph, done before
	//BF128 Digit_cell_Assigned[9];
	__stosq((uint64_t*)zh_g.Digit_cell_Assigned, 0, 18);
	*this = zhou_i;
	for (int ic = 0; ic < n; ic++)   {
		int digit = t[ic].u8[1], cell = t[ic].u8[0], xcell = C_To128[cell];
		if (FD[digit][0].Off(xcell))  return 1;// check not valid entry
		Assign(digit, cell, xcell);
		zh_g.Digit_cell_Assigned[digit].Set(xcell);
	}
	BF128 w = cells_unsolved; w.bf.u32[3] = ~0;
	for (int i = 0; i<9; i++)  FD[i][0] &= w | zh_g.Digit_cell_Assigned[i];
	return 0;
}

int ZHOU::EndInitSudoku(GINT16 * t, int n){// if morph, done before
	BF128 Digit_cell_Assigned[9];
	__movsq((uint64_t*)Digit_cell_Assigned, (uint64_t*)zh_g.Digit_cell_Assigned, 18);
	*this = zhou_ip;
	for (int ic = 0; ic < n; ic++)   {
		int digit = t[ic].u8[1], cell = t[ic].u8[0], xcell = C_To128[cell];
		if (FD[digit][0].Off(xcell))  return 1;// check not valid entry
		Assign(digit, cell, xcell);
		Digit_cell_Assigned[digit].Set(xcell);
	}
	BF128 w = cells_unsolved; w.bf.u32[3] = ~0;
	for (int i = 0; i<9; i++)  FD[i][0] &= w | Digit_cell_Assigned[i];
	return 0;
}


int ZHOU::SolveHiddenPair_Box_Row(){// look for a single or a hidden pair in row or box
	if (zh_g.diag)		cout << "solve hidden pair box row" << endl;
	int vret = 0;
	zh_g.pairs_naked.SetAll_0();
	// process 9 boxes
	for (int iband = 0, dband=0, ibox = 0; iband < 3; iband++,dband+=27){
		int bandpairs = zh_g.pairs.bf.u32[iband];
		//cout << Char27out(bandpairs) << " band pairs" << endl;
		for (int iboxr = 0, mask = 07007007; iboxr < 3; iboxr++, ibox++, mask <<= 3){
			if (box_hidden_pair & (1 << ibox)) continue;

			//cout << "box=" << ibox + 1 << endl;
			int tbiv[9], digbiv[9], nbiv = 0;
			for (int idig = 0; idig < 9; idig++){// collect digits bi values in box
				int hidden = mask & FD[idig][0].bf.u32[iband];
				if (_popcnt32(hidden) == 2){
					digbiv[nbiv] = idig;
					tbiv[nbiv++] = hidden;
					//cout << Char27out(hidden) << " biv dig=" << idig + 1 << endl;
				}
			}

			if (nbiv < 2) continue; // no bi-value in the box
			for (int i = 0; i < nbiv - 1; i++){
				register int R = tbiv[i];
				for (int j = i + 1; j < nbiv; j++) if (R == tbiv[j]){// a pair seen
					if (zh_g.diag)	cout << Char27out(R) << " pair seen ibox=" << ibox << endl;
					if ((bandpairs & R) == R){
						zh_g.pairs_naked.bf.u32[iband]|=R;// give priority to solved naked pairs
						continue; //no effect
					}
					if (zh_g.diag)cout << "pair applied ibox=" << ibox << endl;
					vret = 1;
					R = ~R;// clean other cells
					box_hidden_pair |= 1 << ibox;
					for (int dig = 0; dig < 9; dig++)
						if (dig != digbiv[i] && dig != digbiv[j]) 
							FD[dig][0].bf.u32[iband]&=R;
				}
			}
		}
	}
	zh_g.cpt[6] += vret;
	return vret;
}

int ZHOU::SolveHiddenPair_BoxB(){// hidden pair on digit pairs
	if (zh_g.diag)		cout << "solve hidden pair box row" << endl;
	int vret = 0;
	for (int idig1 = 0; idig1 < 8; idig1++)		for (int idig2 = idig1 + 1; idig2 < 9; idig2++){
		BF128 fdw = FD[idig1][0] | FD[idig2][0];
		fdw &= cells_unsolved;
	}

	return vret;
}

int ZHOU::SolveHiddenTriplet_Box_Row(){// look a hidden triplet in row or box
	if (zh_g.diag)		cout << "solve hidden triplet box row" << endl;
	int vret = 0;
	zh_g.pairs_naked.SetAll_0();
	// process 9 boxes
	for (int iband = 0, dband = 0, ibox = 0; iband < 3; iband++, dband += 27){
		int bandpairs = zh_g.pairs.bf.u32[iband];
		for (int iboxr = 0, mask = 07007007; iboxr < 3; iboxr++, ibox++, mask <<= 3){
			if (box_hidden_pair & (1 << ibox)) continue;
			if (_popcnt32(cells_unsolved.bf.u32[iband] & mask) < 5) continue;
			int t3[9], dig3[9], n3 = 0;
			for (int idig = 0; idig < 9; idig++){// collect digits bi values in box
				int hidden = mask & FD[idig][0].bf.u32[iband], n = _popcnt32(hidden);
				if (n>1 && n<4){
					dig3[n3] = idig;
					t3[n3++] = hidden;
					//cout << Char27out(hidden) << " triplet dig=" << idig + 1 << endl;
				}
			}

			if (n3 < 3) continue; // no triplet in the box
			for (int i = 0; i < n3 - 2; i++){
				register int R = t3[i];
				for (int j = i + 1; j < n3 - 1; j++){
					register int R2 = R | t3[j];
					if (_popcnt32(R2) > 3) continue;
					for (int k = j + 1; k < n3; k++){
						register int R3 = R2 | t3[k];
						if (_popcnt32(R3) != 3) continue;
						int aig = 0;
						for (int dig = 0; dig < 9; dig++)
							if (dig != dig3[i] && dig != dig3[j] && dig != dig3[k]){
								if (FD[dig][0].bf.u32[iband] & R3){
									vret = aig=1;
									box_hidden_pair |= 1 << ibox;
									//cout << "active" << endl;
									FD[dig][0].bf.u32[iband] &= ~R3;
								}
							}
						if (aig){
							//cout << Char27out(R3) << " hidden triplet seen puz=" << zh_g.npuz
							//	<< " box=" << ibox  << endl;
							//ImageCandidats();
						}
					}
				}
			}
		}
	}
	zh_g.cpt[7] += vret;
	return vret;
}
int  ZHOU::UpdateFloor(){
	register uint32_t Shrink = 1, r_free, B, A, S, loop;
loop_upd:
	loop = 0;
	if (FD[0][0].bf.u64[0] != FD[0][1].bf.u64[0]
		|| FD[0][0].bf.u32[2] != FD[9][1].bf.u32[2]){
		r_free = FD[0][0].bf.u32[3];
		UPD_012(0, 0, 1, 2)	if ((r_free & 7) != S){
			r_free &= 0770 | S;	loop = 1;	}	}
	    UPD_012(0, 1, 0, 2)	if (((r_free >> 3) & 7) != S){
		    r_free &= 0707 | (S << 3);	loop = 1;	}}
        UPD_012(0, 2, 0, 1)	if (((r_free >> 6) & 7) != S){
	        r_free &= 077 | (S << 6);	loop = 1;}}
		FD[0][0].bf.u32[3] = r_free;
	}
	if (loop) goto loop_upd;// nothing to do in the next cycle
	return 1;
}
void ZHOU::StartFloor(int digit, ZHOU & o){
	//cout << "start floor for digit " << digit + 1 << endl;
	zh_g.current_digit = digit;
	FD[0][0] = o.FD[digit][0];
	FD[0][1].SetAll_0();
	zh_g.or_floor[digit].SetAll_0();
	zh_g.elim_floor[digit] = FD[0][0];
	zh_g.elim_floor[digit].bf.u32[3] = 0;// be sure to do nothing here
	//DebugDigit(0);
	GuessFloor();
	if (zh_g.or_floor[digit] == zh_g.elim_floor[digit]) return; // no elim
	zh_g.elim_floor[digit] -= zh_g.or_floor[digit];
	zh_g.active_floor |= 1 << digit;
	if (0){
		char ws[82];
		cout << "active floor for digit " << digit + 1 << endl;
		cout << zh_g.elim_floor[digit].String3X(ws) << " elims" << endl;
	}
}

void ZHOU::GuessFloor(){
	BF128 & fd = FD[0][0];

	if (!fd.bf.u32[3]){
		//cout << "une solution" << endl;
		zh_g.or_floor[zh_g.current_digit] |= fd;
		//char ws[82];
		//cout << zh_g.or_floor[zh_g.current_digit].String3X(ws) << " new or" << endl;
		return;
	}
	int hidden, rows = fd.bf.u32[3], dcell = 0, dxcell = 0, bfree = 7, ccmin = 10, min,bmin;
	unsigned long res;
	for (int iband = 0; iband < 3; iband++,dcell+=27,dxcell+=32,bfree<<=3){
		if (!(rows&bfree)) continue; // solved band
		int  band = fd.bf.u32[iband];
		for (int irb = 0; irb < 6; irb++){// boxes and rows of the band
			hidden = band & tband_box[irb];
			int cc = _popcnt32(hidden);
			if (cc < 2)continue;
			if (cc == 2) goto exitok;
			if ( cc < ccmin){ ccmin = cc; min = hidden;	bmin = iband; }
		}
	}
	//cout << "exit no bivalue ccmin=" << ccmin << endl;
	//DebugDigit(0);
	if (ccmin == 10)return;// should never be
	dcell = 27 * bmin;
	dxcell = 32 * bmin;
	while (_BitScanForward(&res, min)){// no bi value, use the smallest set
		min ^= 1 << res;// clear the bit
		ZHOU * mynext = this + 1; // start next guess
		mynext->Copy(*this);
		mynext->SetaCom(0, res + dcell, res + dxcell);
		mynext->Upd1(0);
		mynext->GuessFloor();
	}
	return;

exitok:
	//cout << Char27out(hidden) << " guess band dcell=" << dcell << endl;

	_BitScanForward(&res, hidden);
	ZHOU * mynext = this + 1; // start next guess
	mynext->Copy(*this);
	mynext->SetaCom(0, res + dcell, res + dxcell);
	mynext->Upd1(0);
	mynext->GuessFloor();
	_BitScanReverse(&res, hidden);
	//if (1)cout <<Char27out(1<<res)<< "second digit biv last dcell=" <<dcell<< endl;
	SetaCom(0, res + dcell, res + dxcell);
	Upd1(0);
	GuessFloor();
}

