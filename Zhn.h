
/*
ZhouSolver.h
Based on code posted to <http://forum.enjoysudoku.com/3-77us-solver-2-8g-cpu-testcase-17sodoku-t30470.html>
by user zhouyundong_2012.
The copyright is not specified.

this is a remorphing of the code to work in recursive mode and in a 128 bit field
funtions have been added to work in generation mode
*/

/* the BIT_SET_30 pattern is defined in a band for a digit
   it is a map of candidates plus 3 bits for unsolved rows
   bits 0-8 9-17 18-26 for the rows 27-29 for unsolved rows
*/
//#include "t_128GP.h"
// tables specific to the brute force located in zh4_tables
const extern int TblRowMask[8];// rows where single  found  000 to 111
const extern int  Tblstartblock[27]; // (i%3) * 27 in zhou brute force
const extern int TblShrinkMask[512];// existing minirows 000 to 111
const extern int TblComplexMask[512]; // keep mini rows still valid optimised process
const extern int TblMaskSingle[512]; // kill in other blocks locked column /box
const extern int TblMaskDouble[512];// kill for locked in box / column
const extern int TblColumnSingle[512]; // single in column applied to shrinked bloc
const extern int TblShrinkSingle[512]; // keep only rows with single
const extern int TblRowUniq[512]; // 1 is row not defined in block  mode  to 111
const extern T128 AssignMask_Digit[81];
//const extern T128 AssignMask_OtherDigits[81];
struct ZHOU;

struct ZH_GLOBAL { // global variables for the game table
	uint64_t cpt[10],cptg[10],npuz;
	// morphing a puzzle at the start
	int x3_dmap_inv[9];// digit map 3 bands (start)
	int x3_cmap[81];// cells map at start
	char puz[82]; // the solved puzzle (after morph)
	//given and singles found
	GINT16 tgiven[81];
	int ngiven, tsingles[40], nsingles;

	BF128 digit_sol[9]; // final solution per digit original sort sequence
	int nsol, lim, icount, ntsol, single_applied, new_single_in_Update,
		go_back, goadduas,// usually 0, special treatment for a solution
		rdigit, loop, diag, modeguess,maxindex;
	BF128 Digit_cell_Assigned[9];
	BF128 init_3x, init_digit, pairs,triplets, pairs_naked;
	char * zsol, *puzfinal, *pat, 
		stdfirstsol[82],
		zerobased_sol[81];

	// switching to solver mode
	PM3X pm, pmdiag,pmelims;
	BF128  cells_unsolved_e, cells_unsolved_diag,// pm status direct and diagonal
		cells_assigned;// to avoid redundancy in new assignments 
	// in locked.. column "[2]" is in diagonal mode 
	BF128 locked_nacked_brc_seen[3],// seen nacked in row; column; box (priority box)
		locked_nacked_brc_done[3];// same cleaning done  nacked in row; column; box (priority box)
	int  dig_rows[9][9], dig_cols[9][9];//rows cols in 9 bits mode
	int  dig_boxes[9][9];// box in 9 bits mode
	int dig_cells[81], cells_count[81];//cells digits  in 9 bits mode
	int unsolved_r_count[9], unsolved_c_count[9]; // pm status unsolved rows columns per digit
	int row_col_x2[9][2], dig_unsolved_col[9], oldcount;
	BF128 digits_cells_pair_bf[9]; 
	ZHOU * zhou_current;

    //=================== floor analysis (one digit)
	int current_digit,active_floor;
	BF128  or_floor[9], elim_floor[9];

	// specific to the attempt to optimize the X+Y+27 process
	char *entry_base0, zdebug[82];
	int * grid0; // using a process with known solution grid
	// specific to symmetry of given generation
	USHORT * ptcor;
	// specific to the search 17 process
	int s2_ind, naddtable;// see go_17sol
	BF128 * addtable;
	int band3digits[9], band3nextua;// specific to 17 search
	uint64_t * digsols, b12nextua; // pointer to solution grid per digit
	ZH_GLOBAL();
	inline void InitCount(int elim){
		memset(cpt, 0, sizeof cpt);
		nsol = 0;
		lim = elim;
	}
	void MorphPat(char * ze);
	void Morph_digits(int morph);
	void Morph_digits(GINT16 * td, int nc);
	void Map_Morph_digits(GINT16 * td, int nc);
	void NoMorph();
	int InitSudoku();
	int Go_InitSudoku(char * ze);
	int Go_InitSudoku_NoMorph(char * ze);
	int Go_InitSolve(char * ze);
	int Go_InitSolve(GINT16 * td, int nc);
	void ValidPuzzle(ZHOU * z);
	void Debug();
	// located in go_0xxcpp
	inline void Init_Assign(){ nsingles = 0; cells_assigned.SetAll_0(); }
	void Pm_Status(ZHOU * z);
	void Pm_Status_End(ZHOU * z);// box and cells
	void AddSingle(int band, int vband);
	void AddSingleDiag(int band, int vband);
	void Build_digits_cells_pair_bf();


	// located in solver step
	void DebugNacked();

};
/* 2 BF 128 per digit
	equivalent to F and Fcomp in Zhou
	Last 32 bits in FD[0] is a  bits field for unknown rows
	Last 32 bits in FD[1]] contain  digit mapped
*/
// class encapsulating the brute force 
struct ZHOU{// size 32 bytes 
	BF128 FD[9][2];
	//ZHOUDIG zhd[9];
	BF128 cells_unsolved;
	//ZHOUDIG * dlast;//pointer to last unsolved
	int ndigits,index,unsolved_digits,box_hidden_pair;

	inline void Copy(ZHOU & o);
	inline void Assign(int digit, int cell, int xcell);
	int Upd1(int digit);
	int Update();
	int InitSudoku(GINT16 * t, int n);
	int InitSudoku(char * zpuz,int morph=1);
	char * SetKnown(char * zs);
	int ApplySingleOrEmptyCells_Band3();
	int ApplySingleOrEmptyCells_B12();
	int ApplySingleOrEmptyCells();
	int FullUpdate();
	void SetaCom(int digit, int cell, int xcell);
	inline void SetFloor(int cell){ FD[0][0] &= AssignMask_Digit[cell]; }
	inline void Seta_c(int digit, int cell){ SetaCom(digit, cell, C_To128[cell]); }
	inline void SetG(GINT16 x){ Seta_c(x.u8[1], x.u8[0]); }
	void Setcell(int cell);
	inline void GuessBivalueInCell(BF128 & wc);
	inline int GuessHiddenBivalue();
	int SolveHiddenPair_Box_Row();
	int SolveHiddenPair_BoxB();
	int SolveHiddenTriplet_Box_Row();
	int GuessHiddenTriplet();
	void Guess();
	// located in zh4_doc_debug

	//inline void Copy(ZHOU & o)	{ int rind = index; *this = o; index = rind; }

	inline int Unsolved_Count(){ return cells_unsolved.Count(); }
	// standard calls


	inline void ComputeNext(){ 
		if (FullUpdate()){
			Guess();
		} 
	}


	// other calls and functions in Zh4_calls_variants
	int CheckValidityQuick(char *puzzle);
	int GetSol(char *puzzle);
	int CheckOneStep(char *puzzle, char * zs);
	inline int StatusValid(GINT16 * t, int n){
		InitSudoku(t, n); return Isvalid();
	}

	int StartOneStep(GINT * t, int n);

	int PartialInitSudoku(GINT16 * t, int n);
	int EndInitSudoku( GINT16 * t, int n);
	int EndInitNextUa(ZHOU & o, int bf);// 17 search check know small uas in bloc

	int IsMinimale(GINT16 * t, int n);
	int IsMinimale(char * ze);
	void PatFinal();
	int GetFreeDigits_c(int cell){ return GetFreeDigits(C_To128[cell]); }
	int GetFreeDigits(int xcell);
	int GetSolvedDigitForCell(int cell);

	// inline small functions
	inline int IsOffCandidate_c(int dig, int cell){return FD[dig][0].Off_c(cell); }
	inline int IsOnCandidate_c(int dig, int cell){ return FD[dig][0].On_c(cell); }
	inline void ClearCandidate_x(int dig, int xcell){ FD[dig][0].clearBit(xcell); }
	inline void ClearCandidate_c(int dig, int cell){ FD[dig][0].clearBit(C_To128[cell]); }
	void CleanCellForDigits(int cell, int digits);
	int CleanCellsForDigits(BF128 &  cells, int digits);
	int Clean(PM3X &elims);
	int Isvalid();
	inline int CountDigs(){
		int n =0;
		for (int i = 0; i < 9; i++)
			n += FD[i][0].Count96();
		return n;
	}
	//================== find active floors
	void StartFloor(int digit, ZHOU & o);
	int  UpdateFloor();
	void GuessFloor();

	// debugging code or print code
	void Debug(int all = 0);
	void DebugDigit(int digit);
	void DebugActiveDigits();
	int GetAllDigits(int cell);
	void ImageCandidats();
	void ImageCandidats_b3();
	int FullUpdateAtStart();
	int CheckStatus();// located in solver_step 

	// located in go_0xx.cpp
	int Rate10_LastInUnit();
	int Rate12_SingleBox();
	int Rate15_SingleRow();
	int Rate15_SingleColumn();
	int RateSingle(GINT64 * t, int nt);
	int RateSingleBox(GINT64 * t, int nt);
	int RateSingleRow(GINT64 * t, int nt);
	int RateSingleCol(GINT64 * t, int nt);
	int RateSingleDiag(GINT64 * t, int nt);
	int CollectHiddenPairsBox(GINT64* tp, int & np,int lim=4);
	int CollectHiddenPairsRow(GINT64* tp, int & np, int lim = 4);
	int CollectHiddenPairsCol(GINT64* tp, int & np, int lim = 4);
	int Rate17_lockedBox_Assign();
	void FindNakedPairsTriplets_NoSingle();
	int Rate20_HiddenPair_Assign();
	int Rate23_SingleInCell_Assign();
	int CollectHiddenTripletsBox(GINT64* tp, int & np);
	int CollectHiddenTripletsRow(GINT64* tp, int & np);
	int CollectHiddenTripletsCol(GINT64* tp, int & np);
	int Rate25_HiddenTriplet_Assign();
	int Rate26_lockedBox();
	int Rate28_lockedRowCol();
	int CleanNake(int & naked, int unit, int iband,int modebr);
	int CleanNakeColumn(int & naked, int unit, int iband);
	int Rate30_NakedPair();
	int Rate32_XWing();

	int ApplyHidden(GINT64 * t, int nt,int mode);
	int ApplyHiddenColumn(GINT64 * t, int nt);

	int Rate34_HiddenPair();
	int Rate36_FindClean_NakedTriplet(int unit_triplet_cells,int unit,int iband,int mode);
	int Rate36_FindClean_NakedTripletCol(int unit_triplet_cells, int unit, int iband);
	int Rate36_NakedTriplet();
	int Rate38_SwordFish();
	int Rate40_HiddenTriplet();
	int Rate42_XYWing();
	int Rate44_XYZWing();
	int Rate50_NakedQuad();
	int Rate52_JellyFish();
	int Rate54_HiddenQuad();

	void AssignSolver(int print = 0);
	void XW_template(int idig);
	void Naked_Pairs_Seen();

	/*
    inline void SetPat(char * pat, char * zsol, char * puzfinal){
		zh_g.pat = pat; zh_g.zsol = zsol; zh_g.puzfinal = puzfinal;
	};
    // located in zh_uacollector.cpp  generation of UAs on a multi floor
	void InitGenUas(char * zpuz);
	void GenUas( int floors);
	void GenUasBands12(int floors);
	void GenUasBands12LessMinirow(int floors,int c1,int c2,int c3);
	void GenUas2();
	void GenUas3();
	void GenUas4();
	void GenUas5();
	void GenUas6();
	void GenMinirows2();
	void GenMinirows3();
	int Update2();
	int Update3();
	int Update4();
	int Update5();
	int Update6();
	void Guess2();
	void Guess3();
	void Guess4();
	void Guess6();
	void Guess5();
	void Guess5_4();
	void Guess5_3();
	int CollectFinal(BF128 *td,int &lim10);
   */
 };

 struct ZHBAND{// similar to  ZHOU for one band
	 int FD[9][2],cells_unsolved;
	 int ndigits, index, unsolved_digits;

	 inline void Copy(ZHBAND & o);
	 inline void Assign(int digit, int cell, int xcell);
	 int Update();
	 char * SetKnown(char * zs);
	 int ApplySingleOrEmptyCells_Band3();
	 int ApplySingleOrEmptyCells_B12();
	 int ApplySingleOrEmptyCells();
	 int FullUpdate();
	 void SetaCom(int digit, int cell);
	 void Setcell(int cell);
	 inline void GuessBivalueInCell(BF128 & wc);
	 inline int GuessHiddenBivalue();
	 int GuessHiddenTriplet();
	 void Guess();
	 inline void ComputeNext(){	 if (FullUpdate())	 Guess();	 	 }
	 int GetFreeDigits_c(int cell){ return GetFreeDigits(C_To128[cell]); }
	 int GetFreeDigits(int xcell);
	 int GetSolvedDigitForCell(int cell);
	 // debugging code or print code
	 void Debug(int all = 0);
	 int GetAllDigits(int cell);
	 void ImageCandidats();
	 void ImageCandidats_b3();
 };