
//#include "main.h"

struct GINTBUF{//Storing variable lentgh chains .. start one 64 bit limit
	GINT * buf;
	int n32, ncur,nback;
	inline void Set(GINT * b, int en32){ buf = b; n32=en32; }
	inline void Init(){ ncur = 0; }
	int Get(int n);
	int Store(GINT * t);
	inline int InitGetback(){ nback = 0; return ncur; }
	int GetBack(GINT * d);
	void EraseOlds(int nkeep);
	inline void ResetTo(int x){ ncur = x; }
};
struct EXPLAIN_BUFFER{// kind of stack to store data for back path in a chain
	USHORT t[20000];
	int n, nr;
	inline USHORT * GetBuf(){ return &t[n]; }
	inline void Addn(int ne){ n += ne; }
	void StoreChain(SCAND * tc, int ntc)	{
		t[n++] =(USHORT) ntc;
		for (int i = 0; i<ntc; i++)
			t[n++] = tc[i];
	}
	inline USHORT * StoreStart(){ nr = n; return &t[n]; }
	inline void StoreCancel(){ n = nr; }

};//expbuf;

struct STORE_UL{
	BF128 cells, one_digit_elims;// pattern and elims if one digit
	GINT64 ur2;// ur 2 cells equivalent (see UR)
	int type,digit_one;// keep it nx128 bits

	void Print(const char * lib);
};
struct BUG{// in wpaires, data for bug processing
	int cell,  el_par_ch[27];    // parity of digits for pairs in units
	int or_change, or_plus_tot, wchange;
	BF128 zz, wplus;
	BF16 r_c, c_c, b_c;		// row,col,box having plus
	int tplus[10], tplus_digits[10], change_plus[10],ntplus;  // same in table limited to 8
	USHORT r_cn[9], c_cn[9], b_cn[9], // count plus per unit
		np, aigpmax, aigun;
	int  Init();
	int  DoChange(int iplus,  int el);

};

struct WWUR2{// solving UR/UL 2 cells in a unit
	int cell1, cell2, digs, ul_plus, digitsothers, nothers, digits_ur,
    target, rbase, go_naked, wcell_nserate, degree,unit,
    wcells_nserate,
    dfree,tfree[10],nfree;
	int locdiag;
	char const * det_mess;
	BF128 wcells, cells_ur, cells_others_ur, cells_3
    ,w_b,// cells with extra digits no common digit
    w_c,// cells with only free digits 
    wnaked,wclean ;
	inline void Set(GINT64 & t,const char * lib){
		cell1 = t.u8[0], cell2 = t.u8[1], digs = t.u16[1], ul_plus = t.u8[4], nothers= t.u8[5],
			digitsothers=t.u16[3],rbase = 45 + ul_plus;
		det_mess = lib;
		locdiag=0;
	}
	inline void SetTarget(int etarget, int edegree){
		degree = edegree, target = etarget;
	}
	void Init( int eunit);
	void InitFreeDigits();
	int Hidden_pair();
	int Naked_pair();
	int Hidden_triplet();
	int Naked_triplet();
	int Naked_triplet21();
	int Hidden_quad();
	int Naked_quad();
	int Naked_quad31();
	int Naked_quad22();
	int Go_serate(GINT64 & t,int etarget);
	int Go_serate_unit(int eunit,int hdegree, int ndegree);
	int Go_serate_Fast(GINT64 & t);
	int Go_serate_unit_Fast(int eunit);


};
struct XSTATUS{// fish data for fish processing
	int active, digit,maxpas,npas;// maxpas=20 at start
	int dig_sets[27];
	BF32 bivsets;// unit with bivalue 27 bits field
	BF128 pm, // all candidates of the digit
		pmbiv, // cells belonging to one or more bi values
		elims; // potential eliminations (through brute force
	BF128 expon, expoff;
	void Init(int dig);
	void AddStart(int xcell, int unit, GINT64 * t, int & nt,BF128 * tx,BF128 & x);
	int R65Xexpand(int xcell1, int xcell2, int loop, int * back,BF128 & loopend);
	int XCycle(int fast);
	int XChain(int fast);
	int CleanLoop(int * i,int nt);
	int CleanChain(int * i, int nt,BF128 & clean);
	int Nishio1();
	int Nishio2();
	int XexpandNishio(int cell1);
};

struct STORE_XLC{
	int dig, rating,t[20],nt,loop;
	BF128 pattern;
};
struct YLSEARCH{
	int idig, xcell1, xcell2, maxpas, ncells, 
		ncells1, mode, c0, c1, c2, locdiag, diag;
	uint32_t d2;
	BF128 loop;
	GINT64 tback[30];
	int Search(int fast=0);
	int SearchOut(int fast );
	int Expand();
	int ExpandOut();
	int CleanLoop();
	int CleanLoopOut();
	int Is_To_Clean(int rating);
	void PrintTback();
};
struct XYSEARCH{
	struct PATH{
		GINT64  t[400];// cell,digit,source,step
		int nt, dig,cell;
	}paths[9];
	
	int idig, digit,  cell, ddig,dcell,
		nt,ntd,ntp, maxpas,maxrating, npaths,elim_done,
		nsteps, c1, c2, locdiag, diag, mode,fastmode,opprint;
	int dig_sets[9][27];
	uint32_t d2;
	BF128 pairs, cells_biv_all, cells_all,cells_biv_true,  dig_b_true[9];
	BF128 loop, wb;
	BF32 dig_bivsets[9],dig_sets3[9];
	GINT64 tback[300],t[400];// cell,digit,digit2,source
	PM3X used_off_digits, used_on_digits, cleang, cleanstart,
		active_all,active_unit,dbiv;
	GINT telims[50];	int ntelims ;
	//====== for dynamic and more
	int ind_pm[9][81],nind,is_contradiction,ntcands; // direct index for storing tables
	PM3X off_status[400],contradiction,elim_stored;
	GINT64 off_path[100][400];
	GINT tcands[400]; // candidates tables and index if belong to binary
	int tex[9][81]; // index to first "on" index in the path
	inline void SearchInit(int fast){
		fastmode = fast;
		elim_done = ntelims = 0;
		maxrating = 200;
		elim_stored.SetAll_0();
	}
	inline void Addt(int cell, int di, int source){
		GINT64 & tx = t[nt++];
		tx.u16[0] = (uint16_t)cell; tx.u16[1] = (uint16_t)di;
		tx.u16[2] = (uint16_t)source; tx.u16[3] = (uint16_t)nsteps;
	}
	inline void AddLastInCell(int cell, int di){
		GINT64 & tx = t[nt++];
		tx.u16[0] = (uint16_t)cell; tx.u16[1] = (uint16_t)di;
		tx.u16[2] = (uint16_t)(0x1000|cell); tx.u16[3] = (uint16_t)nsteps;
	}	
	inline void AddLastInUnit(int cell, int di,int unit){
		GINT64 & tx = t[nt++];
		tx.u16[0] = (uint16_t)cell; tx.u16[1] = (uint16_t)di;
		tx.u16[2] = (uint16_t)(0x2000|unit); tx.u16[3] = (uint16_t)nsteps;
	}
	void AddUnit(int unit, int source);
	void Init();
	void Init2();
	void InitCandidatesTable();
	int Search(int fast = 1);
	int SearchMulti(int fast);
	int MultiCell(int c0);
	int MultiCellExpand(int cell, int digit,PM3X & clean);
	int MultiUnit(int udigit, int unit);
	void OffToOn(int i);
	void OffToOn_Dyn(int i);
	void OnToOff(int i);
	void OnToOff_Dyn(int i);
	void StartMulti(int dig, int cell);
	void Expand_Multi(PM3X & cleanstart);
	int Expand_true_false();
	int CleanLoopOut();
	int Is_To_Clean(int rating);
	int CleanXYChain();
	int CleanXYLoop(GINT64 * t, int nt);
	int Do_Clean();
	int AddElim(int d, int c, int rat);
	int SearchDyn(int fast);
	void SearchDynPass(int nmax);
	void SearchDynPassMulti(int nmax);
	void ExpandDynamic(GINT cand);
	int ExpandDynamicToElim(GINT cand,GINT targ);
	int BackDynamic(GINT64 target, GINT64 * tb, int ntb);
	int BackDynamicOff(GINT target);
	void DynamicSolveContradiction(GINT cand,PM3X cont);
	void DynamicSolveContradiction(int dig1, int cell1, int dig2, int cell2, PM3X cont);
	void PrintTback();
	void PrintBackMulti(int elim_dig, int elim_cell);
	void PrintBackCom(const char * lib,GINT64 * ptback, int nback,int mode );
	void DebugT();
};



class BUILDSTRING{
public:
	char * p, *lim, *current_elim_start;
	char *zs;
	int mysize;
	// must be called before any use
	void SetUp(char * zse, int mysizee){ mysize = mysizee; zs = zse; lim = zs + mysize - 100; }
	void Init() { p = current_elim_start = zs; }
	inline void Ach(char x){ if (lim - p)		(*p++) = x; }
	inline void Achint(int x){ if (lim - p)(*p++) = (char)(x + '1'); }
	void Astr(const char * x);
	void Aunit(int x);
	//	inline void Arating(int x){if(lim-p)(*p++)=(char)(x/10+'0');(*p++)=(char)(x%10+'0');}
	//void Aint(int x){ char buf[20]; _itoa_s(x, buf, 20, 10); Astr(buf); } //MD 23.9.2018 commented out
	inline void Aendl(){ if (lim - p) (*p++) = 0; } // put end string end of line
	inline void Aspace(){ if (lim - p)(*p++) = ' '; }
	void AScand(SCAND x);
	void Close();
	void Store(UCAND & w);
	void StoreChain(SCAND * t, int n);
	void StoreX(SCAND & w, int dig);
	void StoreXChain(SCAND * t, int n, int dig);
	inline void SetCurrentElim(){ current_elim_start = p; }
	inline void ClearCurrent(){ p = current_elim_start; }
	void SetCurrentAsFirst(){
		int n = (int)(p - current_elim_start);
		p = zs;
		for (int i = 0; i<n; i++)
			*p++ = *current_elim_start++;
		current_elim_start = zs;
	}
};



class PM_GO;
struct TWO_DIGITS{// one of 36 possible 2 digits 
	BF128 bf2,bf2d;
	int free_units,digit[2],nhp;
	GINT64 thp[20];// band 2 cells ; i band ; 
	int Hiden_Pairs_Box();
	int Hiden_Pairs_Rows();

};
/*
*/

struct PM_DATA{ // keep it open for dynamic expansion process

	/* 
	PMBF pm,pm_diag;
	int free_unit_per_digit [9]; // 9 27 bits fields
	int unknown_count_per_unit [27];
	int locked_count_per_unit [27];
	int free_unit_for_search2,free_unit_for_search3,free_unit_for_search4;
	*/
	PM_GO * parent;
	int unknown_count;
	BF81 unknown;
	PMBF pm;
	BF16 dig_cells[81],emptyrows,emptycolumns;
	RBF27 unassigned_sets, changed_sets;

	BF16 dig_reg[9][27];
	BF81 dig_reg81[9][27];
	char res_zero_based[81], // 0 to 8
         puz_zero_based[81],
		 start_puz[82];
	int	 contradiction_in_set;
	// all above updated live at each assignment
	void AssignStandard(int dig,int cell);
	void CleanAll81(int dig,int cell);
	void CleanDyn(int dig,int cell);
	void Clean(BF16 & digs,USHORT cell);
	void Clean(USHORT dig,BF81 & cells);
	//
	void InitialAndFirstAssignmentx();
	int  Std_SetDigRegsAndSinglesx(int boxonly=0);  //12
	void DoInit_zhou_to_me(); 
	void Do_Zhou_Assign(int dig,int cell);
	void Do_Zhou_Row_Transfer(int dig,int irow,int vrow);
	void DoEndInitAfterZhou();
	int DoZhouCycle();
	int  Std_HiddenSingles(int boxonly);          //15
	int  Std_SinglesInCells();
	int  Std_LastCellInUnit();  //10

	inline int IsAssigned(int digit,int unit)
		{return unassigned_sets.t[digit].Off(unit);}
	int Set(SCAND x); // return 0,1, contradiction_in_set if conflict
	int IsGameLocked(PM_DATA & pmst,PMBF & bf);
	int SetsNowAssigned(RCAND * tass);
	int  FindNewOn(SCAND * tnewon,RCAND * tue,USHORT & nue);

	int  ApplyOn(UCAND x,SCAND * telims);
	void TrackError(char * lib);
	void GetStdCharPuzzleStatus(char * zs);
	void ImageUn_myd(int ch);
	void ImageFloors_myd(BF16 ff);
	void ImageDump(char * lib);
	void ImageCandidats();
	void DynSet(PMBFONOFF & offon,PMBFONOFF & offonold); // create the dyn status for dynamic plus
	void PrepareSubGrids();
};


class PM_GO{
public:

	struct HINT {  // combine all hints having the same "lowest"rating
		PM_GO * parent;
		PM3X pmelims;

		USHORT rating,rating_done;

		inline void Set_Target(int r){ rating = (USHORT)r; rating_done = 0; }
		inline int IsToDo(int r){ return (r <= rating); }
		inline void Done(int r){ if (r > rating_done) rating_done = (USHORT)r; }


		inline void Init() {rating=999;pmelims.SetAll_0();}

		int ChainLengthAdjusted(int base,int length);
		int MaxLengthForRating(int newbase);
		void Add(PM3X & elime, USHORT rate);
		int AddCand(USHORT dig,USHORT cell,USHORT rate);
	}hint;  
	class XYCOM;
	enum Ex_End{
		ex_end_nothing=0,  // process assumes 0 for nothing
		ex_end_contradiction,
		ex_end_loop,
		ex_end_maxpas
	};
	enum Explain_mode{ex_lastinregion,ex_lastincell,ex_pointing,
		ex_NP,ex_HP,ex_xwr,ex_xwc,
		ex_Ntrip,ex_Htrip,ex_swr,ex_swc,ex_chains,ex_UR,ex_exo,
		ex_symg,ex_aahs,ex_kite,ex_fish
	};


	class XYCOM{ // expansion common to base and nested levels
	public:
		PM_GO * parent;
		PM_DATA * tpmd, // pointers to main and dynam expansion areas 
			pmddyncycle; // tpmd[1] at the start of a new cycle

		BF16 * dr,* drcells;// pointers to cells and regions in use
		BF81  *dr81;

		SCAND wc			; //current cand in progress  cell 7 + digit 4 +off 1
		USHORT  wcoff,wcdig,wcell;
		// data used in expansion processes 
		PMBFONOFF onoff,onoff_back,onoffold,lastonoff, usednn;
		PMBF biv_map,ybiv_map,cont, //expon,expoff,cont,
			store_expoff[320]; //keep Init() results for multi chains

		SCAND xtex[2*16*128 ],tcand[640],tback[150],tempty[9],	wtempty[9],
			ttback[9][400]	;  // used inGet MAICS
		UCAND telims[200], store_list[320],store_index[16*128],
			tret1[400],ntret1,tret2[400],ntret2,itret;
		RBF27 biv_sets,touched_sets_assigned,seen_set_empty;
		BF81  biv_cells,touched_cells,touched_assigned,seen_cell_empty;
		USHORT  ispluscycle,plusmode,firstnested,
				nelims,npas,indtcand[50],ntcand,ntcandr,maxpas,
				diag,diagplus,nstore,rating,nempty,wnempty,
				plusnested,myplusnested,minlength;
		Ex_End xend;

		BUILDSTRING buildstring;
		USHORT * texplain [2000],ntexplain,nttback[9];
		EXPLAIN_BUFFER expbuf;

		int ExplainRegion(int dest,int unit ,int dig);
		int ExplainCell(int dest,int cell);

		int LocateStep (SCAND x);
		int LocateStepInd (int x);

		void FishStoreOff(BF81 & wcum,SCAND wsc);
		int FalseAddBack(SCAND sc);
		void NewExpand(SCAND x); // start a new on/off to expand

		int DynLiveAdd(SCAND sc);
		void DynLiveAddBiv(SCAND sc,SCAND xorr);

		void XYexpandDynamicStart(UCAND cand1,PM_DATA & myd);
		void DynLastInSet();
		void DynLastInCell();
		void DynOnToOff(PM_DATA & myd);
		void DiagTcandTexplain(int mainprocess=0);
		void DynCycleInit();
		int DynPlusPointingClaiming();
		void DynPlusCycle();
		int DynBackSkfrCom(UCAND x,USHORT * tret,PM_DATA & myd);

	};
	
/*
	class  NESTED:public XYCOM{ // work area to search nested chains
		//				chaines	" ~x et x"	multi	dynami
		//95	level2		x			
		//105	level 3		x			x??		x	
		//110	level 4		x			x		x		x
		int nnstart,maxdynlength;
	public:
		void InitNested();
		void NewExpand(SCAND x); // start a new on/off to expand
		void All_Biv_AICs(int mode);
		void DynamicChains();
		int XYexpandDynamic(UCAND cand1);
		int XYexpandDynamicToTarget(UCAND cand1,SCAND target);
		void DynMultiChains(REG_LIST &myrl);
		void XYexpandDynamicWhile();
		int DynBackSkfr(UCAND x,USHORT * tret){return DynBackSkfrCom(x, tret,pmddyncycle);}
	}nested;
	struct SCEN_VIRUS{// store a scenario for SK loop or virus chain
//		STEP_STATUS ss;
//		BFCAND cands;
		PMBFONOFF pof;
		USHORT pairs[8] , scen[8];
		void UcandAdd(UCAND cd,BF81 * c);
		void PairAdd(int dig,USHORT * tcells,BF8 cellspat,PM_GO * parent);

		inline void PairAdd(int d1,int d2,USHORT * tcells,BF8 cellspat,PM_GO * parent){
			PairAdd(d1,tcells,cellspat,parent); PairAdd(d2,tcells,cellspat, parent); 
		}
	};
	struct AAHS_AC2{ // 2 cells four digits 
		USHORT cell1,cell2,row_col,tdigits[4];
		// canonical form
		BF8 pdigits[4], //  pattern  01 10 11 exlusively
			pairvalid,scen1valid,scen2valid;//pairs order is given in bf_fix.tp6
		PM_GO * parent;
		void ClearPair(int i);
		void Load( PM_GO * parent, BF81 & cells,USHORT * digits,USHORT unit);
		void ChainVirus(AAHS_AC2 *taahs,int naahs,int iaahs,SCEN_VIRUS & scvold);
//		void AddBelt(GO_SOLVER * gsv,STEP_STATUS *ssw);
		void AddPair(SCEN_VIRUS &scw,int pair);
		void AddOne(SCEN_VIRUS &scw,int pair,int dig);
		int IsBilocal(int pair);
	};
	struct VLOOP{      // collecting the sk loop pattern for later action
		PM_GO * parent;
		BF81 cells[13];      //  8 or 12 pairs of cells
		USHORT digits[13][2]; // linking digits(0;8) start unlinked first
		USHORT mainbeltscleared,seen,
			   units[13],        // and corresponding units
			   units_number;       // 8 or 13 could be less if virus chain
		AAHS_AC2 aahs[12];
		int FirstAction();
		int SecondAction(int do_it);
		int SecondBilocations();
		void GenNotValidEffect();
	}vloop;
	struct TWO_CELLS{  // store 2 cells mini row mini col 3 or 4 digits
		USHORT cell1,cell2,rel_box,row_col,tdigits[4],ndigits;
		BF16 digits16;
		// canonical form
		BF8 pdigits[4], //  pattern  01 10 11 exlusively
		pairvalid,scen1valid,scen2valid;//pairs order is given in bf_fix.tp6
		void LoadInit(BF81 * c,USHORT ecell1,USHORT ecell2,BF16 digits,
			          USHORT endigits,USHORT unit,USHORT relbox);
		void Canonical3(BF81 * c);
		void Canonical4(BF81 * c);
		void ClearPair(int i);
	};
	struct TWO_CELLS_TABLE{
		PM_GO * parent;
		TWO_CELLS			c2[100];
		USHORT nc2;
		int Find_2Cells() ;
	}t2cells;
	struct WDYN_EXO{// working area in dynamic processing for an exocet
	//	BF81 
		struct WDIG{
			USHORT new_in_base,new_in_target,new_on_in_target,
				old_in_base,old_in_target,
				is_in_base,is_in_target,
				on_in_base,on_in_target,
				cell_in_base,cell_in_target,
				cell_in_target1,cell_in_target2;
			BF81 zbase,ztarget,ztarget1,ztarget2,
				base81,t_81,t1_81,t2_81; // zbase is in the box
			BF8 base,target,pdigit ;
			void Init();
		}wdig[4];
	};
	struct EXOCET{ // what is collected in the first step
		PM_GO * parent;
		USHORT cells[8],type; // base and target  (1 or 2 cells per target)
		BF16 digits,abi_digits,exo_digits;
		USHORT units[3],cross_lines[3],ndigits,  bande1,bande2;//,altbox1,altbox2;
		// canonical form
		BF8 pattern,pdigits[4], //  pattern  1111 ou 11111 ou 111111
			pairvalid,scen1valid,scen2valid;//,pairmask[4]; // pair mask is 01 10 11
		USHORT tdigits[4],ncells,npairs,twin1,twin2;

		void Loadinit(PM_GO * parent,TWO_CELLS & source);
		int IsNotCrossPossible(PM_DATA & myd,BF81 & map,USHORT unit,BF16 & w16);
		void Canonical();
		void ClearPair(USHORT d1,USHORT d2);
		int FindEliminations();
		//void DynPrepare(PM_DATA & myd,PM_GO::XYSEARCH & xys,WDYN_EXO & wdex); //MD 23.9.2018: seel next line
		void DynPrepare(PM_DATA & myd,XYSEARCH & xys,WDYN_EXO & wdex);
	};
	struct EXO8: public EXOCET{
	};
	struct R0SEARCH{  
		struct R0S_LOT4{// a 3 rows or 3 colums base for cells
			BF81 c81;
			BF16 rcbm;
			USHORT rc4[5],nrc4;  // keep it open for 5 rows or columns
		}rclot[84+126][2];
		PM_GO * parent;
		BFSETS added_truths;
		BF16 floors;
		BF81 cellsf,cellspuref,cellstruth,cells_sets,cells_link,cells_base,
			 subgrids[235][2];
		USHORT unit_sets[27],ntruth,nlinks,rank,code_ret,opt,
			   tu_sets[10],ntu_sets,tu_links[10],ntu_links,
			   telims[100],nelims,r0count,r1count,debug,ntr1,do_elims,
			   tcovers[50],ncovers,ntruths_cells;
		int diag;
		//BFCAND tr1elims[30],wr1elims;
		PMBF telim_native[30],welim_native;
		void Init(USHORT *tt,USHORT ntt);
		int TryRowCol(USHORT * tu,USHORT ntu,USHORT dtruth);
		int TryRowColBand(USHORT * tu,USHORT ntu,USHORT dtruth);
		int CheckBand(USHORT * tu,USHORT ntu);
		int CheckElims();
		int CheckElims(BF16 ordigs,R0S_LOT4 & lr,R0S_LOT4 & lc);
		int GoForRowCol(USHORT * tu,USHORT ntu);
		int GoForX(USHORT row1,USHORT row2,USHORT col1,USHORT col2);
		int GoForBox4();
		void Add(BFSETS & ss,USHORT unit);
		void GoForCellsPrepare();
		int GoForCells();
		int GoForCells63(); // option no given
		int GoForCells(R0S_LOT4 & lr,R0S_LOT4 & lc);
		void GoForCellsRC(R0S_LOT4 & lx,BF16 digpat,int moderc);
		int TryUsingElims(BF16 floors,BF81 * elims_floors, BF81 & elims_cells);

	}r0search;	
*/
	//          start of PM_GO data and functions
	EXPLAIN_BUFFER expbuf;// locate the source in a path
	USHORT * texplain[1000];// start of each "event"
	//PM_DATA zpmd[3], mydxx, *myd_at_start, dyndx;
	//char	* mybits;	//c = "....|...._....|...._....|....-....|
	//PMBF pmelims,pair_biv;
	BF81 * c, 
		yusedpaires;
	BF16 active_floors,mfloors; 
	ONE_FLOOR one_floor;
	ACTIVERCB activercb;
//	EXOCET	texocet[30],	wexo; // 30 in open band mode
//	EXO8 texo8[5],wexo8;
//	PMBFONOFF pof_store[1000];  // use in vloop expansion

	// settings for builstrings in xysearch and nested
    #define GINTBUFSIZE1 2000
	#define buildstringsize1 200000
	#define buildstringsize2 20000
	char builstr1[buildstringsize1],builstr2[buildstringsize2];
	GINT gintbuffer[GINTBUFSIZE1];
	GINTBUF gintbuf;
	BUILDSTRING * bdsp[2];  //pointers to builstring main and nested
	GG		gg,gr,        // puzzle normalized 
			gsolution;   // final result if valid (one solution)

	int opprint, opprint2, stop_rating,
		ntr0logic,r0logictype,logictype,
		rank0_min, rank0_max,//usually  2/5, can be adjusted
		assigned,cycle,
	   ur_serate_mode,quick_mode,ratlim,
	   find_common_known, // limit of known before start of the search
		rat_er,rat_ep,rat_ed,
		nexocet,npof_store;
	BF32	bits_tasks_done;
	// data for kites
	int nbiv,nempty;
	uint32_t ratfound[17]; // set to 0 in the constructor
	const char * det_mess;

	GINT64 tur[20];	STORE_UL tul[10];	WWUR2 wwur2; BF128 lastul; int ntur, ntul;//==== UR UL handling
	BUG bug; // bug handling
	//===================== fish handling
	int active_digits; //digits with potential eliminations (brute force effect
	XSTATUS xstatus[9];	STORE_XLC store_xlc[20]; int nstore_xlc;
	YLSEARCH ylsearch, store_yl[20];	 int nstore_yl;
	XYSEARCH xysearch;
//============= data for symmetry of given processing
// 0 D1  1 D2  2 stick  3 central  4 R90
	USHORT sgiven_ok ,sgiven_paires[5][9],sgiven_singles[5][3];
	BF16 sgiven_singlesBM[5];


	//========================
	PM_GO();

	//int Assign(int digit,int cell,char * lib);
	//int Assign(int digit,BF81 &cells,char * lib);
	int CleanOr(int d1, int c1, int d2, int c2);
	//void Start();

	//int Solve_All_Singles(char * ze);
	//int  SolveFinalFilter(char * puz);

	//int Solve();
	int SolveGetLow61();// internal call valid puzzle 
	int SolveGetLow44(int pack=0);// internal call valid puzzle 
	void SolveSerate110();
	//void SolveSerate111();
	//void Solve199test();

	void Quickrate(int x) ;
	void Status(const char * lib, int option);
	int Rate10(); int Rate12();	int Rate15(); int Rate17(); 
	int Rate20(); int Rate23(); int Rate25(); int Rate26();
	int Rate28(); int Rate30(); int Rate32(); int Rate34();
	//void SetupActiveDigits();
	int Rate36(); int Rate38(); int Rate40(); int Rate42(); int Rate44();
	int Rate45_52();  int Rate45_52_Fast();
	int Rate45Plus(GINT64 * t, int  nt, int plus);
	int Rate45_el(GINT64 & t, int unit,int plus);// serate mode min 2 cells within unit
	//int Rate_wwur2_serate(int unit);
	int Rate2cellsGo(GINT64 & t);// 2 cells bivalues UR UL
	int Rate45_2cells(GINT64 * t, int & nt);// serate mode min 2 cells bivalues, not diagonal
	int Rate45_URs(GINT64 * t,int & nt);// serate mode min 2 cells bivalues, not diagonal
	int Rate46_Find_ULs();
	int RateUL_base(STORE_UL & s);
	int Rate_ULs(int plus45);
	int Rate52(); int Rate54();
	int Rate56(); int Rate56BUG();
	int Rate62(); int Rate62_APE();
	void SetupActiveDigits();	
	inline void XStatusPrepare(){		for (int i = 0; i < 9; i++)xstatus[i].Init(i);	}
	int Rate65Xcycle(int fast);
	int Rate6xXcycle(int rating);
	int Rate66Xchain(int fast);
	int Rate67_70(); int R67_70(int rating);
	int Rate70_75(int rat_ed);
	int Rate75(); int Rate75_ATE();// aligned triplet exclusion
	int Rate76Nishio(int fast);
	int Rate80Multi(int fast);
	int Rate85Dynamic(int fast);


	int Next10_28();
	int Next28();
	int Next30_44();

	//int NextElim();
	//int NextSolve2();
	//int NextSolve3();
	//int NextSolve3Sym();
	//int NextSolve4();
	//int NextElimQuick90();
	//int NextElimQuick90fast();
	//int NextElim90Solve();
	//void Quick_Split(GG & gg);// called by sk_gsplit

	//void PrepareSet();
	//void PrepareSetDyn(XYCOM * xyc);   // same for dynamic cycle set and fish
	//int DoActiveRCBx(PMBF * doit); // obsolete
	//int WWings();
 

	//int Kites(int dynamic=0);
	//void KitesDynamicGo();
	//void FishDynamicOthers();

	//void TraiteCheck();
	//void LookForBackdoors();
	//void AddSingles();
	//void TraiteCallFinalFilter();
	//void SolveFilter();
	//int CheckIfAssignedR1Logic(PMBF & fn);
	//void GoForRank1Logic(int print);
	//void UseRank1Logic();

	//void Traite_Find_Common();
	//void Traite_Find_Multi_Fish();
	//int FindMultiFish2();
	//void Traite_Find_SKLoop();
	//int Locate_VLoop();
	//int Is_VLoop();
	//int ReductVLoop();
	//void AddChainVirus(SCEN_VIRUS & scv);
	//int CombineStillValidScenario(char * lib);
	//void AAHS_Dynamic(AAHS_AC2 & aahs);
	//void AAHS_Dynamic(AAHS_AC2 & aahs,int ip,int id1,int id2);
	//void Traite_Find_Exocet();
	//void  FindAllJExocets();
	//void  FindAllJExocets2Cells(int option,int ic2);	

	//void  FindBandExocets(int ic2);	
	//void  FindAllExocets(int ic2);	
	//int   BuildPairTable18(USHORT * tt,int cell);
	//void  FindExocets18(int ic2);
	//void  FindExocets183(int ic2);
	//int  Exocet18Det();
	//	void  Exocet183Det();
	//void FindExpandExocet();
	//void ExoAnalysis();
	//void Traite_Find_Conjugated();



	//void Traite_FindR0N();
	//int Traite_FindR0N_Go();
	//void R0Analysis();
	//int R0Cells();
	//int Apply_Rank0Direct();



	//char IsDoubleExocet(EXOCET & wi, EXOCET & wj);
//	int Is_JExocet();
	//int Is_FullJExocet2Cells();
	//int Is_AbiLoop ();
	//int Is_AbiLoop(USHORT * abi_pairs);
	//int Expand_Exocet();
	//void Exocet_Dynamic(); // dynamic cycle exocet effect

	//int Do_Exo_multirank1D(BF81 &x,BF81 &y,USHORT ntruths,USHORT links0);
	//int Do_Exo_multirank1();
 


	//USHORT IsConjugatedAAHS();
	//int IsConjugatedDigit(int dig);
	//int ClearFixGiven(int i);
	//int Do_First_Sym_Given();
	//void SymGiven_Dynamic(); // dynamic cycle sym given effect
	//void Given_DynamicPair(SCAND a,SCAND b);
	//void Dynam_One_Floor(PMBF & mypm);

	//============= debugging
	//void ImagePoints(BF81 & w);
	//void Image(PMBF & pm,char * lib);
	//void Image(PMBFONOFF & pof,char * lib);
	//void Image(UCAND * t,USHORT n, char * lib);
	//void ImagePaires(EXOCET & wexo);
	//void ImageExtended(EXOCET & wexo);
	//void Xdebug(UCAND x,char * lib);
	//void XYdebug(UCAND x,char * lib);
	//void ListImage(SCAND * t,int n);
	//void ListImageX(SCAND * t,int n);
	//inline void Candidats(PM_DATA & pmd){pmd.ImageCandidats();}
	//void ImagePuzzle(char * puz,char * lib);
	//void ImageSolution();
	};
