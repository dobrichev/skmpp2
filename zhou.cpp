#include "sk_t.h"
#include "solver_step.h"
#include "Zhn.h"
#include "Zhtables_cpp.h" // also describes the pattern

extern PM_GO pm_go;
ZHOU zhou[50]; // must host main brute force plus minimality analysis and recursive generation
ZHOU zhou_i,zhou_ip,//zhou_i================== initial for a GAME
    zhou_solve;// basis to solve a puzzle using elimination logic


void ZH_GLOBAL::Pm_Status(ZHOU * z){
	zhou_current = z;
	oldcount = z->CountDigs();
	pmelims.SetAll_0();
	cells_unsolved_e = z->cells_unsolved;
	cells_unsolved_diag.Diag3x27(cells_unsolved_e);    
	for (int idig = 0; idig < 9; idig++){
		BF128 w = z->FD[idig][0];
		pm.pmdig[idig] = w& cells_unsolved_e;
		pmdiag.pmdig[idig].Diag3x27(pm.pmdig[idig]);
		unsolved_r_count[idig] = _popcnt32(w.bf.u32[3]);
		//unsolved_c_count[idig] = __popcnt(pmdiag.pmdig[idig].bf.u32[3]);  idiot
		int row = 0;
		for (int iband = 0; iband < 3; iband++){// fill rows cols 9 bits
			register int band = pm.pmdig[idig].bf.u32[iband],
				bandd = pmdiag.pmdig[idig].bf.u32[iband];
			for (int irow = 0; irow < 3; irow++, row++,band>>=9,bandd>>=9){
				dig_rows[idig][row]= band & 0777;
				dig_cols[idig][row] = bandd & 0777;
			}
		}
	}
}
void ZH_GLOBAL::Pm_Status_End(){// prepare boxes and cells status
	memset(dig_cells, 0, sizeof dig_cells);
	memset(cells_count, 0, sizeof cells_count);
	for (int idig = 0; idig < 9; idig++){
		int box = 0;
		for (int iband = 0; iband < 3; iband++){// fill rows cols 9 bits
			register int band = pm.pmdig[idig].bf.u32[iband];
			for (int ibox = 0; ibox < 3; ibox++,box++ ,band>>=3){
				register int  R = band & 7;// band shifted in the low bits
				R |= (band >> 6) & 070;
				R |= (band >> 12) & 0700;
				dig_boxes[idig][box] = R;
			}
		}
	}
	BF128 wu = cells_unsolved_e;
	int xcell;
	while ((xcell = wu.getFirst128()) >= 0){
		wu.clearBit(xcell);
		int cell = From_128_To_81[xcell];
		int & cc = cells_count[cell], &cdigs = dig_cells[cell];
		for (int idig = 0; idig < 9; idig++){
			if (pm.pmdig[idig].Off(xcell))continue;
			cc++;
			cdigs |= 1 << idig;
		}
	}
}
void ZH_GLOBAL::AddSingle(int band,  int vband){
	if (!vband) return;
	uint32_t cell;
	bitscanforward(cell, vband);
	int xcell = cell + 32 * band;
	if (cells_assigned.On(xcell)) return;
	tsingles[nsingles++] = (int)cell + 27 * band;
	cells_assigned.setBit(xcell);
}
void ZH_GLOBAL::AddSingleDiag(int band, int vband){
	if (!vband) return;
	uint32_t celld;
	bitscanforward(celld, vband);
	int cell = C_transpose_d[celld + 27 * band],
		xcell = C_To128[cell];
	if (cells_assigned.On(xcell)) return;
	tsingles[nsingles++] = cell;
	cells_assigned.setBit(xcell);
}
void ZH_GLOBAL::Build_digits_cells_pair_bf(){
	for (int i = 0; i < 9; i++)		digits_cells_pair_bf[i] = pairs & pm.pmdig[i];
}
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
	memset(count, 0, sizeof count);
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
		__movsd((uint32_t*)count, (uint32_t *) &count[3], 3);
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
	__movsd((uint32_t*)x3_cmap, (uint32_t *)&sourcei[27 * ib1], 27);
	__movsd((uint32_t*)&x3_cmap[27], (uint32_t *)&sourcei[27 * ib2], 27);
	__movsd((uint32_t*)&x3_cmap[54], (uint32_t *)&sourcei[27 * ib3], 27);

}
//td is 8bits cell + 8 bits digits
void ZH_GLOBAL::Morph_digits(int morph){// using the given entry
	ngiven = 0;
	uint32_t count[9];
	memset(count,0,sizeof count);
	for (int i = 0; i < 81; i++){
		register int c = puz[i];
		if (c<'1' || c>'9') continue;
		c -= '1';
		count[c]++;
		tgiven[ngiven++].u16 =(uint16_t)( i | (c << 8));
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
			tgiven[ic].u8[1] = (uint8_t) x3_dmap_inv[tgiven[ic].u8[1]];
	}
}
//void ZH_GLOBAL::Map_Morph_digits(GINT16 * td, int nc){//applying x3_cmap to the cells
//}
void ZH_GLOBAL::Morph_digits(GINT16 * td, int nc){// must be in line with the morphed pattern
	ngiven = 0;
	uint32_t count[9];
	memset(count, 0, sizeof count);
	for (int it = 0; it < nc; it++){
		register int  dig = td[it].u8[1];
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
		tgiven[ic].u8[1] =(uint8_t) x3_dmap_inv[tgiven[ic].u8[1]];
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
		tgiven[ngiven++].u16 =(uint16_t)( i | (c << 8));
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
		tgiven[ngiven++].u16 = (uint16_t)(i | (c << 8));
	}
	if (_popcnt32(digs) < 8) return 1; // don't accept less than 8 digits given
	if( InitSudoku()) return 1;
	zhou_solve = zhou[0];
	zhou[0].ComputeNext();
	if (nsol != 1) return 1;
	for (int i = 0; i < 81; i++)zerobased_sol[i] = (char)(stdfirstsol[i] - '1');
	memset(locked_nacked_brc_done, 0, sizeof locked_nacked_brc_done);
	memset(row_col_x2, 0, sizeof row_col_x2);
	return 0;
}

int ZH_GLOBAL::Go_InitSolve(GINT16 * td, int nc){
	zsol = stdfirstsol;
	InitCount(1);
	strcpy(puz, empty_puzzle);
	NoMorph();
	ngiven = nc;
	for (int i = 0; i < nc; i++){
		puz[td[i].u8[0]] = (char)(td[i].u8[1] + '1');
		tgiven[i] = td[i];
	}
	if (InitSudoku()) return 1;
	zhou_solve = zhou[0];
	zhou[0].ComputeNext();
	if (nsol != 1) return 1;
	for (int i = 0; i < 81; i++)zerobased_sol[i] =(char)( stdfirstsol[i] - '1');
	memset(locked_nacked_brc_done, 0, sizeof locked_nacked_brc_done);
	memset(row_col_x2, 0, sizeof row_col_x2);
	return 0;
}
void ZH_GLOBAL::ValidPuzzle(ZHOU * z){
	if (zsol && (!nsol)){// store the first solution
		z->SetKnown(zsol);
		for (int i = 0; i < 9; i++)digit_sol[i] = z->FD[i][0];
	}
	nsol++;
}

void ZHOU::AssignSolver(int rating){
	if (!zh_g.nsingles) return;
	for (int i = 0; i < zh_g.nsingles; i++){
		if (rating)cout << "assign cell " << cellsFixedData[zh_g.tsingles[i]].pt
			<< " for digit " <<(int)( zh_g.zerobased_sol[zh_g.tsingles[i]]+1) << endl;
		Setcell(zh_g.tsingles[i]);
	}
	if (rating){
		char zs[82];
		SetKnown(zs);
		cout << zs << " new known status after assign rating  "<<rating << endl;
	}

}
void ZHOU::Naked_Pairs_Seen(){
	memset(zh_g.locked_nacked_brc_seen, 0, sizeof zh_g.locked_nacked_brc_seen);
	BF128 & pa = zh_g.pairs, fd[9];
	if (pa.Count() < 2) return;
	int td[9], nd = 0;
	// find digits with >= 2 biv cells
	for (int i = 0; i < 9; i++){
		fd[i] = FD[i][0] & pa;
		if (fd[i].Count()>1) td[nd++] = i;
	}
	if (nd < 2)return;
	// try each stored digits pair
	for (int i1 = 0; i1 < nd - 1; i1++){
		for (int i2 = i1+1; i2 < nd; i2++){
			BF128 w = fd[td[i1]] & fd[td[i2]];
			if (w.Count() < 2)continue;
			int xcell1,xcell2;
			while ((xcell1 = w.getFirst128()) >= 0){
				w.clearBit(xcell1);
				BF128 w2 = w;
				CELL_FIX cf1 = cellsFixedData[From_128_To_81[xcell1]];
				while ((xcell2 = w2.getFirst128()) >= 0){
					w2.clearBit(xcell2);
					CELL_FIX cf2 = cellsFixedData[From_128_To_81[xcell2]];
					BF128 wor; wor.SetAll_0();
					wor.setBit(xcell1); wor.setBit(xcell2);
					//cout << xcell1 << ";" << xcell2 << "\t"
						//<< From_128_To_81[xcell1] << ";" << From_128_To_81[xcell1]<<"\t"
						//<< cf1.eb << ";" << cf2.eb << endl;
					if (cf1.eb == cf2.eb){//same box
						//cout << "same box" << endl;
						zh_g.locked_nacked_brc_seen[0] |= wor;
					}
					if (cf1.el == cf2.el){
						zh_g.locked_nacked_brc_seen[1] |= wor;
					}
					if (cf1.pl == cf2.pl){
						wor.SetAll_0();
						wor.SetDiagX(xcell1); wor.SetDiagX(xcell2);
						zh_g.locked_nacked_brc_seen[2] |= wor;
					}
				}
			}
		}
	}
}
void ZHOU::XW_template(int idig){
	int rx2 = zh_g.row_col_x2[idig][0]/*, cx2 = zh_g.row_col_x2[idig][1]*/;
	int tr[9], nr = 0, nfree = zh_g.unsolved_r_count[idig];
	if (nfree - _popcnt32(rx2) < 3) return;
	int *rows = zh_g.dig_rows[idig]/*, *cols = zh_g.dig_cols[idig]*/;
	// collect all rows/cols 2 cells
	for (int irow = 0; irow < 9; irow++){
		if (rx2 & (1 << irow))continue;
		if (_popcnt32(rows[irow]) == 2)tr[nr++] = irow;
	}
	if (nr < 2) return;
	for (int ir1 = 0; ir1 < nr - 1; ir1++){
		register int R = rows[tr[ir1]];
		for (int ir2 = ir1 + 1; ir2 < nr; ir2++){
			if (R != rows[tr[ir2]]) continue;
			// new XW lock it  and check if active
			//Find column mask
			// aply column and row mask
			//check if active and exit if active
			// exit if now <3
			// next ir1 anyway

		}

	}
}


/*
LastCell=10,				///< last cell in row column box
SingleBox=12,               ///< single in box
Single_R_C=15,              ///< single in row or column
Single_after_Locked=17,		///< locked in box  clearing row/col ?? giving a fix??
PointingClaiming=19,        ///< unknown
HiddenPair_single=20,		///< hidden pair, hidden fix
SingleInCell=23,			///< cell one candidate
HiddenTriplet_single=25,    ///< Hidden triplet, fix
Locked_box=26,				///< locked in box, no fix
Locked_RC=28,				///< locked in row/col  no fix
NakedPair=30,               ///< 2 cells containing 2 digits
XWing=32,                   ///< XWing
HiddenPair=34,              ///< 2 digits locked in 2 cell
Naked_triplet=36,           ///< 3 cells containing 3 digits
swordfish=38,               ///< swordfish
*/


int ZHOU::Rate10_LastInUnit(){// look for a single in box
	if (zh_g.diag)		cout << "last in unit" << endl;
	uint32_t hidden;
	register int R1 = 0, R2 = 0;
	for (int iband = 0; iband < 3; iband++){
		register int  band = cells_unsolved.bf.u32[iband];
		for (int ibox = 0; ibox < 3; ibox++){
			hidden = band & tband_box[ibox];
			if (_popcnt32(hidden) == 1)zh_g.AddSingle(iband,  hidden);
		}
		for (int irow = 0; irow < 3; irow++){
			hidden = band & tband_row[irow];
			if (_popcnt32(hidden) == 1)zh_g.AddSingle(iband, hidden);
		}
		register int R = band & 0777;	R2 |= R&R1; R1 |= R;
		R = (band >> 9) & 0777;	R2 |= R&R1; R1 |= R;
		R = (band >> 18) & 0777;	R2 |= R&R1; R1 |= R;
	}
	R1 &= ~R2;
	if (R1){// last in column low probability (at start)
		for (int icol = 0; icol < 9; icol++) if (R1 & (1 << icol)){
			int rcol = Zhoucol << icol;
			for (int iband = 0; iband < 3; iband++){
				int hidden = cells_unsolved.bf.u32[iband]&rcol;
				if (hidden) {
					zh_g.AddSingle(iband,  hidden);
					break;
				}
			}
		}
	}
	return zh_g.nsingles;
}
int ZHOU::Rate12_SingleBox(){// look for a single in box
	if (zh_g.diag)		cout << "single in box" << endl;
	uint32_t hidden;
	int idig;
	for (idig = 0; idig <9; idig++){// priority to high digits last done
		BF128 & fd = FD[idig][0];
		for (int iband = 0; iband < 3; iband++){
			int  band = fd.bf.u32[iband] & cells_unsolved.bf.u32[iband];
			for (int ibox = 0; ibox < 3; ibox++){
				hidden = band & tband_box[ibox];
				if (hidden && _popcnt32(hidden) == 1)
					zh_g.AddSingle(iband, hidden);
			}
		}
	}
	return zh_g.nsingles;
}
int ZHOU::Rate15_SingleRow(){// look for a single in Rox
	int trow[3] = { 0777, 0777000, 0777000000 };
	if (zh_g.diag)		cout << "single in row" << endl;
	uint32_t hidden;
	int idig;
	for (idig = 0; idig <9; idig++){// priority to high digits last done
		BF128 & fd = FD[idig][0];
		for (int iband = 0; iband < 3; iband++){
			int  band = fd.bf.u32[iband] & cells_unsolved.bf.u32[iband];
			for (int irow = 0; irow < 3; irow++){
				hidden = band & trow[irow];
				if (hidden && _popcnt32(hidden) == 1)
					zh_g.AddSingle(iband,  hidden);
			}
		}
	}
	return zh_g.nsingles;
}
int ZHOU::Rate15_SingleColumn(){// look for a single in column
	if (zh_g.diag)		cout << "single in column" << endl;
	uint32_t hidden;
	for (int idig = 0; idig <9; idig++){
		BF128  fd = FD[idig][0] & cells_unsolved;
		register int col = Zhoucol;
		for (int icol = 0; icol < 9; icol++,col<<=1){
			int nn = 0;
			for (int iband = 0; iband < 3; iband++){
				nn += _popcnt32(fd.bf.u32[iband]&col);
			}
			if (nn!= 1) continue;
			// redo it to catch the band
			for (int iband = 0; iband < 3; iband++){
				hidden = fd.bf.u32[iband] & col;
				if (hidden){
					zh_g.AddSingle(iband,  hidden);
					break;
				}
			}
		}
	}
	return zh_g.nsingles;
}
int ZHOU::Rate17_lockedBox_Assign(){
	if (zh_g.diag)		cout << "rate 17 locked box + assign" << endl;
	int tboxskrink[3] = { 0111, 0222, 0444 };
	int tcol[3];
	zh_g.nsingles = 0;
	for (int idig = 0; idig < 9; idig++){
		BF128  fd = FD[idig][0] & cells_unsolved;
		for (int iband = 0; iband < 3; iband++){
			register int band = fd.bf.u32[iband];
			tcol[iband] = (band | (band >> 9) | (band >> 18)) & 0x1FF;
			if (!band)continue;// solved band
			register int shrink = (TblShrinkMask[band & 0x1FF]) |
				(TblShrinkMask[(band >> 9) & 0x1FF] << 3) |
				(TblShrinkMask[(band >> 18) & 0x1FF] << 6);
			for (int irow1 = 0; irow1 < 3; irow1++){
				register int x = (shrink>>(3*irow1)) & 7 ;
				if (_popcnt32(x) - 1)continue;
				// this is a minirow locked in box
				//cout << "17 row dig=" << idig + 1 << " band=" << iband + 1
				//	<< " rrow=" << irow1+1 << endl;
				uint32_t rbox;// relative box
				bitscanforward(rbox, x);
				if (_popcnt32(shrink & tboxskrink[rbox]) == 1)continue;
				int nbox = ~tband_box[rbox];
				for (int irow2 = 0; irow2 < 3; irow2++) if(irow2-irow1){
					int hidden = band &tband_row[irow2] & nbox;
					if (_popcnt32(hidden) - 1)continue;
					zh_g.AddSingle(iband,  hidden);
					//cout << Char27out(hidden) << " add irow2=" << irow2 + 1 << endl;
				}
			}
		}
		// look now for columns 
		for (int iperm = 0; iperm < 3; iperm++){
			int * p = tperm3[iperm], b1 = tcol[p[0]], b2 = tcol[p[1]], b3 = tcol[p[2]],
				or23=b2|b3,and23=b2&b3,single23=or23^and23;
			for (int ibox = 0; ibox < 3; ibox++){
				int mask = 7 << (3 * ibox), b1m = b1&mask;
				if (_popcnt32(b1m) < 2) continue; // solved or cleaned
				int col23 = or23 & mask, colsingle = b1m & (~col23);
				if (!colsingle)continue;
				int single23box = single23&mask;
				if (!single23box)continue;
				//cout << "17 col dig=" << idig + 1 << " band=" << p[0]+1 << " box=" << ibox + 1 << endl;
				// can be 1 or 2 columns ?? test the 2 bands only one per band can come
				int singleb2 = b2&single23box;
				uint32_t icol;
				if (singleb2){
					bitscanforward(icol, singleb2);
					int colmask = 01001001 << icol;
					int bx = fd.bf.u32[p[1]] & tband_box[ibox] & colmask;
					if (_popcnt32(bx) == 1){
						zh_g.AddSingle(p[1], bx);
						//cout << Char27out(bx) << " add singleb2" << endl;
					}

				}
				int singleb3 = b3&single23box;
				if (singleb3){
					bitscanforward(icol, singleb3);
					int colmask = 01001001 << icol;
					int bx = fd.bf.u32[p[2]] & tband_box[ibox] & colmask;
					if (_popcnt32(bx) == 1){
						zh_g.AddSingle(p[2], bx);
						//cout << Char27out(bx) << " add singleb3" << endl;
					}
				}
			}
		}
	}
	return zh_g.nsingles;
}

void ZHOU::FindNakedPairsTriplets_NoSingle(){
	BF128 R1 = FD[0][0] & cells_unsolved,
		R = FD[1][0] & cells_unsolved,
		R2 = R&R1,R3,R4;
	R1 |= R;
	R = FD[2][0] & cells_unsolved;
	R3 = R2&R; R2 |= R1&R; R1 |= R;
	R = FD[3][0] & cells_unsolved;
	R4 = R3&R; R3 |= R2&R; R2 |= R1&R; R1 |= R;
	for (int i = 4; i < 9; i++){
		R = FD[i][0] & cells_unsolved;
		R4 |= R3&R; R3 |= R2&R; R2 |= R1&R; R1 |= R;
	}
	// if no single R1=R2
	zh_g.pairs = R2-R3;
	zh_g.triplets = R3 - R4;
}

int ZHOU::RateSingle(GINT64 * t, int nt){// look if a single appears after a hidden locked set
	if (!nt) return 0;
	for (int i = 0; i < nt; i++){
		GINT64 w = t[i];
		register int R = w.u32[1]; // cells pattern of the locked set
		int iband = w.u16[0], digs = w.u16[1];
		for (int idig = 0, dig = 1; idig < 9; dig <<= 1, idig++){
			if (dig &digs)continue;// skip digs of the hidden set
			BF128 fd = FD[idig][0] & cells_unsolved;
			int fdb = fd.bf.u32[iband];
			if (!(fdb&R)) continue;// no clearing  
			if (pm_go.opprint2 & 2)
				cout << Char27out(R) << " rate single try for digs 0" << oct << digs <<dec
                << " bande="<<iband+1<< " digit="<<idig+1<< endl;
			for (int iboxr = 0; iboxr < 3; iboxr++){// try a box
				int box = fdb & tband_box[iboxr];		if (!(box&R))continue;
				int hidden = box & (~R);
				if (_popcnt32(hidden) == 1)	zh_g.AddSingle(iband, hidden);
			}	
			if (pm_go.opprint2 & 1)cout << " nsingles after box " << zh_g.nsingles<<endl;
			for (int irow = 0; irow < 3; irow++){// try a row
				int row = fdb & tband_row[irow];		if (!(row&R))continue;
				int hidden = row & (~R);
				if (_popcnt32(hidden) == 1) zh_g.AddSingle(iband, hidden);
 			}
			if (pm_go.opprint2 & 1)cout << " nsingles after row " << zh_g.nsingles << endl;
			for (int icol = 0, col = Zhoucol; icol < 9; icol++, col <<= 1){// try a column
				int Rcol = R & col; if (!Rcol) continue;
				if (!(fdb & Rcol)) continue;
				int * p = tperm3[iband],		b0 = fdb & col & (~R),
                    b1 = fd.bf.u32[p[1]]&col,   b2 = fd.bf.u32[p[2]]&col;
				if ((_popcnt32(b0) + _popcnt32(b1) + _popcnt32(b2)) != 1) continue;
				if (b1)zh_g.AddSingle(p[1], b1);	
				else if (b2)zh_g.AddSingle(p[2], b2);
				else zh_g.AddSingle(iband, b0);
			}
			if (pm_go.opprint2 & 1)cout << " nsingles after col " << zh_g.nsingles << endl;
			//cout << "apr�s col zh_g.nsingles=" << zh_g.nsingles << endl;
		}
	}
	return zh_g.nsingles;
}
int ZHOU::RateSingleDiag(GINT64 * t, int nt){// look if a single appears after a hidden locked set in column
	if (!nt) return 0;
	for (int i = 0; i < nt; i++){
		GINT64 w = t[i];
		register int R = w.u32[1]; // cells pattern of the locked set
		int iband = w.u16[0], digs = w.u16[1];
		for (int idig = 0, dig = 1; idig < 9; dig <<= 1, idig++){
			if (dig &digs)continue;// skip digs of the hidden set
			BF128 fd = zh_g.pmdiag.pmdig[idig];
			int fdb = fd.bf.u32[iband];
			if (!(fdb&R)) continue;// no clearing    
			for (int iboxr = 0; iboxr < 3; iboxr++){// try a box
				int box = fdb & tband_box[iboxr];		if (!(box&R))continue;
				int hidden = box & (~R);
				if (_popcnt32(hidden) == 1)	zh_g.AddSingleDiag(iband, hidden);
			}
			for (int irow = 0; irow < 3; irow++){// try a row
				int row = fdb & tband_row[irow];		if (!(row&R))continue;
				int hidden = row & (~R);
				if (_popcnt32(hidden) == 1) zh_g.AddSingleDiag(iband, hidden);
			}
			for (int icol = 0, col = Zhoucol; icol < 9; icol++, col <<= 1){// try a column
				int Rcol = R & col; if (!Rcol) continue;
				if (!(fdb & Rcol)) continue;
				int * p = tperm3[iband],
					b0 = fdb & col & (~R),
					b1 = fd.bf.u32[p[1]] & col,
					b2 = fd.bf.u32[p[2]] & col;
				if ((_popcnt32(b0) + _popcnt32(b1) + _popcnt32(b2)) != 1) continue;
				if (b1)zh_g.AddSingleDiag(p[1], b1);
				else if (b2)zh_g.AddSingleDiag(p[2], b2);
				else zh_g.AddSingleDiag(iband, b0);
			}
			//cout << "apr�s col zh_g.nsingles=" << zh_g.nsingles << endl;
		}
	}
	return zh_g.nsingles;
}
int ZHOU::RateSingleBox(GINT64 * t, int nt){// look if a single appears after a hidden locked set
	if (!nt) return 0;
	for (int i = 0; i < nt; i++){
		GINT64 w = t[i];
		register int R = w.u32[1]; // cells pattern of the locked set
		int iband = w.u16[0], digs = w.u16[1];
		for (int idig = 0, dig = 1; idig < 9; dig <<= 1, idig++){
			if (dig &digs)continue;// skip digs of the hidden set
			BF128 fd = FD[idig][0] & cells_unsolved;
			int fdb = fd.bf.u32[iband];
			if (!(fdb&R)) continue;// no clearing  
			if (pm_go.opprint2 & 1)
				cout << Char27out(R) << "RateSingleBox  try for digs 0" << oct << digs << dec
				<< " bande=" << iband + 1 << " digit=" << idig + 1 << endl;
			for (int iboxr = 0; iboxr < 3; iboxr++){// try a box
				int box = fdb & tband_box[iboxr];		if (!(box&R))continue;
				int hidden = box & (~R);
				if (_popcnt32(hidden) == 1)	zh_g.AddSingle(iband, hidden);
			}
			if (pm_go.opprint2 & 1)cout << " nsingles after box " << zh_g.nsingles << endl;
		}
	}
	return zh_g.nsingles;
}
int ZHOU::RateSingleRow(GINT64 * t, int nt){// look if a single appears after a hidden locked set
	if (!nt) return 0;
	for (int i = 0; i < nt; i++){
		GINT64 w = t[i];
		register int R = w.u32[1]; // cells pattern of the locked set
		int iband = w.u16[0], digs = w.u16[1];
		for (int idig = 0, dig = 1; idig < 9; dig <<= 1, idig++){
			if (dig &digs)continue;// skip digs of the hidden set
			BF128 fd = FD[idig][0] & cells_unsolved;
			int fdb = fd.bf.u32[iband];
			if (!(fdb&R)) continue;// no clearing  
			if (pm_go.opprint2 & 1)
				cout << Char27out(R) << " RateSingleRow try for digs 0" << oct << digs << dec
				<< " bande=" << iband + 1 << " digit=" << idig + 1 << endl;
			for (int irow = 0; irow < 3; irow++){// try a row
				int row = fdb & tband_row[irow];		if (!(row&R))continue;
				int hidden = row & (~R);
				if (_popcnt32(hidden) == 1) zh_g.AddSingle(iband, hidden);
			}
			if (pm_go.opprint2 & 1)cout << " nsingles after row " << zh_g.nsingles << endl;
		}
	}
	return zh_g.nsingles;
}
int ZHOU::RateSingleCol(GINT64 * t, int nt){// look if a single appears after a hidden locked set in column
	if (!nt) return 0;
	for (int i = 0; i < nt; i++){
		GINT64 w = t[i];
		register int R = w.u32[1]; // cells pattern of the locked set
		int iband = w.u16[0], digs = w.u16[1];
		for (int idig = 0, dig = 1; idig < 9; dig <<= 1, idig++){
			if (dig &digs)continue;// skip digs of the hidden set
			BF128 fd = zh_g.pmdiag.pmdig[idig];
			int fdb = fd.bf.u32[iband];
			if (!(fdb&R)) continue;// no clearing    
			if (pm_go.opprint2 & 1)
				cout << Char27out(R) << " RateSingleCol try for digs 0" << oct << digs << dec
				<< " bande=" << iband + 1 << " digit=" << idig + 1 << endl;
			for (int irow = 0; irow < 3; irow++){// try a row
				int row = fdb & tband_row[irow];		if (!(row&R))continue;
				int hidden = row & (~R);
				if (_popcnt32(hidden) == 1) zh_g.AddSingleDiag(iband, hidden);
			}
		}
		if (pm_go.opprint2 & 1)cout << "RateSingleCol nsingles after cols " << zh_g.nsingles << endl;
	}
	return zh_g.nsingles;
}


int ZHOU::CollectHiddenPairsBox(GINT64* tp, int & np, int lim){
	np = 0;
	BF128 wfree = cells_unsolved - zh_g.locked_nacked_brc_done[0];
	for (int iband = 0, dband = 0; iband < 3; iband++, dband += 27){
		int bandpairs = zh_g.locked_nacked_brc_seen[0].bf.u32[iband],
			bandfree = wfree.bf.u32[iband];
		for (int iboxr = 0, mask = 07007007; iboxr < 3; iboxr++,  mask <<= 3){
			if ((int)_popcnt32(bandfree & mask) < lim) continue;
			int tbiv[9], digbiv[9], nbiv = 0;
			for (int idig = 0; idig < 9; idig++){// collect digits bi values in box
				int hidden = mask & FD[idig][0].bf.u32[iband];
				if (_popcnt32(hidden) == 2){
					digbiv[nbiv] = idig;
					tbiv[nbiv++] = hidden;
				}
			}
			if (nbiv < 2) continue; // no bi-value in the box
			for (int i = 0; i < nbiv - 1; i++){
				register int R = tbiv[i];
				for (int j = i + 1; j < nbiv; j++) if (R == tbiv[j]){// a pair seen
					if ((bandpairs & R) == R)						continue; //no effect
					tp[np].u32[1] = R;
					tp[np].u16[0] = (uint16_t)iband;
					tp[np++].u16[1] = (uint16_t)((1 << digbiv[i]) | (1 << digbiv[j]));
				}
			}
		}
	}
	return np;
}
int ZHOU::CollectHiddenPairsRow(GINT64* tp, int & np, int lim){
	np = 0;
	BF128 wfree = cells_unsolved - zh_g.locked_nacked_brc_done[1];
	for (int iband = 0, dband = 0/*, ibox = 0*/; iband < 3; iband++, dband += 27){
		int bandpairs = zh_g.locked_nacked_brc_seen[1].bf.u32[iband],
			bandfree = wfree.bf.u32[iband];
		for (int irowr = 0, mask = 0777; irowr < 3; irowr++, mask <<= 9){
			if ((int)_popcnt32(bandfree & mask) < lim) continue;
			int tbiv[9], digbiv[9], nbiv = 0;
			for (int idig = 0; idig < 9; idig++){// collect digits bi values in row
				int hidden = mask & FD[idig][0].bf.u32[iband];
				if (_popcnt32(hidden) == 2){
					digbiv[nbiv] = idig;
					tbiv[nbiv++] = hidden;
				}
			}
			if (nbiv < 2) continue; // no bi-value in the row
			for (int i = 0; i < nbiv - 1; i++){
				register int R = tbiv[i];
				for (int j = i + 1; j < nbiv; j++) if (R == tbiv[j]){// a pair seen
					if ((bandpairs & R) == R)						continue; //no effect
					tp[np].u32[1] = R;
					tp[np].u16[0] = (uint16_t)iband;
					tp[np++].u16[1] = (uint16_t)((1 << digbiv[i]) | (1 << digbiv[j]));
				}
			}
		}
	}
	return np;
}
int ZHOU::CollectHiddenPairsCol(GINT64* tp, int & np, int lim){
	np = 0;
	//int tboxskrink[3] = { 0111, 0222, 0444 };
	BF128 wfree = zh_g.cells_unsolved_diag - zh_g.locked_nacked_brc_done[2];
	for (int iband = 0, dband = 0/*, ibox = 0*/; iband < 3; iband++, dband += 27){
        int bandpairs = zh_g.locked_nacked_brc_seen[2].bf.u32[iband], //digonal mode
			bandfree = wfree.bf.u32[iband];
		for (int irowr = 0, mask = 0777; irowr < 3; irowr++, mask <<= 9){
			if ((int)_popcnt32(bandfree & mask) < lim) continue;
			int tbiv[9], digbiv[9], nbiv = 0;
			for (int idig = 0; idig < 9; idig++){// collect digits bi values in row
				int hidden = mask & zh_g.pmdiag.pmdig[idig].bf.u32[iband];
				if (_popcnt32(hidden) == 2){
					digbiv[nbiv] = idig;
					tbiv[nbiv++] = hidden;
				}
			}
			if (nbiv < 2) continue; // no bi-value in the row
			for (int i = 0; i < nbiv - 1; i++){
				register int R = tbiv[i];
				for (int j = i + 1; j < nbiv; j++) if (R == tbiv[j]){// a pair seen
					if ((bandpairs & R) == R)						continue; //no effect
					tp[np].u32[1] = R;
					tp[np].u16[0] = (uint16_t)iband;
					tp[np++].u16[1] = (uint16_t)((1 << digbiv[i]) | (1 << digbiv[j]));
				}
			}
		}
	}
	return np;
}

int ZHOU::Rate20_HiddenPair_Assign(){
	if (zh_g.diag)		cout << "rate 20_HiddenPair_Assign" << endl;
	zh_g.pairs_naked.SetAll_0();
	GINT64 tp[50];
	int np;
	if (CollectHiddenPairsBox(tp, np,3))		if (RateSingleBox(tp, np)) return zh_g.nsingles;
	if (CollectHiddenPairsRow(tp, np,3))		if (RateSingleRow(tp, np)) return zh_g.nsingles;
	if (CollectHiddenPairsCol(tp, np, 3))	RateSingleCol(tp, np);
	return zh_g.nsingles;
}


int ZHOU::Rate23_SingleInCell_Assign(){
	if (zh_g.diag)	cout << "rating 23 single in cell" << endl;
	BF128 R1 = FD[0][0] & cells_unsolved,
		R = FD[1][0] & cells_unsolved,	R2 = R&R1;
	R1 |= R;
	for (int i = 2; i < 9; i++){
		R = FD[i][0] & cells_unsolved;
		R2 |= R1&R; R1 |= R;
	}	
	//char ws[82];
	//cout << cells_unsolved.String3X(ws) << " unsolved" << endl;
	//cout << R1.String3X(ws) << " R1" << endl;
	//cout << R2.String3X(ws) << " R2" << endl;
	R1 -= R2;// if no single R1=R2
	//cout << R1.String3X(ws) << " R1 actif" << endl;
	if (R1.isEmpty()) return 0;
	int cell;
	while ((cell = R1.getFirst128()) >= 0){
		R1.clearBit(cell);
		zh_g.tsingles[zh_g.nsingles++] = From_128_To_81[cell];
	}
	//cout << "zh_g.nsingles= "<< zh_g.nsingles << endl;
	return 	zh_g.nsingles;
}

int ZHOU::CollectHiddenTripletsBox(GINT64* tp, int & np){
	np = 0;
	BF128 wfree = cells_unsolved - zh_g.locked_nacked_brc_done[0];
	for (int iband = 0, dband = 0; iband < 3; iband++, dband += 27){
		int /*bandpairs = zh_g.locked_nacked_brc_seen[0].bf.u32[iband],*/
			bandfree = wfree.bf.u32[iband];
		for (int iboxr = 0, mask = 07007007; iboxr < 3; iboxr++, mask <<= 3){
			if ((int)_popcnt32(bandfree & mask) < 6) continue;
			int ttrip[9], digtrip[9], ntrip = 0;
			for (int idig = 0; idig < 9; idig++){// collect digits bi values in box
				int hidden = mask & zh_g.pm.pmdig[idig].bf.u32[iband],
					nn = _popcnt32(hidden);
				if (nn && nn<4){// no single so n>1
					digtrip[ntrip] = idig;
					ttrip[ntrip++] = hidden;
				}
			}
			if (ntrip < 3) continue; // no bi-value in the box
			//cout << "box band=" << iband + 1 << " box=" << iboxr + 1 << " ntrip=" << ntrip << endl;
			for (int i = 0; i < ntrip - 2; i++){
				register int R = ttrip[i];
				for (int j = i + 1; j < ntrip-1; j++) {// a pair seen
					register int R2 = R | ttrip[j];
					if (_popcnt32(R2) > 3) continue;
					for (int k = j + 1; k < ntrip; k++) {
						register int R3 = R2 | ttrip[k];
						if (_popcnt32(R3)!= 3) continue;
						tp[np].u32[1] = R3;
						tp[np].u16[0] = (uint16_t)iband;
						tp[np++].u16[1] = (uint16_t)((1 << digtrip[i]) | (1 << digtrip[j]) | (1 << digtrip[k]));
						if (pm_go.opprint&2)cout << Char27out(R3) << " band=" << iband + 1
							<< " digs 0" << oct << tp[np - 1].u16[1]<<dec << endl;
					}
				}
			}
		}
	}
	//cout << "box np=" << np << endl;
	return np;
}
int ZHOU::CollectHiddenTripletsRow(GINT64* tp, int & np){
	//cout << "collect hiden triplet row" << endl;
	np = 0;
	//int tboxskrink[3] = { 0111, 0222, 0444 };
	BF128 wfree = cells_unsolved - zh_g.locked_nacked_brc_done[1];
	for (int iband = 0, dband = 0; iband < 3; iband++, dband += 27){
		int /*bandpairs = zh_g.locked_nacked_brc_seen[1].bf.u32[iband],*/
			bandfree = wfree.bf.u32[iband];
		for (int irowr = 0, mask = 0777; irowr < 3; irowr++, mask <<= 9){
			if ((int)_popcnt32(bandfree & mask) < 6) continue;
			int ttrip[9], digtrip[9], ntrip = 0;
			for (int idig = 0; idig < 9; idig++){// collect digits bi values in row
				int hidden = mask & zh_g.pm.pmdig[idig].bf.u32[iband],
					nn = _popcnt32(hidden);
				if (nn && nn<4){// no single so n>1
					digtrip[ntrip] = idig;
					ttrip[ntrip++] = hidden;
				}
			}
			if (ntrip < 3) continue; // no bi-value in the box
			//cout << "row band=" << iband + 1 << " row=" << irowr + 1 << " ntrip=" << ntrip << endl;
			for (int i = 0; i < ntrip - 2; i++){
				register int R = ttrip[i];
				for (int j = i + 1; j < ntrip - 1; j++) {// a pair seen
					register int R2 = R | ttrip[j];
					if (_popcnt32(R2) > 3) continue;
					for (int k = j + 1; k < ntrip; k++) {
						register int R3 = R2 | ttrip[k];
						if (_popcnt32(R3) != 3) continue;
						tp[np].u32[1] = R3;
						tp[np].u16[0] = (uint16_t)iband;
						tp[np++].u16[1] = (uint16_t)((1 << digtrip[i]) | (1 << digtrip[j]) | (1 << digtrip[k]));
						if (pm_go.opprint & 2)cout << Char27out(R3) << " band=" << iband + 1
							<< " digs 0" << oct << tp[np - 1].u16[1]<<dec << endl;
					}
				}
			}
		}
	}
	//cout << "row np=" << np << endl;
	return np;
}
int ZHOU::CollectHiddenTripletsCol(GINT64* tp, int & np){
	if (pm_go.opprint & 2)cout << "start collect hidden triplet column" << endl;
    np = 0;
	BF128 wfree = zh_g.cells_unsolved_diag - zh_g.locked_nacked_brc_done[2];
	for (int iband = 0, dband = 0/*, ibox = 0*/; iband < 3; iband++, dband += 27){
		int /*bandpairs = zh_g.locked_nacked_brc_seen[2].bf.u32[iband],*/ //digonal mode
			bandfree = wfree.bf.u32[iband];
		for (int irowr = 0, mask = 0777; irowr < 3; irowr++, mask <<= 9){
			if ((int)_popcnt32(bandfree & mask) < 6) continue;
			for (int idig = 0; idig < 9; idig++){// collect digits bi values in row
				int ttrip[9], digtrip[9], ntrip = 0;
				for (int idig = 0; idig < 9; idig++){// collect digits bi values in row
					int hidden = mask & zh_g.pmdiag.pmdig[idig].bf.u32[iband],
						nn = _popcnt32(hidden);
					if (nn && nn<4){// no single so n>1
						digtrip[ntrip] = idig;
						ttrip[ntrip++] = hidden;
					}
				}
				if (ntrip < 3) continue; // no bi-value in therow
				//cout << "col stack=" << iband + 1 << " col=" << irowr + 1 << " ntrip=" << ntrip << endl;
				for (int i = 0; i < ntrip - 2; i++){
					register int R = ttrip[i];
					for (int j = i + 1; j < ntrip - 1; j++) {// a pair seen
						register int R2 = R | ttrip[j];
						if (_popcnt32(R2) > 3) continue;
						for (int k = j + 1; k < ntrip; k++) {
							register int R3 = R2 | ttrip[k];
							if (_popcnt32(R3) != 3) continue;
							tp[np].u32[1] = R3;
							tp[np].u16[0] = (uint16_t)iband;
							tp[np++].u16[1] = (uint16_t)((1 << digtrip[i]) | (1 << digtrip[j]) | (1 << digtrip[k]));
							if (pm_go.opprint & 2)cout << Char27out(R3) << " stack=" << iband + 1
								<< " digs 0" << oct << tp[np - 1].u16[1]<<dec << endl;
						}
					}
				}
			}
		}
	}
	//cout << "col np=" << np << endl;
	return np;
}


int ZHOU::Rate25_HiddenTriplet_Assign(){
	if (pm_go.opprint & 2)		cout << "rate 25_Hiddentriplet_Assign" << endl;
	GINT64 tp[100];
	int np;
	if (CollectHiddenTripletsBox(tp, np))		if (RateSingleBox(tp, np)) return zh_g.nsingles;
	if (CollectHiddenTripletsRow(tp, np))		if (RateSingleRow(tp, np)) return zh_g.nsingles;
	if (CollectHiddenTripletsCol(tp, np))	RateSingleCol(tp, np);
	return zh_g.nsingles;
}

int ZHOU::Rate26_lockedBox(){
	int iret = 0;
	for (int idig = 0; idig < 9; idig++){
		BF128  fd = FD[idig][0] & cells_unsolved;
		for (int iband = 0; iband < 3; iband++){
			register int band = fd.bf.u32[iband];
			if (!band)continue;// solved band
			register int shrink = (TblShrinkMask[band & 0x1FF]) |
				(TblShrinkMask[(band >> 9) & 0x1FF] << 3) |
				(TblShrinkMask[(band >> 18) & 0x1FF] << 6);
			for (int ibox = 0; ibox < 3; ibox++,shrink>>=1){
				register int x = shrink & 0111;
				if (_popcnt32(x) - 1)continue;// not locked in row
				uint32_t irow;	bitscanforward(irow, x);// catch the row
				irow /= 3;
				int to_clean = tband_row[irow] & (~tband_box[ibox]) & band;
				if (to_clean){
					//if (pm_go.opprint2 & 1){
					//	cout << Char27out(to_clean) << "clean 26 row digit=" << idig + 1 << " box=" << ibox + 1
					//		<< " band=" << iband + 1 << " row="<<irow+1<< endl;
					//	cout << hex << "bande=0x" << band << dec << endl;
					//}
					int clean = FD[idig][0].bf.u32[iband] & ~to_clean;
					iret = 1;
					FD[idig][0].bf.u32[iband]  = clean;
				}
			}
		// look now for columns 
			int colband = (band | (band >> 9) | (band >> 18)) & 0x1FF;
			for (int ibox = 0; ibox < 3; ibox++){
				int mask = 7 << (3 * ibox), bm = colband&mask;
				if (_popcnt32(bm) -1) continue; // solved or cleaned
				uint32_t icol;
				bitscanforward(icol, bm);// column to clean
				int to_clean = Zhoucol << icol, *p = tperm3[iband];
				if (to_clean){
					for (int j = 1; j < 3; j++){
						int jband = p[j];
						int jclean = FD[idig][0].bf.u32[jband] & to_clean;
						if (jclean){
							//if (pm_go.opprint2 & 1){
							//	cout <<Char27out(jclean)<< "clean 26 col digit=" << idig + 1 << " box=" << ibox + 1
							//		<< " jband=" << jband + 1 << endl;
							//	cout << hex << "bande=0x" << band << dec << endl;
							//}
							iret = 1;
							FD[idig][0].bf.u32[jband] &= ~jclean;
						}
					}
				}
			}
		}
	}

	return iret;
}
int ZHOU::Rate28_lockedRowCol(){// 26 (box) done 
	int iret = 0;
	for (int idig = 0; idig < 9; idig++){
		BF128  &fd = FD[idig][0];
		for (int iband = 0; iband < 3; iband++){
			register int band = fd.bf.u32[iband];
			if (!(band& cells_unsolved.bf.u32[iband]))continue;// solved band
			register int shrink = (TblShrinkMask[band & 0x1FF]) |
				(TblShrinkMask[(band >> 9) & 0x1FF] << 3) |
				(TblShrinkMask[(band >> 18) & 0x1FF] << 6);
			register int b = band &TblComplexMask[shrink];
			if (b != band){
				//cout << "28 lockedrow dig=" << idig + 1 << " band=" << iband + 1 << endl;
				iret = 1;
				fd.bf.u32[iband] = b;
			}
		}
		// look for locked in column
		int tcol[3];
		for (int iband = 0; iband < 3; iband++){
			register int A = fd.bf.u32[iband] & cells_unsolved.bf.u32[iband];
			tcol[iband] = ((A | (A >> 9) | (A >> 18)) & 0x1FF);// active mini columns
		}
		static int tbox_to_clean[8] = { 0, 01001001, 02002002, 03003003, 04004004, 05005005, 06006006, 07007007 };
		for (int iband = 0; iband < 3; iband++){
			int * p = tperm3[iband];
			register int R = tcol[iband],
				W = R & ~(tcol[p[1]] | tcol[p[2]]);// one minicol
			for (int iboxr = 0; iboxr < 3; iboxr++, R >>= 3, W >>= 3){// check the thre boxes
				int active_col = R & 7;
				if (_popcnt32(active_col) < 2) continue;
				if (!(W & 7))continue;
				int columns_to_clean_in_box = active_col &  ~W;
				int clean = ~(fd.bf.u32[iband] & (tbox_to_clean[columns_to_clean_in_box] << (3 * iboxr)));
				// remark clean is not~ 0
				//cout << "28 locked col dig=" << idig + 1 << " band=" << iband + 1 << " box=" << iboxr + 1
				//	<< "active_cols=" << active_col
				//	<< " columns_to_clean_in_box=" << columns_to_clean_in_box << endl;
				iret = 1;
				clean &= fd.bf.u32[iband];
				fd.bf.u32[iband] = clean;
			}
		}
	}
	return iret;
}

int ZHOU::CleanNake(int & naked, int unit, int iband, int modebr){
	int iret = 0;
	while (naked){
		int cell1bf = _blsi_u32(naked); // get lower bit as cell 1
		int digits[5], nd = 0, nakp = 0;// accept naked 5
		for (int i = 0; i < 9; i++){
			register uint32_t & w = zh_g.pm.pmdig[i].bf.u32[iband];
			if (w & cell1bf){
				digits[nd++] = i;
				nakp |= w&unit&naked;
			}
		}
		if (0){
			cout << Char27out(nakp) << " br pair nd=" << nd << endl;
			if (nd != 2) return 0;
		}
		naked ^= nakp; // erase the digits for several pairs
		// clean the 2 digits
		int zclean = unit & ~nakp/*, go_clean = 0*/;
		for (int i = 0; i < nd; i++){
			int idig = digits[i];
			register uint32_t & w = FD[idig][0].bf.u32[iband];
			int clean = w&zclean;
			if (clean){
				//cout << Char27out(clean) << " br pair clean digit=" << idig+1 << endl;
				iret = 1;
				w &= ~clean;
			}
		}
		zh_g.locked_nacked_brc_done[modebr].bf.u32[iband] |= nakp;
	}
	return iret;
}
int ZHOU::CleanNakeColumn(int & naked, int unit, int iband){
	int iret = 0;
	while (naked){
		//cout << Char27out(naked) << " column status" << endl;
		int cell1bf = _blsi_u32(naked); // get lower bit as cell 1
		int digits[5], nd = 0, nakp = 0;// accept naked 5
		for (int i = 0; i < 9; i++){
			register uint32_t & w = zh_g.pmdiag.pmdig[i].bf.u32[iband];
			if (w & cell1bf){
				digits[nd++] = i;
				nakp |= w&naked&unit;
			}
		}
		if (0){
			cout << Char27out(nakp) << " col pair nd=" << nd << endl;
			if (nd != 2) return 0;
		}
		naked ^= nakp; // erase the digits for several pairs
		// clean the 2 digits
		int zclean = unit & ~nakp/*, go_clean = 0*/;
		for (int i = 0; i < nd; i++){
			int idig = digits[i];
			register uint32_t  w = zh_g.pmdiag.pmdig[idig].bf.u32[iband];
			int clean = w&zclean;
			if (clean){// clean the source 
				//cout << Char27out(clean) << " col clean digit=" << idig + 1 << endl;
				FD[idig][0].ClearDiag(clean, iband);
				iret = 1;
			}
		}
		zh_g.locked_nacked_brc_done[2].bf.u32[iband] |= nakp;
    }
	return iret;
}
int ZHOU::Rate30_NakedPair(){
	int iret = 0;
    // first band => box row
	BF128 lb = zh_g.locked_nacked_brc_seen[0] - zh_g.locked_nacked_brc_done[0],
	     lr = zh_g.locked_nacked_brc_seen[1] - zh_g.locked_nacked_brc_done[1],
		 lc = zh_g.locked_nacked_brc_seen[2] - zh_g.locked_nacked_brc_done[2];
	//char ws[82];
	//cout << lb.String3X(ws) << " lb" << endl;
	//cout << lr.String3X(ws) << " lr" << endl;
	//cout << lc.String3X(ws) << " lc" << endl;

	for (int iband = 0; iband < 3; iband++){
		int boxband = lb.bf.u32[iband];
		if (_popcnt32(boxband) >= 2){
			for (int ibox = 0; ibox < 3; ibox++){
				int box = boxband & tband_box[ibox];
				if (_popcnt32(box) >= 2){// one or more naked pair in the box
					//cout << Char27out(box) << " box naked pair band=" << iband + 1 << endl;
					iret+=CleanNake(box, tband_box[ibox], iband,0);
				}
			}
		}
		int rowband = lr.bf.u32[iband];
		if (_popcnt32(rowband) >= 2){
			for (int irow = 0; irow < 3; irow++){
				int row = rowband & tband_row[irow];
				if (_popcnt32(row) >= 2){// one or more naked pair in the row
					//cout << Char27out(row) << " row naked pair band=" << iband + 1 << endl;
					iret += CleanNake(row, tband_row[irow], iband,1);
				}
			}
		}
		int colband = lc.bf.u32[iband];
		if (_popcnt32(colband) >= 2){
			for (int icol = 0; icol < 3; icol++){
				int row = colband & tband_row[icol];// diagonal mode, the column looks like a row
				if (_popcnt32(row) >= 2){// one or more naked pair in the row
					//cout << Char27out(row) << " col naked pair band=" << iband + 1 << endl;
					iret += CleanNakeColumn(row, tband_row[icol], iband);
				}
			}
		}
	}
	//ImageCandidats();
    return iret;
}

int ZHOU::Rate32_XWing(){
	int iret = 0;
	for (int idig = 0; idig < 9; idig++){
		int locdiag = 0;
		int & rx2 = zh_g.row_col_x2[idig][0], 
			& cx2 = zh_g.row_col_x2[idig][1],
			& rfree = (int &)FD[idig][0].bf.u32[3],
			& cfree = zh_g.dig_unsolved_col[idig],
			* rows = zh_g.dig_rows[idig], 
			* cols = zh_g.dig_cols[idig],
			nfree = zh_g.unsolved_r_count[idig],
			tr[9], nr = 0;
		cfree = 0;
		for (int i = 0; i < 9; i++)if (cols[i]) cfree |= 1 << i;
		if (0){
			cout << "analyse bug digit 8" << endl;
			cout << Char9out(rx2) << " rx2" << endl;
			cout << Char9out(cx2) << " cx2" << endl;
			cout << Char9out(rfree) << " rfree" << endl;
			cout << Char9out(cfree) << " cfree" << endl;
			cout << "nfree=" << nfree << endl;
			for (int i = 0; i < 9; i++){
				cout << Char9out(rows[i]) << "\t";
				cout << Char9out(cols[i]) << " r/c i=" <<i<< endl;				
			}
			locdiag = 1;
		}
		//============================ XWing in rows
		int rw = rfree & ~rx2;
		if (_popcnt32(rw) < 3) continue;//not  goto columns; same count in both ways
		for (int irow = 0; irow < 9; irow++)
			if(rw & (1<<irow))if (_popcnt32(rows[irow]) == 2)tr[nr++] = irow;
		if (nr < 2) goto columns;
		for (int ir1 = 0; ir1 < nr - 1; ir1++){
			int r1 = tr[ir1];
			register int R = rows[r1];
			for (int ir2 = ir1 + 1; ir2 < nr; ir2++){
				int r2 = tr[ir2];
				if (R != rows[r2]) continue;
				uint32_t c1, c2;
				bitscanforward(c1, R);	bitscanreverse(c2, R);
				int mask = (1 << r1) | (1 << r2);
				rx2 |= mask; 
				cx2 |= R;
				int clean1 = cols[c1] & ~mask, clean2 = cols[c2] & ~mask;
				if (clean1 | clean2){//some clening appears
					iret = 1;
					if (pm_go.opprint2 & 2)cout << "XWING row digit=" << idig + 1 << " " << r1 + 1 << r2 + 1 << " " << c1 + 1 << c2 + 1 << endl;
					if (clean1)FD[idig][0].ClearCol(clean1, c1);
					if (clean2)FD[idig][0].ClearCol(clean2, c2);
					break;
				}
			}
		}
		//============================ XWing in columns
		columns:
		if (locdiag) cout << "see columns" << endl;
		nr = 0;
		int cw = cfree & ~cx2;
		if (_popcnt32(cw) < 3) continue;
		for (int icol = 0; icol < 9; icol++)
			if (cw & (1 << icol))if (_popcnt32(cols[icol]) == 2)tr[nr++] = icol;
		if (nr < 2) continue;
		if (locdiag) cout << "nr=" << nr << endl;
		for (int ic1 = 0; ic1 < nr - 1; ic1++){
			int c1 = tr[ic1];
			register int R = cols[c1];
			for (int ic2 = ic1 + 1; ic2 < nr; ic2++){
				int c2 = tr[ic2];
				if (locdiag) cout << "c1;c2=" <<c1<<c2 << endl;

				if (R != cols[c2]) continue;
				uint32_t r1, r2;
				bitscanforward(r1, R);	bitscanreverse(r2, R);
				if (locdiag) cout << "r1;r2=" << r1 << r2 << endl;
				int mask = (1 << c1) | (1 << c2);
				cx2 |= mask;
				rx2 |= R;
				int clean1 = rows[r1] & ~mask, clean2 = rows[r2] & ~mask;
				if (clean1 | clean2){//some cleaning appears
					iret = 1;
					if (pm_go.opprint2 & 2)cout << "XWING col digit=" << idig + 1 << " " << c1 + 1 << c2 + 1 << " " << r1 + 1 << r2 + 1 << endl;
					if (clean1)FD[idig][0].ClearRow(clean1, r1);
					if (clean2)FD[idig][0].ClearRow(clean2, r2);
					break;
				}
			}
		}
		// exit if now <3
	}
	return iret;
}

int ZHOU::ApplyHidden(GINT64 * t, int nt,int mode){
	if (!nt) return 0;
	int iret = 0;
	for (int i = 0; i < nt; i++){
		GINT64 w = t[i];
		register int R = w.u32[1]; // cells pattern of the locked set
		int iband = w.u16[0], digs = w.u16[1];
		for (int idig = 0, dig = 1; idig < 9; dig <<= 1, idig++){
			if(dig &digs)continue;// skip digs of the hidden set
			uint32_t  & band = FD[idig][0].bf.u32[iband] ,
				clean=band & R;
			if (clean){
				iret = 1;
				band &= ~clean;
			}
		}
		zh_g.locked_nacked_brc_done[mode].bf.u32[iband] |= R;
	}
	return iret;
}
int ZHOU::ApplyHiddenColumn(GINT64 * t, int nt){
	if (!nt) return 0;
	int iret = 0;
	for (int i = 0; i < nt; i++){
		GINT64 w = t[i];
		register int R = w.u32[1]; // cells pattern of the locked set
		int iband = w.u16[0], digs = w.u16[1];
		for (int idig = 0, dig = 1; idig < 9; dig <<= 1, idig++){
			if (dig &digs)continue;// skip digs of the hidden set
			int   band = zh_g.pmdiag.pmdig[idig].bf.u32[iband],
				clean = band & R;
			if (clean){
				iret = 1;
				FD[idig][0].ClearDiag(clean, iband);
			}
		}
		zh_g.locked_nacked_brc_done[2].bf.u32[iband] |= R;
	}
	return iret;
}
int ZHOU::Rate34_HiddenPair(){
	GINT64 tp[50];
	int np,iret=0;
	if (CollectHiddenPairsBox(tp, np, 3))iret += ApplyHidden(tp, np,0);
	if (CollectHiddenPairsRow(tp, np, 3))iret += ApplyHidden(tp, np,1);
	if (CollectHiddenPairsCol(tp, np, 3))iret += ApplyHiddenColumn(tp, np);
	return iret;
}

int ZHOU::Rate36_FindClean_NakedTriplet(int unit_triplet_cells, int unit, int iband, int mode){
	int iret = 0, t3[10], t3b[10], n3 = 0;
	BitsInTable32(t3, n3, unit_triplet_cells);
	if (n3 < 3) return 0;
	for (int i = 0; i < n3; i++) {
		t3b[i] = 1 << t3[i];
		t3[i] += 27 * iband;
	}
	for (int i1 = 0; i1 < n3 - 2; i1++){
		int digs1 = zh_g.dig_cells[t3[i1]];
		for (int i2 = i1+1; i2 < n3 - 1; i2++){
			int digs2 = digs1 |zh_g.dig_cells[t3[i2]];
			if (_popcnt32(digs2) > 3) continue;
			for (int i3 = i2 + 1; i3 < n3; i3++){
				int digs3 = digs2 | zh_g.dig_cells[t3[i3]];
				if (_popcnt32(digs3) > 3) continue;
				// this is a naked triplet 
				int btrip = t3b[i1] | t3b[i2] | t3b[i3];
				//cout << Char27out(btrip) << " triplet � voir band=" << iband + 1 << " mode=" << mode << endl;
				zh_g.locked_nacked_brc_done[mode].bf.u32[iband] |= btrip;
				int zclean = unit & ~btrip;
				for (int idig = 0; idig < 9; idig++) {// box or row
					register uint32_t & w = FD[idig][0].bf.u32[iband];
					if (!(w&btrip)) continue;
					int clean = w&zclean;
					if (clean){
						iret = 1;
						w &= ~clean;
					}
				}
				if (iret) return 1;// nearly no chance to have more see later
				unit_triplet_cells ^= btrip; // loop if room for another triplet
				if (_popcnt32(unit_triplet_cells) >= 3)
					return Rate36_FindClean_NakedTriplet(unit_triplet_cells, unit, iband, mode);
			}
		}
	}
	return iret;
}

int ZHOU::Rate36_FindClean_NakedTripletCol(int unit_triplet_cells, int unit, int iband){
	int iret = 0, t3[10], t3b[10], n3 = 0;
	BitsInTable32(t3, n3, unit_triplet_cells);
	if (n3 < 3) return 0;
	for (int i = 0; i < n3; i++) {
		int cell = t3[i];
		t3b[i] = 1 << cell;
		t3[i]= C_transpose_d[cell + 27 * iband];;
	}
	for (int i1 = 0; i1 < n3 - 2; i1++){
		int digs1 = zh_g.dig_cells[t3[i1]];
		for (int i2 = i1 + 1; i2 < n3 - 1; i2++){
			int digs2 = digs1 | zh_g.dig_cells[t3[i2]];
			if (_popcnt32(digs2) > 3) continue;
			for (int i3 = i2 + 1; i3 < n3; i3++){
				int digs3 = digs2 | zh_g.dig_cells[t3[i3]];
				if (_popcnt32(digs3) > 3) continue;
				// this is a naked triplet 
				int btrip = t3b[i1] | t3b[i2] | t3b[i3];
				//cout << Char27out(btrip) << "col  triplet � voir band=" << iband + 1  << endl;
				zh_g.locked_nacked_brc_done[2].bf.u32[iband] |= btrip;
				int zclean = unit & ~btrip;
				for (int idig = 0; idig < 9; idig++){// box or row
					register uint32_t & w = zh_g.pmdiag.pmdig[idig].bf.u32[iband];
					if (!(w&btrip)) continue;
					int clean = w&zclean;
					if (clean){
						iret = 1;
						FD[idig][0].ClearDiag(clean, iband);
					}
				}
				if (iret) return 1;// nearly no chance to have more see later
				unit_triplet_cells ^= btrip; // loop if room for another triplet
				if (_popcnt32(unit_triplet_cells) >= 3)
					return Rate36_FindClean_NakedTripletCol(unit_triplet_cells, unit, iband);
			}
		}
	}
	return iret;
}
int ZHOU::Rate36_NakedTriplet(){
	int iret = 0;
	// unsolved cells with 2 or 3 candidates not yet processed
	BF128 wu = zh_g.triplets | zh_g.pairs,
		wudiag;
	wudiag.Diag3x27(wu);
	BF128 lb = wu - zh_g.locked_nacked_brc_done[0],
		lr = wu - zh_g.locked_nacked_brc_done[1],
		lc = wudiag - zh_g.locked_nacked_brc_done[2],
		lbu = cells_unsolved - zh_g.locked_nacked_brc_done[0],
		lru = cells_unsolved - zh_g.locked_nacked_brc_done[1],
		lcu = zh_g.cells_unsolved_diag - zh_g.locked_nacked_brc_done[2];

	for (int iband = 0; iband < 3; iband++){
		int boxband = lb.bf.u32[iband],
			buns = lbu.bf.u32[iband];
		if (_popcnt32(buns) >= 6){// no naked pair no hidden pair seen
			for (int ibox = 0; ibox < 3; ibox++){
				int box = boxband & tband_box[ibox],
					boxuns = buns & tband_box[ibox];
				if (_popcnt32(box) >= 3 && _popcnt32(boxuns)>=6){
					// one or more naked triplets possibles in the box
					//cout << Char27out(box) << " box naked pair band=" << iband + 1 << endl;
					iret += Rate36_FindClean_NakedTriplet(box, tband_box[ibox], iband, 0);
				}
			}
		}
		int rowband = lr.bf.u32[iband],
			runs = lru.bf.u32[iband];
		if (_popcnt32(runs) >= 6){
			for (int irow = 0; irow < 3; irow++){
				int row = rowband & tband_row[irow],
					rowuns = runs & tband_row[irow];
				if (_popcnt32(row) >= 3 && _popcnt32(rowuns) >= 6){
					// one or more naked triplet possible  in the row
					//cout << Char27out(row) << " row naked pair band=" << iband + 1 << endl;
					iret += Rate36_FindClean_NakedTriplet(row, tband_row[irow], iband, 1);
				}
			}
		}
		int colband = lc.bf.u32[iband],
			cuns = lcu.bf.u32[iband];;
		if (_popcnt32(cuns) >= 6){
			for (int icol = 0; icol < 3; icol++){
				int col = colband & tband_row[icol],
					coluns = cuns & tband_row[icol];// diagonal mode, the column looks like a row
				if (_popcnt32(col) >= 3 && _popcnt32(coluns) >= 6){// one or more naked pair in the row
					//cout << Char27out(col) << " col possible naked triplet band=" << iband + 1 << endl;
					iret += Rate36_FindClean_NakedTripletCol(col, tband_row[icol], iband);
				}
			}
		}
	}
	//ImageCandidats();
	return iret;
}


int ZHOU::Rate38_SwordFish(){
	int iret = 0;
	for (int idig = 0; idig < 9; idig++){
		int & rx2 = zh_g.row_col_x2[idig][0],
			&cx2 = zh_g.row_col_x2[idig][1],
			&rfree = (int &)FD[idig][0].bf.u32[3],
			&cfree = zh_g.dig_unsolved_col[idig],
			*rows = zh_g.dig_rows[idig],
			*cols = zh_g.dig_cols[idig],
			nfree = zh_g.unsolved_r_count[idig],
			tr[9], nr = 0;
		cfree = 0;
		for (int i = 0; i < 9; i++)if (cols[i]) cfree |= 1 << i;
		if (0){//pm_go.opprint2 & 1){
			cout << "analyse swordfish digit "<<idig+1 << endl;
			cout << Char9out(rx2) << " rx2" << endl;
			cout << Char9out(cx2) << " cx2" << endl;
			cout << Char9out(rfree) << " rfree" << endl;
			cout << Char9out(cfree) << " cfree" << endl;
			cout << "nfree=" << nfree << endl;
			for (int i = 0; i < 9; i++){
				cout << Char9out(rows[i]) << "\t";
				cout << Char9out(cols[i]) << " r/c i=" << i << endl;
			}
		}
		//============================ SwordFish in rows
		int rw = rfree & ~rx2;
		if (_popcnt32(rw) < 4) continue;// must be at minimum 4 "free" rows/cols
		for (int irow = 0; irow < 9; irow++)
			if (rw & (1 << irow))if (_popcnt32(rows[irow]) <4)tr[nr++] = irow;
		if (nr < 3) goto columns;
		for (int ir1 = 0; ir1 < nr - 2; ir1++){
			int r1 = tr[ir1];
			register int R1 = rows[r1];
			for (int ir2 = ir1 + 1; ir2 < nr-1; ir2++){
				int r2 = tr[ir2];
				int R2 = R1 | rows[r2];
				if (_popcnt32(R2) > 3) continue;
				for (int ir3 = ir2 + 1; ir3 < nr; ir3++){
					int r3 = tr[ir3];
					int R3 = R2 | rows[r3];
					if (_popcnt32(R3) > 3) continue;
					uint32_t c1, c2, c3;
					bitscanforward(c1, R3);	bitscanreverse(c2, R3);
					int w =R3 ^ ((1 << c1) | (1 <<c2));//last bit
					bitscanforward(c3, w);
					int mask = (1 << r1) | (1 << r2) | (1 << r3);
					rx2 |= mask;
					cx2 |= R3;
					int clean1 = cols[c1] & ~mask, clean2 = cols[c2] & ~mask, clean3 = cols[c3] & ~mask;
					if (clean1 | clean2 | clean3){//some cleaning appears
						iret = 1;
						if (pm_go.opprint2 & 2)cout << "swordfish row digit=" << idig + 1 << " " << r1 + 1 << r2 + 1 << r3 + 1
							<< " " << c1 + 1 << c2 + 1 << c3 + 1 << endl;
						if (clean1)FD[idig][0].ClearCol(clean1, c1);
						if (clean2)FD[idig][0].ClearCol(clean2, c2);
						if (clean3)FD[idig][0].ClearCol(clean3, c3);
						break;
					}
				}
			}
		}
		//============================ SwordFish in columns
	columns:
		nr = 0;
		int cw = cfree & ~cx2;
		for (int icol = 0; icol < 9; icol++)
			if (cw & (1 << icol))if (_popcnt32(cols[icol]) <4)tr[nr++] = icol;
		if (nr < 3) continue;
		for (int ic1 = 0; ic1 < nr - 2; ic1++){
			int c1 = tr[ic1];
			int R1 = cols[c1];
			for (int ic2 = ic1 + 1; ic2 < nr-1; ic2++){
				int c2 = tr[ic2];
				int R2 = R1 | cols[c2];
				if (_popcnt32(R2) > 3) continue;
				for (int ic3 = ic2 + 1; ic3 < nr; ic3++){
					int c3 = tr[ic3];
					int R3 = R2 | cols[c3];
					if (_popcnt32(R3) > 3) continue;
					uint32_t r1, r2, r3;
					bitscanforward(r1, R3);	bitscanreverse(r2, R3);
					int w = R3 ^ ((1 << r1) | (1 << r2));//last bit
					bitscanforward(r3, w);
					int mask = (1 << c1) | (1 << c2) | (1 << c3);
					cx2 |= mask;
					rx2 |= R3;
					int clean1 = rows[r1] & ~mask, clean2 = rows[r2] & ~mask, clean3 = rows[r3] & ~mask;
					if (clean1 | clean2 | clean3){//some cleaning appears
						iret = 1;
						if (pm_go.opprint2 & 2)cout << "swordfish col digit=" << idig + 1 << " " << c1 + 1 << c2 + 1 << c3 + 1
							<< " " << r1 + 1 << r2 + 1 << r3 + 1 << endl;
						if (clean1)FD[idig][0].ClearRow(clean1, r1);
						if (clean2)FD[idig][0].ClearRow(clean2, r2);
						if (clean3)FD[idig][0].ClearRow(clean3, r3);
						break;
					}
				}
			}
		}
		// exit if now <3
	}
	return iret;
}
int ZHOU::Rate40_HiddenTriplet(){
	GINT64 tp[50];
	int np, iret = 0;
	if (CollectHiddenTripletsBox(tp, np))	iret += ApplyHidden(tp, np, 0);
	if (CollectHiddenTripletsRow(tp, np))	iret += ApplyHidden(tp, np, 1);
	if (CollectHiddenTripletsCol(tp, np))	iret += ApplyHiddenColumn(tp, np);
	return iret;
}

int ZHOU::Rate42_XYWing(){
	zh_g.Build_digits_cells_pair_bf();
	BF128 wp = zh_g.pairs;
	if (wp.Count() < 3) return 0;
	int xcell,xcell1,xcell2,iret=0;
	while ((xcell = wp.getFirst128()) >= 0){
		wp.clearBit(xcell);
		int cell = From_128_To_81[xcell],d1,d2,d3,nd=0;
		for (int i = 0; i < 9; i++){// find the 2 digits
			if (zh_g.digits_cells_pair_bf[i].On(xcell)){
				if (nd){ d2 = i; break; }
				else{ d1 = i; nd = 1; }
			}
		}
		BF128 zd1 = zh_g.digits_cells_pair_bf[d1];
		zd1 &= cell_z3x[cell];// seen by cell
		zd1 -= zh_g.digits_cells_pair_bf[d2];// not a naked pair
		if (zd1.isEmpty()) continue;
		while ((xcell1 = zd1.getFirst128()) >= 0){// try each cell zd1
			zd1.clearBit(xcell1);
			int cell1 = From_128_To_81[xcell1];
			for (int i = 0; i < 9; i++)if((i-d1)&& (i-d2))// find  digit3
				if (zh_g.digits_cells_pair_bf[i].On(xcell1)){	 d3 = i; break; }
			BF128 zd2 = zh_g.digits_cells_pair_bf[d2];
			zd2 &= cell_z3x[cell];// seen by cell
			zd2 &= zh_g.digits_cells_pair_bf[d3];// second digit is d3
			if (zd2.isEmpty()) continue;// usually empty
			while ((xcell2 = zd2.getFirst128()) >= 0){// try each cell zd1
				zd2.clearBit(xcell2);
				int cell2 = From_128_To_81[xcell2];
				BF128 clean = FD[d3][0]; 
				clean &= cell_z3x[cell1]; 		clean &= cell_z3x[cell2];
				if (clean.isNotEmpty()){                    
					FD[d3][0]-= clean;
					iret = 1;
				}
			}
		}
	}
return iret;
}
int ZHOU::Rate44_XYZWing(){
	int iret = 0;
	// build zh_g.digits_cells_pair_bf[9]; 	zh_g.pairs   	zh_g.triplets 
	//zh_g.Build_digits_cells_pair_bf();
	BF128 wp = zh_g.pairs, wp3 = zh_g.triplets;
	if (wp.Count() < 2 || wp3.isEmpty()) return 0;
	int xcell;
	while ((xcell = wp3.getFirst128()) >= 0){// try each triplet as pivot
		wp3.clearBit(xcell);
		int cell = From_128_To_81[xcell], dx[3], nd = 0;
		BF128 wpw = wp; wpw &= cell_z3x[cell];
		if (wpw.Count() < 2) continue; // must see 2 pairs
		int digits = zh_g.dig_cells[cell];
		for (int i = 0, id = 1; i < 9; i++, id <<= 1)// find the 2 digits
			if (digits & id)dx[nd++] = i;
		if (nd != 3) return 0; // safety should never be

		for (int ind = 0; ind < 3; ind++){// try each digit as "common to 3 cells"
			BF128 zd1 = zh_g.digits_cells_pair_bf[dx[ind]] & wpw;
			if (zd1.Count() < 2) continue; // must see 2 pairs or more
			int tpz1[20], ntpz1 = zd1.Table3X27(tpz1);// collect pairs seen by pivot
			for (int i1 = 0; i1 < ntpz1 - 1; i1++){
				int cell1 = tpz1[i1], digs1 = zh_g.dig_cells[cell1];
				if (digs1 & ~digits) continue; // not fitting XYZWing
				for (int i2 = i1 + 1; i2 < ntpz1; i2++){
					int cell2 = tpz1[i2], digs2 = zh_g.dig_cells[cell2];
					if (digs2 & ~digits) continue; // not fitting XYZWing
					if (digs1 == digs2) continue;// not fitting XYZWing
					BF128 clean = FD[dx[ind]][0]; clean.Clear(xcell);
					clean &= cell_z3x[cell]; 		
					clean &= cell_z3x[cell1]; 		clean &= cell_z3x[cell2];
					if (clean.isNotEmpty()){
						if (pm_go.opprint2 & 2)cout << "XYZWing " << cell << " " << cell1 << " " << cell2 << endl;
						FD[dx[ind]][0] -= clean;
						iret = 1;
					}
				}
			}
		}
	}
	return iret;
}
int ZHOU::Rate50_NakedQuad(){// pair/triplets done
	int iret = 0, cell, tmode[3] = { 1, 2, 0 };
   
	for (int iu = 0; iu < 27; iu++){// look in unit
		BF128 wu = units3xBM[iu]; wu &= cells_unsolved;
 		if (wu.Count() < 8) continue;
		//cout << "try unit=" << iu << endl;
		BF128 ww = wu, wuu = wu;
		int mode = tmode[iu/9],tcells[10],tdigs[10],ncells=0;
		if (mode == 2)wuu.Diag3x27(wu);// locked done[2] is diagonal
		if ((wuu&zh_g.locked_nacked_brc_done[mode]).isNotEmpty()) continue;
		//cout << "try go unit=" << iu << endl;
		while ((cell = ww.getFirsCell()) >= 0){// collect cells with max 4 digits
			ww.Clear_c(cell);
			register int d = zh_g.dig_cells[cell];
			if (_popcnt32(d)>4)continue;
			tdigs[ncells] = d;
			tcells[ncells++] = cell; 
		}
		if (ncells < 4)continue;
		//cout << "try go unit=" << iu << "ncells4 = "<< ncells << endl;
		for (int i1 = 0; i1 < ncells - 3; i1++){
			int dig1 = tdigs[i1];
			for (int i2 = i1+1; i2 < ncells - 2; i2++){
				int dig2 = tdigs[i2]|dig1;
				if (_popcnt32(dig2) > 4)continue;
				for (int i3 = i2 + 1; i3 < ncells - 1; i3++){
					int dig3 = tdigs[i3] | dig2;
					if (_popcnt32(dig3) > 4)continue;
					for (int i4 = i3 + 1; i4 < ncells; i4++){
						int dig4 = tdigs[i4] | dig3;
						if (_popcnt32(dig4) != 4)continue;
						//cout << "seen a nacked quas" << endl;
						ww = wu; ww.Clear_c(tcells[i1]); ww.Clear_c(tcells[i2]);
						ww.Clear_c(tcells[i3]); ww.Clear_c(tcells[i4]);
						for (int idig = 0; idig < 9; idig++)if (dig4 & (1 << idig)){
							BF128 clean = FD[idig][0]&ww;
							if (clean.isNotEmpty()){
								FD[idig][0] -= clean;
								iret = 1;
							}
						}
					}
				}
			}
		}
	}

	return iret;
}
int ZHOU::Rate52_JellyFish(){
	if (pm_go.opprint2 & 2) cout << "look for Jelly fish rating 52" << endl;
	int iret = 0;
	for (int idig = 0; idig < 9; idig++){
		int & rx2 = zh_g.row_col_x2[idig][0],
			&cx2 = zh_g.row_col_x2[idig][1],
			&rfree = (int &)FD[idig][0].bf.u32[3],
			&cfree = zh_g.dig_unsolved_col[idig],
			*rows = zh_g.dig_rows[idig],
			*cols = zh_g.dig_cols[idig],
			//nfree = zh_g.unsolved_r_count[idig],
			tr[9], nr = 0;
		cfree = 0;
		for (int i = 0; i < 9; i++)if (cols[i]) cfree |= 1 << i;
		//============================ SwordFish in rows
		int rw = rfree & ~rx2;
		if (_popcnt32(rw) < 5) continue;// must be at minimum 5 "free" rows/cols
		for (int irow = 0; irow < 9; irow++)
			if (rw & (1 << irow))if (_popcnt32(rows[irow]) <5)tr[nr++] = irow;
		if (nr < 4) goto columns;
		for (int ir1 = 0; ir1 < nr - 3; ir1++){
			int r1 = tr[ir1];
			register int R1 = rows[r1];
			for (int ir2 = ir1 + 1; ir2 < nr - 2; ir2++){
				int r2 = tr[ir2];
				int R2 = R1 | rows[r2];
				if (_popcnt32(R2) >4 ) continue;
				for (int ir3 = ir2 + 1; ir3 < nr-1; ir3++){
					int r3 = tr[ir3];
					int R3 = R2 | rows[r3];
					if (_popcnt32(R3) > 4) continue;
					for (int ir4 = ir3 + 1; ir4 < nr; ir4++){
						int r4 = tr[ir4];
						int R4 = R3 | rows[r4];
						if (_popcnt32(R4) > 4) continue;
						uint32_t c1, c2, c3, c4;
						bitscanforward(c1, R3);	bitscanreverse(c2, R3);
						int w = R3 ^ ((1 << c1) | (1 << c2));//last bit
						bitscanforward(c3, w); bitscanreverse(c4, w);
						int mask = (1 << r1) | (1 << r2) | (1 << r3) | (1 << r4);
						//cx2 |= mask;//?? risk to have it downgraded later to Xwing
						//rx2 |= R4;
						int clean1 = cols[c1] & ~mask, clean2 = cols[c2] & ~mask,
							clean3 = cols[c3] & ~mask, clean4 = cols[c4] & ~mask;
						if (clean1 | clean2 | clean3 | clean4){//some cleaning appears
							iret = 1;
							if (pm_go.opprint2 & 2)cout << "Jellyfish row digit=" << idig + 1 << " rows " << r1 + 1 << r2 + 1 << r3 + 1 << r4 + 1
								<< " cols " << c1 + 1 << c3 + 1 << c4 + 1 << c2 + 1 << endl;
							if (clean1)FD[idig][0].ClearCol(clean1, c1);
							if (clean2)FD[idig][0].ClearCol(clean2, c2);
							if (clean3)FD[idig][0].ClearCol(clean3, c3);
							if (clean4)FD[idig][0].ClearCol(clean4, c4);
							break;
						}
					}
				}
			}
		}
		//============================ SwordFish in columns
	columns:
		nr = 0;
		int cw = cfree & ~cx2;
		for (int icol = 0; icol < 9; icol++)
			if (cw & (1 << icol))if (_popcnt32(cols[icol]) <5)tr[nr++] = icol;
		if (pm_go.opprint2 & 2) cout << "Jelly fish columns nr=" <<nr
        << endl;


		if (nr < 4) continue;
		for (int ic1 = 0; ic1 < nr - 3; ic1++){
			int c1 = tr[ic1];
			int R1 = cols[c1];
			for (int ic2 = ic1 + 1; ic2 < nr - 2; ic2++){
				int c2 = tr[ic2];
				int R2 = R1 | cols[c2];
				if (_popcnt32(R2) > 4) continue;
				for (int ic3 = ic2 + 1; ic3 < nr-1; ic3++){
					int c3 = tr[ic3];
					int R3 = R2 | cols[c3];
					if (_popcnt32(R3) > 4) continue;
					for (int ic4 = ic3 + 1; ic4 < nr; ic4++){
						int c4 = tr[ic4];
						int R4 = R3 | cols[c4];
						if (_popcnt32(R4) > 4) continue;
						uint32_t r1, r2, r3, r4;
						bitscanforward(r1, R3);	bitscanreverse(r2, R3);
						int w = R3 ^ ((1 << r1) | (1 << r2));//last bit
						bitscanforward(r3, w); bitscanreverse(r4, w);
						int mask = (1 << c1) | (1 << c2) | (1 << c3) | (1 << c4);
						//rx2 |= mask;//?? risk to have it downgraded later to Xwing
						//cx2 |= R3;
						int clean1 = rows[r1] & ~mask, clean2 = rows[r2] & ~mask, 
							clean3 = rows[r3] & ~mask, clean4 = rows[r4] & ~mask;
						if (clean1 | clean2 | clean3 | clean4){//some cleaning appears
							iret = 1;
							if (pm_go.opprint2 & 2)cout << "swordfish col digit=" << idig + 1 << " cols " << c1 + 1 << c2 + 1 << c3 + 1 << c4 + 1
								<< " rows " << r1 + 1 << r3 + 1 << r4 + 1 << r2 + 1 << endl;
							if (clean1)FD[idig][0].ClearRow(clean1, r1);
							if (clean2)FD[idig][0].ClearRow(clean2, r2);
							if (clean3)FD[idig][0].ClearRow(clean3, r3);
							if (clean4)FD[idig][0].ClearRow(clean4, r4);
							break;
						}
					}
				}
			}
		}
		// exit if now <3
	}
	return iret;
}
int ZHOU::Rate54_HiddenQuad(){// no naked quad in fact coded as naked 5
	int iret = 0;
	if (pm_go.opprint2 & 2){
		cout << "look for hidden quad rating 54 ( in fact naked 5)" << endl;
	}

	for (int iu = 0; iu < 27; iu++){// look in box
		BF128 wu = units3xBM[iu]; wu &= cells_unsolved;
		if (wu.Count() < 9) continue;
		byte *cig = cellsInGroup[iu];
		for (int i1 = 0; i1 < 5; i1++){// cells 0_5
			int c1=cig[i1],dig1 = zh_g.dig_cells[c1];
			for (int i2 = i1 + 1; i2 < 6; i2++){
				int c2 = cig[i2], dig2 = zh_g.dig_cells[c2] | dig1;
				if (_popcnt32(dig2) > 5)continue;
				for (int i3 = i2 + 1; i3 <7; i3++){
					int c3 = cig[i3], dig3 = zh_g.dig_cells[c3] | dig2;
					if (_popcnt32(dig3) > 5)continue;
					for (int i4 = i3 + 1; i4 < 8; i4++){
						int c4 = cig[i4], dig4 = zh_g.dig_cells[c4]  | dig3;
						if (_popcnt32(dig4) > 5)continue;
						for (int i5 = i4 + 1; i5 <9; i5++){
							int c5 = cig[i5], dig5 = zh_g.dig_cells[c5] | dig4;
							if (_popcnt32(dig5) != 5)continue;
							BF128 ww = wu; ww.Clear_c(c1); ww.Clear_c(c2);
							ww.Clear_c(c3); ww.Clear_c(c4); ww.Clear_c(c5);
							for (int idig = 0; idig < 9; idig++)if (dig5 & (1 << idig)){
								BF128 clean = FD[idig][0] & ww;
								if (clean.isNotEmpty()){
									FD[idig][0] -= clean;
									iret = 1;
								}
							}
							goto nextiu;// nothing more in this iu
						}
					}
				}
			}
		}
	nextiu:;
	}
	return iret;
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
	for (int i = 0, bit = 1; i < 9; i++, bit <<= 1) {
		if(bit & digits)ClearCandidate_c(i, cell);
	}
}
int ZHOU::CleanCellsForDigits(BF128 &  cells, int digits){
	int iret = 0;
	for (int i = 0, bit = 1; i < 9; i++, bit <<= 1) {
		if (bit & digits) {
			BF128 clean = FD[i][0] & cells;
			if (clean.isNotEmpty()) {
				iret = 1;
				FD[i][0] -= clean;
			}
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
	memset(Digit_cell_Assigned, 0, sizeof Digit_cell_Assigned);
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
				uint32_t  irow;
				bitscanforward(irow, row);
				int	cell = Tblstartblock[ib] + TblMult9[j] + irow;
				zs[cell] = (char)(digit + '1');
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
	if (zh_g.diag ) {
		cout << "opening situation after init"  << endl;
		Debug(1);
		ImageCandidats();
		if (zh_g.diag>2) return 0;
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
			Debug(1);
			ImageCandidats();
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
	int tp[30], np = 0;
	BitsInTable32(tp, np,R1);
	for(int i=0;i<np;i++){
		register int res= tp[i],R2= 1<< res ;
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
	int tp[60], np=0;
	BitsInTable64(tp, np, R1);
	for (int i = 0; i < np; i++) {
		register int res = tp[i], r3 = From_128_To_81[res];
		register uint64_t R2 = 1; R2 <<= res;
		if (R2 & FD[0][0].bf.u64[0]){ Assign(0, r3, res); goto nextr1; }
		if (R2 & FD[1][0].bf.u64[0]){ Assign(1, r3, res); goto nextr1; }
		if (R2 & FD[2][0].bf.u64[0]){ Assign(2, r3, res); goto nextr1; }
		if (R2 & FD[3][0].bf.u64[0]){ Assign(3, r3, res); goto nextr1; }
		if (R2 & FD[4][0].bf.u64[0]){ Assign(4, r3, res); goto nextr1; }
		if (R2 & FD[5][0].bf.u64[0]){ Assign(5, r3, res); goto nextr1; }
		if (R2 & FD[6][0].bf.u64[0]){ Assign(6, r3, res); goto nextr1; }
		if (R2 & FD[7][0].bf.u64[0]){ Assign(7, r3, res); goto nextr1; }
		if (R2 & FD[8][0].bf.u64[0]){ Assign(8, r3, res); goto nextr1; }
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
	if (wc.isEmpty())return;//
	int res = wc.getFirst128();
	// do digit update for the 2 digits
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
	uint32_t res;
	bitscanforward(res, hidden);
	ZHOU * mynext = this + 1; // start next guess
	mynext->Copy(*this);
	mynext->SetaCom(idig, res + dcell, res + dxcell);
	mynext->Upd1(idig);
	mynext->ComputeNext();
	bitscanreverse(res, hidden);
	SetaCom(idig, res + dcell, res + dxcell);
	Upd1(idig);
	ComputeNext();
	return 1;
}
int ZHOU::GuessHiddenTriplet(){// look for a triplet in row or box in band 3
	if (zh_g.diag)		cout << "guess hidden triplet in band 3" << endl;
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
	int tp[5], ntp=0;
	BitsInTable32(tp, ntp, hidden);
	for (int i = 0; i < ntp; i++) {
		int res = tp[i];
		if (i!=2) {//not last cell
			ZHOU * mynext = this + 1; // start next guess
			mynext->Copy(*this);
			mynext->SetaCom(idig, res + 54, res + 64);
			mynext->ComputeNext();
		}
		else {// this is the last do it in the same place
			SetaCom(idig, res + 54, res + 64);
			ComputeNext();
			return 1;
		}
	}
	return 1;// should never be reached
}


int ZHOU::PartialInitSudoku(GINT16 * t, int n){// if morph, done before
	memset(zh_g.Digit_cell_Assigned, 0, sizeof zh_g.Digit_cell_Assigned);
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
		/*int bandpairs = zh_g.pairs.bf.u32[iband];*/
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
		|| FD[0][0].bf.u32[2] != FD[0][1].bf.u32[2]){
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
	uint32_t res;
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
	{
	int tp[10], ntp=0;
	BitsInTable32(tp, ntp, min);
	for (int i = 0; i < ntp; i++) {// no bi value, use the smallest set
		int res = tp[i];
		ZHOU * mynext = this + 1; // start next guess
		mynext->Copy(*this);
		mynext->SetaCom(0, res + dcell, res + dxcell);
		mynext->Upd1(0);
		mynext->GuessFloor();
	}
	} // end of int tp[10], ntp=0 scope
	return;

exitok:
	//cout << Char27out(hidden) << " guess band dcell=" << dcell << endl;

	bitscanforward(res, hidden);
	ZHOU * mynext = this + 1; // start next guess
	mynext->Copy(*this);
	mynext->SetaCom(0, res + dcell, res + dxcell);
	mynext->Upd1(0);
	mynext->GuessFloor();
	bitscanreverse(res, hidden);
	//if (1)cout <<Char27out(1<<res)<< "second digit biv last dcell=" <<dcell<< endl;
	SetaCom(0, res + dcell, res + dxcell);
	Upd1(0);
	GuessFloor();
}
