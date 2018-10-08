extern PM_GO pm_go;

void Go_c0(){
	if (!sgo.finput_name) return;
	long long cptg[20];
	memset(cptg, 0, sizeof cptg);
	cout << "Go_0 entry " << sgo.finput_name << " input" << endl;
	// temporary code extract solution with a known 17 6 6 5
	finput.open(sgo.finput_name);
	if (!finput.is_open()){
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	char ze[82]; ze[81] = 0;
	uint32_t npuz = 0;
	int nn[3] = { 0, 0, 0 };
	long ttdeb = GetTimeMillis();
	while (finput.GetPuzzle(ze)){
		npuz++;
		if (npuz < sgo.vx[0]) continue;
		if (sgo.vx[8])cout << ze << "to process npuz=" << npuz << endl;
		long tdeb = GetTimeMillis();
#ifdef ZHOU_OLD
		zhou[0].glb.diag = sgo.vx[9];
		for (uint32_t i = 0; i< sgo.vx[2]; i++){
			int ir = (zhou[0].CheckValidityQuick(ze));
			nn[ir]++;
			//			if (!ir && !i) cout << finput.ze << " invalide" << endl;
		}
		if (sgo.vx[8]){
			cout << "loop sommaire nsol=" << zhou[0].glb.nsol << endl;
			for (int i = 0; i < 10; i++) if (zhou[0].glb.cpt[i])
				cout << zh_g_cpt[i] << "\t" << zhou[0].glb.cpt[i] << endl;
			long tfin = GetTimeMillis();
			cout << "puz loop t=" << tfin - tdeb << endl;

		}

		for (int i = 0; i < 10; i++) cptg[i] += zhou[0].glb.cpt[i];
		if (npuz >= sgo.vx[1]) break;

#else
		zh_g.diag = sgo.vx[9];
		for (uint32_t loop = 0; loop < sgo.vx[2]; loop++){
			zh_g.InitCount(1);
			if (0){
				zh_g.Go_InitSudoku(ze);
				if (sgo.vx[8] && zh_g.diag)cout << zh_g.puz << "morphed" << endl << endl;
				//zh_g.Debug();
				zhou[0].ComputeNext();
			}
			else{
				zhou[0].CheckValidityQuick(ze);
			}
			if (!sgo.vx[8]) nn[zh_g.nsol]++;

			if (1){
				GINT16 tgiven[40];
				int ngiven = 0;
				for (int i = 0; i < 81; i++){
					register int c = ze[i];
					if (c<'1' || c>'9') continue;
					c -= '1';
					tgiven[ngiven++].u16 =(uint16_t)( i | (c << 8));
				}
				cout << "check minimale=" << zhou[0].IsMinimale(tgiven,ngiven) << endl;
				if (1) return;
			}
		}
		if (sgo.vx[8]){
			cout << "loop sommaire nsol=" << zh_g.nsol << endl;
			for (int i = 0; i < 10; i++) if (zh_g.cpt[i])
				cout << zh_g_cpt[i] << "\t" << zh_g.cpt[i] << endl;
			long tfin = GetTimeMillis();
			cout << "puz loop t=" << tfin - tdeb << endl;

		}
		for (int i = 0; i < 10; i++) cptg[i] += zh_g.cpt[i];
		if (npuz >= sgo.vx[1]) break;

#endif
		if (1) return;
	}
	if (1) {
		cout << "compte nsol 0;1;2 \t"
			<< nn[0] << "\t" << nn[1] << "\t" << nn[2] << endl;
		long ttfin = GetTimeMillis();
		cout << "t=" << ttfin - ttdeb << endl;
		for (int i = 0; i < 10; i++) if (cptg[i])
			cout << zh_g_cpt[i] << "\t" << cptg[i] << endl;

	}

}

int TWO_DIGITS::Hiden_Pairs_Box(){// etude recherche hidden biv
	nhp = 0;
	for (int iband = 0; iband < 3; iband++){
		int band = bf2.bf.u32[iband];
		if (!iband) continue;
		for (int ibox = 0; ibox < 3; ibox++){
			int box = band & tband_box[ibox];
			if (_popcnt32(box) - 2) continue;
			// if ! naked pair store it  set up naked pair 
			int locked = box & zh_g.pairs.bf.u32[ibox];
			if (box & (~locked)){// some eliminations 
				zh_g.pairs.bf.u32[ibox] |= box;
			}
		}
	}
	return nhp;
}

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
