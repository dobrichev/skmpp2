

struct GENSTEP{
	struct CLUE{
		int free,er,ec,eb,unlocked;	
		int level_dig;
	}tc[40];
	int row_free[9], col_free[9], box_free[9];	
	GINT16 tclues[40];
	int nclues,iclue;
	int tfree[40], ntcf;
	USHORT * tcor;
	char puz[82];
	void Init(){ 
		nclues = iclue= 0; 
		tc[0].free = 0777;
		for (int i = 0; i < 9; i++)row_free[i] = col_free[i] = box_free[i] = 0777;
	}
	void ReInit(int nfix){
		for (int i = 0; i < 9; i++)row_free[i] = col_free[i] = box_free[i] = 0777;
	}
	void SetClues(int n){
		for (int i = 0; i < nclues; i++){
			CLUE & c = tc[i];
			int cell = tclues[i].u8[0], digit = tclues[i].u8[1];
			CELL_FIX &cfd = cellsFixedData[cell];
			c.er = cfd.el;
			c.ec = cfd.pl;
			c.eb = cfd.eb;
			c.unlocked = 0777;
			if (i < n){
				int bit = 1 << digit;
				row_free[c.er] ^= bit;
				col_free[c.ec] ^= bit;
				box_free[c.eb] ^= bit;
			}
		}
	}
	void Addclue(int cell){
		CLUE & c = tc[nclues];
		tclues[nclues++].u16 = cell;
		CELL_FIX &cfd = cellsFixedData[cell];
		c.er = cfd.el;
		c.ec = cfd.pl;
		c.eb = cfd.eb;
		c.unlocked = 0777;
	}
	void AddClueClear123(int cell, int * boxes){
		CLUE & c = tc[nclues];
		Addclue(cell);
		int box = C_box[cell]; 
		for (int i = 0; i < 3;i++)
			if (boxes[i] == box){
				c.unlocked = 0770;
				return;
			}
	}
	inline void Assign(int bit){
		register CLUE & c = tc[iclue];
		c.free ^= bit;
		row_free[c.er] ^= bit;
		col_free[c.ec] ^= bit;
		box_free[c.eb] ^= bit;
	}
	inline void Restore(int bit){
		register CLUE & c = tc[iclue];
		row_free[c.er] ^= bit;
		col_free[c.ec] ^= bit;
		box_free[c.eb] ^= bit;
	}	
	inline void Restore(){
		int bit = 1 << tclues[iclue].u8[1];
		Restore(bit);
	}
	void Assign2(int icl, int d1, int d2){
		iclue = icl;
		Assign(1 << d1);
		tclues[iclue].u8[1] = (uint8_t)d1;
		iclue++;
		Assign(1 << d2);
		tclues[iclue].u8[1] = (uint8_t)d2;
	}
	void Restore2(int icl,int &d1,int & d2){// icl is first of 2
		iclue = icl + 1;
		Restore();
		d2 = tclues[iclue].u8[1];
		iclue = icl;
		d1 = tclues[iclue].u8[1];
		Restore();
	}
	int PuzzleToTest();
	void PuzzleToSplit();
	inline int Find_free_minimal_change(int lim=2,int diag=0){
		register CLUE & c = tc[iclue];
		int free = row_free[c.er] & col_free[c.ec] & box_free[c.eb] & c.unlocked;
		tc[iclue].free = free;
		if (diag){
			cout << "find free iclue=" << iclue << " er=" << c.er << "  ec=" << c.ec << " eb=" << c.eb
			 << " cell=" << (int)tclues[iclue].u8[0] << " digit=" << (int)tclues[iclue].u8[1] + 1 
			 <<"bfree=0"<<oct << free << dec << endl;
		}
		if ((int)_popcnt32(free) < lim) return 0;// locked or not minimal
		return 1;
	}
	void Gengo(int istart);
	void Debug(int n){
		cout << "gscom status for istart=" << n << endl;
		for (int i = 0; i < 9; i++) cout << oct << " 0" << row_free[i] << dec;
		cout << " row status"<<endl;
		for (int i = 0; i < 9; i++) cout << oct << " 0" << col_free[i] << dec;
		cout << " col status" << endl;
		for (int i = 0; i < 9; i++) cout << oct << " 0" << box_free[i] << dec;
		cout << " box status" << endl;
	}
	void PrintPartial(int n){
		char zs[82];
		strcpy(zs, empty_puzzle);
		for (int i = 0; i <= n; i++)
			zs[tclues[i].u8[0]] = tclues[i].u8[1] + '1';
		cout << zs << " n= "<<n << endl;
	}
	void PrintClues(){
		cout << "clues ";
		for (int i = 0; i < nclues; i++){
			cout << cellsFixedData[tclues[i].u8[0]].pt << " ";
		}
		cout << endl;
	}
	void GenSym36(int ntcf);
	void GenSym_Loop2(int sym36=1);
}gscom;


int GENSTEP::PuzzleToTest(){
	int digits = 0;
	for (int i = 0; i < nclues; i++) digits |= 1<<tclues[i].u8[1];
	if (_popcnt32(digits) < 8) return 0;// minimum 8 digits given to have a sudoku
	if (zh_g.Go_InitSolve(tclues, nclues))goto no;
	{
	int nguess = (int)zh_g.cpt[1];
	strcpy(puz, zh_g.puz);
	zh_g.zsol = 0; // be sure to keep the solution 
	if (!zhou[0].IsMinimale(tclues, nclues)) goto no;
	int ir = pm_go.SolveGetLow44(1);// pack the low ratings
	if (ir < 0) return 0; 
	if (ir) {
		fout1 << puz
			<< ";" << pm_go.rat_er << ";" << pm_go.rat_ep << ";" << pm_go.rat_ed
		<< endl; return 1;
	}
	if (pm_go.rat_ed>23){ fout3 << puz << ";" << nguess << endl; return 3; }
	fout2 << puz << ";" << nguess << endl;
	return 2;
	} //end of int nguess scope
no:
	if (sgo.command == 202)	{
		char puz[82];
		strcpy(puz, empty_puzzle);
		for (int i = 0; i < nclues; i++) puz[tclues[i].u8[0]] = tclues[i].u8[1] + '1';
		fout4 << puz << "seed" << endl;
	}
	return 0;
}
void GENSTEP::PuzzleToSplit(){
	char puz[82];
	strcpy(puz, empty_puzzle);
	for (int i = 0; i < nclues; i++) puz[tclues[i].u8[0]] = tclues[i].u8[1] + '1';
	zh_g.InitCount(0);
	if (zhou[0].InitSudoku(tclues, nclues))return;
	if (zhou[0].Isvalid() != 1)return;
	uint64_t nguess = zh_g.cpt[1];
	if (!nguess) { fout1 << puz << ";" << nguess << endl; return; }
	fout3 << puz << ";" << nguess << endl;
}

void GENSTEP::Gengo(int istart){
	int ilim = nclues - 1;
	iclue = istart;
	Find_free_minimal_change();
next:
	unsigned long iw;
	while ( tc[iclue].free){
		_BitScanForward(&iw, tc[iclue].free);
		//cout << "next iclue=" << iclue << " digit=" << iw+1 << endl;
		Assign(1 << iw);
		tclues[iclue].u8[1] = (uint8_t)iw;
		if (iclue < ilim){// next step find free
			iclue++;
			if(Find_free_minimal_change())			goto next;
			iclue--;
			goto loop;
		}
		PuzzleToTest();// this is one puzzle to check
		loop:
		Restore(1 << iw); // reset the digit as free in row column box
	}
	// back one clue if not closed
	if (iclue > istart){
		iclue--;
		iw = tclues[iclue].u8[1];
		Restore(1 << iw);
		goto next;
	}
}

void GENSTEP::GenSym36(int lim1){
	ntcf = lim1;
	if (!lim1){
		tc[0].unlocked = 9;// first cell digits 1;4
		GenSym_Loop2();
		return;
	}
	int tleveldig[4] = { 1, 3, 7, 7 };
	unsigned long digit;
	int  i1 = -1; // loop1 on fix clues
	tc[0].unlocked = 1;// first cell only digit 1
	tc[0].level_dig = 0;
nextindf:
	i1++;
	{	register CLUE & c = tc[i1];
	int free = row_free[c.er] & col_free[c.ec] & box_free[c.eb] & c.unlocked;
	tc[i1].free = tfree[i1] = free;
	}
	int digits = 0, level = 0;
	for (int ic = 0; ic < i1; ic++) digits |= 1 << tclues[ic].u8[1];
	if (i1) level = tc[i1 - 1].level_dig;
	if (digits & ~tleveldig[level])level++;
	tc[i1].level_dig = level;
	tfree[i1] &= tleveldig[level + 1];
nextf:
	while (tfree[i1]){
		_BitScanForward(&digit, tfree[i1]);
		tfree[i1] ^= 1 << digit;
		iclue = i1;
		tclues[i1].u8[1] = (uint8_t)digit;
		Assign(1 << digit);
		if (i1 < ntcf - 1) goto nextindf;
		PrintPartial(i1);
		GenSym_Loop2();

		iclue = i1;//be sure to have the good one
		Restore();// back while next loop1
	} 
// back one step in loop 1
	if (i1--){
		iclue = i1;  
		Restore();  
		goto nextf;
	}
}

void GENSTEP::GenSym_Loop2(int sym36){// loop on pairs of clues
	//if (1) return;
	// start with cell ntcf
	int tleveldig[6] = { 1,7, 0x1f, 0x7f, 0x1ff, 0x1ff };
	unsigned long digit;
	int  i = ntcf - 2, d1, d2;
	tc[0].level_dig = sym36;
	if (ntcf)tc[ntcf-1].level_dig = sym36;
nextind2:
	i += 2;
	{	register CLUE & c = tc[i];
		int free = row_free[c.er] & col_free[c.ec] & box_free[c.eb] & c.unlocked;
		if (_popcnt32(free) <2)  goto back2; // locked or not minimal
		tc[i].free = tfree[i]=free;
	}
	{	register CLUE & c = tc[i+1];
		int free = row_free[c.er] & col_free[c.ec] & box_free[c.eb] & c.unlocked;
		tc[i + 1].free = tfree[i+1] = free;
	}
	{
	int digits = 0,level=sym36;
	for (int ic = 0; ic < i; ic++) digits |= 1 << tclues[ic].u8[1];
	if (i) level = tc[i - 2].level_dig;
	if (digits & ~tleveldig[level])level++;
	if (i == nclues - 2){// must have reached level 2
		if (level < 2+sym36) goto back2;
	}
	tc[i].level_dig = level;
	tfree[i] &= tleveldig[level + 1];
	if (i == nclues - 2 && level==2+sym36){// must be 89
		tfree[i] &= 0x180;
	}
	} //end of int digits scope
next2:
	while ( tfree[i]){
		_BitScanForward(&digit, tfree[i]);
		//if (i == 1)cout << "loop1 0x" << hex << tfree[i]<<dec << endl;
		tfree[i] ^= 1 << digit;
		d2 = tcor[digit];
		d1 = digit;
		tfree[i] &=~( 1 << d2);// erase also pair digit

		if ((tc[i].free &(1 << d1)) && (tc[i+1].free &(1 << d2))){
			//cout << "go1 i="
			Assign2(i, d1, d2);
			if (i < (nclues - 2)) goto nextind2;
			else PuzzleToTest();
			Restore2(i, d1, d2);
		}
		if ((tc[i].free &(1 << d2)) && (tc[i+1].free &(1 << d1))
			&& d1!=d2){
			Assign2(i, d2, d1);
			if (i < (nclues - 2)) goto nextind2;
			else PuzzleToTest();
			Restore2(i, d1, d2);
		}
	}
back2:
	if (i <= ntcf) return;
	i -= 2;
	Restore2(i, d1, d2);
	int lvl = tc[i].level_dig, bit = 1 << tclues[i].u8[1];
	if (~tleveldig[lvl] & bit) goto back2;// redundancy filter
	goto next2;
}

/*
case 200: Go_c200(); break;// split the entry file in files 1;2;3
case 201: Go_c201(); break;// change n clues or 1 to n clues
case 202: Go_c202(); break;// gen symmetry of given
case 210:  Go_c210(); break;// create a seed on a pattern
case 211:  Go_c211(); break;// create a seed gen interim file
case 212:  Go_c212(); break;// create a seed gen interim file
case 217:  Go_c217(); break;// restore part of a data base
case 218:  Go_c218(); break;// extract played seeds
case 219:  Go_c219(); break;// restore a data base
*/

void Go_c200(){// just split the entry file 
	if (!sgo.foutput_name) return;	
	memset(zh_g.cptg, 0, sizeof zh_g.cptg);
	zh_g.npuz = 0;
	zh_g.diag = (int)sgo.vx[9];
	if (!sgo.finput_name) return;
	finput.open(sgo.finput_name);
	if (!finput.is_open()){ cerr << "error open " << sgo.finput_name << endl; return; }
	char ze[82]; ze[81] = 0;
	while (finput.GetPuzzle(ze)){
		if (zh_g.diag)cout << ze << "trait�" << endl;
		zh_g.npuz++;
		gscom.Init();
		for (int i = 0; i < 81; i++)if (ze[i] != '.'){// catch given
			int c = ze[i] - '1';
			if (c < 0 || c>9){ cerr << "invalid file" << endl; return; }
			gscom.tclues[gscom.nclues++].u16 = (c << 8) | i;
		}
		gscom.PuzzleToSplit();
		for (int i = 0; i < 10; i++) zh_g.cptg[i] += zh_g.cpt[i];
		//if (zh_g.cptg[7]>20)break;
		//if (zh_g.npuz>1000)break;
	}
	cout << "summary npuz=" << zh_g.npuz << endl;
	for (int i = 0; i < 10; i++) if (zh_g.cptg[i])
		cout << zh_g_cpt[i] << "\t" << zh_g.cptg[i] << endl;
}
void Go_c201(){
	if (!sgo.foutput_name){
		cerr << "missing output name" << endl; return;
	}
	int cpt = 0;
	if (!sgo.finput_name) return;
	cout << "c201 entry for input " << sgo.finput_name << " mode upto=" << sgo.vx[1] << endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()){ cerr << "error open " << sgo.finput_name << endl; return; }
	char ze[82]; ze[81] = 0;
	GINT16 myclues[40];
	int myn;
	uint32_t change = sgo.vx[0], idep = change, nfix = sgo.vx[2];
	if (sgo.vx[1])idep = 1;

	while (finput.GetPuzzle(ze)){
		cout << ze << " to process" << endl;
		myn = 0;
		gscom.Init();
		for (int i = 0; i < 81; i++)if (ze[i] != '.'){// catch given
			int c = ze[i] - '1';
			if (c < 0 || c>9){ cerr << "invalid file" << endl; return; }
			myclues[myn++].u16 = (c << 8) | i;
		}
		if (myn <=(int) (change + nfix)){
			cout << "cancelled due to nfix+change too high" << endl;
			return;
		}

		gscom.nclues = myn;
		COMBINE combi;
		USHORT tclues[40], tclues_s[40];
//		int istart = myn - change;
		if (nfix)	cout << "start combi nfix=" << nfix << " change="<<change<< endl;
		
		// skip the fix clues here and load the nfix in gscom
		for (int i = 0; i <(int) nfix; i++)gscom.tclues[i] = myclues[i];
		for (int i = 0; i < (int)(myn-nfix); i++)tclues[i] = i;
		for (uint32_t ichange = idep; ichange <= change; ichange++){
			int istart = myn - ichange;
			combi.First(myn - nfix, ichange, tclues, tclues_s);
			while (1){
				for (int i = nfix; i < myn; i++)
					gscom.tclues[i] = myclues[tclues_s[i-nfix]+nfix];
				//  process the 
				if (nfix)gscom.ReInit(nfix);
				gscom.SetClues(istart);
				// lock known digit in cells to change
				for (int i = istart; i < myn; i++){
					int digit = gscom.tclues[i].u8[1];
					gscom.tc[i].unlocked = 0777 ^ (1 << digit);
				}
				if (0 &&nfix){
					char zs[82];
					strcpy(zs, empty_puzzle);
					for (int i = 0; i < istart; i++)
						zs[gscom.tclues[i].u8[0]] = gscom.tclues[i].u8[1] + '1';
					for (int i = istart; i < myn; i++)
						zs[gscom.tclues[i].u8[0]] = 'x';


					cout <<zs<< "call gengo " << endl;
				}
				gscom.Gengo(istart);
				if (!combi.Next()) break;
			}
		}
	}
}


/* tables used in symmetry processing
[36][2] sd1_36  diagonal  sd2_36 diag 2   sst_36 stick
sym_81[5][81] direct sym   d1  d2 st cc r90
sym_f9[3][9] fix positions  d1 d2 stick
sym_tcor[1][9]  pairing 34 56 78  1;2;3 self
*/

void Go_c202(){
	if (!sgo.foutput_name){	cerr << "missing output name" << endl; return;	}
	if (!sgo.finput_name) return;
	cout << "c202 symmetry entry for input " << sgo.finput_name << endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()){ cerr << "error open " << sgo.finput_name << endl; return; }
	char ze[82]; ze[81] = 0;
	if (sgo.foutput_name){
		char zn[200];
		strcpy(zn, sgo.foutput_name);
		int ll = (int)strlen(zn);
		strcpy(&zn[ll], "_file4.txt");
		fout4.open(zn);
	}
	if(finput.GetPuzzle(ze)){
		int 	nclues = 0, cclue = 0;
		cout << ze << " pattern to process" << endl;
		// is it central
		gscom.Init();
		gscom.tcor = sym_tcor[0];
		if (ze[40] - '.'){// central cell s a given
			gscom.Addclue(40);
			cclue = 1;
			gscom.Assign(1);// assign digit 0 to central cell as given
		}
		for (int i = 0; i<40; i++) {
			USHORT * vv = sc_40[i];
			int ij = vv[0], ji = vv[1];
			if (ze[ij] - '.') {// un clue
				if (ze[ji] == '.') goto notcentral; // not the expected pattern
				gscom.Addclue(ji);
				gscom.Addclue(ij);
			}
			else if (ze[ji] - '.') goto notcentral; //not the expected pattern
		}
		gscom.ntcf=cclue;
		cout << "study central nclues=" << nclues << " cclue=" << cclue << endl;
		{// process central symmetry of given 
			gscom.GenSym_Loop2(0);
			cout << "exit central" << endl;
		}
	notcentral:
		//return;
		//        diagonal or stick 9 self symmetric + 36 x 2
		int boxes_with_fix[3][3] = { 0, 4, 8, 2, 4, 6, 3, 4, 5 };
		char * ttyp[3] = { "diag1", "diag2", "stick" };
		gscom.tcor = sym_tcor[1];
		for (int itype = 0; itype < 3; itype++){// diagonal or stock
			cout << " entry sym_36 type= " << itype << endl;
			nclues = 0;
			gscom.Init();
			int  ntcf = 0;  // find clues on fix cells
			USHORT *vf = sym_f9[itype];
			for (int i = 0; i<9; i++)   {
				int ii = vf[i]; if (ze[ii] - '.') gscom.Addclue(ii);
			}
			ntcf = gscom.nclues;
			USHORT(*p_36_2)[2] = (!itype) ? sd1_36 : (itype<2) ? sd2_36 : sst_36;
			for (int i = 0; i<36; i++) {
				USHORT * vv = p_36_2[i];
				int ij = vv[0], ji = vv[1];
				if (ze[ij] - '.') {// un clue
					if (ze[ji] == '.') goto exit_itype; // not the expected pattern
					gscom.AddClueClear123(ij, boxes_with_fix[itype]);
					gscom.AddClueClear123(ji, boxes_with_fix[itype]);
					// lock digits 0;1;2 if necessary (relevant boxes)

				}
				else if (ze[ji] - '.') goto exit_itype; //not the expected pattern
			}
			nclues = gscom.nclues;
			cout << "seen symmetry " << ttyp[itype] << " ndiag="<< ntcf << endl;
			gscom.GenSym36(ntcf);
		exit_itype:;
		}// endl loop diagonal; stick
	}
}


// 210;211;212 seed on a pattern
// 210 in once skip then one ever xx; one seed per 'n' first given
// 211 create a primary file size 'n' other set to 1
// 212 create one seed per entry if possible given 'n' 
// 217 restore part of a database given start string
// 218 extract top clues and count print small number of clues in file 1 
// 219 restore all the data base (usually after 218)
void Go_c210(){// create a seed file on a pattern
	if (!sgo.foutput_name){
		cerr << "missing output name" << endl; return;
	}
	int cpt = 0, cpt2 = 0;
	if (!sgo.finput_name) return;
	cout << "c203 entry create a seed file for input " << sgo.finput_name << endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()){ cerr << "error open " << sgo.finput_name << endl; return; }
	char ze[82]; ze[81] = 0;
	int nclues, cclue = 0, n1 = sgo.vx[0];
	if (n1 < 4 || n1>14)n1 = 9;
	n1--; //set n1 to test limit
	unsigned long digit;
	int tfree[40], tmaxdig[40];

	if (finput.GetPuzzle(ze)){
		cout << ze << " pattern to process" << endl;
		gscom.Init();
		for (int i = 0; i<81; i++) if (ze[i] - '.') gscom.Addclue(i);
		nclues = gscom.nclues;
		cout << "nclues=" << nclues << endl;
		if (nclues < 17 || nclues>30) return;// not a good range
		// set a digit limit for the first 8 clues
		int bit = 1, n2 = nclues - 1;
		for (int i = 0; i < 8; i++){
			gscom.tc[i].unlocked = bit;
			bit = (bit << 1) | 1;
		}
		{
			int i = -1;
			tmaxdig[0] = 0;
		nextindc:
			i++;
			gscom.iclue = i;
			gscom.Find_free_minimal_change();
			if (i)tmaxdig[i] = tmaxdig[i - 1];
			tfree[i] = gscom.tc[i].free;
		nextc:
			while ( tfree[i]) {
				_BitScanForward(&digit, tfree[i]);
				tfree[i] ^= 1 << digit;
				gscom.iclue = i;
				gscom.tclues[i].u8[1] = (uint8_t)digit;
				gscom.Assign(1 << digit);
				if ((int)digit>tmaxdig[i])tmaxdig[i] = digit;
				if (i < n1) goto nextindc;
				char zs[82];
				if (sgo.vx[1]){ sgo.vx[1] --; goto exitloop2; }// skip
				sgo.vx[1] = sgo.vx[3];
				if (!sgo.vx[2]--) return;
				if (1){// print to control first loop
					strcpy(zs, empty_puzzle);
					for (int j = 0; j <= n1; j++){
						GINT16 w = gscom.tclues[j];
						int cell = w.u8[0], digit = w.u8[1] + '1';
						zs[cell] = digit;
					}
					cout << zs << "go  reste " << sgo.vx[2] << endl;
				}
				//if (1) goto exitloop2;
				// start loop 2 one valid gscom.PuzzleToTest();
				// expand other digit till a valid grid is seen
			nextind2:
				i++;
				gscom.iclue = i;
				gscom.Find_free_minimal_change();
				tmaxdig[i] = tmaxdig[i - 1];
				tfree[i] = gscom.tc[i].free;
				if (!tfree[i]) goto back2; // locked
			next2:
				int tp2[10], ntp2;
				BitsInTable32(tp2, ntp2, tfree[i]);
				while (tfree[i]){
					_BitScanForward(&digit, tfree[i]);
					tfree[i] ^= 1 << digit;
					gscom.iclue = i;
					gscom.tclues[i].u8[1] = (uint8_t)digit;
					gscom.Assign(1 << digit);
					if ((int)digit>tmaxdig[i])tmaxdig[i] = digit;
					if (i < n2) goto nextind2;
					if (gscom.PuzzleToTest()){// good puzzle back to next start
						while (1){//back to next 1
							gscom.Restore();
							gscom.iclue = --i;
							if (i == n1) goto exitloop2;
						}// end of a good puzzle

					}// end of while loop2 not a good puzzle
					if (0){
						for (int ic = n1; ic < nclues; ic++){
							GINT16 w = gscom.tclues[ic];
							int cell = w.u8[0], digit = w.u8[1] + '1';
							zs[cell] = digit;
						}
						if (!(++cpt & 1023)){
							cout << zs << " pas ok cpt/1024=" << (cpt >> 10) << endl;
						}

					}
					gscom.Restore();
					if ((int)digit > tmaxdig[i - 1]) break;
				}
			back2:
				if (--i>n1){
					gscom.iclue = i;
					gscom.Restore();
					digit = gscom.tclues[i].u8[1];
					if ((int)digit > tmaxdig[i - 1]) goto back2;
					goto next2;
				}
				// unlikely no good puzzle found
			exitloop2:
				//cout << "exitloop2 i=" << i << " free=0" << oct << tfree[i] <<dec << endl;
				gscom.iclue = i;//be sure to have the good one
				gscom.Restore();// back while next loop1 
				digit = gscom.tclues[i].u8[1];
				if ((int)digit > tmaxdig[i - 1]) break;

			}// end while next loop1
		backc: // back one step in loop 1
			if (i--){
				gscom.iclue = i;
				gscom.Restore();
				digit = gscom.tclues[i].u8[1];
				if ((int)digit > tmaxdig[i - 1]) goto backc;
				goto nextc;
			}
		}// end loop1

	} // end if get puzzle
}

void Go_c211(){// create a seed file on a pattern
	if (!sgo.foutput_name){
		cerr << "missing output name" << endl; return;
	}
	int cpt = 0, cpt2 = 0;
	if (!sgo.finput_name) return;
	cout << "c203 entry create a seed file for input " << sgo.finput_name << endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()){ cerr << "error open " << sgo.finput_name << endl; return; }
	char ze[82]; ze[81] = 0;
	int nclues, cclue = 0, n1 = sgo.vx[0];
	if (n1 < 4 || n1>14)n1 = 9;
	n1--; //set n1 to test limit
	unsigned long digit;
	int tfree[40], tmaxdig[40];

	if (finput.GetPuzzle(ze)){
		cout << ze << " pattern to process" << endl;
		gscom.Init();
		for (int i = 0; i<81; i++) if (ze[i] - '.') gscom.Addclue(i);
		nclues = gscom.nclues;
		cout << "nclues=" << nclues << endl;
		if (nclues < 17 || nclues>30) return;// not a good range
		// set a digit limit for the first 8 clues
		int bit = 1, n2 = nclues - 1;
		for (int i = 0; i < 8; i++){
			gscom.tc[i].unlocked = bit;
			bit = (bit << 1) | 1;
		}
		{
			int i = -1;
			tmaxdig[0] = 0;
		nextindc:
			i++;
			gscom.iclue = i;
			gscom.Find_free_minimal_change();
			if (i)tmaxdig[i] = tmaxdig[i - 1];
			tfree[i] = gscom.tc[i].free;
		nextc:
			while ( tfree[i]){
				_BitScanForward(&digit, tfree[i]);
				tfree[i] ^= 1 << digit;
				gscom.iclue = i;
				gscom.tclues[i].u8[1] = (uint8_t)digit;
				gscom.Assign(1 << digit);
				if ((int)digit>tmaxdig[i])tmaxdig[i] = digit;
				if (i < n1) goto nextindc;
				char zs[82];
				strcpy(zs, empty_puzzle);
				for (int j = 0; j <nclues; j++){// digit 0 if not assigned
					GINT16 w = gscom.tclues[j];
					int cell = w.u8[0], digit = w.u8[1] + '1';
					zs[cell] = digit;
				}
				fout1 << zs << ";0;0"<< endl;
				gscom.Restore();// back while next loop1 
				digit = gscom.tclues[i].u8[1];
				if ((int)digit > tmaxdig[i - 1]) break;

			}// end while next loop1
		backc: // back one step in loop 1
			if (i--){
				gscom.iclue = i;
				gscom.Restore();
				digit = gscom.tclues[i].u8[1];
				if ((int)digit > tmaxdig[i - 1]) goto backc;
				goto nextc;
			}
		}// end loop1

	} // end if get puzzle
}
void Go_c212(){// create a seed file on a pattern based on primary status
	if (!sgo.foutput_name){
		cerr << "missing output name" << endl; return;
	}
	int cpt = 0, cpt2 = 0;
	if (!sgo.finput_name) return;
	cout << "c203 entry create a seed file for input " << sgo.finput_name << endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()){ cerr << "error open " << sgo.finput_name << endl; return; }
	char ze[82], zs[82]; ze[81] = 0;
	int tclues[40],nclues, cclue = 0, n1 = sgo.vx[0];
	if (n1 < 4 || n1>14){
		cout << n1 << " given clues not an accepted value" << endl;
		return;
	}
	unsigned long digit;
	int tfree[40];

	while (finput.GetPuzzle(ze)){
		cout << ze << " pattern to process" << endl;
		gscom.Init();
		nclues = 0;
		for (int i = 0; i < 81; i++) if (ze[i] - '.'){
			gscom.Addclue(i);
			tclues[nclues++] = ze[i] - '1';
		}
		if (nclues < 17 || nclues>30) continue;// not a good range
		// assign fixclues
		for (int i = 0; i < n1; i++){
			gscom.iclue = i;
			digit = tclues[i];
			gscom.tclues[i].u8[1] = (uint8_t)digit;
			gscom.Assign(1 << digit);
		}
		if (1){// print to control first loop
			strcpy(zs, empty_puzzle);
			for (int j = 0; j < nclues; j++){
				GINT16 w = gscom.tclues[j];
				int cell = w.u8[0], digit = w.u8[1] + '1';
				zs[cell] = digit;
			}
			cout << zs << "go  reste " << nclues-n1 << endl;
		}
		{
			int i = n1 - 1;
			// expand other digit till a valid grid is seen
		nextind2:
			i++;
			gscom.iclue = i;
			gscom.Find_free_minimal_change();
			tfree[i] = gscom.tc[i].free;
			if (!tfree[i]) goto back2; // locked
		next2:
			while (tfree[i]){
				_BitScanForward(&digit, tfree[i]);
				tfree[i] ^= 1 << digit;
				gscom.iclue = i;
				gscom.tclues[i].u8[1] = (uint8_t)digit;
				gscom.Assign(1 << digit);
				if (i < nclues-1) goto nextind2;
				if (gscom.PuzzleToTest()) goto nextentry; // good puzzle back to next start
				gscom.Restore();
			}// end of while loop2 not a good puzzle
			
		back2:
			if (--i >=n1){
				gscom.iclue = i;
				gscom.Restore();
				goto next2;
			}
		nextentry:;
		}

	} // end if get puzzle
}

void Go_c217(){// restore part of a database on  a pattern
	if (!sgo.s_strings[0]){
		cerr << "missins gbase input name" << endl; return;
	}
	if (!sgo.finput_name) return;
	cout << "c219 restore a data base file input " << sgo.finput_name
		<< " file to restore " << sgo.s_strings[0] << endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()){ cerr << "error open " << sgo.finput_name << endl; return; }
	FINPUT input2;
	input2.open(sgo.s_strings[0]);
	if (!input2.is_open()){ cerr << "error open 2" << sgo.foutput_name << endl; return; }

	char ze[82], *ze2 = input2.ze;
	ze[81] = 0;

	int tclues[40], nclues = 0;
	if (finput.GetPuzzle(ze)){
		cout << ze << " pattern to process row 1 known" << endl;
		for (int i = 9; i < 81; i++) if (ze[i] - '.') 		tclues[nclues++] = i;
		char * scomp = sgo.s_strings[1];
		if (!scomp) {
			cout << "missing seclection string" << endl;
			return;
		}
		int lcomp = (int)strlen(scomp);
		if (lcomp > nclues){
			cout << "seclection string too long" << endl;
			return;
		}
		uint64_t nlim = sgo.vx[1], npuz = 0;;
		while (input2.GetLigne()){
			if (strlen(ze2) != nclues)continue;// not the right file
			int ir = strncmp(ze2, scomp, lcomp);
			if (ir > 0) return;
			npuz++;
			if (!ir){
				for (int j = 0; j < nclues; j++)	ze[tclues[j]] = ze2[j];
				fout1 << ze << endl;
			}
			if (nlim && npuz == nlim) return;
		}

	} // end if get puzzle
}

void Go_c218(){// get known seeds on a pattern
	if (!sgo.s_strings[0]){
		cerr << "missins gbase input name" << endl; return;
	}
	int cpt = 0, cpt2=0;
	if (!sgo.finput_name) return;
	cout << "c218 find known seeds input " << sgo.finput_name << endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()){ cerr << "error open " << sgo.finput_name << endl; return; }
	FINPUT input2;
	input2.open(sgo.s_strings[0]);
	if (!input2.is_open()){ cerr << "error open 2" << sgo.foutput_name << endl; return; }

	char ze[82], *ze2 = input2.ze, ze2comp[40];
	char store[10][40],wstore[40];
	ze[81] = 0; ze2comp[0] = 0;

	int tclues[40],nclues=0,  n1 = sgo.vx[0];
	if (n1 < 3 || n1>10){	cout << "nb clues invalide" << endl; return;}
	if (finput.GetPuzzle(ze)){
		cout << ze << " pattern to process row 1 known" << endl;
		for (int i = 9; i < 81; i++) if (ze[i] - '.') {			
			if (nclues<n1)tclues[nclues++] = i;
			else ze[i] = '1';// unused clues set to '1'
		}
		cout << ze << " pattern unknown set to 1" << endl;
		while (input2.GetLigne()){
			//if (++cpt2 < 10) cout << ze2 << " base" << endl;
			if (strlen(ze2) < n1)return;// not the right file
			// store the 10 first in case
			strcpy(wstore, ze2);
			ze2[n1]	 = 0;// cut entry to length used
			if (!strcmp(ze2, ze2comp)){// =  store the 10 first in cse
				if (cpt < 10)strcpy(store[cpt], wstore);
				cpt++; continue;
			}
			if (cpt){// not the first set output to ze2
				for (int j = 0; j < nclues; j++)	ze[tclues[j]] = ze2comp[j];
				cout << ze << ";" << sgo.vx[1] << ";" << cpt << endl;
				if (cpt <= (int)sgo.vx[2]){
					//cout << "write " << cpt << " records in fout1" << endl;
					for (int j = 0; j < cpt; j++)
						fout1 << store[j] << endl;
				}
			}
			strcpy(ze2comp, ze2);
			strcpy(store[0], wstore);
			cpt = 1;
		}
		if (cpt){// output last
			for (int j = 0; j < nclues; j++)	ze[tclues[j]] = ze2comp[j];
			cout << ze << ";" << sgo.vx[1] << ";" << cpt << endl;
			if (cpt <= (int)sgo.vx[2]){
				for (int j = 0; j < cpt; j++)
					fout1 << store[j] << endl;
			}
		}

	} // end if get puzzle
}

void Go_c219(){// restore a database on  a pattern
	if (!sgo.s_strings[0]){
		cerr << "missins gbase input name" << endl; return;
	}
	if (!sgo.finput_name) return;
	cout << "c219 restore a data base file input " << sgo.finput_name
		<< " file to restore " << sgo.s_strings[0] << endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()){ cerr << "error open " << sgo.finput_name << endl; return; }
	FINPUT input2;
	input2.open(sgo.s_strings[0]);
	if (!input2.is_open()){ cerr << "error open 2" << sgo.s_strings[0] << endl; return; }

	char ze[82], *ze2 = input2.ze;
	ze[81] = 0;

	int tclues[40], nclues = 0;
	if (finput.GetPuzzle(ze)){
		for (int i = 9; i < 81; i++) if (ze[i] - '.') 		tclues[nclues++] = i;
		uint64_t nlim = sgo.vx[1], npuz = 0;;
		int gfilter = sgo.vx[0];
		cout << ze << " pattern to process row 1 known" << endl;
		if (gfilter) cout << "only for guesses over " << gfilter << endl;
		if (nlim)cout << "number of entries processed " <<nlim << endl;
		while (input2.GetLigne()){
			if (strlen(ze2) != nclues)continue;// not the right file
			npuz++;
			for (int j = 0; j < nclues; j++)	ze[tclues[j]] = ze2[j];
			if (gfilter){
				if (!zh_g.Go_InitSolve(ze)){
					int nguess = (int)zh_g.cpt[1];
					if (nguess>=gfilter)
						fout1 << ze << ";" << nguess << endl;
					else cout << ze << "ignored" << endl;
				}
				else cout << ze << "not valid" << endl;

			}
			else fout1 << ze  << endl;
			if (nlim && npuz == nlim) return;
		}

	} // end if get puzzle
}
