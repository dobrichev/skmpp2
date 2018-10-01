
/*
442 dat file to txt file

..\release\skmpp    -i"extractclues"  -c"446"   -v24         "//split =idpos clues"
 #  19    1.2/1.2/1.2 - gsf
 1.2/1.2/1.2 - gsf
*/

//____________________________________________________________

/*	400 small tasks
		0x add 0,sequence 1,string0 2 nclues 
			9 std puzzle
			1x replace 0,empty'.' 1,erase'"'  2 cut to vx[1]
	401 .dat to .txt
	402 morph a puzzle

*/
void Go_c400_40p(char *ze, int task) {// work on bitfields std puz
	BF128 pattern; pattern.SetAll_0();
	int d_bands[6], d_units[27],digits=0;
	int c_bands[6], c_units[27];
	memset(d_bands, 0, sizeof d_bands);
	memset(d_units, 0, sizeof d_units);
	memset(c_bands, 0, sizeof c_bands);
	memset(c_units, 0, sizeof c_units);
	const char* rcb  = "RCB";
	for (int i = 0; i < 81; i++)if (ze[i] != '.') {
		pattern.Set_c(i);
		int d = ze[i] - '1',bit=1<<d;
		if (d < 0 || d>8) return;// invalid puzzle
		CELL_FIX & cf = cellsFixedData[i];
		int band = cf.el / 3, stack = cf.pl / 3 + 3;
		d_bands[band] |= bit; d_bands[stack] |= bit; digits |= bit;
		d_units[cf.el] |= bit; d_units[cf.plu] |= bit; d_units[cf.ebu] |= bit;
		c_bands[band]++; c_bands[stack]++;
		c_units[cf.el]++; c_units[cf.plu]++; c_units[cf.ebu]++;
	}
	ze[81] = 0;
	fout1 << ze;
	switch (task) {
	case 40: // count digits
		fout1 << ";" << _popcnt32(digits);
		break;
	case 41: // count given per band
		for (int i = 0; i < 6; i++) fout1 << ";" << c_bands[i];
		break;
	case 42: // count digits per band
		for (int i = 0; i < 6; i++) fout1 << ";" << _popcnt32(d_bands[i]);
		break;
	case 43: // count given per unit
		for (int i1 = 0, iu = 0; i1 < 3; i1++) {
			fout1 << ";" << rcb[i1];
			for (int i = 0; i < 9; i++,iu++) fout1 << ";" << c_units[iu];
		}
		break;
	case 44: // count digits per unit
		for (int i1 = 0, iu = 0; i1 < 3; i1++) {
			fout1 << ";" << rcb[i1];
			for (int i = 0; i < 9; i++, iu++) fout1 << ";" << _popcnt32(d_units[iu]);
		}
		break;
	}
	fout1 << endl;
}
void Go_c400() {// small tasks on entry -v0- is the task
	int task = sgo.vx[0];// parameter defining the small task
	cout << "Go_400 entry " << sgo.finput_name 
		<< " small tasks on entry subtask "<<task << endl;
	uint64_t npuz = 0;
	char * ze = finput.ze;
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	while (finput.GetLigne()) {
		int ll = (int)strlen(ze);
		if (ll < 81 || ll > 250) continue;
		if (task > 40 && task<50) {	Go_c400_40p(ze, task); continue;}
		switch (task) {

		//=========== add to ze

		case 0:// add sequence
			fout1 << ze << ";" << npuz++ << endl;			break;
		case 1:// add string 0
			fout1 << ze << ";" << sgo.s_strings[0] << endl; break;
		case 2:// add nclues 
			{	int nc = 0;
				for (int i = 0; i < 81; i++)
					if (ze[i]>= '1' && ze[i]<= '9') nc++;
				fout1 << ze << ";" << nc << endl; break;
			}
		case 9:// add std puzzle to entry out puz + std puz
			{
				char zs[82];	zs[81] = 0;
				if (ll < 81) break;
				ze[81] = 0;
				for (int i = 0; i < 81; i++) {
					char c = ze[i];
					if (c > '0' && c <= '9') zs[i] = c;
					else zs[i] = '.';
				}
				fout1 << ze << ";" << zs << endl;
				break;
		}

		//===================== replace 

		case 10://'.' if not '0' '9' in puzzle area 
			{	if (ll < 81) break;
				for (int i = 0; i < 81; i++) {
					int c = ze[i];
					if (c < '1' || c > '9')ze[i] = '.';
				}
				fout1 << ze  << endl; break;
			}
		case 11://erase " in the entry
			for (int i = 0; i < ll; i++) if (ze[i] - '"') fout1 << ze[i];
			fout1 << endl; break;
		case 12:// cut entry to sgo.vx[1]
			if (ll < (int)sgo.vx[1]) break;
			ze[sgo.vx[1]] = 0;
			fout1 << ze << endl; break;
		case 15:// mantext in output
			{	int known = 0;
				char vec[9], vi = '9';
				for (int i = 0; i < 81; i++)if (ze[i] - '.') {
					int c = ze[i] - '1', bit = 1 << c;;
					if (! (known&bit)) {
						vec[c] = vi--;
						known|=bit;
					}
					ze[i] = vec[c];
				}
				break;
			}
		case 16:// mintext in output
		{	int known = 0;
			char vec[9], vi = '1';
			for (int i = 0; i < 81; i++)if (ze[i] - '.') {
				int c = ze[i] - '1', bit = 1 << c;;
				if (!(known&bit)) {
					vec[c] = vi++;
					known |= bit;
				}
				ze[i] = vec[c];
			}
			break;
		}


		//================= extract  
		case 21:// extract 81 start sgo.vx[1]
			if (ll >= (int)sgo.vx[1]) {
				ze[sgo.vx[1] + 81] = 0;
				fout1 << &ze [sgo.vx[1]] << endl;
			}
			break;
		case 22:// extract first sgo.vx[1] puzzles
			if (++npuz < sgo.vx[1]) fout1 << ze << endl;
			else fout2 << ze << endl;
			break;
		case 23://sampling start sgo.vx[1] one every sgo.vx[2]
			if (++npuz < sgo.vx[1]) break;
			{	int rn = (npuz - sgo.vx[1]) % sgo.vx[2];
				if (!rn)fout1 << ze << endl;
			}
			break;
		}
	}
}
void Go_c401() {// .dat to .txt
	cout << "Go_401 entry " << sgo.finput_name << "dat file to txt file" << endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	char zs[200];  int n = 0;
	while (!finput.eof()) {
		char c;
		finput.get(c); if (finput.eof()) break; // fin ligne normale
		if (c < 20) {//fin ligne 0a
			if (n) {// no line length 0
				zs[n] = 0;
				fout1 << zs << endl;
				n = 0;
			}
		}
		else {
			if (n < 190) zs[n++] = c; // cut line more than 190
		}
	}
	if (n) {// last line if exist
		zs[n] = 0;
		fout1 << zs << endl;
	}
}
void Go_c402() {// morph row s[1] cols s[2] v[1] diag
	cout << "Go_c402 entry " << sgo.finput_name << " morph a puzzle" << endl;
	char ze[82],zout[82]; ze[81] = zout[81] = 0;
	char * sr = sgo.s_strings[1], *sc = sgo.s_strings[2];
	if (!sr || !sc) {
		cout << "check missing -s1- for rows and/or -s2- for cols " << endl;
		return;
	}
	int rows = 0, cols = 0, tr[9], tc[9];
	int lrows = (int)strlen(sr), lcols = (int)strlen(sc);
	if (lrows != 9 || lcols != 9) {
		cout << "check -s1- for rows -s2- for cols must digits 1-9" << endl;
		return;
	}
	for (int i = 0; i < 9; i++) {
		int r = sr[i], c = sc[i];
		if (r<'1' || r>'9' || c<'1' || c>'9') {
			cout << "check -s1- -s2- for  digits 1-9" << endl;
			return;
		}
		r -= '1'; c -= '1';
		rows |= 1 << r; cols |= 1 << c;
		tr[i] = r; tc[i] = c;
	}
	if(rows!=0x1ff || cols!=0x1ff) {
		cout << "check -s1- -s2- for  digits 1-9 must all be there" << endl;
		return;
	}

	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	// build the morphing table
	int tmorph[81];
	for (int ir = 0; ir < 9; ir++) for(int ic=0;ic<9;ic++){
		tmorph[9 * ir + ic] = 9 * tr[ir] + tc[ic];
	}
	if (sgo.vx[1]) for(int i=0;i<81;i++)// then make it diagonal 
		tmorph[i]= C_transpose_d[tmorph[i]];
	while (finput.GetPuzzle(ze)) {
		for (int i = 0; i < 81; i++) zout[i] = ze[tmorph[i]];
		fout1 << zout << endl;
	}
}


void Go_cxx() {
	cout << "Go_xx entry " << sgo.finput_name << " xxxxx" << endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	char * ze = finput.ze;
	while (finput.GetLigne()) {
	}

	char zep[82]; zep[81] = 0;
	while (finput.GetPuzzle(zep)) {
	}
}

//_________________________________________________________________________

void Go_c440(){
	cout << "Go_440 entry " << sgo.finput_name  <<" game results to parse"<< endl;
	//int v = sgo.vx[0];
	char * ze = finput.ze;
	finput.open(sgo.finput_name);
	if (!finput.is_open()){
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	while (finput.GetLigne()){
		int ll = (int)strlen(ze);
		if (ll<106 || ll > 200) continue;
		char zout[200];
		for (int i = 0; i < 81; i++)
			if (ze[i] - '0')zout[i] = ze[i];
			else zout[i] = '.';
		zout[81] = ';';
		int n = 82;
		for (int i = 88; i < 102; i++){
			if (ze[i] == ' ')continue;
			if (ze[i] == '/'){
				zout[n++] = ';'; continue;
			}
			zout[n++] = ze[i];
		}
		zout[n++] = ';';
		if (ze[103] == '+')zout[n++] = '+'; else zout[n++] = ' ';
		zout[n++] = ';';
		strcpy(&zout[n], &ze[104]);
		fout1 << zout << endl;
		//break;
	}

}
void Go_c445(){
	cout << "Go_445 entry " << sgo.finput_name << " param=" << sgo.bfx[0] << endl;
	if (_popcnt32(sgo.bfx[0]) != 1) return;//pointer to the  parameter to consider
	uint32_t ipar; bitscanforward(ipar, sgo.bfx[0]);
	cout << "split o, parameter rank=" << ipar << " file 1 <=" << sgo.vx[0] << endl;
	int v = sgo.vx[0];
	//ipar--;// switch to index;
	char * ze = finput.ze;
	finput.open(sgo.finput_name);
	if (!finput.is_open()){
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	while (finput.GetLigne()){
		int ll = (int)strlen(ze);
		if (ll<82 || ll > 200) continue;
		sgo.ParseInt(ze, ';');
		//cout << ze << " check v=" << sgo.tparse[ipar] << endl;
		if (sgo.tparse[ipar] <= v) fout1 << ze << endl;
		else fout2 << ze << endl;
		//break;
	}

}

void Go_c480() {//add  compressed clues to entry
	cout << "Go_481 entry " << sgo.finput_name << " base check" << endl;
	if (!sgo.foutput_name) {
		cerr << "missing output root" << endl;
		return;
	}
	if (!sgo.s_strings[0]) {
		cerr << "missing base in" << endl;
		return;
	}
	if (!sgo.s_strings[1]) {
		cerr << "missing base out" << endl;
		return;
	}
	//int update = sgo.vx[0];
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	char ze[200]; ze[81] = ';';
	int lcmp = 0;
	while (finput.GetPuzzle(ze)) {
		int /*n = 0,*/ p = 82;
		for (int i = 9; i < 81; i++)if (ze[i] - '.')
			ze[p++] = ze[i];

		int ll = p-82;
		if (lcmp && lcmp != ll) {
			cerr << "stop wrong input file" << endl;
		}
		else lcmp = ll;
		ze[p] = 0;
		fout1 << ze << endl;
	}

}
void Go_c481() {//base check -i ads -s1- base -o add root
	FINPUT fin2;
	cout << "Go_481 entry " << sgo.finput_name << " base check" << endl;
	if (!sgo.foutput_name) {
		cerr << "missing output root" << endl;
		return;
	}
	if (!sgo.s_strings[0]) {
		cerr << "missing base in" << endl;
		return;
	}
	int update  = sgo.vx[0];
	char * ze = finput.ze;
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	fin2.open(sgo.s_strings[0]);
	if (!fin2.is_open()) {
		cerr << "error open base " << sgo.s_strings[0] << endl;
		return;
	}
	

	//__________________________________
	int lcmp = 0;
	char scomp[10]; scomp[0] = 0;
	char * zcompress = &ze[81],*zeb= fin2.ze;;
	zeb[0] = 0;
	while (finput.GetLigne()) {
		int ll = (int)strlen(ze);
		if (ll < 91 || ll > 120) {
			cerr << "wrong add file cancelled " << endl;
			return;
		}
		if(lcmp){
			if (ll != lcmp) {
				cerr << "stop wrong add file" << endl;
				return;
			}
		}
		else lcmp = ll;
		//int lbase = (int)strlen(zcompress);
		// read/ out if update base as long as below
	loop_base:
		int icomp = strcmp(zcompress, zeb);
		if (!icomp)goto next_add;
		if (icomp > 0) {// next in base
			if (zeb[0] && update) fout2 << zeb << endl; //write old if update
			if(fin2.GetLigne())	goto loop_base;
		}
		// new or end of file in old base
		if (update) {
			fout2 << zcompress << endl;// new in base
			ze[81] = 0;
			fout1 << ze << endl;// just the puzzle in output
		}
		else fout1 << ze << endl;// just recopy in output
			//char zout[200];
	next_add:;
	}
	if (update && zeb[0] < 255) {// copy the rest of the base
		fout2 << zeb << endl;
		while(fin2.GetLigne())fout2 << zeb << endl;
	}
}


void Go_c484() {
	cout << "Go_484 entry " << sgo.finput_name << " data base to restore" << endl;
	if (!sgo.s_strings[0]) {
		cerr << "missing -s1- pattern " << endl;
		return;
	}
	if (strlen(sgo.s_strings[0]) != 81) {
		cerr << "-s1- not length 81 " << endl;
		return;
	}
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	char * ze = finput.ze;
	char zout[82];
	strcpy(zout, empty_puzzle);
	int tclues[50], nclues = 0,puz_int[81];
	memset(puz_int, 0, sizeof puz_int);
	char *w = sgo.s_strings[0];
	for (int i = 0; i < 81; i++)
		if (w[i] >= '1' && w[i] <= '9')puz_int[i] = w[i];
	for (int i = 0; i < 9; i++) if (puz_int[i]) zout[i] = puz_int[i];
	for (int i = 9; i < 81; i++) 
		if (puz_int[i]&& nclues<40) tclues[nclues++]=i;
	if (nclues > 30) {
		cerr << "-s1- too many clues cancel " << endl;
		return;
	}
	uint64_t ne = 0,nd= sgo.vx[1],nf= sgo.vx[2];
	while (finput.GetLigne()) {
		if (++ne < nd)continue;
		if (nf && ne >= nf)break;
	}

}

/*


void GO_MISC::Do_70(){  // start a game from gremlin pattenr
	char * ze=myin->ze,zs[82];
	zs[81]=0;
	for(int i=0;i<9;i++){
		myin->GetLigne();
		int length=strlen(ze);
		(*myout1)<<ze<<endl;
		for(int j=0;j<9;j++)
			zs[9*i+j]=ze[2*j];
	}
	(*myout1)<<endl<<endl<<zs<<endl;
	for(int i=0;i<81;i++)
		if(zs[i]-'.') zs[i]='1';
	(*myout1)<<endl<<zs<<endl<<endl;
	char *a="1__2__3__",*b="4__5__6__",*c="7__8__9__";
	(*myout1)<<a<<a<<a<<b<<b<<b<<c<<c<<c<<endl;

}

void Go_c484(){  // restore game database for a pattern
	char * ze=myin->ze,zs[82];
	zs[81]=0;
	PUZ0 puzd;
	puzd.Empty();
	if((!options.first) ||	(strlen(options.first)-81) ){
		(*myout2)<< " pattern missing or length not correct "<<endl;
		return;
	}
	(*myout2)<<options.first<< " pattern "<<endl;
	int nd=0;
	for(int i=9;i<81;i++)
		if(options.first[i]-'.') nd++;

	while(myin->GetLigne()){
		if(strlen(ze)-nd) continue;
		int nx=0,lig1='1';
		for(int i=0;i<81;i++)
			if(options.first[i]=='.')
				zs[i]='.';
			else
				if(i<9)
					zs[i]=lig1++;
				else
					zs[i]=ze[nx++];

		(*myout1)<<zs<<endl;
	}

}

void GO_MISC::Do_72(){  // find puzzles close to a pattern
	char * ze=myin->ze,ztarget[82],maxf=11;
	PUZ0 puzd;
	puzd.Empty();
	if((!options.first) ||	(strlen(options.first)-81) ){
		(*myout2)<< " pattern missing or length not correct "<<endl;
		return;
	}
	(*myout2)<<options.first<< " target "<<endl;
	int nd=0;
	for(int i=9;i<81;i++)
		if(options.first[i]-'.') ztarget[nd++]=options.first[i];

	while(myin->GetLigne()){
		if(strlen(ze)-nd) {
			(*myout2)<< " length not correct in the file"<<endl;
			return;
		}
		int nx=0;
		for(int i=0;i<nd;i++)
			if(ztarget[i]-ze[i]) nx++;
		if(nx>=maxf) continue;
		if(nx>=6)maxf=nx;
		(*myout2)<<ze<<";"<<nx<<endl;
		if(nx<5)
			(*myout1)<<ze<<endl;

	}

}

*/
