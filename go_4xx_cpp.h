
/*
..\release\skmpp    -i"d6_6bb92"  -c"444"
440 parse games results 
441 just add subgroup in ourput
442 dat file to txt file
443 puzzle + std puzzle empty '.'
444 switch from  ed=  /  /  to ;  ;  ;
445 select on quick ed -v is ed mini to keep

..\release\skmpp    -i"extractclues"  -c"446"   -v24         "//split =idpos clues"
 #  19    1.2/1.2/1.2 - gsf
 1.2/1.2/1.2 - gsf
*/


/*

void GO_MISC::Do_71(){  // restore game database for a pattern
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

}*/


void Go_c440(){
	cout << "Go_440 entry " << sgo.finput_name  <<" game results to parse"<< endl;
	int v = sgo.vx[0];
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
	if (__popcnt(sgo.bfx[0]) != 1) return;//pointer to the  parameter to consider
	unsigned long ipar; _BitScanForward(&ipar, sgo.bfx[0]);
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
