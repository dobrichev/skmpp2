
// working in solver mode

/*
DynamicForcingChain=85,
DynamicForcingChainPlus=90,
// 90 is DynamicForcingChainPlus
// all "events"  claiming, pointing, pair, hidden pair, XWing
// follows an empty step 85
// consider each false candidate as start
// search for new bi values. If none, skip it
// look for new false thru basic sets
NestedForcingChain=95,
Nestedmultiple chains=100,
NestedForcingChain dynamic subchains dynamic=105,
NesttedLevel5=110
*/


extern PM_GO pm_go;
void Go_c110(){// template serate mode
    if (!sgo.finput_name) return;
    int64_t cptg[20];
	memset(cptg, 0, sizeof cptg);
	std::cout << "Go_110 entry " << sgo.finput_name << " input" << std::endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()){
		std::cerr << "error open file " << sgo.finput_name << std::endl;
		return;
	}
	char ze[200]; ze[81] = 0;
	uint32_t npuz = 0;
	int32_t ttdeb = GetTimeMillis();
	while (finput.GetPuzzle(ze)){
		npuz++;
		if (npuz < sgo.vx[0]) continue;
		if(pm_go.opprint2)std::cout << finput.ze << "to process npuz=" << npuz << std::endl;
		zh_g.npuz = npuz;
		int32_t tdeb = GetTimeMillis();
		if (zh_g.Go_InitSolve(ze)) {
			std::cout << finput.ze << "invalid or multiple solutions npuz=" << npuz << std::endl;
			continue;
		}
		pm_go.SolveSerate110();
		if (sgo.vx[1] && npuz >= sgo.vx[1]) break;
	}
}
//fout_diam, fout_pearl, fout_l45, fout_l65, fout_solved, fout_unsolved;
//fout1=l45;  fout2 solved;   fout3 unsolved
void Go_c111(){// fast serate mode 
	if (!sgo.finput_name) return;
	int64_t cptg[20];
	memset(cptg, 0, sizeof cptg);
	std::cout << "Go_111 entry " << sgo.finput_name << " input" << std::endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()){
		std::cerr << "error open file " << sgo.finput_name << std::endl;
		return;
	}
	if (!sgo.foutput_name)return;
    {// create output files to split rated puzzles
		char zn[200];
		strcpy(zn, sgo.foutput_name);
		int ll = (int)strlen(zn);
		strcpy(&zn[ll], "_0diam.txt");		fout_diam.open(zn);// quasi diam 45_,,
		strcpy(&zn[ll], "_1pearl.txt");		fout_pearl.open(zn);// quasi pearl 45_,,
		strcpy(&zn[ll], "_2l65.txt");		fout_l65.open(zn);// rating 45 62
  	}

	char ze[200]; ze[81] = 0;
	uint32_t npuz = 0;
	int32_t ttdeb = GetTimeMillis();
	while (finput.GetPuzzle(ze)){
		npuz++;
		if (0)std::cout << finput.ze << "to process npuz=" << npuz << std::endl;
		zh_g.npuz = npuz;
		int32_t tdeb = GetTimeMillis();
		if (zh_g.Go_InitSolve(ze)) {
			std::cout << finput.ze << "invalid or multiple solutions npuz=" << npuz << std::endl;
			continue;
		}
		pm_go.SolveSerate111();
	}



}
void Go_c199(){// test on demand
	if (!sgo.finput_name) return;
	std::cout << "Go_199 entry " << sgo.finput_name << " input" << std::endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()){
		std::cerr << "error open file " << sgo.finput_name << std::endl;
		return;
	}
	char ze[200]; ze[81] = 0;
	uint32_t npuz = 0;
	while (finput.GetPuzzle(ze)){
		npuz++;
		if (npuz < sgo.vx[0]) continue;
		std::cout << finput.ze << "to process npuz=" << npuz << std::endl;
		zh_g.npuz = npuz;
		if (zh_g.Go_InitSolve(ze)) {
			std::cout << finput.ze << "invalid or multiple solutions npuz=" << npuz << std::endl;
			continue;
		}
		pm_go.Solve199test();
		if (npuz >= sgo.vx[1]) break;
	}
}
