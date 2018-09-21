

#define _CRT_SECURE_NO_DEPRECATE
#include <sys/timeb.h>
#include "sk_t.h"  // tab0  
#include "Zhn.h"

struct SGO{
	int command;
	uint32_t * vx, *bfx;
	char ** s_strings, *finput_name, *foutput_name;
	int tparse[8], nparse;

	void ParseInt(char * ze, int  delimiter);
};

class FINPUT: public ifstream{
	word 		endf; //! Indicate if end of file has been found 1=EOF, 0=not EOF
public:	
	char * namex, //Input file name (without file type). default "puzzles"
	      ze[600]; // big buffer to avoid problems linked to getline

	FINPUT(){	endf = 0;	namex = "puzzles";	}
	void SetName(char * nn){	namex = nn;	}

	int OpenFI(){
		char zn[300];
		strcpy(zn, namex);
		int ll = (int)strlen(zn);
		//	strcpy(&zn[ll],".txt"); // add file type
		ifstream::open(zn);
		if (is_open()) return 0;
		endf = 1;
		return 1;
	}

	int GetLigne(){
		if (endf) return 0;
		if (eof() || (!getline(ze, 550))){ // assumed end of file
			endf = 1;
			return 0;
		}
		ze[551] = 0;
		int i = (int)strlen(ze);
		if (i>500)		{ // not correct treat as end of file
			endf = 1; return 0;
		}
		return 1;
	}

	int GetPuzzle(char * zp){
		while (GetLigne())	{
			if (strlen(ze) < 81)continue;
			zp[81] = 0; // file truncated to 81
			// standard file  clues or "." minimum 17 clues
			int n = 0;
			for (int i = 0; i<81; i++) 		{
				char c = ze[i];
				if (c<'1' || c>'9') zp[i] = '.';
				else { n++; zp[i] = ze[i]; }
			}
			if (n<17) continue;
			return 1;
		}// end while
		return 0;
	}

};

extern long GetTimeMillis();
extern  void PrintTime(long ts, long te);
extern void PrintTimeCout(long ts, long te);