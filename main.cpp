
#define _CRT_SECURE_NO_DEPRECATE
#include "main.h"
#include "Zhn_cpp.h"


// catching time as seconds+millis  (seconds since year 1970)
int32_t GetTimeMillis() {
	struct _timeb tbuf;
	_ftime64_s(&tbuf); 
	return ((int32_t)(1000 * tbuf.time) + tbuf.millitm);
}
// builing an appropriate message depending on the elapsed time te-ts
void PrintTime(int32_t ts,int32_t te){
	UINT dt=te-ts,dtmil=dt%1000,dts=dt/1000,dth=dts/3600;   dth=dth%1000;
	cerr << endl<<"total elapsed time "; 
    UINT dtm=dts/60; dts=dts%60 ,   dth=dtm/60, dtm=dtm%60;
    if(dth) cerr <<dth<<"h "; 
	if(dth || dtm) cerr <<dtm<<"m "; 
	cerr	<<dts <<"s ";
	if(dtmil<10) cerr << "00"; else  if(dtmil<100) cerr << '0';
	cerr <<dtmil<<"ms "<<endl;   return;
}

void PrintTimeCout(int32_t ts, int32_t te){
	UINT dt = te - ts, dtmil = dt % 1000, dts = dt / 1000, dth = dts / 3600;   dth = dth % 1000;
	cout << endl << "total elapsed time ";
	UINT dtm = dts / 60; dts = dts % 60, dth = dtm / 60, dtm = dtm % 60;
	if (dth) cout << dth << "h ";
	if (dth || dtm) cout << dtm << "m ";
	cout << dts << "s ";
	if (dtmil<10) cout << "00"; else  if (dtmil<100) cout << '0';
	cout << dtmil << "ms " << endl;   return;
}


int Search_ccd(const char * ww)
{	// List of 2 char command, 9 commands
	const char * ccd[]={"-i" ,    // input name including extension
				  "-o" ,	//	output or second filename
				  "-c",  // main command option + check boxes
				  "-v" ,  // value   0 to 9 default 0
				  "-b" ,  // bit field 0 to 9 default is 0
				  "-s",  // strings 0 to 9
	};   // puzzle processing specificities
	char wt[4]; 
	strncpy_s(wt,4,ww,2);
	wt[2]=0;
	for(int i=0;i<6;i++)
		if(!strcmp(wt,ccd[i])) 
			return i;
	return -1;
}

//#include "Go_0.cpp"

SGO sgo;
extern void Go_0();
int main(int narg, char *argv[]) {
	cerr << "mainstart" << endl;
	int32_t tdeb=GetTimeMillis();
	char * finput_name=0,*foutput_name=0,* ww;
	char * s_strings[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };//optionnal 10 strings

	uint32_t command = 0, 
		vx[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, //integers 0 to 9
		bfx[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };// bit fields 

	for(int i=0;i<narg;i++)	{
		ww = argv[i];
		int ir=Search_ccd(ww);
		if(ir<0) continue;
		if (ir == 3){// -vn-xxxx
			if (ww[3] - '-') continue; //must be -vn-  with n 0_9
			int ind=ww[2] - '0';
			if (ind < 0 || ind>9)continue;
			vx[ind] = atoi(&ww[4]);
			continue;
		}	
		else if (ir == 4){//  -bn- followed by a bit field bits rigth to left max 8 bits
			if (ww[3] - '-') continue; 
			int ind = ww[2] - '0';
			if (ind < 0 || ind>9)continue;
			int length = (int)strlen(&ww[4]);
			if (length > 8)continue; 
			for (int i = 0; i < length; i++) if (ww[4 + i] == '1') bfx[ind] |= (1 << i);
		}
		else if (ir == 5){//  -sn- followed by a string
			if (ww[3] - '-') continue; //must be -vn-  with n 0_9
			int ind = ww[2] - '0';
			if (ind < 0 || ind>9)continue;
			s_strings[ind] = &ww[4];
		}
		else{
			switch (ir)	{
			case 0: finput_name = &ww[2]; 	break;  // -i
			case 1: foutput_name = &ww[2]; 	break;  // -o
			case 2: command = atoi(&ww[2]);	break; //-c
			}// end command  
		}
	}// end loop on options
	if(finput_name) cerr <<" file1 (input) " << finput_name<<endl;
	if(foutput_name) cerr <<" file2 (output) " << foutput_name<<endl;
	cerr << "command " << command<<endl;
	// set zh_g general commands
	zh_g.modeguess =(int) bfx[9];
	if (!(zh_g.modeguess & 1))zh_g.modeguess = 0;
	zh_g.maxindex = (int)vx[7];
	// store command line parameters 
	sgo.command = command;
	sgo.bfx = bfx;
	sgo.finput_name = finput_name;
	sgo.foutput_name = foutput_name;
	sgo.s_strings = s_strings;
	sgo.vx = vx;
	
	Go_0();
	cerr << " print cout time "  << endl;

	int32_t tfin=GetTimeMillis();
    PrintTimeCout(tdeb,tfin);
	PrintTime(tdeb, tfin);
	return 0;
}

void SGO::ParseInt(char * ze, int  delimiter){
// bfx[0] 1 to 8 parameters
	__stosd((uint32_t *)tparse, 0, 8); nparse = 0;
	if (!bfx[0]) return;
	//cout << ze << "go parse  delimiter "<<(char) delimiter << endl;
	char * w = ze,temp[20];
	int pos = 0, i = -1;
	while (1){// scan the entry
		++i;
		if (i > 200) break; // safety code
		if (w[i] == delimiter || w[i]==0){
			//cout << &w[i] << " seen param" << endl;
			int n = i - pos; // length of the parameter
			if (n > 15 || n < 1)goto next;
			if (!(bfx[0] & (1 << nparse))) goto next;
			strncpy(temp, &w[pos], n); temp[n] = 0;
			//cout << temp << "parse" << endl;
			tparse[nparse] = atoi(temp);
		next:
			if (!w[i]) return;
			if (++nparse > 7)return;
			pos = i + 1;
		}
	}
	return;
}
