
#include "sk_t.h"
//#include "sk_bitfields.h"
USHORT BF16::CountEtString(char *s) {
	USHORT n = 0;
	for (int i = 0; i < 9; i++) {
		if (On(i))
			s[n++] = (char)(i + '1');
	}
	s[n] = 0;
	return n;
}
char * BF16::String(char * ws, int lettre ) {
	int n = 0;
	for (int i = 0; i < 9; i++) {
		if (On(i))
			ws[n++] = (char)(lettre ? i + 'A' : i + '1');
	}
	ws[n] = 0;
	return ws;
}
USHORT BF16::String(USHORT * r){ 
	if (!f)	return 0;
	unsigned long cc;
	int n = 0;
	register  unsigned int Rw = f;
	while (_BitScanForward(&cc, Rw)){// look for  possible cells
		Rw ^= 1 << cc;// clear bit
		r[n++] = (USHORT)cc;
	}
	return n;
}

USHORT BF32::String(USHORT * r){ 
	if (!f)	return 0;
	unsigned long cc;
	int n = 0;
	register  unsigned int Rw = f;
	while (_BitScanForward(&cc, Rw)){// look for  possible cells
		Rw ^= 1 << cc;// clear bit
		r[n++] = (USHORT)cc;
	}
	return n;
}
USHORT BF32::String(USHORT * r, int i1, int i2){
	UINT  n = 0, x = 1 << i1;
	for (int i = i1; i<i2; i++){
		if (f&x) r[n++] = i;
		x = x << 1;
	}
	return n;
}

void  BF128::ClearDiag(int band, int stack){
// stack appears here as a band
	register int R=band;
	unsigned long cell;
	while (_BitScanForward(&cell, R)){
		R ^= (1 << cell);// clear bit
		int celld = C_transpose_d[cell + 27 * stack];
		Clear_c(celld);
	}
}
void BF128::ClearRow(int clear, int row){
	register int R = clear;
	unsigned long col;
	while (_BitScanForward(&col, R)){
		R ^= (1 << col);// clear bit
		Clear_c(9 * row + col);
	}
}
void BF128::ClearCol(int clear, int col){
	register int R = clear;
	unsigned long row;
	while (_BitScanForward(&row,R)){		
		R ^= (1 << row);// clear bit
		Clear_c(9 * row + col);
	}
}


void BF128::Diag3x27(BF128 & r){
	__stosq(bf.u64, 0, 2);
	if (r.isEmpty()) return;
	for (int iband = 0; iband < 3; iband++){// build diagonal sym
		register int band = r.bf.u32[iband];
		unsigned long cell;
		while (_BitScanForward(&cell, band)){
			band ^= (1 << cell);// clear bit
			int celld = C_transpose_d[cell + 27 * iband];
			Set_c(celld);
		}
	}

}

int BF128::String3X27_to_gint_c(GINT * t, int digit){// active cells 0-80
	unsigned long cc;
	int n = 0,drb=0;
	for (int i32 = 0; i32 < 3; i32++,drb+=27){// process 32 bits units
		register  unsigned int Rw = bf.u32[i32];
		while (_BitScanForward(&cc, Rw)){// look for  possible cells
			Rw ^= 1 << cc;;// clear bit
			t[n++].u32 =cc + drb + (digit << 8);
		}
	}
	return n;
}
int BF128::String3X27(int * r){
	unsigned long cc;
	int n = 0, drb = 0;
	for (int i32 = 0; i32 < 3; i32++, drb += 27){// process 32 bits units
		register  unsigned int Rw = bf.u32[i32];
		while (_BitScanForward(&cc, Rw)){// look for  possible cells
			Rw ^= 1 << cc;;// clear bit
			r[n++] = cc + drb ;
		}
	}
	return n;	
}
int BF128::String81GP(int * r){// derived from bitfield.h supposed  best process
	unsigned long cc;
	int n = 0, drb = 0;
	for (int i32 = 0; i32 < 3; i32++, drb += 32){// process 32 bits units
		register  unsigned int Rw = bf.u32[i32];
		while (_BitScanForward(&cc, Rw)){// look for  possible cells
			Rw ^= 1 << cc;;// clear bit
			r[n++] = cc + drb;
		}
	}
	return n;

}
int BF128::StringGP(int * t, int bloc){// derived from bitfield.h supposed  best process
	unsigned long cc;
	int vbloc = bloc << 7, n = 0; // 128 bits per bloc start index bloc * 128
	for (int i64 = 0; i64 < 2; i64++, vbloc += 64){// process 64 bits units
		register  uint64_t Rw = bf.u64[i64];
		while (_BitScanForward64(&cc, Rw)){// look for  possible cells
			register uint64_t bit = 1; bit <<= cc;
			Rw ^= bit;// clear bit
			t[n++] = cc + vbloc;
		}
	}
	return n;
}
char * BF128::String3X_Rev(char * ws){
	strcpy(ws, empty_puzzle);
	uint32_t *bfw = bf.u32;
	unsigned int w = bfw[0];
	for (int j = 26; j >= 0; j--) if (w & (1 << j))
		ws[j] = '1';
	w = bfw[1];
	for (int j = 26; j >= 0; j--) if (w & (1 << j))
		ws[j + 27] = '1';
	w = bfw[2];
	for (int j = 26; j >= 0; j--) if (w & (1 << j))
		ws[j + 54] = '1';
	return ws;

}
char * BF128::String81(char * ws){
	strcpy(ws, empty_puzzle);
	for (int j = 0; j<81; j++) if (On(j))
		ws[j] = '1';
	return ws;
}
char * BF128::String3X(char * ws){
	strcpy(ws, empty_puzzle);
	for (int j = 0; j<81; j++) if (On(C_To128[j]))
		ws[j] = '1';
	return ws;
}
char * BF128::String128(char * ws){
	ws[128] = 0;
	for (int j = 0; j<128; j++)
		if (On(j))		ws[j] = '1';		else ws[j] = '.';
	return ws;

}

/* pas encore fait dans PM3X
void Set(BF16 & digs, int cell);
void SetRegion(int dig, int unit, BF16 & pdigs);
void operator &= (const BF128 * bf128);
void Image(BUILDSTRING & zs);
USHORT String(UCAND * t);
*/

void PM3X::operator &= (const PM3X &z2){
	for (register int i = 0; i < 9; i++) pmdig[i] &= z2.pmdig[i];
}
void PM3X::operator |= (const PM3X &z2){
	for (register int i = 0; i < 9; i++) pmdig[i] |= z2.pmdig[i];
}
void PM3X::operator -= (const PM3X &z2){
	for (register int i = 0; i < 9; i++) pmdig[i] -= z2.pmdig[i];
}
int PM3X::IsEmpty(){
	for (register int i = 0; i<9; i++)
		if (pmdig[i].isNotEmpty()) return 0;
	return 1;
}

int PM3X::Count(){
	register uint64_t *w = pmdig[0].bf.u64,
    n=0;
	for (register int i = 0; i < 18; i++)	n += _popcnt64(w[i]);
	return (int)n;
}

void PM3X::Print(char * lib){
	cout << "pm3x status for " << lib << endl;
	for (int i = 0; i < 3;i++)	cout << ".........+++++++++---------";
	cout << endl;
	char ws[82];
	for (int i = 0; i < 9; i++)if (pmdig[i].isNotEmpty())
		cout << pmdig[i].String3X(ws) << " " << i + 1<<endl;
	for (int i = 0; i < 3; i++)	cout << ".........+++++++++---------";
	cout << endl << endl;

}



BF81::BF81(int i1, int i2) {
	SetToBit(i1);
	setBit(i2);
}
BF81::BF81(char * mode, int i1, int i2) {
	clear();
	if ((*mode) == 'z'){
		bf = cellsFixedData[i1].z;
		*this &= cellsFixedData[i2].z;
	}
}
BF81::BF81(char * mode, int i1) {
	clear();
	if ((*mode) == 'z')
		bf = cellsFixedData[i1].z;
}
BF81::BF81(const T128 &r, const BF81 & r2) {
	bf = r;
	*this &= r2;
}
void BF81::PackRows(BF16 * rows){
	bf.u32[3] = 0;
	bf.u32[2] = (rows[7].f >> 1) | (rows[8].f << 8);
	bf.u32[1] = (rows[7].f << 31) | (rows[6].f << 22) | (rows[5].f << 13)
		| (rows[4].f << 4) | (rows[3].f >> 5);
	bf.u32[0] = (rows[3].f << 27) | (rows[2].f << 18) | (rows[1].f << 9) | rows[0].f;
}
void BF81::OrBand(int F, int iband){
	if (!iband) { bf.u32[0] |= F; return; }
	if (iband == 1){
		bf.u32[0] |= ((F & 0x1f) << 27);
		bf.u32[1] |= (F >> 5);
		return;
	}
	bf.u32[1] |= ((F & 0x3ff) << 22);
	bf.u32[2] |= (F >> 10);
}
void BF81::LoadZhou(int * F){ //
	SetAll_0(); // could be  u32[3]=0 
	bf.u32[0] =	(F[0] & BIT_SET_27) | ((F[1] & 0x1f) << 27);
	bf.u32[1] |=	((F[1] & BIT_SET_27) >> 5) | ((F[2] & 0x3ff) << 22);
	bf.u32[2] =	((F[2] & BIT_SET_27) >> 10);
}
void BF81::BackZhou(int  * Fx){
	Fx[0] = bf.u32[0] & BIT_SET_27;
	Fx[1] = ((bf.u32[0] >> 27) |	(bf.u32[1] << 5))		& BIT_SET_27;
	Fx[2] = ((bf.u32[1] >> 22) |	(bf.u32[2] << 10))		& BIT_SET_27;
}
USHORT BF81::String( USHORT *r) {
	unsigned long cc;
	int n = 0, drb = 0;
	for (int i32 = 0; i32 < 3; i32++, drb += 32){// process 32 bits units
		register  unsigned int Rw = bf.u32[i32];
		while (_BitScanForward(&cc, Rw)){// look for  possible cells
			Rw ^= 1 << cc;;// clear bit
			r[n++] = (USHORT)(cc + drb);
		}
	}
	return n;
}
USHORT BF81::String(USHORT *r, USHORT digit){
	USHORT n = String(r), x = digit << 7;
	for (int i = 0; i<n; i++)
		r[i] |= x; // insert digit to get a UCAND
	return n;
}
USHORT BF81::StringUnit(int unit, USHORT *r) {
	USHORT n = 0;
	byte *t = cellsInGroup[unit];
	for (int i = 0; i<9; i++, t++)
		if (On(*t)) r[n++] = *t;
	return n;
}
USHORT BF81::StringUnit(int unit, USHORT *r, USHORT digit) {
	USHORT n = 0;
	byte *t = cellsInGroup[unit];
	for (int i = 0; i<9; i++, t++)
		if (On(*t)) r[n++] = *t | (digit << 7);
	return n;
}
USHORT BF81::GetRegion(int unit){
	byte *t = cellsInGroup[unit];
	BF16 rr;
	for (int i = 0; i<9; i++, t++)
		if (On(*t))rr.Set(i);
	return rr.f;
}
BF16 BF81::GetBFRegion(int unit){
	byte *t = cellsInGroup[unit];
	BF16 rr;
	for (int i = 0; i<9; i++, t++)
		if (this->On(*t))rr.Set(i);
	return rr;
}
/*
void BF81::Image(BUILDSTRING & zs, int digit, char * lib, int doinit){
if (doinit)
zs.Init();
if (lib)
zs.Astr(lib);
USHORT tp[100], ntp = String(tp);
for (int i = 0; i<ntp; i++){
if (digit >= 0 && digit <9)
zs.Achint(digit);
zs.Astr(cellsFixedData[tp[i]].pt);
zs.Aspace();
}

zs.Close();
}
*/

void BF81::Store(USHORT * tstore){
	for (int i = 0; i<8; i++) tstore[i] = (USHORT)bf.u16[i];
}
void BF81::Re_Load(USHORT * tstore){
	for (int i = 0; i<8; i++)  bf.u16[i] = (uint16_t)tstore[i];
}
