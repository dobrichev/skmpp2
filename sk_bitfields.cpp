<<<<<<< HEAD

#include "sk_t.h"

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

int BF32::Table(int * r) {
	if (!f)	return 0;
	int n = 0;
	BitsInTable32(r, n, f, 0);
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

void  BF128::ClearDiag(int band, int stack) {
	// stack appears here as a band
	int tp[32], ntp=0;
	BitsInTable32(tp, ntp, band);
	for (int i = 0; i < ntp; i++) {
		int cell=tp[i], celld = C_transpose_d[cell + 27 * stack];
		Clear_c(celld);
	}
}
void BF128::ClearRow(int clear, int row) {
	for (int i = 0,bit=1; i < 9; i++,bit<<=1) {
		if(bit&clear)Clear_c(9 * row + i);
	}
}
void BF128::ClearCol(int clear, int col) {
	for (int i = 0, bit = 1; i < 9; i++, bit <<= 1) {
		if (bit&clear)Clear_c(9 * i + col);
	}
}
void BF128::Diag3x27(BF128 & r){
	__stosq(bf.u64, 0, 2);
	if (r.isEmpty()) return;
	for (int iband = 0; iband < 3; iband++){// build diagonal sym
		int tc[27], ntc=0;
		BitsInTable32(tc, ntc, r.bf.u32[iband], 27 * iband);
		for (int i = 0; i < ntc; i++) {
			int celld = C_transpose_d[tc[i]];
			Set_c(celld);
		}
	}
}

int BF128::Table3X27(int * r) {	
	int n = 0;
	BitsInTable32(r, n, bf.u32[0]);
	BitsInTable32(r, n, bf.u32[1],27);
	BitsInTable32(r, n, bf.u32[2], 54);
	return n;
}
int BF128::Table128(int * r) {
	int n = 0;
	BitsInTable64(r, n, bf.u64[0]);
	BitsInTable64(r, n, bf.u64[1], 64);
	return n;
}
int BF128::Table64_0(int * r) {
	int n = 0;
	BitsInTable64(r, n, bf.u64[0]);
	return n;
}
int BF128::Table64_1(int * r) {
	int n = 0;
	BitsInTable64(r, n, bf.u64[1], 64);
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

void PM3X::Print(const char * lib){
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

=======

#include "sk_t.h"

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

int BF32::Table(int * r) {
	if (!f)	return 0;
	int n = 0;
	BitsInTable32(r, n, f, 0);
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

void  BF128::ClearDiag(int band, int stack) {
	// stack appears here as a band
	int tp[32], ntp=0;
	BitsInTable32(tp, ntp, band);
	for (int i = 0; i < ntp; i++) {
		int cell=tp[i], celld = C_transpose_d[cell + 27 * stack];
		Clear_c(celld);
	}
}
void BF128::ClearRow(int clear, int row) {
	for (int i = 0,bit=1; i < 9; i++,bit<<=1) {
		if(bit&clear)Clear_c(9 * row + i);
	}
}
void BF128::ClearCol(int clear, int col) {
	for (int i = 0, bit = 1; i < 9; i++, bit <<= 1) {
		if (bit&clear)Clear_c(9 * i + col);
	}
}
void BF128::Diag3x27(BF128 & r){
	bf.u64[0] = bf.u64[1] = 0;
	if (r.isEmpty()) return;
	for (int iband = 0; iband < 3; iband++){// build diagonal sym
		int tc[27], ntc=0;
		BitsInTable32(tc, ntc, r.bf.u32[iband], 27 * iband);
		for (int i = 0; i < ntc; i++) {
			int celld = C_transpose_d[tc[i]];
			Set_c(celld);
		}
	}
}

int BF128::Table3X27(int * r) {	
	int n = 0;
	BitsInTable32(r, n, bf.u32[0]);
	BitsInTable32(r, n, bf.u32[1],27);
	BitsInTable32(r, n, bf.u32[2], 54);
	return n;
}
int BF128::Table128(int * r) {
	int n = 0;
	BitsInTable64(r, n, bf.u64[0]);
	BitsInTable64(r, n, bf.u64[1], 64);
	return n;
}
int BF128::Table64_0(int * r) {
	int n = 0;
	BitsInTable64(r, n, bf.u64[0]);
	return n;
}
int BF128::Table64_1(int * r) {
	int n = 0;
	BitsInTable64(r, n, bf.u64[1], 64);
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

>>>>>>> refs/remotes/GPenet/skmpp2/master
