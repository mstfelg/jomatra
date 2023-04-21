int nthBit(int a, int n) {
	return AND(a, 2^n)#(2**n);
}
int nthSgn(int a, int n) { return 2*nthBit(a,n)-1; }
