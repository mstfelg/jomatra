Circ PPP(pair A, pair B, pair C) { return Circ(A,B,C); }

Circ LLL(pair A, pair B, pair C, int n=0) {
	pair bis_B = bisect(A,B,C)[nthBit(n,0)];
	pair bis_C = bisect(A,C,B)[nthBit(n,1)];
	pair I = cut(B,bis_B, C,bis_C);
	real r = dis(I,B,C);
	return Circ(I,r);
}

Circ LPP(pair lp, pair lq, pair A, pair B, int n=0) {
	if (lp == lq)
		return Circ(lp,A,B);
	if (are_parallel(lp,lq,A,B)) {
		pair[] bis = bisect(B,A);
		pair T = cut(bis[0], bis[1], lp, lq);
		return Circ(T,A,B);
	}
	pair P = cut(lp,lq,A,B);
	real ir = sqrt(dis(P,A)*dis(P,B));
	pair T = P+nthSgn(n,0)*ir*unit(P-lp);
	return Circ(T,A,B);
}
