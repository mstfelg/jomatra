circle PPP(point A, point B, point C) { return circle(A,B,C); }

circle LLL(point A, point B, point C, int n=0) {
	point bis_B = bisect(A,B,C)[nthBit(n,0)];
	point bis_C = bisect(A,C,B)[nthBit(n,1)];
	point I = cut(B,bis_B, C,bis_C);
	real r = dis(I,B,C);
	return circle(I,r);
}

circle LPP(point lp, point lq, point A, point B, int n=0) {
	if (lp == lq)
		return circle(lp,A,B);
	if (are_parallel(lp,lq,A,B)) {
		point[] bis = bisect(B,A);
		point T = cut(bis[0], bis[1], lp, lq);
		return circle(T,A,B);
	}
	point P = cut(lp,lq,A,B);
	real ir = sqrt(dis(P,A)*dis(P,B));
	point T = P+nthSgn(n,0)*ir*unit(P-lp);
	return circle(T,A,B);
}
