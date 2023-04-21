real circumradius(real a, real b, real c) {
	return a*b*c/sqrt((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c));
}
real circumradius(pair A, pair B, pair C) {
	return circumradius(abs(B-C), abs(C-A), abs(A-B));
}
pair circumcenter(pair A, pair B, pair C) {
	real f(real, real, real) = new real(real a, real b, real c) {
		return sin(2*a);
	};
	return bary(A,B,C,f);
}
real exradius(real a, real b, real c) {
	real s=(a+b+c)/2;
	return sqrt(s*(s-b)*(s-c)/(s-a));
}
real exradius(pair A, pair B, pair C) {
	return exradius(abs(B-C), abs(C-A), abs(A-B));
}
real inradius(real a, real b, real c) {
	real s=(a+b+c)/2;
	return sqrt(s*(s-a)*(s-b)*(s-c))/s;
}
real inradius(pair A, pair B, pair C) {
	return inradius(abs(B-C), abs(C-A), abs(A-B));
}

pair centroid(pair A, pair B, pair C) {
	return (A+B+C)/3;
}

pair centroid(pair[] polygon) {
	pair G = (0,0);
	for (int i = 0; i < polygon.length; ++i) {
		G += polygon[i];
	}
	return G/polygon.length;
}

pair orthocenter(pair A, pair B, pair C) {
	real f(real, real, real) = new real(real a, real b, real c) {
		return tan(a);
	};
	return bary(A,B,C,f);
}

pair excenter(pair A, pair B, pair C) {
	real[] ang = angles_sss(A,B,C);
	return bary(A,B,C, -sin(ang[0]), sin(ang[1]), sin(ang[2]));
}

pair incenter(pair A, pair B, pair C) {
	real f(real, real, real) = new real(real a, real b, real c) {
		return sin(a);
	};
	return bary(A,B,C,f);
}

pair incenter(pair[] tABC) { return incenter(tABC[0], tABC[1], tABC[2]); }

pair symmedian(pair A, pair B, pair C) {
	real[] ang = angles_sss(A,B,C);
	return bary(A,B,C, -sin(ang[0])**2, sin(ang[1])**2, sin(ang[2])**2);
}
pair lemoine(pair A, pair B, pair C) {
	real f(real, real, real) = new real(real a, real b, real c) {
		return sin(a)**2;
	};
	return bary(A,B,C,f);
}
