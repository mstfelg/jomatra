real exradius(real a, real b, real c) {
	real s=(a+b+c)/2;
	return sqrt(s*(s-b)*(s-c)/(s-a));
}
real exradius(point A, point B, point C) {
	return exradius(abs(B-C), abs(C-A), abs(A-B));
}
real inradius(real a, real b, real c) {
	real s=(a+b+c)/2;
	return sqrt(s*(s-a)*(s-b)*(s-c))/s;
}
real inradius(point A, point B, point C) {
	return inradius(abs(B-C), abs(C-A), abs(A-B));
}

point centroid(point A, point B, point C) {
	return (A+B+C)/3;
}

point centroid(point[] polygon) {
	point G = (0,0);
	for (int i = 0; i < polygon.length; ++i) {
		G += polygon[i];
	}
	return G/polygon.length;
}

point orthocenter(point A, point B, point C) {
	real f(real, real, real) = new real(real a, real b, real c) {
		return tan(a);
	};
	return bary(A,B,C,f);
}

point excenter(point A, point B, point C) {
	real[] ang = angles_sss(A,B,C);
	return bary(A,B,C, -sin(ang[0]), sin(ang[1]), sin(ang[2]));
}

point incenter(point A, point B, point C) {
	real f(real, real, real) = new real(real a, real b, real c) {
		return sin(a);
	};
	return bary(A,B,C,f);
}

point incenter(point[] tABC) { return incenter(tABC[0], tABC[1], tABC[2]); }

point symmedian(point A, point B, point C) {
	real[] ang = angles_sss(A,B,C);
	return bary(A,B,C, -sin(ang[0])**2, sin(ang[1])**2, sin(ang[2])**2);
}
point lemoine(point A, point B, point C) {
	real f(real, real, real) = new real(real a, real b, real c) {
		return sin(a)**2;
	};
	return bary(A,B,C,f);
}
