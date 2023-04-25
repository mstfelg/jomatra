import graph;
// Radian
path carc(pair O, real r, real alpha, real beta) {
	return arc(O, r, degrees(alpha), degrees(beta));
}
path carc(pair O, real r, real alpha) {
	return carc(O, r, 0, alpha);
}
path carc(circle c, real alpha, real beta) {
	return arc(c.C, c.r, degrees(alpha), degrees(beta));
}
path carc(circle c, real alpha) {
	return carc(c, 0, alpha);
}

bool are_cyclic(pair A, pair B, pair C, pair D) {
	pair O1, O2;
	O1 = circumcenter(A,B,C);
	O2 = circumcenter(A,B,D);
	return abs(O1.x-O2.x) < 1/10^(5)
		&& abs(O1.y-O2.y) < 1/10^(5);
}

pair[] circumscribe(pair O, real r, pair[] pts) {
	if (pts.length == 2) {
		return null;
	}
	real R = circumradius(pts[0], pts[1], pts[2]);
	pair O2 = circumcenter(pts[0], pts[1], pts[2]);
	return shift(O-O2) * scale(r/R) * pts;
}

path circ_arc(pair A, pair B, pair C) {
	return arc(circumcenter(A,B,C), A, C);
}
