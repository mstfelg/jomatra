import graph;
// Barycentric coordinates
pair bary(pair A, pair B, pair C, real x, real y, real z) {
	real k = x+y+z;
	x /= k;
	y /= k;
	z /= k;
	return x*A + y*B + z*C;
}

real[] angles_sss(real a, real b, real c) {
	real alpha = acos((b*b+c*c-a*a)/(2*b*c));
	real beta = acos((c*c+a*a-b*b)/(2*c*a));
	real gamma = pi-alpha-beta;
	return new real[] {alpha, beta, gamma};
}

real[] angles_sss(pair A, pair B, pair C) {
	return angles_sss(abs(B-C), abs(C-A), abs(A-B));
}

pair bary(pair A, pair B, pair C, real f(real, real, real)) {
	real[] angles = angles_sss(A,B,C);
	real aa = angles[0];
	real bb = angles[1];
	real cc = angles[2];
	real x = f(aa, bb, cc);
	real y = f(bb, cc, aa);
	real z = f(cc, aa, bb);
	real k = x+y+z;
	x /= k;
	y /= k;
	z /= k;
	return x*A + y*B + z*C;
}

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

struct Circ {
	pair O;
	real r;
	void operator init(pair O, real r) {
		this.O=O;
		this.r=r;
	}
	Circ operator init() { return Circ.Circ((0,0), 1.0); }
}
Circ Circ(pair O, real r) {
	Circ c = new Circ;
	c.O = O;
	c.r = r;
	return c;
}
Circ operator cast(pair P) { return Circ(P, 0); }
pair operator cast(Circ c) { return c.O; }
path operator cast(Circ c) { return Circle(c.O, c.r); }
Circ Circ(pair P) { return Circ(P, 1.0); }
Circ Circ(pair P, pair Q) { return Circ(P, abs(Q-P)); }
Circ Circ() { return Circ((0,0)); }

// fixes complaints about ambiguity
Circ Circ(pair P, int r) { return Circ(P, (real)r); }

// Radian
path carc(pair O, real r, real alpha, real beta) {
	return arc(O, r, degrees(alpha), degrees(beta));
}
path carc(pair O, real r, real alpha) {
	return carc(O, r, 0, alpha);
}
path carc(Circ c, real alpha, real beta) {
	return arc(c.O, c.r, degrees(alpha), degrees(beta));
}
path carc(Circ c, real alpha) {
	return carc(c, 0, alpha);
}
pair invert(pair P, Circ c) {
	return c.O + unit(P-c.O)*c.r**2/(abs(P-c.O)**2);
}

Circ Circ(pair A, pair B) { return Circ(A, abs(B-A)); }
Circ Circ(pair A, pair B, pair C) {
	return Circ(circumcenter(A,B,C),circumradius(A,B,C));
}
Circ incircle(pair A, pair B, pair C) {
	return Circ(incenter(A,B,C),inradius(A,B,C));
}
Circ excircle(pair A, pair B, pair C) {
	return Circ(excenter(A,B,C),exradius(A,B,C));
}
Circ Circ(pair A, pair B, pair C) {
	return Circ(circumcenter(A,B,C),circumradius(A,B,C));
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
