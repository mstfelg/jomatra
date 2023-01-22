import graph;
import math;

include style;
include "./modules/bitmanip";

// Polar coordinates
pair pol(real r, real theta) { return r*dir(theta); }

transform rot(real angle, pair z=(0,0)) {
	return rotate(degrees(angle), z);
}
// basic geometry
pair foot(pair P, pair A, pair B) {
	real s;
	s=dot(P-A,unit(B-A));
	return (scale(s)*unit(B-A)+A);
}

pair perp(pair A, pair B, pair O=(0,0)) {
	return (0,1)*(B-A)+O;
}

// Paths
pair waypoint(path p, real r) {
	return point(p,reltime(p,r));
}
pair midpoint(pair A, pair B){ return (A+B)/2;}
pair midpoint(path p){ return waypoint(p,.5);}


// TODO
pair bisect(pair A, pair O, pair B) {
	return (0,0);
}
// Perpendicular bisector
pair bisect(pair A, pair B) { return bisect(A, (A+B)/2, B); }

// Calculations
real[] angles_sss(real a, real b, real c) {
	real alpha = acos((b*b+c*c-a*a)/(2*b*c));
	real beta = acos((c*c+a*a-b*b)/(2*c*a));
	real gamma = pi-alpha-beta;
	return new real[] {alpha, beta, gamma};
}
real[] angles_sss(pair A, pair B, pair C) {
	return angles_sss(abs(B-C), abs(C-A), abs(A-B));
}

real circumradius(real a, real b, real c) {
	return a*b*c/sqrt((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c));
}
real circumradius(pair A, pair B, pair C) {
	return circumradius(abs(B-C), abs(C-A), abs(A-B));
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

// Polygons
pair[] polygon(int n) {
	pair[] gon = new pair[n];
	for (int i=0; i < n; ++i) gon[i]=expi(2pi*(i+0.5)/n-0.5*pi);
	return gon;
}

// TODO
pair[] polygon(int n, pair A, pair B) {
	real s = abs(B-A);
	real apothem = s/2*cot(pi/n);
	pair O = (A+B)/2+(0,apothem)*unit(B-A);
	real R = s/2*csc(2*pi/n);
	pair[] gon = new pair[n];
		for (int i=0; i < n; ++i) gon[i]= O + expi(2pi*(i+0.5)/n-0.5*pi);
	return gon;
}

pair[] square(pair A, pair B) { return polygon(4, A, B); }

// Square by diagonal
pair[] square_d(pair B, pair D) {
	pair M = midpoint(B,D);
	return new pair[] {rotate(90, M)*D, B, rotate(90, M)*B, D};
}

// TODO make into transform
pair mirror(pair A, pair B) { return 2B-A; }

pair bary(pair A, pair B, real x, real y) {
	return x*A + y*B;
}
// basic centers
pair bary(pair A, pair B, pair C, real x, real y, real z) {
	real k = x+y+z;
	x /= k;
	y /= k;
	z /= k;
	return x*A + y*B + z*C;
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

pair centroid(pair A, pair B, pair C) {
	return (A+B+C)/3;
}

pair orthocenter(pair A, pair B, pair C) {
	real f(real, real, real) = new real(real a, real b, real c) {
		return tan(a);
	};
	return bary(A,B,C,f);
}

pair circumcenter(pair A, pair B, pair C) {
	real f(real, real, real) = new real(real a, real b, real c) {
		return sin(2*a);
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

// Circles
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

pair[] tangent(pair P, Circ c) {
	pair O = c.O;
	real r = c.r;
	real d = abs(P-O);
	if (d == r)
		return new pair[] {P+(0,1)*unit(O-P), P+(0,-1)*unit(O-P)};
	if (d < r) {
		P = r**2/d*unit(P-O);
		d = abs(P-O);
	}

	pair A=O+r*expi(acos(r/d))*unit(P-O);
	pair B=O+r*expi(-acos(r/d))*unit(P-O);
	return new pair[] {A,B};
}

// Directed tangent
pair[] dirtangent(Circ c, Circ d, real sn=1) {
	pair A = c.O;
	pair B = d.O;
	real r = c.r;
	real s = d.r;
	real d = abs(B-A);
	if (d + s < r) {
		if (d == s)
			return new pair[] {(0,1)*unit(B-A), (0,-1)*unit(B-A)};
		B = r**2/(d**2-d*s**2)*unit(B-A);
		s *= r**2/(d**4-d**2*s**2);
		d = abs(B-A);
	}

	pair P=A+r*expi(acos((r-sn*s)/d))*unit(B-A);
	pair Q=B+sn*s*expi(acos((r-sn*s)/d))*unit(B-A);
	return new pair[] {P,Q};
}

// triangle solving

// SSS construction around circle
pair[] tri_sss(real a, real b, real c, pair O=(0,0)) {
	real[] angles = angles_sss(a,b,c);

	real alpha = angles[0];
	real beta = angles[1];
	real gamma = angles[2];
	real R = circumradius(a,b,c);

	pair A = O + expi(beta-gamma)*(0,R);
	pair B = O + expi(-pi+beta+gamma)*(0,-R);
	pair C = O + expi(pi-beta-gamma)*(0,-R);
	return new pair[] {A,B,C};
}

// SSS construction around vertex
pair[] tri_sss(pair A, real a, real b, real c, pair O=(0,0)) {
	real alpha = acos((b*b+c*c-a*a)/(2*b*c));
	real beta = acos((c*c+a*a-b*b)/(2*c*a));
	real gamma = pi-alpha-beta;

	pair B = O + rot(2*beta, O)*A;
	pair C = O + rot(-2*gamma, O)*A;
	return new pair[] {A,B,C};
}

pair[] tri_sas(real b, real alpha, real c, pair O=(0,0)) {
	return tri_sss(a=sqrt(b*b+c*c-2*b*c*cos(alpha)), b=b, c=c, O=O);
}

pair[] tri_sas(pair A, real b, real alpha, real c, pair O=(0,0)) {
	real a = sqrt(b*b+c*c-2*b*c*cos(alpha));
	real beta = acos((c*c+a*a-b*b)/(2*c*a));
	real gamma = pi - alpha - beta;
	real R = a/(2*sin(alpha));

	pair B = rot(2*beta, O) * A;
	pair C = rot(-2*gamma, O) * A;
	return new pair[] {A,B,C};
}

pair[] tri_aa(real beta, real gamma, Circ c=Circ()) {
	pair O = c.O;
	real R = c.r;
	pair A = O + expi(beta-gamma)*(0,R);
	pair B = O + expi(-pi+beta+gamma)*(0,-R);
	pair C = O + expi(pi-beta-gamma)*(0,-R);
	return new pair[] {A,B,C};
}

pair[] tri_aa(pair A, real beta, real gamma, pair O=(0,0)) {
	pair B = rot(2*beta, O) * A;
	pair C = rot(-2*gamma, O) * A;
	return new pair[] {A,B,C};
}

// Bipolar coordinates
pair bipolar(pair O1, real r1, pair O2, real r2) {
	real a = abs(O2-O1);
	pair D = ((r2**2+a**2-r1**2)*O1 + (a**2+r1**2-r2**2)*O2)/(2*a**2);
	real R = circumradius(a, r2, r1);
	return D + (0,r2*r1/2/R)*unit(O2-O1);
}

pair bipolar_aa(pair O1, real beta, pair O2, real gamma) {
	real a = abs(O1-O2);
	real alpha = pi-beta-gamma;
	real R = a/(2*sin(alpha));
	pair D = (sin(beta)*cos(gamma)*O1 + sin(gamma)*cos(beta)*O2)/sin(alpha);
	return D + (0,sin(beta)*sin(gamma)/sin(alpha))*(O2-O1);
}

pair bipolar_sss(pair O1, pair O2, real a, real b, real c) {
	return bipolar(O1, b/a*abs(O2-O1), O2, c/a*abs(O2-O1));
}
// aliases
pair join(pair A, pair B, pair C, pair D) { return extension(A,B,C,D); }
pair[] join(pair A, pair B, Circ c) { return intersectionpoints(A--B,c); }
pair[] join(path a, path b) { return intersectionpoints(a, b); }

include "./modules/decor.asy";
include "./modules/shapes.asy";
include "./modules/symbolic.asy";
