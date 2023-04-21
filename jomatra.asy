import graph;
import math;

// Style
size(4cm);
real markscalefactor=0.03;
pair O = (0,0);
pair[] acute_t = {dir(110), dir(215), dir(325)};
pair[] obtuse_t = {dir(110), dir(145), dir(35)};
pair[] right_t = {dir(110), dir(180), dir(0)};

// Helper functions
int nthBit(int a, int n) {
	return AND(a, 2^n)#(2**n);
}

int nthSgn(int a, int n)
{
	return 2*nthBit(a,n)-1;
}

// Point in infinity
struct InftyPoint {
	pair base;
	real dir;
	void operator init(pair base=(0,0), pair dest) {
		this.base = base;
		this.dir = angle(dest-base);
	}
}
InftyPoint InftyPoint(pair base=(0,0), real theta) {
	InftyPoint ip = new InftyPoint;
	ip.base = base;
	ip.dir = theta;
	return ip;
}
InftyPoint InftyPoint(pair base=(0,0), pair dest) {
	InftyPoint ip = new InftyPoint;
	ip.base = base;
	ip.dir = angle(dest-base);
	return ip;
}

InftyPoint operator init() { return InftyPoint((0,0),0.0); }
InftyPoint operator cast(pair P) { return InftyPoint((0,0),angle(P)); }
pair operator cast(InftyPoint ip) {
	return ip.base + (cos(ip.dir), sin(ip.dir));
}

// Points
// Polar coordinates
pair pol(real r, real theta) { return r*dir(theta); }
real dis(pair A, pair B) {
	return sqrt((A.x-B.x)**2 + (A.y-B.y)**2);
}


// Weights & Parametrization
pair point(pair A, pair B, real t) {
	return (1-t)*A + t*B;
}

pair waypoint(path p, real r) {
	return point(p,reltime(p,r));
}

pair midpoint(path p) { return waypoint(p,.5); }

pair waypoint(pair A, pair B, real r) {
	return point(A--B,reltime(A--B,r));
}

pair midpoint(pair A, pair B) { return (A+B)/2; }

path operator ++(path p, pair pt) {
	return p--(waypoint(p,1)+pt);
}

path operator ++(... pair[] pts) {
	path p = pts[0];
	for (int i = 1; i < pts.length; ++i)
		p = p--(pts[i-1]+pts[i]);
	return p;
}

real ratio(pair P, pair A, pair B) {
	return dis(P,A)/dis(P,B);
}

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

pair bary(pair A, pair B, real x, real y) {
	return x*A + y*B;
}

// Orthogonality
pair perp(pair A, pair B, pair O=(0,0)) {
	return (0,1)*(B-A)+O;
}

pair foot(pair P, pair A, pair B) {
	real s;
	s=dot(P-A,unit(B-A));
	return (scale(s)*unit(B-A)+A);
}

real dis(pair A, pair B, pair C) {
	return dis(A,foot(A,B,C));
}

// Lines
bool are_concurrent(pair A, pair B, pair C, pair D, pair E, pair F) {
	return (extension(A,B,C,D) == (infinity,infinity) // parallel case
		 && (infinity,infinity)==extension(C,D,E,F))
		 || (abs(extension(A,B,C,D).x-extension(C,D,E,F).x) < 1/10^(5)
		 && abs(extension(A,B,C,D).y-extension(C,D,E,F).y) < 1/10^(5)
		 );
}

bool are_collinear(pair A, pair B, pair C) {
	return A == B || B == C || C == A
			|| abs(unit(C-A)-unit(A-B)) < 1/10^5
			|| abs(unit(B-A)+unit(C-A)) < 1/10^5;
}

bool are_parallel(pair A, pair B, pair C, pair D) {
	return extension(A,B,C,D).x == infinity;
}

pair cut(pair A, pair B, pair C, pair D) {
	if (are_parallel(A,B,C,D))
		return InftyPoint((A+B+C+D)/4, angle(B-A));
	return extension(A,B,C,D);
}
pair[] cut(path a, path b) { return intersectionpoints(a, b); }

pair midline(pair A, pair B, pair C, pair D) {
	if (are_parallel(A,B,C,D))
		return (A+B+C+D)/4;
	pair P = (A+B)/2;
	pair Q = (C+D)/2;
	pair R = cut(A,B, C,D);
	real a = dis(R,P);
	real b = dis(R,Q);
	return a*Q/(a+b) + b*P/(a+b);
}

// Angles
real angle(pair A, pair B, pair O=(0,0)) { return angle(A-O)-angle(B-O); }
real angle_d(pair A, pair B, pair O=(0,0)) { return degrees(angle(A,B,O)); }

transform rot(real angle, pair z=(0,0)) {
	return rotate(degrees(angle), z);
}
pair mirror(pair A, pair B) { return 2B-A; }

// Parametrized directed arc
pair arc_param(pair A, pair B, real t=0, pair O=(0,0)) {
	pair BB = unit(B-O);
	pair AA = unit(A-O);
	pair C = unit(BB-AA);
	real ang = angle(BB,AA);
	return rot(-pi/2+.5*ang*t,O) * C;
}

// Parametrized symmetric (non-directed) arc
pair arc_symparam(pair A, pair B, real t=0, pair O=(0,0)) {
	pair BB = unit(B-O);
	pair AA = unit(A-O);
	pair C = O+unit(BB+AA);
	real ang = angle(BB,AA);
	return rot(ang*t,O) * C;
}

pair mid_arc(pair A, pair B, pair O=(0,0)) {
	return arc_param(A,B,O);
}

pair[] bisect(pair A, pair O, pair B) {
	pair inBisector = arc_symparam(A,B,O);
	pair outBisector = rot(pi/2, O)*inBisector;
	return new pair[] {inBisector, outBisector};
}

pair[] bisect(pair A, pair B) { return bisect(A, B, (A+B)/2); }

// Triangles
include "./geo-modules/centers.asy";
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
pair[] cut(pair A, pair B, Circ c) { return intersectionpoints(A--B,c); }

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
// Used for constructing triangles on segments.
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

pair isosceles(pair A, pair B, real theta=pi/3) {
	return bipolar_aa(A, .5*pi-.5*theta, B, .5*pi-.5*theta);
}

pair bipolar_sss(pair O1, pair O2, real a, real b, real c) {
	return bipolar(O1, b/a*abs(O2-O1), O2, c/a*abs(O2-O1));
}

// Tangency
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


// Polygons
pair[] polygon(int n) {
	pair[] gon = new pair[n];
	for (int i=0; i < n; ++i) gon[i]=expi(2pi*(i+0.5)/n-0.5*pi);
	return gon;
}

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

include "./geo-modules/decor.asy";
include "./geo-modules/apollonius.asy";
