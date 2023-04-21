// DECORATION

/* void dirty_label(pair[] gon) { */
/* 	for (int i = 0; i < gon.length; ++i) */
/* 		dot() */
/* } */

/* void draw(... pair[] gon) { */
/* 	path p; */
/* 	for (int i = 0; i < gon.length; ++i) */
/* 		p = p--gon[i]; */
/* 	p = p--cycle; */
/* 	draw(p); */
/* } */

void draw(Circ c) {
	draw(circle(c.O, c.r));
}

void draw(pair[] gon ... Label[] labels) {
	path p;
	for (int i = 0; i < gon.length; ++i)
		p = p--gon[i];
	p = p--cycle;
	draw(p);
	dot(labels, p);
}

void draw(pair O, pair[] Ps) {
	for (int i = 0; i < Ps.length; ++i)
		draw(O--Ps[i]);
}

// returns interval [a, b] on number line A=(-1,0),B=(1,0)
path segment(pair A, pair B, real a, real b=-a) {
	pair AA = (1-a)*A/2 + (1+a)*B/2;
	pair BB = (1-b)*A/2 + (1+b)*B/2;
	return AA--BB;
}

path rightanglemark(pair A, pair B, pair C, real s=8) {
	pair P,Q,R;
	P=s*markscalefactor*unit(A-B)+B;
	R=s*markscalefactor*unit(C-B)+B;
	Q=P+R-B;
	return P--Q--R;
}

path anglemark(pair A, pair B, pair C, real t=8 ... real[] s) {
	pair M,N,P[],Q[];
	path mark;
	int n=s.length;
	M=t*markscalefactor*unit(A-B)+B;
	N=t*markscalefactor*unit(C-B)+B;
	for (int i=0; i<n; ++i)
	{
		P[i]=s[i]*markscalefactor*unit(A-B)+B;
		Q[i]=s[i]*markscalefactor*unit(C-B)+B;
	}
	mark=arc(B,M,N);
	for (int i=0; i<n; ++i)
	{
		if (i%2==0)
		{
			mark=mark--reverse(arc(B,P[i],Q[i]));
		}
		else
		{
			mark=mark--arc(B,P[i],Q[i]);
		}
	}
	if (n%2==0 && n!=0)
	mark=(mark--B--P[n-1]);
	else if (n!=0)
	mark=(mark--B--Q[n-1]);
	else mark=(mark--B--cycle);
	return mark;
}

picture pathticks(path g, int n=1, real r=.5, real spacing=6, real s=8, pen p=currentpen) {
	picture pict;
	pair A,B,C,direct;
	real t,l=arclength(g), space=spacing*markscalefactor, halftick=s*markscalefactor/2, startpt;
	if (n>0)
	{
		direct=unit(dir(g,arctime(g,r*l)));
		startpt=r*l-(n-1)/2*space;
		for (int i=0; i<n; ++i)
		{
			t=startpt+i*space;
			B=point(g,arctime(g,t))+(0,1)*halftick*direct;
			C=B+2*(0,-1)*halftick*direct;
			draw(pict,B--C,p);
		}
	}
	return pict;
}

// Ticks a point P with length d in dir theta.
path tick(pair P, real d, real theta) {
	return (P+pol(d,theta))--(P-pol(d,theta));
}

// Endpoint ticks
path[] endtick(pair A, pair B, real d, real theta=pi/2) {
	theta += angle(B-A);
	return new path[] {tick(A, d, theta), tick(B, d, theta)};
}

path[] arcendtick(pair O, pair A, pair B, real d, real theta=0) {
	return new path[] {tick(A, d, theta+angle(A-O)), tick(B, d, theta+angle(B-O))};
}

// Ticks around a point
// creates n ticks over the interval of intrv units with direction and tick
// length of d
path[] pointticks(pair P, real d, real theta, int n, real intrv) {
	path[] ticks = new path[n];
	for (int i = -n#2; i < n#2; ++i)
		ticks[i+n#2] = tick(P + pol(i*intrv/n, theta+pi/2), d, theta);
	return ticks;
}

// Endpoints are not ticked
path[] nticks(pair A, pair B, real d, int n) {
	path[] ticks = new path[n];
	real theta = angle(B-A) + pi/2;
	for (int i = 1; i < n; ++i)
		ticks[i-1] = tick(A + (B-A)*i/n, d, theta);
	return ticks;
}
