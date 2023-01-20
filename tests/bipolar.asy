import jomatra;

dot("$O$", O, S);

pair O1 = pol(2, 180); dot("$O1$", O1, dir(O1));
pair O2 = pol(2, 0); dot("$O2$", O2, dir(O2));

real r1=5, r2=3;
real beta=pi/2, gamma=pi/3;

pair P = bipolar(O1,r1,O2,r2); dot(P);
pair Q = bipolar_aa(O1,beta,O2,gamma); dot(Q);

draw(new pair[] {P, O1, O2});
draw(new pair[] {Q, O1, O2});
