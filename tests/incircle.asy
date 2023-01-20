import jomatra;
size(4cm);

pair[] tABC = acute;
pair A,B,C;
A = tABC[0]; B = tABC[1]; C = tABC[2];
pair I = incenter(tABC);

pair D = foot(I,B,C);
pair E = foot(I,C,A);
pair F = foot(I,A,B);

draw(new pair[] {D,E,F}, "$D$", "$E$", "$F$");
dot(new pair[] {I, D,E,F});
draw(tABC, "$A$", "$B$", "$C$");
label(I, "$I$", NW);
draw(incircle(A,B,C));
