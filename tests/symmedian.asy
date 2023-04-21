import jomatra;

pair[] tABC = acute_t; draw(tABC);
pair A = tABC[0]; dot("$A$", A, dir(A));
pair B = tABC[1]; dot("$B$", B, W);
pair C = tABC[2]; dot("$C$", C, E);

Circ c = Circ(O,1); draw(c);

pair[] tb = tangent(B,c);
pair[] tc = tangent(C,c);

pair S = symmedian(B,C,A); dot("$S$", S, dir(S));
draw(S--B);
draw(S--C);

pair L = lemoine(A,B,C); dot("$L$", L, dir(L));
draw(segment(S,A, -1, 1.3));
draw(segment(S,B, -1, 1.3));
draw(segment(S,C, -1, 1.3));
