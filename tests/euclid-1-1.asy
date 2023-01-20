import jomatra;

dot("$O$", O, S);

pair B = (1,0); dot("$B$", B, dir(B));
pair C = (-1,0); dot("$C$", C, dir(C));

path c = Circ(B,C); draw(c);
path d = Circ(C,B); draw(d);

pair A = intersectionpoints(c,d)[0]; dot("$A$", A, dir(A));

draw(new pair[] {A,B,C});
