import jomatra;

pair[] tABC = acute_t; draw(tABC);
pair A = tABC[0]; dot("$A$", A, dir(A));
pair B = tABC[1]; dot("$B$", B, dir(B));
pair C = tABC[2]; dot("$C$", C, dir(C));

Circ c = LLL(A,B,C,3); draw(c);
