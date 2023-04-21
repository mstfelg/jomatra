import jomatra;

pair O = (0,0); dot("$O$", O, S);

pair A = pol(2,10); dot("$A$", A, dir(A));
pair B = pol(1,20); dot("$B$", B, dir(B));

/* pair C = arc_param(B,A); dot("$C$", C, dir(C)); */
/* label(angle_d(C,B)); */

pair[] bis_pts = bisect(A,O,B);
pair C = bis_pts[0]; dot("$C$", C, dir(C));
pair D = bis_pts[1]; dot("$D$", D, dir(D));

draw(A--O--B^^O--C^^O--D);
