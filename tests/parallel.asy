import jomatra;

pair A = (0,0); dot("$A$", A, N);
pair B = A+(2,2); dot("$B$", B, dir(B));
pair C = A+(4,-1); dot("$C$", C, dir(C));
pair D = C+B-A; dot("$D$", D, dir(D));

draw(A--B^^C--D);

/* InftyPoint P = InftyPoint(A,B); dot("$P$", P); */
pair PP = cut(A,B,C,D); dot("$PP$", PP, dir(PP));

//pair Q = (A+B+C+D)/4; dot("$Q$", Q, dir(Q));
/* dot(P.base); */
