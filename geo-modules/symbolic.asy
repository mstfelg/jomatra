// Symbolic calculations
bool are_concurrent(pair A, pair B, pair C, pair D, pair E, pair F) {
	return (extension(A,B,C,D) == (infinity,infinity) // parallel case
		 && (infinity,infinity)==extension(C,D,E,F))
		 || (abs(extension(A,B,C,D).x-extension(C,D,E,F).x) < 1/10^(5)
		 && abs(extension(A,B,C,D).y-extension(C,D,E,F).y) < 1/10^(5)
		 );

}

bool are_parallel(pair A, pair B, pair C, pair D) {
	return extension(A,B,C,D).x == infinity;
}

// TODO: performance improved by comparing line slopes
bool are_collinear(pair A, pair B, pair C) {
	return A == B || B == C || C == A
			|| abs(unit(C-A)-unit(A-B)) < 1/10^5
			|| abs(unit(B-A)+unit(C-A)) < 1/10^5;
}

// TODO: performance improved by Ptolomy's identity
bool are_cyclic(pair A, pair B, pair C, pair D) {
	pair O1, O2;
	O1 = circumcenter(A,B,C);
	O2 = circumcenter(A,B,D);
	return abs(O1.x-O2.x) < 1/10^(5)
		&& abs(O1.y-O2.y) < 1/10^(5);
}
