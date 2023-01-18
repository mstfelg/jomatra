pair[] rectangle(pair A, pair B) {
	return new pair[] {A, (xpart(B), ypart(A)), B, (xpart(A), ypart(B))};
}

pair[] rectangle(pair A, pair B, real w) {
	pair u = (0,w)*unit(B-A);
	return new pair[] { A+u, B+u, B-u, A-u };
}

pair[] rod(pair A, pair B, real w) { return rectangle(A, B, w); }
