import jomatra;
size(4cm);

pair B = (-1,0); dot("$B$", B, dir(B));
pair C = (1,0); dot("$C$", C, dir(C));

Circ b = Circ(B,C); draw(b);
Circ c = Circ(C,B); draw(c);

picture circ_or(Circ a, Circ b, pen p=currentpen) {
	picture ineq_a;
	picture ineq_b;
	fill(ineq_a, (path)a, p);
	fill(ineq_b, (path)b, p);
	add(ineq_a, ineq_b);
	return ineq_a;
}

picture circ_and(Circ a, Circ b, pen p=currentpen) {
	picture ineq_ab;
	fill(ineq_ab, (path)a, p);
	clip(ineq_ab, (path)b, p);
	return ineq_ab;
}

picture circ_xor(Circ a, Circ b, pen p=currentpen) {
	picture ineq_ab;
	filloutside(ineq_ab, buildcycle((path)a,(path)b), p);
	return ineq_ab;
}

/* add(circ_and(b,c, blue)); */
/* add(circ_or(b,c, blue)); */
add(circ_xor(b,c, blue));
