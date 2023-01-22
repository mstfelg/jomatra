import jomatra;
size(4cm);

pair B = (-1,0); dot("$B$", B, dir(B));
pair C = (1,0); dot("$C$", C, dir(C));

Circ b = Circ(B,C);
Circ c = Circ(C,B);

// Paths
path lens(Circ a, Circ b) {
	pair[] ip = intersectionpoints(a,b);
	pair A = a.O;
	pair B = b.O;
	return arc(A, ip[1], ip[0])&arc(B, ip[0], ip[1])&cycle;
}

// Shading
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
	fill(ineq_ab, lens(a,b), p);
	return ineq_ab;
}

picture circ_nand(Circ a, Circ b, pen p=currentpen) {
	picture ineq_ab;
	filloutside(ineq_ab, lens(a, b), p);
	return ineq_ab;
}

picture circ_xnor(Circ a, Circ b, pen p=currentpen) {
	picture ineq_ab;
	fill(ineq_ab, lens(a,b), p);
	filloutside(ineq_ab, buildcycle((path)a,(path)b), p);
	return ineq_ab;
}

picture circ_xor(Circ a, Circ b, pen p=currentpen) {
	picture ineq_a, ineq_b;
	fill(ineq_a, buildcycle((path)a,reverse(b)), p);
	fill(ineq_b, buildcycle(reverse(a),(path)b), p);
	add(ineq_a, ineq_b);
	return ineq_a;
}

/* add(circ_and(b,c, blue)); */
/* add(circ_or(b,c, blue)); */
/* add(circ_xnor(b,c, blue)); */

//TODO
/* add(circ_xor(b,c, blue)); */
/* add(circ_nand(b,c, blue)); */
