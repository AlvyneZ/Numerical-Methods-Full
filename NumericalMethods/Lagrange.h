
#ifndef LAGRANGE
#define LAGRANGE

#include "Polynomials.h"
#include "DifferenceTable.h"
#include <string>
#include <sstream>


template<class S>
polynomial<S> getLagrangePolynomial(diff_tbl<S> & dt) {

	std::cout << "The set of values given:\n" << dt;

	std::cout << "\nn = " << dt.n << "\n";

	std::vector<S> xCoeff = { 0, 1 };
	polynomial<S> x(xCoeff);


	std::string Formula = "(l_0 * y_0) + (l_1 * y_1) + (l_2 * y_2) + ... + (l_n * y_n)\n";
	Formula = Formula + "\tl_0 = {(x - x_n)...(x - x_2)(x - x_1)} / {(x_0 - x_n)...(x_0 - x_2)(x_0 - x_1)}\n";
	Formula = Formula + "\tl_1 = {(x - x_n)...(x - x_2)(x - x_0)} / {(x_1 - x_n)...(x_1 - x_2)(x_1 - x_0)}\n";
	Formula = Formula + "\t...\n";
	Formula = Formula + "\tl_n = {(x - x_(n-1))...(x - x_1)(x - x_0)} / {(x_n - x_(n-1))...(x_n - x_1)(x_n - x_0)}";

	std::cout << "\nNewton's Divided Diff Interpolating Polynomial Formula:\n" << Formula << "\n\n";

	std::vector< polynomial<S> > lagranges(dt.n + 1);
	std::vector<S> lagrangesDenominators(dt.n + 1);

	std::vector<S> start = { 1 };
	std::vector<S> partOfLCoeff = { 0, 1 };
	polynomial<S> partOfL;
	std::stringstream strstream;
	for (unsigned int i = 0; i < dt.x.size(); i++) {
		lagranges[i] = polynomial<S>(start);
		strstream.clear();
		std::cout << "l_" << i << " = {(";
		for (unsigned int j = 0; j < dt.x.size(); j++) {
			if (i != j) {
				partOfLCoeff[0] = 0 - dt.x[j];
				partOfL = polynomial<S>(partOfLCoeff);
				std::cout << partOfL << ")";
				lagranges[i] = lagranges[i] * partOfL;
				strstream << "(" << dt.x[i] << " - " << dt.x[j] << ")";
			}
		}
		std::cout << "} / {" << strstream.str() << "}\n";
		lagrangesDenominators[i] = lagranges[i].evaluate(dt.x[i]);
		std::cout << "l_" << i << " = {" << lagranges[i] << "} / " << lagrangesDenominators[i] << "\n\n";
	}

	start[0] = 0;
	polynomial<S> lagrangePoly(start);

	std::cout << "\nP_n(x) =";
	for (unsigned int i = 0; i < lagranges.size(); i++) {
		std::cout << ((i == 0) ? " " : " + ") << "{" << lagranges[i] << "} * {";
		std::cout << dt.y[i] << " / " << lagrangesDenominators[i] << "}";
		lagrangePoly = lagrangePoly + (lagranges[i] * (dt.y[i] / lagrangesDenominators[i]));
	}

	std::cout << "\n\nThus the polynomial becomes:\n";

	std::cout << "P_n(x)=";
	std::cout << lagrangePoly << "\n";
	return lagrangePoly;
}


#endif  // !LAGRANGE