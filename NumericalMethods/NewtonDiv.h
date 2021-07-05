
#ifndef NEWTON_DIV
#define NEWTON_DIV

#include "Polynomials.h"
#include "DifferenceTable.h"
#include <string>

namespace Newton {
	std::string Divided = "y_0 + (x - x_0)*f[x_0, x_1] + [(x - x_0)(x - x_1)]*f[x_0, x_1, x_2] + ... + [(x - x_0)(x - x_1)...(x - x_(n-1))]*f[x_0, x_1, ..., x_n]";
}

template<class S>
polynomial<S> getNewtonDivPolynomial(diff_tbl<S> & dt) {

	std::cout << "The set of values given:\n" << dt;

	dt.populate(true);
	std::cout << "\nThe difference table:\n" << dt;

	std::cout << "\nn = " << dt.n << "\n";

	std::vector<S> xCoeff = { 0, 1 };
	polynomial<S> x(xCoeff);

	std::vector<S> diffToUse(dt.n + 1);
	diffToUse[0] = dt.y[0];
	std::cout << "y_0 = y_" << diffToUse[0] << "\n";
	for (unsigned int i = 0; i < dt.diff.size(); i++) {
		diffToUse[i + 1] = dt.diff[i][0];
		std::cout << "f[x_0, x_1";
		for (int j = 2; j < (i + 2); j++) {
			std::cout << ", x_" << j;
		}
		std::cout << "] = " << diffToUse[i + 1] << "\n";
	}


	std::cout << "\nNewton's Divided Diff Interpolating Polynomial Formula:\n" << Newton::Divided << "\n";

	std::cout << "\nP_n(x) = " << diffToUse[0];
	for (unsigned int i = 1; i < diffToUse.size(); i++) {
		std::cout << " + {";
		for (int j = 0; j <= (i - 1); j++) {
			std::cout << "[" << (x - dt.x[j]) << "]";
		}
		std::cout << " * " << diffToUse[i] << "}";
	}

	std::cout << "\n\nThus the polynomial becomes:\n";

	std::vector<S> start = { diffToUse[0] };
	polynomial<S> newton(start);
	start[0] = 1;
	polynomial<S> calc;

	std::cout << "\nP_n(x) = " << newton;
	for (unsigned int i = 1; i < diffToUse.size(); i++) {
		calc = polynomial<S>(start);
		for (int j = 0; j <= (i - 1); j++) {
			calc = calc * (x - dt.x[j]);
		}
		calc = calc * diffToUse[i];
		std::cout << " + " << calc;
		newton = newton + calc;
	}
	std::cout << "\n\nP_n(x)=";
	std::cout << newton << "\n";
	return newton;
}


#endif  // !NEWTON_DIV