
#ifndef NEWTON
#define NEWTON

#include "Polynomials.h"
#include "DifferenceTable.h"
#include "Factorial.h"
#include <string>

namespace Newton {
	std::string Forward = "y_0 + p*[del y_0] + [p(p - 1) / 2!]*[del^2 y_0] + [p(p - 1)(p - 2) / 3!]*[del^3 y_0] + ... + [p(p - 1)(p - 2)...(p - (n-1)) / 4!]*[del^n y_0]";
	std::string Backward = "y_n + p*[nabla y_n] + [p(p + 1) / 2!]*[nabla^2 y_n] + [p(p + 1)(p + 2) / 3!]*[nabla^3 y_n] + ... + [p(p + 1)(p + 2)...(p + (n-1)) / 4!]*[nabla^n y_n]";
}

template<class S>
polynomial<S> getNewtonPolynomial(bool forward, diff_tbl<S> & dt) {

	std::cout << "The set of values given:\n" << dt;

	dt.populate(false);
	std::cout << "\nThe difference table:\n" << dt;

	std::cout << "\nn = " << dt.n << "\n";
	std::cout << "h = " << dt.h << "\n";
	int ind = (forward) ? 0 : dt.n;
	std::cout << "x_" << ((forward) ? "0" : "n") << " = x_" << ind << " = " << dt.x[ind] << "\n";

	std::vector<S> pCoeff = { (0 - dt.x[ind] / dt.h), (1 / dt.h) };
	polynomial<S> p(pCoeff);
	std::cout << "p = " << p << "\n";

	std::vector<S> diffToUse(dt.n + 1);
	diffToUse[0] = dt.y[ind];
	std::cout << "y_" << ((forward) ? "0" : "n") << " = y_" << ind << " = " << diffToUse[0] << "\n";
	for (unsigned int i = 0; i < dt.diff.size(); i++) {
		if (!forward)
			ind = dt.diff[i].size() - 1;
		diffToUse[i + 1] = dt.diff[i][ind];
		std::cout << ((forward) ? "del" : "nabla") << "^" << (i + 1) << " y_" << ((forward) ? "0" : "n") << " = " << diffToUse[i + 1] << "\n";
	}


	std::cout << "\nNewton's " << ((forward) ? "Forward" : "Backward") << " Diff Interpolating Polynomial Formula:\n" << ((forward) ? Newton::Forward : Newton::Backward) << "\n";

	std::cout << "\nP_n(x) = " << diffToUse[0];
	for (unsigned int i = 1; i < diffToUse.size(); i++) {
		std::cout << " + {";
		for (int j = 0; j <= (i - 1); j++) {
			if (forward) {
				std::cout << "[" << (p - j) << "]";
			}
			else {
				std::cout << "[" << (p + j) << "]";
			}
		}
		std::cout << " / " << factorial(i) << "} * " << diffToUse[i];
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
			if (forward) {
				calc = calc * (p - j);
			}
			else {
				calc = calc * (p + j);
			}
		}
		calc = calc * (1 / S(factorial(i)));
		calc = calc * diffToUse[i];
		std::cout << " + " << calc;
		newton = newton + calc;
	}
	std::cout << "\n\nP_n(x)=";
	std::cout << newton << "\n";
	return newton;
}


#endif  // !NEWTON