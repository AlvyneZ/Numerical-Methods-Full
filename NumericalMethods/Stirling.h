
#ifndef STIRLING
#define STIRLING

#include "Polynomials.h"
#include "DifferenceTable.h"
#include "Factorial.h"
#include <string>

namespace Stirling {
	std::string Formula = "y_(n/2) + p*[mu delta y_(n/2)] + [p^2 / 2!]*[delta^2 y_(n/2)] + [p(p^2 - 1) / 3!]*[mu delta^3 y_(n/2)] + [(p^2)(p^2 - 1) / 4!]*[delta^4 y_(n/2)] + ...";
}

template<class S>
void makeDiffTableHaveEvenN(diff_tbl<S> & dt) {
	if (dt.n % 2) {
		dt.n--;
		dt.x.resize(dt.x.size() - 1);
		dt.y.resize(dt.y.size() - 1);
		if (dt.diff.size() > 0) {
			dt.diff.resize(dt.diff.size() - 1);
			for (int i = 0; i < dt.diff.size(); i++) {
				dt.diff[i].resize(dt.diff[i].size() - 1);
			}
		}
	}
}

template<class S>
polynomial<S> getStirlingPolynomial(diff_tbl<S> & dt) {

	makeDiffTableHaveEvenN(dt);
	std::cout << "The set of values given:\n" << dt;

	dt.populate(false);
	std::cout << "\nThe difference table:\n" << dt;

	std::cout << "\nn = " << dt.n << "\n";
	std::cout << "h = " << dt.h << "\n";
	int ind = dt.n / 2;
	std::cout << "x_(n/2) = x_" << ind << " = " << dt.x[ind] << "\n";

	std::vector<S> pCoeff = { (0-dt.x[ind] / dt.h), (1 / dt.h) };
	polynomial<S> p(pCoeff);
	std::cout << "p = " << p << "\n";

	std::vector<S> diffToUse(dt.n + 1);
	diffToUse[0] = dt.y[ind];
	std::cout << "y_(n/2) = y_" << ind << " = " << diffToUse[0] << "\n";
	for (unsigned int i = 0; i < dt.diff.size(); i++) {
		//Notice: order of differences starts from 1, but i starts from 0
		//Thus even order is odd i
		if (i % 2) {
			//Even order differences (mu not needed)
			ind = (dt.diff[i].size() - 1) / 2;
			diffToUse[i + 1] = dt.diff[i][ind];
			std::cout << "delta^" << (i + 1) << " y_(n/2) = " << diffToUse[i + 1] << "\n";
		}
		else {
			//Odd order differences (mean[mu] needed)
			ind = dt.diff[i].size() / 2;
			diffToUse[i + 1] = (dt.diff[i][ind] + dt.diff[i][ind - 1]) / 2;
			std::cout << "mu delta^" << (i + 1) << " y_(n/2) = " << diffToUse[i + 1] << "\n";
		}
	}

	std::cout << "\nStirling's Central Diff Interpolating Polynomial Formula:\n" << Stirling::Formula << "\n";

	std::cout << "\nP_n(x) = " << diffToUse[0];
	for (unsigned int i = 1; i < diffToUse.size(); i++) {
		std::cout << " + {";
		for (int j = 0; j <= ((i - 1) / 2); j++) {
			if (j == 0) {
				if (i % 2) std::cout << "[" << p << "]";
				else std::cout << "[" << (p * p) << "]";
			}
			else {
				std::cout << "[" << ((p * p) - (j * j)) << "]";
			}
		}
		std::cout << " / " << factorial(i) << "} * " << diffToUse[i];
	}

	std::cout << "\n\nThus the polynomial becomes:\n";

	std::vector<S> start = {diffToUse[0]};
	polynomial<S> stirling(start);
	start[0] = 1;
	polynomial<S> calc;

	std::cout << "\nP_n(x) = " << stirling;
	for (unsigned int i = 1; i < diffToUse.size(); i++) {
		calc = polynomial<S>(start);
		for (int j = 0; j <= ((i - 1) / 2); j++) {
			if (j == 0) {
				if (i % 2) calc = calc * p;
				else calc = calc * (p * p);
			}
			else {
				calc = calc * ((p * p) - (j * j));
			}
		}
		calc = calc * (1 / S(factorial(i)));
		calc = calc * diffToUse[i];
		std::cout << " + " << calc;
		stirling = stirling + calc;
	}
	std::cout << "\n\nP_n(x)=";
	std::cout << stirling << "\n";
	return stirling;
}


#endif  // !STIRLING