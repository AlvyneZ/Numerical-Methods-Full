
#ifndef POLYNOMIAL
#define POLYNOMIAL

#include <vector>
#include <math.h>
#include <iostream>
#include <algorithm>

template<class S>struct polynomial {
	std::vector<S>coeffecients;
	polynomial(std::vector<S> arr) :coeffecients(arr) {};
	polynomial() {};
	void operator = (polynomial<S> rhs) {
		coeffecients = rhs.coeffecients;
	}
	polynomial<S> operator * (polynomial<S> rhs);
	polynomial<S> operator * (S num);
	polynomial<S> operator +(polynomial<S> rhs);
	polynomial<S> operator + (S num);
	polynomial<S> operator - (S num);
	friend std::ostream& operator << (std::ostream & os, polynomial p) {
		int my_len = p.coeffecients.size();
		for (int i = (my_len - 1); i >= 0; --i) {
			os << p.coeffecients[i] << " x^" << i << (i == 0 ? " . " : " + ");
		}
		//os << "\n";
		return os;
	}
	S evaluate(S valueOfX);
	polynomial<S> differentiate();
};


/*
	Given a list of coeffecients for a polynomial , a_0 X^0 + a_1 x^1 + a_2 x ^ 2 +.., in the order a_1, a_2,a_3 stores the polynomial with the properties * which is polynomial multiplication with another polynomial and + which is polynomial addition
	Example Usage :
	vector<double>num = {1.2,3.6,4.5, 8.9,9.1};
	polynomial<double>nums(num);
	auto nums2 = nums * nums;
	auto 2nums2 = nums2 + nums2;
	cout << twonums2;

	outputs 2.88 x ^ 0 + 17.28 x ^ 1 + 47.52 x ^ 2 + 107.52 x ^ 3 + 212.34 x ^ 4 + 291.24 x ^ 5 + 322.22 x ^ 6 + 323.96 x ^ 7 + 165.62 x ^ 8 + 0 x ^ 9
	which is correct by wolfram alpha
*/




template<class S>
polynomial<S> polynomial<S>:: operator*(polynomial<S>rhs) {
	int her_len = rhs.coeffecients.size();
	int my_len = this->coeffecients.size();
	std::vector<S>sol(my_len + her_len - 1);
	for (int i = 0; i < my_len; ++i) {
		for (int j = 0; j < her_len; ++j) {
			sol[i + j] += this->coeffecients[i] * rhs.coeffecients[j];
		}
	}
	polynomial<S>result(sol);
	return result;
}
template<class S>
polynomial<S> polynomial<S>:: operator + (polynomial<S>rhs) {
	int her_len = rhs.coeffecients.size();
	int my_len = this->coeffecients.size();
	std::vector<S> sol(std::max(my_len, her_len), 0);
	for (int i = 0; i < her_len; ++i) {
		sol[i] += rhs.coeffecients[i];
	}
	for (int i = 0; i < my_len; ++i) {
		sol[i] += this->coeffecients[i];
	}
	polynomial<S>result(sol);
	return result;
}
template<class S>
polynomial<S> polynomial<S>:: operator * (S num) {
	std::vector<S> sol(this->coeffecients.size());
	for (int i = 0; i < (int)this->coeffecients.size(); ++i) {
		sol[i] = this->coeffecients[i] * num;
	}
	polynomial<S> result(sol);
	return result;
}
template<class S>
polynomial<S> polynomial<S>:: operator + (S num) {
	std::vector<S> sol(this->coeffecients.size());
	for (int i = 0; i < (int)this->coeffecients.size(); ++i) {
		sol[i] = this->coeffecients[i];
	}
	sol[0] += num;
	polynomial<S> result(sol);
	return result;
}
template<class S>
polynomial<S> polynomial<S>:: operator - (S num) {
	std::vector<S> sol(this->coeffecients.size());
	for (int i = 0; i < (int)this->coeffecients.size(); ++i) {
		sol[i] = this->coeffecients[i];
	}
	sol[0] -= num;
	polynomial<S> result(sol);
	return result;
}

template<class S>
S polynomial<S>::evaluate(S valueOfX) {
	S OGValueOfX = valueOfX;
	S valueOfY = this->coeffecients[0];
	for (unsigned int i = 1; i < this->coeffecients.size(); i++) {
		valueOfY += this->coeffecients[i] * valueOfX;
		valueOfX *= OGValueOfX;
	}
	return valueOfY;
}

template<class S>
polynomial<S> polynomial<S> ::differentiate() {
	std::vector<S> newpoly;
	for (int i = 1; i < (int)this->coeffecients.size(); ++i) {
		newpoly.push_back(this->coeffecients[i] * S(i));
	}
	polynomial<S>derivative(newpoly);
	return derivative;
}




#endif // !POLYNOMIAL