
#ifndef DIFF_TBL
#define DIFF_TBL

#include <vector>
#include <iostream>
#include <algorithm>

template<class S>
struct diff_tbl {
	std::vector<S> x;
	std::vector<S> y;
	std::vector< std::vector<S> > diff;
	int n;
	bool equispaced, typeErr, sizeErr;
	S h;

	diff_tbl(std::vector<S> xin, std::vector<S> yin);

	friend std::ostream& operator << (std::ostream & os, diff_tbl dt) {
		bool populated = (dt.diff.size() == dt.n);

		if (populated)
			os << "i\t|\tx\t|\ty\t|\t1st\t|\t2nd\t|\t3rd\t|\t...\n";
		else
			os << "i\t|\tx\t|\ty\n";

		if (dt.typeErr) return os;

		for (int i = 0; i <= dt.n; i++) {
			os << i << "\t|\t" << dt.x[i] << "\t|\t" << dt.y[i];
			if (populated) {
				for (int j = 0; j < (dt.n - i); j++) {
					if ((j < dt.diff.size()) && (i < dt.diff[j].size()))
						os << "\t|\t" << dt.diff[j][i];
					else
						os << "\t|\t_____";
				}
			}
			os << "\n";
		}
		
		return os;
	}

	void populate(bool div);
};



template<class S>
diff_tbl<S>::diff_tbl(std::vector<S> xin, std::vector<S> yin) {
	this->x = xin;
	this->y = yin;
	this->sizeErr = (xin.size() != yin.size());
	this->n = std::min(xin.size(), yin.size()) - 1;
	try {
		bool eq = true;
		S xDiff = xin[1] - xin[0], xd;
		for (int i = 1; i < n; i++) {
			xd = xin[i + 1] - xin[i];
			if (std::fabs(xDiff - xd) > (1e-13)) {
				eq = false;
				break;
			}
		}
		this->equispaced = eq;
		if (eq)
			this->h = xDiff;
		this->typeErr = false;
	}
	catch (...) {
		this->typeErr = true;
		std::cout << "Type Error occured when attempting to get difference between x's\n";
	}
}


template<class S>
void diff_tbl<S>::populate(bool div) {
	///Stirling's central, Newton's Forward and Backward differences table have the same values
	///	Calculations will be done as per Newton's Forward Difference table
	this->diff.resize(n);
	if (!typeErr) {
		std::vector<S> prevVect;
		prevVect = this->y;
		for (int i = 0; i < n; i++) {
			this->diff[i].resize(prevVect.size() - 1);
			for (int j = 0; j < diff[i].size(); j++) {
				this->diff[i][j] = prevVect[j + 1] - prevVect[j];

				if (div) {
					//Newton's Divided Difference Table
					this->diff[i][j] /= (this->x[i + j + 1] - this->x[j]);
				}

				if (std::fabs(this->diff[i][j]) < 1e-13)
					this->diff[i][j] = 0;
			}
			prevVect = this->diff[i];
		}
	}
}

#endif // !DIFF_TBL
