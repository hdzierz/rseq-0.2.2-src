#ifndef ISOFORM_OPT_H
///Define this macro to prevent from including this header file more than once.
#define ISOFORM_OPT_H

#include "math_utils.h"
//#include "exon_junction_extractor.h"
//#include "gsnake.h"

#define MAX_ITR 10000
bool do_EM = false;
double lambda = 0;

/*double sum_matrix(const vector<vector<double> > &matrix) {
	double result = 0;
	for (int i = 0; i < (int)matrix.size(); i++)
		for (int j = 0; j < (int)matrix[i].size(); j++)
			result += matrix[i][j];
	return result;
}*/

/*vector<vector<double> > inv_matrix(const vector<vector<double> > &matrix) {
	int m = (int)matrix.size(), n = (int)matrix[0].size();
	MATRIX input, *output;
	input.init((short)m, (short)n);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			input.put(i, j, matrix[i][j]);
	output = input.inverse();
	vector<vector<double> > result;
	result.resize(m);
	for (int i = 0; i < m; i++) {
		result[i].resize(m);
		for (int j = 0; j < m; j++)
			result[i][j] = output->get(i,j);
	}
	delete output;
	return result;
}*/

class isoform_opt{
public:
	bool verbose;
	int M; //number of categories
	int K; //number of isoforms
// for parameters
//	vector<double> N; //category read counts
//	vector<vector<double> > A; //category sampling rate matrix; A[i][j] is the sampling rate for isoform i and category j

	isoform_opt() {verbose = false;}

	double log_likelihood(const vector<vector<double> > &A, const vector<double> &N, const vector<double> &X, vector<double> *W = NULL) {
		double result = 0;
		double temp1 = 0;
		if (W == NULL) {
			for (int i = 0; i < M; i++) {
				double temp2 = 0;
				for (int j = 0; j < K; j++) {
					temp2 += X[j]*A[j][i];
				}
				temp1 += temp2;
			}
		} else temp1 = inner_prod(*W, X);
		result -= temp1;
		temp1 = 0;
		for (int i = 0; i < M; i++) {
			double temp2 = 0;
			for (int j = 0; j < K; j++) {
				temp2 += X[j] * A[j][i];
			}
			temp1 += log(temp2 + epsilon) * N[i];
		}
		result += temp1;
		return result;
	}

	double grad_log_likelihood(const vector<vector<double> > &A, const vector<double> &N, const vector<double> &X, int k, vector<double> *W = NULL) {
		double result = 0;
		double temp1 = 0;
		if (W == NULL) {
			for (int i = 0; i < M; i++) {
				temp1 += A[k][i];
			}
		} else temp1 = (*W)[k];
		result -= temp1;
		temp1 = 0;
		for (int i = 0; i < M; i++) {
			double temp2 = 0;
			for (int j = 0; j < K; j++) {
				temp2 += X[j] * A[j][i];
			}
			temp1 += N[i] * A[k][i] / (temp2 + epsilon);
		}
		result += temp1;
		return result;
	}

	vector<double> grad_log_likelihood(const vector<vector<double> > &A, const vector<double> &N, const vector<double> &X, vector<double> *W = NULL) {
		vector<double> grad;
		grad.resize(K);
		for (int i = 0; i < K; i++) {
			grad[i] = grad_log_likelihood(A, N, X, i, W);
		}
		return grad;
	}

	vector<vector<double> > hessian_log_likelihood(const vector<vector<double> > &A, const vector<double> &N, const vector<double> &X) {
		vector<vector<double> > hessian;
		hessian.resize(K);
		for (int i = 0; i < K; i++) hessian[i].resize(K);
		vector<double> temp1;
		temp1.resize(M);
		for (int i = 0; i < M; i++) {
			double temp2 = 0;
			for (int j = 0; j < K; j++) {
				temp2 += X[j] * A[j][i];
			}
			temp1[i] = N[i] / (temp2 + epsilon) / (temp2 + epsilon);
		}
		for (int k = 0; k < K; k++) {
			for (int l = 0; l < K; l++) {
				hessian[k][l] = 0;
				for (int i = 0; i < M; i++) {
					hessian[k][l] -= temp1[i] * A[k][i] * A[l][i];
				}
				if (k == l) hessian[k][l] -= 1e-6;
			}
		}
		return hessian;
	}

	vector<vector<double> > obs_fisher(const vector<vector<double> > &A, const vector<double> &N, const vector<double> &X) {
		vector<vector<double> > temp = hessian_log_likelihood(A, N, X);
		for (int i = 0; i < K; i++)
			for (int j = 0; j < K; j++)
				temp[i][j] = -temp[i][j];
		return temp;
	}

/*	vector<vector<double> > inv_obs_fisher(const vector<double> &X) {
		return inv_matrix(obs_fisher(X));
	}*/
	
	vector<vector<double> > fisher(const vector<vector<double> > &A, const vector<double> &X) {
		vector<vector<double> > fisher;
		fisher.resize(K);
		for (int i = 0; i < K; i++) fisher[i].resize(K);
		vector<double> temp1;
		temp1.resize(M);
		for (int i = 0; i < M; i++) {
			double temp2 = 0;
			for (int j = 0; j < K; j++) {
				temp2 += X[j] * A[j][i];
			}
			temp1[i] = 1 / (temp2 + epsilon);
		}
		for (int k = 0; k < K; k++) {
			for (int l = 0; l < K; l++) {
				fisher[k][l] = 0;
				for (int i = 0; i < M; i++) {
					fisher[k][l] += temp1[i] * A[k][i] * A[l][i];
				}
				if (k == l) fisher[k][l] += 1e-6;
			}
		}
		return fisher;
	}

/*	vector<vector<double> > inv_fisher(const vector<double> &X) {
		return inv_matrix(fisher(X));
	}*/

	double fmax_coord(const vector<vector<double> > &A, const vector<double> &N, vector<double> &X, int k, vector<double> *W = NULL) {
		double f_value = log_likelihood(A, N, X, W);
		double grad = grad_log_likelihood(A, N, X, k, W);
		if (grad > 0) grad = 1; else if (grad < 0) grad = -1; else return f_value;
		double x_old = X[k];
		double step = 1e6;
		int iteration = 0;
		while (true) {
			X[k] = x_old + step * grad;
			if (X[k] < 0) X[k] = 0;
			double new_f_value = log_likelihood(A, N, X, W);
			if (new_f_value > f_value) return new_f_value;
			if (step > epsilon) step /= 2; else break;
			if (iteration < MAX_ITR) iteration++; else break;
		}
		//if (verbose && step <= epsilon) cout << "warning: step <= epsilon\n";
		if (iteration >= MAX_ITR) cout << "warning: iteration >= max_iteration\n";
		X[k] = x_old;
		return f_value;
	}

	double fmax_coord(const vector<vector<double> > &A, const vector<double> &N, vector<double> &X, vector<double> *W = NULL, bool use_EM = do_EM, double lam = lambda) {
		double f_value = log_likelihood(A, N, X, W);
		int iteration = 0;
		vector<double> bias(M, 1.0);
		vector<double> W_bias(K, 1.0);
		double max_N = *max_element(N.begin(), N.end());
		if (lam < 0) lam = sqrt(max_N);
		while (true) {
			if (verbose) {
				cout << "iteration: " << iteration << " value: " << f_value << " X: ";
				for (int i = 0; i < K; i++) cout << X[i] << ",";
				cout << endl;
			}
			vector<vector<double> > A1;
			A1.resize(K);
			for (int i = 0; i < K; i++) {
				A1[i].resize(M);
				for (int j = 0; j < M; j++) {
					A1[i][j] = A[i][j] * bias[j];
				}
			}

			vector<double> *W1 = NULL;
			if (W != NULL) {
				W1 = new vector<double>;
				W1->resize(K);
				for (int i = 0; i < K; i++) {
					(*W1)[i] = (*W)[i];
					for (int j = 0; j < M; j++) {
						(*W1)[i] -= A[i][j];
					}
					(*W1)[i] *= W_bias[i];
					for (int j = 0; j < M; j++) {
						(*W1)[i] += A1[i][j];
					}
				}
			}

			if (use_EM) {
				vector<double> temp1(M, 0);
				for (int j = 0; j < M; j++) {
					for (int i = 0; i < K; i++) {
						temp1[j] += A1[i][j] * X[i];
					}
					if (temp1[j] < epsilon) temp1[j] = epsilon;
				}
				vector<double> X_old = X;
				for (int i = 0; i < K; i++) {
					double temp2 = 0, temp3 = 0;
					for (int j = 0; j < M; j++) {
						temp2 += A1[i][j] * N[j] / temp1[j];
						temp3 += A1[i][j];
					}
					if (W1 != NULL) temp3 = (*W1)[i];
					if (temp3 < epsilon) temp3 = epsilon;
					X[i] = X_old[i] * temp2 / temp3;
					if (X[i] < 1e-4) X[i] = 0;
				}
			} else { //coordinate wise descent
				for (int i = 0; i < K; i++) {
					fmax_coord(A1, N, X, i, W1);
				}
			}

			if (W1 != NULL) delete W1;

			if (lam > 0) {
				//solve for bias term
				for (int j = 0; j < M; j++) {
					double temp = 0;
					for (int i = 0; i < K; i++) {
						temp += A[i][j] * X[i];
					}
/*
					double b1 = (N[j] - lam) / temp;
					double b2 = (N[j] + lam) / temp;
					if (b1 < 1) b1 = 1;
					if (b2 > 1) b2 = 1;
					double f1 = N[j] * log(b1) - b1 * temp - lam * fabs(b1);
					double f2 = N[j] * log(b2) - b2 * temp - lam * fabs(b2);

					if (f1 > f2) bias[j] = b1; else bias[j] = b2;
*/
					if (fabs(temp) < epsilon) bias[j] = 1;
					else if (temp <= N[j] - lam) bias[j] = (N[j] - lam) / temp;
					else if (temp >= N[j] + lam) bias[j] = (N[j] + lam) / temp;
					else bias[j] = 1.0;
				}
				if (W != NULL) {
					for (int i = 0; i < K; i++) {
						double temp = (*W)[i];
						for (int j = 0; j < M; j++) {
							temp -= A[i][j];
						}
						temp *= X[i];
						if (temp > lam) W_bias[i] = lam / temp;
						else W_bias[i] = 1.0;
					}
				}
				vector<double> b;
				for (int i = 0; i < M; i++) b.push_back(log(bias[i]));
				if (W != NULL) for (int i = 0; i < K; i++) b.push_back(log(W_bias[i]));
				double med = median(b);
				for (int i = 0; i < M; i++) bias[i] /= exp(med);
				if (W != NULL) for (int i = 0; i < K; i++) W_bias[i] /= exp(med);
			}
			double new_f_value = log_likelihood(A, N, X, W);
			if (fabs(f_value - new_f_value) < epsilon) {
				if (verbose) cout << "number of iteration: " << iteration << endl;
				return new_f_value;
			}
			f_value = new_f_value;
			if (iteration < MAX_ITR) iteration++; else break;
		}
		if (iteration >= MAX_ITR) cout << "warning: iteration >= max_iteration\n";
		return f_value;
	}
};

inline void test_solve_likelihood() {
	isoform_opt iso_opt;
	int M = 4;
	int K = 2;
	iso_opt.M = M;
	iso_opt.K = K;
	vector<double> N;
	N.resize(M);
	N[0] = 300;
	N[1] = 5000;
	N[2] = 1000;
	N[3] = 100;
	vector<vector<double> > A;
	A.resize(K);
	A[0].resize(M);
	A[0][0] = 1;
	A[0][1] = 1;
	A[0][2] = 0;
	A[0][3] = 0;
	A[1].resize(M);
	A[1][0] = 0;
	A[1][1] = 1;
	A[1][2] = 1;
	A[1][3] = 1;
	iso_opt.verbose = true;
	vector<double> X;
	X.resize(2);
	X[0] = 1;
	X[1] = 1;

	iso_opt.fmax_coord(A, N, X, do_EM);
	vector<double>* W = NULL;
/*	vector<double>* W = new vector<double>;
	W->resize(K);
	(*W)[0] = 3;
	(*W)[1] = 4;*/
	do_EM = true;
	lambda = -1;
	iso_opt.fmax_coord(A, N, X, W, do_EM, lambda);
	cout << "X: " << X[0] << "\t" << X[1] << endl;
	vector<vector<double> > H = iso_opt.hessian_log_likelihood(A, N, X);
	cout << "H: ";
	cout << "\t" << H[0][0] << "\t" << H[0][1] << endl;
	cout << "\t" << H[1][0] << "\t" << H[1][1] << endl;	
	vector<vector<double> > OF = iso_opt.obs_fisher(A, N, X);
	cout << "OF: ";
	cout << "\t" << OF[0][0] << "\t" << OF[0][1] << endl;
	cout << "\t" << OF[1][0] << "\t" << OF[1][1] << endl;	
	vector<vector<double> > F = iso_opt.fisher(A, X);
	cout << "F: ";
	cout << "\t" << F[0][0] << "\t" << F[0][1] << endl;
	cout << "\t" << F[1][0] << "\t" << F[1][1] << endl;	
/*		vector<vector<double> > IOF = inv_obs_fisher(X);
	cout << "IOF: ";
	cout << "\t" << IOF[0][0] << "\t" << IOF[0][1] << endl;
	cout << "\t" << IOF[1][0] << "\t" << IOF[1][1] << endl;	
	vector<vector<double> > IF = inv_fisher(X);
	cout << "IF: ";
	cout << "\t" << IF[0][0] << "\t" << IF[0][1] << endl;
	cout << "\t" << IF[1][0] << "\t" << IF[1][1] << endl;	*/
}

inline void solve_likelihood(const vector<vector<double> > &rates, const vector<double> &counts, vector<double> &exp, vector<double> *W = NULL, bool use_EM = do_EM, double lam = lambda) {
	isoform_opt iso_opt;
	iso_opt.M = (int)counts.size();
	iso_opt.K = (int)rates.size();
	exp.resize(iso_opt.K);
	for (int i = 0; i < iso_opt.K; i++) {
		exp[i] = 1;
	}
	iso_opt.fmax_coord(rates, counts, exp, W, use_EM, lam);
}

inline void solve_likelihood_t(const vector<vector<double> > &rates, const vector<double> &counts, vector<double> &exp, vector<double> *W = NULL, bool use_EM = do_EM, double lam = lambda) {
	isoform_opt iso_opt;
	iso_opt.M = (int)counts.size();
	__ASSERT(iso_opt.M > 0, "internal error: empty rates.\n");
	iso_opt.K = (int)rates[0].size();
	vector<vector<double> > A;
	A.resize(iso_opt.K);
	for (int i = 0; i < iso_opt.K; i++) A[i].resize(iso_opt.M);
	for (int i = 0; i < iso_opt.M; i++) {
		for (int j = 0; j < iso_opt.K; j++) {
			A[j][i] = rates[i][j];
		}
	}
	exp.resize(iso_opt.K);
	for (int i = 0; i < iso_opt.K; i++) {
		exp[i] = 1;
	}
	iso_opt.fmax_coord(A, counts, exp, W, use_EM, lam);
}

#endif // ISOFORM_OPT_H
