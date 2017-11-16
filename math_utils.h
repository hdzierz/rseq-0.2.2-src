#ifndef MATH_UTILS_H
///Define this macro to prevent from including this header file more than once.
#define MATH_UTILS_H

#include "stl.h"
#include "type.h"

const double epsilon = 1e-10;

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

inline vector<vector<double> > transpose(const vector<vector<double> > &X) {
	vector<vector<double> > results;
	int m = (int)X.size();
	if (m == 0) return results;
	int n = (int)X[0].size();
	results.resize(n);
	for (int i = 0; i < n; i++) {
		results[i].resize(m);
	}
	for (int i = 0; i < m; i++) {
		__ASSERT((int)X[i].size() == n, "internal error: wrong matrix dimensions.\n");
		for (int j = 0; j <n; j++) {
			results[j][i] = X[i][j];
		}
	}
	return results;
}

inline int sum(const vector<int> &data) {
	int result = 0;
	for (int i = 0; i < (int)data.size(); i++) {
		result += data[i];
	}
	return result;
}

inline double sum(const vector<double> &data) {
	double result = 0;
	for (int i = 0; i < (int)data.size(); i++) {
		result += data[i];
	}
	return result;
}

inline double mean(const vector<double> &data) {
	return sum(data) / (int)data.size();
}

inline double median(const vector<double> &data) {
	vector<double> temp = data;
	sort(temp.begin(), temp.end());
	int len = (int)temp.size();
	if (len % 2 == 0) return ((temp[len/2-1] + temp[len/2])/2);
	else return temp[len / 2];
}

inline double var(const vector<double> &data) {
	double sum1 = 0, sum2 = 0;
	for (int i = 0; i < (int)data.size(); i++) {
		sum1 += data[i];
		sum2 += data[i] * data[i];
	}
	return (sum2 - sum1 * sum1 / (int)data.size()) / (int)data.size();
}

inline double L1_norm(const vector<double> &v1) {
	double result = 0;
	for (int i = 0; i < (int)v1.size(); i++) {
		result += fabs(v1[i]);
	}
	return result;
}

inline double L1_dist(const vector<double> &v1, const vector<double> &v2) {
	__ASSERT(v1.size() == v2.size(), "fatal error: non-equal vector length.\n");
	double result = 0;
	for (int i = 0; i < (int)v1.size(); i++) {
		result += fabs(v1[i] - v2[i]);
	}
	return result;
}

inline double inner_prod(const vector<double> &v1, const vector<double> &v2) {
	__ASSERT(v1.size() == v2.size(), "fatal error: non-equal vector length.\n");
	double result = 0;
	for (int i = 0; i < (int)v1.size(); i++) {
		result += v1[i] * v2[i];
	}
	return result;
}

inline double L2_norm(const vector<double> &v1) {
	return sqrt(inner_prod(v1, v1));
}

inline double cos_angle(const vector<double> &v1, const vector<double> &v2) {
	return inner_prod(v1, v2) / L2_norm(v1) / L2_norm(v2);
}

inline bool vector_parallel(const vector<double> &v1, const vector<double> &v2, double e = epsilon) {
	return (fabs(cos_angle(v1, v2)) > 1 - e);
}

inline bool nearly_equal(const vector<double> &v1, const vector<double> &v2, double e = epsilon) {
	__ASSERT(v1.size() == v2.size(), "fatal error: non-equal vector length.\n");
	for (int i = 0; i < (int)v1.size(); i++) {
		if (fabs(v1[i] - v2[i]) > e) return false;
	}
	return true;
}

inline void L1_normalize(vector<double> &rate) {
	double rate_sum = L1_norm(rate);
	__ASSERT(fabs(rate_sum) > epsilon, "fatal error: sum too small.\n");
	for (int i = 0; i < (int)rate.size(); i++) {
		rate[i] /= rate_sum;
	}
}

inline double stddev(const vector<double> &data) {
	return sqrt(var(data));
}

inline int round_double(double x){
	if (x > 0) return int(x + 0.5);
	else return int(x-0.5);
}

inline bool is_int(double d) {
	return (fabs(d - round_double(d)) < epsilon);
}

inline double rand_double(){
	return ((double)(rand()))/RAND_MAX;
}

inline int rand_int(int low, int high){
	if (low == high) return low;
	__ASSERT(low <= high, "internal error: low > high");
//	return (rand() % (high - low)) + low;
	int temp = int((rand_double() * (high - low + 1))) + low;
	if (temp > high) temp = high;
	if (temp < low) temp = low;
	return temp;
}

inline int rand_int(int high) {
	return rand_int(1, high);
}

inline int rand_int(const vector<double> &p, bool sorted = false) {
	int n = (int)p.size();
	double total_p = (sorted ? p[n-1] : sum(p));
	double temp = rand_double() * total_p;
	int i;
	if (sorted) {
		i = (int)(lower_bound(p.begin(), p.end(), temp) - p.begin());
	} else {
		for (i = 0; i < n; i++) {
			if (temp <= p[i] + epsilon) break;
			temp -= p[i];
		}
	}
	__ASSERT(0 <= i && i < n, "internal error: rand_int.\n");
	return i;
}

inline vector<int> mnrnd(int n, const vector<double> &p) {
	int k = (int)p.size();
	__ASSERT(k > 0, "internal error: empty p.\n");
	vector<int> results(k, 0);
	vector<double> sort_p = p;
	for (int i = 1; i < (int)sort_p.size(); i++) {
		sort_p[i] = sort_p[i] + sort_p[i-1];
	}
	sort(sort_p.begin(), sort_p.end());
	for (int i = 0; i < n; i++) {
		results[rand_int(sort_p, true)]++;
	}
	return results;
}

class interval{
public:
	double start, end;
	interval(const interval &temp) {
		start = temp.start;
		end = temp.end;
	}
	interval(const double start, const double end){
		this->start = start;
		this->end = end;
	}
	inline double length() const {
		return end - start;
	}
	inline bool is_empty() const{
		return (length() < 0);
	}
	inline void outer_union_with(const interval &temp) { //not real union, which may generate a interval_set
		if (temp.is_empty()) return;
		if (is_empty()) {
			start = temp.start;
			end = temp.end;
		}
		start = min(start, temp.start);
		end = max(end, temp.end);
	}
	inline void intersect_with(const interval &temp){
		if (is_empty()) return;
		start = max(start, temp.start);
		end = min(end, temp.end);
	}
};

class interval_set{
public:
	vector<double> end_points;
	interval_set(){
		end_points.clear();
	}

	interval_set(const interval_set &temp) {
		end_points = temp.end_points;
	}

	interval_set(const interval &temp) {
		end_points.clear();
		end_points.push_back(temp.start);
		end_points.push_back(temp.end);
	}

	interval_set(const double start, const double end) {
		end_points.clear();
		end_points.push_back(start);
		end_points.push_back(end);
	}

	interval_set(const vector<pair<int, int> > pairs) {
		for (int i = 0; i < (int)pairs.size(); i++) {
			end_points.push_back(pairs[i].first);
			end_points.push_back(pairs[i].second);
		}
	}

	inline bool convert_to_int_pairs(vector<pair<int, int> > &pairs) {
		if (!check_int()) return false;
		pairs.clear(); 
		for (int i = 0; i < (int)end_points.size(); i += 2) {
			pair<int, int> temp_pair(round_double(end_points[i]), round_double(end_points[i+1]));
			pairs.push_back(temp_pair);
		}
		return true;
	}

	inline bool check_valid() {
		if ((int)end_points.size() %2 != 0) return false;
		for (int i = 1; i < (int)end_points.size(); i++) {
			if (end_points[i] + epsilon < end_points[i-1]) return false;
		}
		return true;
	}

	inline bool check_int() {
		if (!check_valid()) return false;
		for (int i = 0; i < (int)end_points.size(); i++) {
			if (!is_int(end_points[i])) return false;
		}
		return true;
	}

	inline double length() {
		double total = 0;
		for (int i = 0; i < (int)end_points.size(); i += 2) {
			total += end_points[i+1] - end_points[i];
		}
		return total;
	}

	inline double get_point(double coord) {
		double result = 0;
		__ASSERT(coord >= 0, "internal error: coord < 0");
		int i;
		for (i = 0; i < (int)end_points.size(); i += 2) {
			if (coord <= end_points[i+1] - end_points[i]) {
				result = end_points[i] + coord;
				break;
			} else {
				coord -= end_points[i+1] - end_points[i];
			}
		}
		__ASSERT(i < (int)end_points.size(), "internal error: i == (int)end_points.size()");
		return result;
	}

	inline void get_fragment(double begin, double end) {
		intersect_with(interval(get_point(begin), get_point(end)));
	}

	inline bool operator == (const interval_set &is) {
		interval_set is1, is2 = is;
		is1.end_points = end_points;
		is1.remove_empty_intervals();
		is1.remove_redundant_inner_points();
		is2.remove_empty_intervals();
		is2.remove_redundant_inner_points();
		return is1.end_points == is2.end_points;
	}

	inline bool operator <= (const interval_set &is) {
		interval_set is1, is2 = is;
		is1.end_points = end_points;
		is2.intersect_with(is1);
		return is2 == is1;
	}

	inline interval bracket() const{
		if (end_points.size() > 1) return interval(end_points[0], end_points[end_points.size() - 1]);
		else return interval(0, -1);
	}

	inline bool has_empty_intervals() {
		for (int i = (int)end_points.size() - 1; i > 0; i -= 2) {
			if (end_points[i] <= end_points[i - 1] + epsilon) {
				return true;
			}
		}
		return false;
	}

	inline void remove_empty_intervals() {
		for (int i = (int)end_points.size() - 1; i > 0; i -= 2) {
			if (end_points[i] <= end_points[i - 1] + epsilon) {
				end_points.erase(end_points.begin() + i - 1, end_points.begin() + i + 1);
			}
		}
	}

	inline void remove_redundant_inner_points(){
		for (int i = (int)end_points.size() - 1; i > 0; i--) {
			if (end_points[i] <= end_points[i - 1] + epsilon) {
				end_points.erase(end_points.begin() + i - 1, end_points.begin() + i + 1);
				i--;
			}
		}
	}

	inline void complement_with(const interval &temp) {
		interval bracket1 = bracket();
		if (end_points.size() > 0 && (temp.start > bracket1.start || temp.end < bracket1.end)) return;
		end_points.insert(end_points.begin(), temp.start);
		end_points.push_back(temp.end);
		remove_redundant_inner_points();
	}

	inline void intersect_with(const interval_set &temp) {
		interval bracket1 = bracket();
		bracket1.outer_union_with(temp.bracket());
		if (bracket1.is_empty()) return;
		complement_with(bracket1);
		interval_set temp1 = temp;
		temp1.complement_with(bracket1);
		union_with(temp1);
		complement_with(bracket1);
	}

	inline void union_with(const interval_set &temp) {
		vector<pair<double, bool> > points;
		int i;
		for (i = 0; i < (int)end_points.size(); i++) {
			points.push_back(pair<double, bool>(end_points[i], i%2==0));
		}
		for (i = 0; i < (int)temp.end_points.size(); i++) {
			points.push_back(pair<double, bool>(temp.end_points[i], i%2==0));
		}
		sort(points.begin(), points.end());
		end_points.clear();
		int layer = 0;
		for (i = 0; i < (int)points.size(); i++) {
			if (points[i].second) {
				layer++;
				if (layer == 1) end_points.push_back(points[i].first);
			} else {
				layer--;
				if (layer == 0) end_points.push_back(points[i].first);
			}
		}
		remove_redundant_inner_points();
	}

	inline void break_union_with(const interval &temp) {
		int i = 0;
		for (i = 0; i < (int)end_points.size(); i++) {
			if (temp.start < end_points[i] - epsilon) break;
		}
		bool overlap = (i > 0 && fabs(temp.start - end_points[i-1]) < epsilon);
		bool inside = (i % 2 == 1);
		if (!inside) {
			end_points.insert(end_points.begin() + i, temp.start);
			i += 1;
		} else if (!overlap) {
			end_points.insert(end_points.begin() + i, 2, temp.start);
			i += 2;
		}
		while (true) {
			if (i == (int)end_points.size()) break;
			if (temp.end <= end_points[i] + epsilon) break;
			double temp_pos = end_points[i];
			if (i == (int)end_points.size() - 1 || fabs(end_points[i+1] - end_points[i]) > epsilon) end_points.insert(end_points.begin() + i, temp_pos);
			i += 2;
		}
		overlap = (i < (int)end_points.size() && fabs(temp.end - end_points[i]) < epsilon);
		inside = (((int)end_points.size() - i) % 2 == 1);
		if (!inside) {
			end_points.insert(end_points.begin() + i, temp.end);
			i += 1;
		} else if (!overlap) {
			end_points.insert(end_points.begin() + i, 2, temp.end);
			i += 2;
		}
		__ASSERT(!has_empty_intervals(), "internal error: has_empty_intervals().\n");
	}

	inline void break_union_with(const interval_set &temp) {
		for (int i = 0; i < (int)temp.end_points.size(); i += 2) {
			break_union_with(interval(temp.end_points[i], temp.end_points[i+1]));
		}
	}
};

class folding{
private:
	double equal_gap_length;
	double gap_coef, interval_coef;
public:
	vector<double> points;
	vector<double> mapped_points;
	interval_set intervals;
	bool begin_with_gap;
	double start, end;
	double mapped_start, mapped_end;
	bool equal_gap; //false: ratio_gap
	double gap_ratio;

	inline folding() {
	}
	
	inline folding(const interval_set &intervals, const double start, const double end, const double mapped_start, const double mapped_end, const bool equal_gap = true, const double gap_ratio = 0.1){
		configure(intervals, start, end, mapped_start, mapped_end, equal_gap, gap_ratio);
	}

	inline void configure(const interval_set &intervals, const double start, const double end, const double mapped_start, const double mapped_end, const bool equal_gap = true, const double gap_ratio = 0.1){
		if (start >= end - epsilon) panic("internal error: start > end.");

		this->intervals = intervals;
		this->start = start;
		this->end = end;
		this->mapped_start = mapped_start;
		this->mapped_end = mapped_end;
		this->equal_gap = equal_gap;
		this->gap_ratio = gap_ratio;
		
		this->intervals.intersect_with(interval_set(start, end));
		points = this->intervals.end_points;
		begin_with_gap = false;
		if (points.size() == 0 || start < points[0] - epsilon) {
			begin_with_gap = true;
			points.insert(points.begin(), start);
		}
		if (end > points[points.size() - 1] + epsilon) {
			points.push_back(end);
		}

		double total_interval = 0;
		double total_gap = 0;
		for (int i = 0; i < (int)points.size() - 1; i++) {
			if (begin_with_gap == (i%2==0)) total_gap += points[i+1] - points[i];
			else total_interval += points[i+1] - points[i];
		}

		double mapped_length = mapped_end - mapped_start;
		double mapped_interval_length = mapped_length * (1 - gap_ratio);
		double mapped_gap_length = mapped_length - mapped_interval_length;
		int num_gaps = (int)points.size()/2;
		if (!begin_with_gap && (points.size()%2==0)) num_gaps--;

		if (num_gaps == 0) {
			mapped_gap_length = 0;
			mapped_interval_length = mapped_length;
		}

		if (points.size() == 2 && num_gaps == 1) {
			mapped_gap_length = mapped_length;
			mapped_interval_length = 0;
		}

		if (num_gaps == 0) {
			equal_gap_length = 0;
			gap_coef = 0;
		} else {
			equal_gap_length = mapped_gap_length / num_gaps;
			gap_coef = mapped_gap_length / total_gap;
		}

		if (points.size() == 2 && num_gaps == 1) {
			interval_coef = 0;
		} else {
			interval_coef = mapped_interval_length / total_interval;
		}

		mapped_points.resize(points.size());
		for (int i = 0; i < (int)points.size(); i++) {
			if (i == 0) { 
				mapped_points[i] = mapped_start;
			} else if (begin_with_gap == (i%2==1)) {
				if (equal_gap) {
					mapped_points[i] = mapped_points[i-1] + equal_gap_length;
				} else {
					mapped_points[i] = mapped_points[i-1] + (points[i] - points[i-1]) * gap_coef;
				}
			} else {
				mapped_points[i] = mapped_points[i-1] + (points[i] - points[i-1]) * interval_coef;
			}
		}
		if (fabs(mapped_points[mapped_points.size() - 1] - mapped_end) > epsilon) panic("internal error, mapping length inconsistent.");
	}

	double map(double point) {
       int low = 0, high = (int)points.size();
       while (low < high) {
           int mid = (low + high)/2;
           if (points[mid] < point)
               low = mid + 1;
           else
               high = mid; 
       }
	   if (low == (int)points.size()) {
		   if (point < end) panic("internal error, bad mapping.");
		   return mapped_end + (point - end) * interval_coef;
	   } else if (low == 0) {
		   if (point > start) panic("internal error, bad mapping.");
		   return mapped_start - (start - point) * interval_coef;
	   } else {
		   return ((point - points[low - 1]) * mapped_points[low] + (points[low] - point) * mapped_points[low - 1]) / (points[low] - points[low - 1]);
	   }
	}

	double inverse_map(double mapped_point) {
       int low = 0, high = (int)mapped_points.size();
       while (low < high) {
           int mid = (low + high)/2;
           if (mapped_points[mid] < mapped_point)
               low = mid + 1;
           else
               high = mid; 
       }
	   if (low == (int)mapped_points.size()) {
		   if (mapped_point < mapped_end) panic("internal error, bad mapping.");
		   return end;
	   } else if (low == 0) {
		   if (mapped_point > mapped_start) panic("internal error, bad mapping.");
		   return start;
	   } else {
		   return ((mapped_point - mapped_points[low - 1]) * points[low] + (mapped_points[low] - mapped_point) * points[low - 1]) / (mapped_points[low] - mapped_points[low - 1]);
	   }
	}
};

inline double KolmogorovProb(double z)
{
	// Calculates the Kolmogorov distribution function,
	//Begin_Html
	/*
	<img src="gif/kolmogorov.gif">
	*/
	//End_Html
	// which gives the probability that Kolmogorov's test statistic will exceed
	// the value z assuming the null hypothesis. This gives a very powerful
	// test for comparing two one-dimensional distributions.
	// see, for example, Eadie et al, "statistocal Methods in Experimental
	// Physics', pp 269-270).
	//
	// This function returns the confidence level for the null hypothesis, where:
	//   z = dn*sqrt(n), and
	//   dn  is the maximum deviation between a hypothetical distribution
	//       function and an experimental distribution with
	//   n    events
	//
	// NOTE: To compare two experimental distributions with m and n events,
	//       use z = sqrt(m*n/(m+n))*dn
	//
	// Accuracy: The function is far too accurate for any imaginable application.
	//           Probabilities less than 10^-15 are returned as zero.
	//           However, remember that the formula is only valid for "large" n.
	// Theta function inversion formula is used for z <= 1
	//
	// This function was translated by Rene Brun from PROBKL in CERNLIB.

	double fj[4] = {-2,-8,-18,-32}, r[4];
	const double w = 2.50662827;
	// c1 - -pi**2/8, c2 = 9*c1, c3 = 25*c1
	const double c1 = -1.2337005501361697;
	const double c2 = -11.103304951225528;
	const double c3 = -30.842513753404244;

	double u = fabs(z);
	double p;
	if (u < 0.2) {
		p = 1;
	} else if (u < 0.755) {
		double v = 1./(u*u);
		p = 1 - w*(exp(c1*v) + exp(c2*v) + exp(c3*v))/u;
	} else if (u < 6.8116) {
		r[1] = 0;
		r[2] = 0;
		r[3] = 0;
		double v = u*u;
		int maxj = max(1, round_double(3./u));
		for (int j=0; j<maxj;j++) {
			r[j] = exp(fj[j]*v);
		}
		p = 2*(r[0] - r[1] +r[2] - r[3]);
	} else {
		p = 0;
	}
	return p;
}

inline double K_S_Test(vector<double> &a, vector<double> &b) {
	int na = (int)a.size(), nb = (int)b.size();

	double prob = -1;
	//      Require at least two points in each graph
	if (na <= 2 || nb <= 2) {
//		printf("Error: Sets must have more than 2 points");
		return prob;
	}
	//     Constants needed
	double rna = na;
	double rnb = nb;
	double sa  = 1./rna;
	double sb  = 1./rnb;
	double rdiff;
	int ia,ib;
	//     Starting values for main loop
	if (a[0] < b[0]) {
		rdiff = -sa;
		ia = 2;
		ib = 1;
	} else {
		rdiff = sb;
		ib = 2;
		ia = 1;
	}
	double rdmax = fabs(rdiff);

	//    Main loop over point sets to find max distance
	//    rdiff is the running difference, and rdmax the max.
	bool ok = false;
	for (int i=0;i<na+nb;i++) {
		if (a[ia-1] < b[ib-1]) {
			rdiff -= sa;
			ia++;
			if (ia > na) {ok = true; break;}
		} else if (a[ia-1] > b[ib-1]) {
			rdiff += sb;
			ib++;
			if (ib > nb) {ok = true; break;}
		} else {
			double x = a[ia-1];
			while(a[ia-1] == x && ia <= na) {
				rdiff -= sa;
				ia++;
			}
			while(b[ib-1] == x && ib <= nb) {
				rdiff += sb;
				ib++;
			}
			if (ia > na) {ok = true; break;}
			if (ib > nb) {ok = true; break;}
		}
		rdmax = max(rdmax, fabs(rdiff));
	}
	//    Should never terminate this loop with ok = kFALSE!

	if (ok) {
		rdmax = max(rdmax, fabs(rdiff));
		double z = rdmax * sqrt(rna*rnb/(rna+rnb));
		prob = KolmogorovProb(z);
	}

	return prob;
}

inline int rand_poisson(double lambda) {
	double L = -lambda;
	int k = 0;
	double p = 1;
	do{
		k++;
		double u = rand_double();
		if (u == 0) break; else p += log(u);
	} while (p>=L);
	return k - 1;
}

inline vector<double> poisson_break(double lambda, double len) {
	int k = rand_poisson(lambda*len);
	vector<double> breaks;
	for (int i = 0; i < k; i++) {
		breaks.push_back(rand_double() * len);
	}
	sort(breaks.begin(), breaks.end());
	return breaks;
}


inline bool in_range(int x, int a, int b = -1) {
	if ((b == -1 && x == a) || (x >= a && x <= b)) {
		return true;
	}
	return false;
}

template<typename T> 
class interval_list {
public:
	vector<pair<pair<double, double>, pair<int, T> > > intervals;
	vector<double> starts;
	
	void add_interval(const double start, const double end, const T& t) {
		intervals.push_back(pair<pair<double, double>, pair<int, T> >(pair<double, double>(start, end), pair<int, T>(0, t)));
	}
	
	void prepare() {
		sort(intervals.begin(), intervals.end());
		for (int i = 0; i <(int)intervals.size(); i++) {
			starts.push_back(intervals[i].first.first);
		}
		for (int i = (int)intervals.size() - 1; i >= 0; i--) {
			int j = i;
			while (j < (int)intervals.size() && intervals[i].first.second >= intervals[j].first.first) {
				intervals[j].second.first = j - i;
				j++;
			}
		}
	}
	
	void search_interval(const double start, const double end, vector<T> &t) {
		int i = (int)(lower_bound(starts.begin(), starts.end(), start) - starts.begin());
		if (i > 0) i--;
		i -= intervals[i].second.first;
		while (i < (int)intervals.size() && intervals[i].first.first <= end) {
			if (intervals[i].first.second >= start) t.push_back(intervals[i].second.second);
			i++;
		}
	}
};

inline double get_quantile(const vector<double> &data, double quantile) {
	int index = (int)((data.size() - 1) * quantile);
	return data[index];
}

inline double NormalProb(double x) {

//****************************************************************************80
//
//  Purpose:
//
//    ALNORM computes the cumulative density of the standard normal distribution.
//
//  Modified:
//
//    17 January 2008
//
//  Author:
//
//    David Hill
//    C++ version by John Burkardt
//
//  Reference:
//
//    David Hill,
//    Algorithm AS 66:
//    The Normal Integral,
//    Applied Statistics,
//    Volume 22, Number 3, 1973, pages 424-427.
//
//  Parameters:
//
//    Input, double X, is one endpoint of the semi-infinite interval
//    over which the integration takes place.
//
//    Input, bool UPPER, determines whether the upper or lower
//    interval is to be integrated:
//    .TRUE.  => integrate from X to + Infinity;
//    .FALSE. => integrate from - Infinity to X.
//
//    Output, double ALNORM, the integral of the standard normal
//    distribution over the desired interval.
//
  double a1 = 5.75885480458;
  double a2 = 2.62433121679;
  double a3 = 5.92885724438;
  double b1 = -29.8213557807;
  double b2 = 48.6959930692;
  double c1 = -0.000000038052;
  double c2 = 0.000398064794;
  double c3 = -0.151679116635;
  double c4 = 4.8385912808;
  double c5 = 0.742380924027;
  double c6 = 3.99019417011;
  double con = 1.28;
  double d1 = 1.00000615302;
  double d2 = 1.98615381364;
  double d3 = 5.29330324926;
  double d4 = -15.1508972451;
  double d5 = 30.789933034;
  double ltone = 7.0;
  double p = 0.398942280444;
  double q = 0.39990348504;
  double r = 0.398942280385;
  bool up;
  double utzero = 18.66;
  double value;
  double y;
  double z;

  up = true;
  z = x;

  if ( z < 0.0 )
  {
    up = !up;
    z = - z;
  }

  if ( ltone < z && ( ( !up ) || utzero < z ) )
  {
    if ( up )
    {
      value = 0.0;
    }
    else
    {
      value = 1.0;
    }
    return value;
  }

  y = 0.5 * z * z;

  if ( z <= con )
  {
    value = 0.5 - z * ( p - q * y 
      / ( y + a1 + b1 
      / ( y + a2 + b2 
      / ( y + a3 ))));
  }
  else
  {
    value = r * exp ( - y ) 
      / ( z + c1 + d1 
      / ( z + c2 + d2 
      / ( z + c3 + d3 
      / ( z + c4 + d4 
      / ( z + c5 + d5 
      / ( z + c6 ))))));
  }

  if ( !up )
  {
    value = 1.0 - value;
  }

  return value;
}

// lgamma.cpp -- log gamma function of real argument.
//      Algorithms and coefficient values from "Computation of Special
//      Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
//  Returns log(gamma) of real argument.
//  NOTE: Returns 1e308 if argument is 0 or negative.
//#include <math.h>
inline double lgamma(double x)
{
    double x0,x2,xp,gl,gl0;
    int n=0,k;
    static double a[] = {
        8.333333333333333e-02,
       -2.777777777777778e-03,
        7.936507936507937e-04,
       -5.952380952380952e-04,
        8.417508417508418e-04,
       -1.917526917526918e-03,
        6.410256410256410e-03,
       -2.955065359477124e-02,
        1.796443723688307e-01,
       -1.39243221690590};
    
    x0 = x;
    if (x <= 0.0) return 1e308;
    else if ((x == 1.0) || (x == 2.0)) return 0.0;
    else if (x <= 7.0) {
        n = (int)(7-x);
        x0 = x+n;
    }
    x2 = 1.0/(x0*x0);
    xp = 2.0*M_PI;
    gl0 = a[9];
    for (k=8;k>=0;k--) {
        gl0 = gl0*x2 + a[k];
    }
    gl = gl0/x0+0.5*log(xp)+(x0-0.5)*log(x0)-x0;
    if (x <= 7.0) {
        for (k=1;k<=n;k++) {
            gl -= log(x0-1.0);
            x0 -= 1.0;
        }
    }
    return gl;
}

static std::mt19937 engine; // Mersenne twister MT19937

/*
inline double rand_gamma(double mean = 1, double std = 1) {
		std::gamma_distribution<double> gamma(mean*mean/std/std, mean/std/std);
		return gamma(engine);
}
*/

inline double rand_normal(double mean = 0, double std = 1) {
		std::normal_distribution<double> normal(mean, std);
		return normal(engine);
}

inline double rand_exponential(double mean = 1) {
		std::exponential_distribution<double> exponential(1/mean);
		return exponential(engine);
}

#endif // MATH_UTILS_H
