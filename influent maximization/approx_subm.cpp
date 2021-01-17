// Copyright 2020, Grigorios Loukides
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Takuya Akiba nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <cfloat> // DBL_MAX
#include <cmath> // std::nextafter
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <time.h> // clock
#include <math.h>
#include <random>
#include <memory> // auto_ptr
#include "jlog.h"
#include "tools.h"

using namespace std;
using namespace jlog_internal;
typedef long long LL;
typedef double(*Get_xi)(const vector<double> & ,const set<pair<int, int> >&, const pair<int, int> *);
double epsilon;
default_random_engine* generator;
uniform_real_distribution<double>* distribution;
typedef unsigned long long ULL;
map<pair<set<pair<int, int> >, pair<int, int> >, double> records;

double entropy(vector<vector<vector<int> > > &vals, set<pair<int, int> > &S,
		int n_sensors, int n_ticks, int n_pos) {
	map<vector<int>, int> C;
	for (int t = 0; t < n_ticks; t++) {
		vector<int> x;
		for (auto item:S) {
			int p = item.first, s = item.second;
			x.push_back(vals[s][t][p]);
		}
		C[x]++;
	}
	double H = 0;
	for (auto it = C.begin(); it != C.end(); it++) {
		double pr = 1.0 * it->second / n_ticks;
		H += -pr * log(pr);
	}
	return H;
}

// [s][t][p] = val
// <pos,sensor>
double mutual_information(vector<vector<vector<int> > > &vals,
		set<pair<int, int> > &S, int n_sensors, int n_ticks, int n_pos) {
	// <p,s>
	// U = S \cup O
	set<pair<int, int> > O, U;
	set<pair<int, int> > Sset(S.begin(), S.end());
	for (int p = 0; p < n_pos; p++) {
		for (int s = 0; s < n_sensors; s++) {
			pair<int, int> ps = make_pair(p, s);
			if (!Sset.count(ps)) {
				O.insert(ps);
			}
			U.insert(ps);
		}
	}

	return entropy(vals, S, n_sensors, n_ticks, n_pos);
}

void generate_xi(vector<double> &xi, int n_sensors, default_random_engine &generator, uniform_real_distribution<double> &distribution){
	xi.resize(n_sensors);
	for(int i=0;i<n_sensors;++i)
		xi[i]=distribution(generator);
}
// < <item, node> >
double max_xi(const vector<double> &xi, const set<pair<int, int> > &S, const pair<int, int> *add){
	if(S.size()==0 && add==0)return 1.0;
	double ret = 0.0;
	for(auto item:S)
		ret = max(ret, xi[item.first]);
	if(add!=0) ret = max(ret, xi[add->first]);
	return ret;
}

double sum_xi(const vector<double> &xi, const set<pair<int, int> > &S, const pair<int, int> *add){
	if(S.size()==0&&add==0)return 1.0;
	if(S.size()==0)return xi[add->first];
	double ret = 0.0;
	for(auto item:S)
		ret += xi[item.first];
	if(add!=0) return (ret+xi[add->first])/(S.size()+1.0);
	return ret/S.size();
}

double new_xi(const vector<double> &xi, const set<pair<int, int> > &S, const pair<int, int> *add){
	auto tmp = make_pair(S, *add);
	auto rec = records.find(tmp);
	if(rec!=records.end())return rec->second;
	double ret = (*distribution)(*generator);
	records.insert(make_pair(tmp, ret));
	return ret;
}


void lazy_greedy(vector<vector<vector<int> > > &vals, int k, int tick, int n,
                vector<int> budgets, vector<double> &xi, Get_xi get_xi, bool TS = false, bool f_known = false, bool adr=false) {

	priority_queue<pair<double, pair<pair<int, int>, int> > > que;
	for (int v = 0; v < n; v++) {
		for (int z = 0; z < k; z++) {
			que.push(make_pair(1e12, make_pair(make_pair(v, z), -1)));
		}
	}

	int B = 0;
	for (int j = 0; j < budgets.size(); j++) {
		B += budgets[j];
	}

	vector<bool> used(n);
	double val = 0;
	set<pair<int, int> > S;
	int sum_num = 0;
	for (int j = 0; j < B; j++) {
		double gain = 0;
		pair<int, int> next;
		for (int num = 0;; num++) {
			pair<double, pair<pair<int, int>, int> > pp = que.top();
			que.pop();
			pair<int, int> s = pp.second.first;
			int last = pp.second.second;
			if (used[s.first] || (!TS && budgets[s.second] == 0)) 
				continue;
			if (last == j) {
				next = s;
				gain = pp.first;
				sum_num += num;
				if (adr) {
					if(get_xi == new_xi) val+=(1.0+epsilon)*gain;
					else val += (f_known?get_xi(xi, S, &s):1.0)*gain;
				} else val += gain;
				break;
			}
			set<pair<int, int> > T(S);
			T.insert(s);
			double delta;
			if(adr) {
				double H1 = mutual_information(vals, T, k, tick, n);
				double H2 = mutual_information(vals, S, k, tick, n);
				delta = (f_known?1.0:get_xi(xi, S, &s)) * (H1-H2);
			} else {
				double H1 = (f_known?1.0:get_xi(xi, S, &s)) * mutual_information(vals, T, k, tick, n);
				double H2 = (f_known?1.0:get_xi(xi, S, &s)) * mutual_information(vals, S, k, tick, n);
				delta =(H1-H2);
			}
			que.push(make_pair(delta, make_pair(s, j)));
		}
		S.insert(next);
		used[next.first] = true;
		if(!TS) budgets[next.second]--;
		if (adr)
			cout<<"B = "<<S.size()<<" | Influenced number: "<<val<<endl;
		else 
			cout<<"B = "<<S.size()<<" | Influenced number: "<<(f_known?get_xi==new_xi?(1+epsilon)*val:get_xi(xi, S, 0)*val:val)<<endl;
	}
}
void random(vector<vector<vector<int> > > &vals, int k, int tick, int n,
		vector<int> budgets, vector<double>& xi, Get_xi get_xi, bool TS = false, bool adr=false) {
	int B = 0;
	for (int i = 0; i < budgets.size(); i++) {
		B += budgets[i];
	}
	Xorshift xs(0);
	vector<bool> used(n);
	set<pair<int, int> > S;
	double val = 0, pre = 0;
	for (int j = 0; j < B; j++) {
		int e = xs.nextInt(n);
		int z = xs.nextInt(k);
		if (used[e] || (!TS && budgets[z] == 0)) {
			j--;
			continue;
		}
		used[e] = true;
		if(!TS)budgets[z]--;
		pair<int, int> s = make_pair(e, z);
		double tmp = get_xi(xi, S, &s);
		S.insert(s);
		if(adr) {
			double cur = mutual_information(vals, S, k, tick, n);
			val+=tmp*(cur-pre);
			pre = cur;
		} else 
			val = get_xi(xi, S, 0) * mutual_information(vals, S, k, tick, n);
		cout<<"B = "<<S.size()<<" | Influenced number: "<<val<<endl;
	}
}

int main(int argc, char *argv[]) {
	JLOG_INIT(&argc, argv);

	map<string, string> args;
	init_args(argc, argv, args);

	string algo = get_or_die(args, "algo");
	int n = atoi(get_or_die(args, "n").c_str()); // n_pos
	int k = atoi(get_or_die(args, "k").c_str()); // n_sensors
	int tick = atoi(get_or_die(args, "tick").c_str()); // n_ticks
	int budget = atoi(get_or_die(args, "budget").c_str());
	string Seed, Eps, method;
	Eps = get_or_die(args, "epsilon");
	Seed = get_or_die(args, "seed");
	double seed = atof(Seed.c_str());
	epsilon = atof(Eps.c_str());
	method  = get_or_die(args, "method");
	bool TS = atoi(get_or_die(args, "TS").c_str());
	JLOG_PUT("code", algo);
	Get_xi get_xi = method=="max"?max_xi:(method=="sum"?sum_xi:new_xi);
	// [s][t][p] = val
	vector<vector<vector<int> > > vals(k);
	for (int s = 0; s < k; s++) {
		char cs[256];
		sprintf(cs, "%d", s);
		ifstream is(get_or_die(args, cs));
		vals[s].resize(tick);
		for (int t = 0; t < tick; t++) {
			vals[s][t].resize(n);
			for (int p = 0; p < n; p++) {
				int v = 0;
				is >> v;
				vals[s][t][p] = v;
			}
		}
		is.close();
	}

	generator = new default_random_engine(seed);
	distribution = new uniform_real_distribution<double>(1.0-epsilon, nextafter(1.0, DBL_MAX));

	vector<double> xi;
	generate_xi(xi, n, *generator, *distribution);

	double user_start = get_current_time_sec();
	clock_t cpu_start = clock();
	cout<<"B = "<<budget<<", k = "<<k<<", eps = "<<epsilon<<", seed="<<seed<<endl;
	vector<int> budgets(k, budget);
	if(algo=="AS") {
		cout<<"lazy H: "<<endl;
		lazy_greedy(vals, k, tick, n, budgets, xi, get_xi, TS, true);
		cout<<"lazy tilde H: "<<endl;
		lazy_greedy(vals,k, tick, n, budgets, xi, get_xi, TS);
		cout<<"Random: "<<endl;
		random(vals, k, tick, n, budgets, xi, get_xi, TS);
	}else if(algo=="ADR") {
		cout<<"lazy H: "<<endl;
		lazy_greedy(vals, k, tick, n, budgets, xi, get_xi, TS, true, true);
		cout<<"lazy tilde H: "<<endl;
		lazy_greedy(vals,k, tick, n, budgets, xi, get_xi, TS, true);
		cout<<"Random: "<<endl;
		random(vals, k, tick, n, budgets, xi, get_xi, TS, true);
	}
	double user_end = get_current_time_sec();
	clock_t cpu_end = clock();
	JLOG_PUT("time.user", user_end - user_start);
	JLOG_PUT("time.cpu", (double) (cpu_end - cpu_start) / CLOCKS_PER_SEC );
}

