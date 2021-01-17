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

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>	
#include <cfloat>
#include <random>
#include <time.h> // clock
#include <math.h>
#include <memory> // auto_ptr
#include "jlog.h"
#include "tools.h"
#include "mt19937ar.c"

using namespace std;
using namespace jlog_internal;
typedef long long LL;
typedef unsigned long long ULL;
typedef double(*Get_xi)(const vector<double> &, const set<pair<int, int> >&, const pair<int, int> *add);
vector<bool> X;
vector<bool> active;
double epsilon;
default_random_engine* generator;
uniform_real_distribution<double>* distribution;
map<set<pair<int, int> >, double> records;
map<pair<set<pair<int, int> >, pair<int, int> >, double> xi_records;
// [z][u][i] = <v,p>
double simulate(int n, vector<vector<vector<pair<int, double> > > > &es, int k,
		set<pair<int, int> > &S, int R) {
	Xorshift xs(0);
	LL sum = 0;
	for (int t = 0; t < R; t++) {
		for (int z = 0; z < k; z++) {
			vector<int> tmp;
			queue<int> Q;
			for (auto si: S) {
				if (si.second == z) {
					Q.push(si.first);
					X[si.first] = true;
				}
			}
			for (; !Q.empty();) {
				int u = Q.front();
				Q.pop();
				active[u] = true;
				tmp.push_back(u);
				for (int i = 0; i < es[z][u].size(); i++) {
					int v = es[z][u][i].first;
					double p = es[z][u][i].second;
					if (!X[v] && xs.nextDouble() < p) {
						X[v] = true;
						Q.push(v);
					}
				}
			}
			for (int i = 0; i < tmp.size(); i++) {
				X[tmp[i]] = false;
			}
		}
		int n1 = 0;
		for (int v = 0; v < n; v++) {
			if (active[v]) {
				n1++;
				active[v] = false;
			}
		}
		sum += n1;
	}
	return 1.0 * sum / R;
}

void generate_xi(vector<double> &xi, int n, default_random_engine &generator, uniform_real_distribution<double> &distribution){
	xi.resize(n);
	for(int i=0;i<n;++i)
		xi[i]=distribution(generator);
}

inline double max_xi(const vector<double> &xi, const set<pair<int, int> >& S, const pair<int, int> *add = NULL){
	if(S.size()==0&&add==0)return 1.0;
	if(S.size()==0)return xi[add->first];
	double ret = 0.0;
	for(auto item:S)
		ret = max(ret, xi[item.first]);
	if(add) return max(ret,xi[add->first]);
	return ret;
}

inline double average_xi(const vector<double> &xi, const set<pair<int, int> >& S, const pair<int, int> *add = NULL) {
	if(S.size()==0&&add==0)return 1;
	if(S.size()==0)return xi[add->first]; 
	double ret = 0.0;
	for(auto item:S)
		ret += xi[item.first];
	if(add) return(ret+xi[add->first])/(S.size()+1.0);
	return ret/S.size();
}

inline double new_xi(const vector<double> &xi, const set<pair<int, int> >& S, const pair<int, int> *add){
	auto tmp = make_pair(S, *add);
	auto rec = xi_records.find(tmp);
	if(rec!=xi_records.end())return rec->second;
	double ret = (*distribution)(*generator);
	xi_records.insert(make_pair(tmp, ret));
	return ret;
}


inline double get_simulate(int n, vector<vector<vector<pair<int, double> > > > &es, 
			int k, set<pair<int, int> > &S, int R){
	auto t = records.find(S);
	if(t!=records.end())return t->second;
	double ret = simulate(n, es, k, S, R);
	records.insert(make_pair(S, ret));
	return ret;
}
// [u][i] = <v, p[1..k]>
void greedy(int n, vector<vector<vector<pair<int, double> > > > &es, int k,
		int budget, int beta, vector<double>& xi, Get_xi get_xi, bool f_known = false, bool adr=false) {
	double val = 0;
	// < gain, < <node,item>, tick > >
	priority_queue<pair<double, pair<pair<int, int>, int> > > que;
	for (int v = 0; v < n; v++) {
		for (int z = 0; z < k; z++) {
			que.push(make_pair(1e12, make_pair(make_pair(v, z), -1)));
		}
	}

	vector<bool> used(n);
	set<pair<int, int> > S;
	for (int j = 0; j < budget; j++) {
		//printf(INFO "|S| = %d" DEF, j);
		double gain;
		pair<int, int> next;
		double start = get_current_time_sec();
		for (int num = 0;; num++) {
			double now = get_current_time_sec();

			pair<double, pair<pair<int, int>, int> > pp = que.top();
			que.pop();
			pair<int, int> s = pp.second.first;
			int last = pp.second.second;
			if (used[s.first]) {
				continue;
			}
			if (last == j) {
				next = s;
				gain = pp.first;
				if (adr) {
					if(get_xi == new_xi) val+=(1+epsilon)*gain;
					else val += (f_known?get_xi(xi, S, &s):1)*gain;
				} else val += gain;
				break;
			}
			set<pair<int, int> > SS(S);
			SS.insert(s);
			double sigma, psigma, delta;
			if (adr) {
				sigma = get_simulate(n, es, k, SS, beta);
				psigma = get_simulate(n, es, k, S, beta);
				delta = (f_known?1:get_xi(xi, S, &s)) * (sigma-psigma);
			} else {
				sigma = (f_known?1:get_xi(xi, S, &s)) * get_simulate(n, es, k, SS, beta);
				psigma = (f_known?1:get_xi(xi, S, &s)) * get_simulate(n, es, k, S, beta);
				delta = (sigma-psigma);
			}
			que.push(make_pair(delta, make_pair(s, j)));
		}
		S.insert(next);
		used[next.first] = true;
		if(S.size()%5==0){
			if (adr)
				cout<<"B = "<<S.size()<<" | Influenced number: " << val <<endl;
			else 
				cout<<"B = "<<S.size()<<" | Influenced number: " << (f_known?get_xi==new_xi?(1+epsilon)*val:get_xi(xi, S, 0)*val:val) <<endl;
		}
		// JLOG_ADD("seed", next.first);
		// JLOG_ADD("item", next.second);
		// JLOG_ADD("gain", gain);
	}

}

void Random(int n, vector<vector<vector<pair<int, double> > > > &es, int k,
		int budget, int beta, vector<double>& xi, Get_xi get_xi, bool adr = false) {
	vector<bool> used(n);
	set<pair<int, int> > S;
	Xorshift xs(0);
	double val = 0, pre = 0;
	for (int i = 0; i < budget; i++) {
		int v = xs.nextInt(n);
		if(used[v]){
			i--;
			continue;
		}
		int z = xs.nextInt(k);
		pair<int, int> s = make_pair(v, z);
		double tmp = get_xi(xi, S, &s);
		S.insert(s);
		if (adr) {
			double cur = get_simulate(n, es, k, S, beta); 
			val += tmp * (cur - pre);
			pre = cur;
		} else {
			double cur = get_xi(xi, S, 0) * get_simulate(n, es, k, S, beta);
			val = cur;
		}
		used[v]=true;
		if(S.size()%5==0){
			cout<<"B = "<<S.size()<<" | Influenced number: " << val <<endl;
		}
		// JLOG_ADD("seed", v);
		// JLOG_ADD("item", z);
	}

}

void degree(int n, vector<vector<vector<pair<int, double> > > > &es, int k,
		int budget, int beta, vector<double>& xi, Get_xi get_xi, bool adr = false) {
	vector<pair<int, int> > odeg;
	Xorshift xs(0);
	set<pair<int, int> > S;
	double val = 0, pre = 0;
	for (int v = 0; v < n; v++) {
		int o = es[v].size();
		odeg.push_back(make_pair(-o, v));
	}
	sort(odeg.begin(), odeg.end());
	for (int i = 0; i < budget; i++) {
		int v = odeg[i].second;
		int z = xs.nextInt(k);
		pair<int, int> s = make_pair(v, z);
		double tmp = get_xi(xi, S, &s);
		S.insert(s);
		if (adr) {
			double cur = get_simulate(n, es, k, S, beta); 
			val += tmp * (cur - pre);
			pre = cur;
		} else val = get_xi(xi, S, 0) * get_simulate(n, es, k, S, beta);
		// JLOG_ADD("seed", v);
		// JLOG_ADD("item", z);
		if(S.size()%5==0){
			cout<<"B = "<<S.size()<<" | Influenced number: " << val <<endl;
		}
	}
}

int main(int argc, char *argv[]) {
	JLOG_INIT(&argc, argv);		

	map<string, string> args;
	init_args(argc, argv, args);

	string algo = get_or_die(args, "algo");
	string edge = get_or_die(args, "edge");
	string prob = get_or_die(args, "prob");
	int k = atoi(get_or_die(args, "k").c_str());
	int budget = atoi(get_or_die(args, "budget").c_str());
	string Seeds, Eps;
	Eps = get_or_die(args, "epsilon");
	Seeds = get_or_die(args, "seed");
	double seed = atof(Seeds.c_str());
	epsilon = atof(Eps.c_str());
	int beta = atoi(get_or_die(args, "beta").c_str());
	string method = get_or_die(args, "method");
	Get_xi get_xi = method=="max"?max_xi:(method=="sum"?average_xi:new_xi);
	JLOG_PUT("code", algo);
	vector<vector<pair<int, vector<double> > > > es0 = load_graph(edge, prob,
			k);
	int n = es0.size();
	vector<vector<vector<pair<int, double> > > > es(k);
	for (int z = 0; z < k; z++) {
		es[z].resize(n);
		for (int u = 0; u < n; u++) {
			for (int i = 0; i < es0[u].size(); i++) {
				int v = es0[u][i].first;
				double p = es0[u][i].second[z];
				es[z][u].push_back(make_pair(v, p));
			}
		}
	}
	generator = new default_random_engine(seed);
	distribution = new uniform_real_distribution<double>(1.0-epsilon, std::nextafter(1.0, DBL_MAX));
	vector<double> xi;
	generate_xi(xi, n, *generator, *distribution);
	double user_start = get_current_time_sec();
	clock_t cpu_start = clock();

	X.resize(n);
	X.assign(n, false);
	active.resize(n);
	active.assign(n, false);
	cout<<"seed="<<seed<<", k = "<<k<<", eps = "<<epsilon<<endl;
	if (algo == "AS") {
		cout<<"tilde I AS: "<<endl;
		greedy(n, es, k, budget, beta, xi, get_xi);
		cout<<"I AS: "<<endl;
		greedy(n, es, k, budget, beta, xi, get_xi, true);
		cout<<"Random AS:"<<endl;
		Random(n, es, k, budget, beta, xi, get_xi);
		cout<<"Degree AS:"<<endl;
		degree(n, es, k, budget, beta, xi, get_xi);
	} else
	if (algo == "ADR") {
		cout<<"tilde I ADR: "<<endl;
		greedy(n, es, k, budget, beta, xi, get_xi, false, true);
		cout<<"I ADR: "<<endl;
		greedy(n, es, k, budget, beta, xi, get_xi, true, true);
		cout<<"Random ADR:"<<endl;
		Random(n, es, k, budget, beta, xi, get_xi, true);
		cout<<"Degree ADR:"<<endl;
		degree(n, es, k, budget, beta, xi, get_xi, true);
	}
	double user_end = get_current_time_sec();
	clock_t cpu_end = clock();
	JLOG_PUT("time.user", user_end - user_start);
	JLOG_PUT("time.cpu", (double) (cpu_end - cpu_start) / CLOCKS_PER_SEC );
}

