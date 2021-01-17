//Provided by the authors of https://papers.nips.cc/paper/5709-monotone-k-submodular-function-maximization-with-size-constraints.pdf#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <bitset>
#include <algorithm>
#include <map>
#include <set>
//#include <unordered_map>
//#include <unordered_set>
#include <ctime>
#include <cmath>

#include <sys/time.h>
#include <unistd.h> // getpid
#include <memory> // auto_ptr
#include "jlog.h"

using namespace std;
using namespace jlog_internal;
typedef long long LL;
typedef unsigned long long ULL;

#define ERROR "\x1b[41mERROR: "
#define WARN "\x1b[45mWARN: "
#define INFO "\x1b[32m\x1b[1mINFO: "
#define CALL "\x1b[44mCALL: "
#define QUERY "\x1b[45m\x1b[1mQUERY: "
#define EMPH "\x1b[37m\x1b[1m"
#define DEF "\x1b[0m\n"

class Xorshift {
public:
	Xorshift(int seed) {
		x = _(seed, 0);
		y = _(x, 1);
		z = _(y, 2);
		w = _(z, 3);
	}

	int _(int s, int i) {
		return 1812433253 * (s ^ (s >> 30)) + i + 1;
	}

	// 32bit signed
	inline int gen_int() {
		unsigned int t = x ^ (x << 11);
		x = y;
		y = z;
		z = w;
		return w = w ^ (w >> 19) ^ t ^ (t >> 8);
	}

	// error = O(n*2^-32)
	inline int nextInt(int n) {
		return (int) (n * nextDouble());
	}

	// [0, 1) (53bit)
	inline double nextDouble() {
		unsigned int a = ((unsigned int) gen_int()) >> 5, b =
				((unsigned int) gen_int()) >> 6;
		return (a * 67108864.0 + b) * (1.0 / (1LL << 53));
	}

private:
	unsigned int x, y, z, w;
};

LL mem_usage() {
	auto_ptr<istream> is(new ifstream("/proc/self/stat"));
	for (int i = 0; i < 23; i++) {
		string s;
		*is >> s;
		if (i == 22) {
			return atol(s.c_str());
		}
	}
	return -1;
}

vector<vector<pair<int, vector<double> > > > load_graph(string input,
		string prob, int k) {
	FILE *fp = fopen(input.c_str(), "r");
	vector<pair<pair<int, int>, vector<double> > > ps;
	int n = 0, m = 0;
	for (int u, v; fscanf(fp, "%d\t%d", &u, &v) != EOF;) {
		vector<double> vec;
		ps.push_back(make_pair(make_pair(u, v), vec));
		n = max(n, max(u, v) + 1);
		m++;
	}
	fclose(fp);
	JLOG_PUT("graph.n", n);
	JLOG_PUT("graph.m", m);

	fp = fopen(prob.c_str(), "r");
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < k; j++) {
			double p;
			fscanf(fp, "%lf", &p);
			ps[i].second.push_back(p);
		}
	}
	fclose(fp);
	vector<vector<pair<int, vector<double> > > > es(n);
	for (int i = 0; i < m; i++) {
		int u = ps[i].first.first, v = ps[i].first.second;
		es[u].push_back(make_pair(v, ps[i].second));
	}

	return es;
}

int get_n(vector<pair<int, int> > &ps) {
	int n = 0;
	for (LL i = 0; i < ps.size(); i++) {
		n = max(n, max(ps[i].first, ps[i].second) + 1);
	}
	return n;
}

void init_args(int argc, char *argv[], map<string, string> &args) {
	for (int i = 1; i < argc; i++) {
		string a(argv[i]);
		if (a[0] != '-') {
			continue;
		}
		int at = a.find("=");
		string key = a.substr(1, at - 1);
		string val = a.substr(at + 1);
		args[key] = val;
		cerr << key << "=" << val << endl;
		JLOG_PUT(("params." + key).c_str(), val);
	}
}

string get_or_die(map<string, string> &argv, string key) {
	if (argv.count(key) == 0) {
		printf(ERROR "Illegal key %s (get_or_die) " DEF, key.c_str());
		exit(1);
	}
	return argv[key];
}

double average(vector<LL> &vec) {
	LL sum = 0;
	for (int i = 0; i < vec.size(); i++) {
		sum += vec[i];
	}
	return 1.0 * sum / vec.size();
}
