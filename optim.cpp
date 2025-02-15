#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include <array>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <map>

#include "weights.h"

using namespace std;

const int LVL = 165;
const int PA = 12;
const int PM = 5;
static_assert(PM > 4);

constexpr int EFFECTS = 52;
constexpr int SLOTS = 16;
constexpr int NDOFUS = 6;
constexpr int SLOTS_WD = SLOTS - NDOFUS;

unordered_set<string> blacklist = {
	"Dofus_Ocre"
};

struct Effect {
	int id, min, max;
};
array<string, EFFECTS> effects_name;
unordered_map<string, int> s2i;

struct Condition {
	enum TYPE { AND, OR, INF, SUP };
	TYPE type;
	union { int a, effect_id; };
	union { int b, value; };
};
vector<Condition> conditions;

struct Item {
	enum TYPE {
		Amulette,
		Anneau,
		Arme,
		Bottes,
		Bouclier,
		Cape,
		Ceinture,
		Chapeau,
		Dofus,
		Familier,
		SIZE
	};

	string name;
	int level;
	TYPE type;
	vector<Effect> effects; 
	int condition;
	int pano = 0;
};
vector<Item> items;

struct Pano {
	string name;
	vector<int> items;
	vector<vector<Effect>> effects;
};
vector<Pano> panos;

int read_condition(istream &stream) {
	string cond_type;
	stream >> cond_type;
	if(cond_type == "none") return -1;
	const int ans = conditions.size();
	if(cond_type.size() == 1) {
		conditions.push_back({cond_type[0] == '<' ? Condition::TYPE::INF : Condition::TYPE::SUP});
		stream >> conditions.back().effect_id >> conditions.back().value;
		return ans;
	}
	conditions.push_back({cond_type == "and" ? Condition::TYPE::AND : Condition::TYPE::OR});
	conditions[ans].a = read_condition(stream);
	conditions[ans].b = read_condition(stream);
	return ans;
}

double score(const array<double, EFFECTS> &W, const Item &it) {
	double sc = 0.;
	for(const Effect &e : it.effects) sc += W[e.id] * e.max;
	return sc;
}

bool check_cond(const array<int, EFFECTS> &stats, int c) {
	const Condition &cond = conditions[c];
	if(cond.type == Condition::AND || cond.type == Condition::OR) {
		cout << '[';
		const bool a = check_cond(stats, cond.a);
		cout << (cond.type == Condition::AND ? " AND " : " OR ");
		const bool b = check_cond(stats, cond.b);
		cout << ']';
		return cond.type == Condition::AND ? (a&&b) : (a||b);
	}
	cout << effects_name[cond.effect_id] << '(' << stats[cond.effect_id] << ") " << "<>"[cond.type-2] << ' ' << cond.value;
	return cond.type == Condition::INF ? stats[cond.effect_id] < cond.value : stats[cond.effect_id] > cond.value;
}

bool check_cond_PA_PM(int c) {
	if(c == -1) return true;
	const Condition &cond = conditions[c];
	if(cond.type == Condition::AND || cond.type == Condition::OR) {
		const bool a = check_cond_PA_PM(cond.a);
		const bool b = check_cond_PA_PM(cond.b);
		return cond.type == Condition::AND ? (a&&b) : (a||b);
	}
	if(cond.effect_id == s2i["PA"]) return cond.type == Condition::INF ? PA < cond.value : PA > cond.value;
	if(cond.effect_id == s2i["PM"]) return cond.type == Condition::INF ? PM < cond.value : PM > cond.value;
	return true;
}

bool check_item(const Item &it) {
	return it.level <= LVL && check_cond_PA_PM(it.condition) && !blacklist.count(it.name);
}

array<int, SLOTS> optimize(const array<double, EFFECTS> &W) {
	constexpr int PMA = PM+2;
	const int PMID = s2i["PM"];
	array<int, SLOTS> sol;
	vector<pair<double, int>> dofus[3];

	// Init DP
	array<
		array<array<
			pair<double, array<int, SLOTS_WD>>,
			PMA>, 4>,
		(1<<SLOTS_WD)
	> A; // A[m][bonus][pm]
	for(int m = 0; m < (int) A.size(); ++m)
		for(int b = 0; b < 4; ++b)
			for(int pm = 0; pm < PMA; ++pm) {
				A[m][b][pm].first = numeric_limits<double>::min();
				A[m][b][pm].second.fill(-2);
			}
	A[0][0][2].first = 0.;
	A[0][0][2].second.fill(-1);

	// Update for each item not int pano
	map<pair<int, int>, pair<int, int>> usefull_items;
	for(int i = 1; i < (int) items.size(); ++i) {
		const Item &it = items[i];
		if(it.pano) continue;
		if(!check_item(it)) continue;;
		int pm = 0;
		for(const Effect &e : it.effects) if(e.id == PMID) pm += e.max;
		if(it.type == Item::Dofus) {
			assert(-1 <= pm && pm <= 1);
			dofus[pm+1].emplace_back(score(W, it), i);
			continue;
		}
		pair<int, int> tpm {it.type, pm};
		if(usefull_items.count(tpm)) {
			const double sc = score(W, it);
			auto &best = usefull_items[tpm];
			if(it.type == Item::Anneau) {
				if(best.second == -1 || score(W, items[best.second]) < sc) {
					best.second = i;
					if(score(W, items[best.first]) < sc)
						swap(best.first, best.second);
				}
			} else if(score(W, items[best.first]) < sc)
				best.first = i;
		} else usefull_items[tpm] = make_pair(i, -1);
	}
	for(const auto &[k, pi] : usefull_items) {
		const int i = pi.first;
		const auto &[t, pmi] = k;
		const double sci = score(W, items[i]);
		const int m2 = 1<<t;
		const int PM_0 = max(-pmi, 0);
		const int PM_1 = PMA - max(pmi, 0);
		const auto update = [&](const int m, const int m2, const int PM_0, const int PM_1, const int pm2, const double sc2, const int j=-1) {
			if(m&m2) return;
			const int m3 = m|m2;
			for(int pm = PM_0; pm < PM_1; ++pm) {
				if(A[m][0][pm].first == numeric_limits<double>::min()) continue;
				const double sc3 = A[m][0][pm].first + sc2;
				const int pm3 = pm + pm2;
				if(A[m3][0][pm3].first >= sc3) continue;
				A[m3][0][pm3].first = sc3;
				A[m3][0][pm3].second = A[m][0][pm].second;
				A[m3][0][pm3].second[t] = i;
				if(j != -1) A[m3][0][pm3].second[Item::Dofus] = j;
			}
		};
		if(pi.second != -1) {
			const int j = pi.second;
			const int m2b = m2 | (1<<Item::Dofus);
			const int PM_0b = max(-2*pmi, 0);
			const int PM_1b = PMA - max(2*pmi, 0);
			const int pmib = 2*pmi;
			const int scib = sci + score(W, items[j]);
			for(int m = A.size()-1; m >= 0; --m) {
				if((m>>Item::Dofus)&1 && !((m>>Item::Anneau)&1)) continue;
				update(m, m2, PM_0, PM_1, pmi, sci);
				update(m, m2b, PM_0b, PM_1b, pmib, scib, j);
			}
		} else {
			for(int m = A.size()-1; m >= 0; --m) {
				if((m>>Item::Dofus)&1 && !((m>>Item::Anneau)&1)) continue;
				update(m, m2, PM_0, PM_1, pmi, sci);
			}
		}
	}

	// Update for each pano
	for(const Pano &p : panos) {
		struct ItemStats {
			double score;
			int pm;
			int i;
		};
		vector<ItemStats> its;
		vector<double> bonus;
		vector<int> bonus_pm;
		for(int i : p.items) if(check_item(items[i])) {
			int pm = 0;
			for(const Effect &e : items[i].effects) if(e.id == PMID) pm += e.max;
			its.push_back({score(W, items[i]), pm, i});
		}
		for(const vector<Effect> &es : p.effects) {
			double sc = 0;
			int pm = 0;
			for(const Effect &e : es)
				if(e.id == PMID) pm += e.max;
				else sc += W[e.id] * e.max;
			bonus.push_back(sc);
			bonus_pm.push_back(pm);
		}

		struct Subset {
			vector<pair<int, int>> updt;
			double score = 0.;
			int m = 0;
			int bonus, pm;
		};
		vector<Subset> subsets;
		const int nits = its.size();
		const int MP = 1<<nits;
		for(int mp = 1; mp < MP; ++mp) {
			Subset &ss = subsets.emplace_back();
			for(int u = 0; u < nits; ++u) if((mp>>u)&1) {
				const auto &[sc, pm, i] = its[u];
				int t = items[i].type;
				if(ss.m&(1<<t)) {
					if(t == Item::Arme) goto bad_sub;
					if(t == Item::Anneau) {
						t = Item::Dofus;
						if(ss.m&(1<<t)) goto bad_sub;
					} else assert(false);
				}
				ss.score += sc;
				ss.pm += pm;
				ss.m |= 1<<t;
				ss.updt.emplace_back(t, i);
			}
			ss.bonus = min(ss.updt.size()-1, bonus.size()-1);
			ss.score += bonus[ss.bonus];
			ss.pm += bonus_pm[ss.bonus];
			if(ss.pm >= PMA) goto bad_sub;
			if((ss.m>>Item::Anneau)&1 && !((ss.m>>Item::Dofus)&1)) {
				subsets.push_back(ss);
				subsets.back().m ^= (1<<Item::Anneau) | (1<<Item::Dofus);
				for(auto &[t, i] : subsets.back().updt) if(t == Item::Anneau) {
					t = Item::Dofus;
					break;
				}
			}
			continue;
			bad_sub:
			subsets.pop_back();
			continue;
		}

		for(int m = A.size()-1; m >= 0; --m) {
			if((m>>Item::Dofus)&1 && !((m>>Item::Anneau)&1)) continue;
			for(const auto &[updt, sc2, m2, bonus2, pm2] : subsets) {
				if(m&m2) continue;
				const int m3 = m|m2;
				const int PMA_0 = max(-pm2, 0);
				const int PMA_1 = PMA - max(pm2, 0);
				for(int bonus = 0; bonus < 4; ++bonus) for(int pm = PMA_0; pm < PMA_1; ++pm) {
					if(A[m][bonus][pm].first == numeric_limits<double>::min()) continue;
					const double sc3 = A[m][bonus][pm].first + sc2;
					const int bonus3 = min(bonus + bonus2, 3);
					const int pm3 = pm + pm2;
					if(A[m3][bonus3][pm3].first >= sc3) continue;
					A[m3][bonus3][pm3].first = sc3;
					A[m3][bonus3][pm3].second = A[m][bonus][pm].second;
					for(const auto &[t, i] : updt) A[m3][bonus3][pm3].second[t] = i;
				}
			}
		}
	}

	const auto get_dofus_sets = [&]() {
		assert(dofus[0].size() + dofus[2].size() <= NDOFUS);
		assert(dofus[0].size() + dofus[1].size() + dofus[2].size() >= NDOFUS);
		array<pair<double, array<int, NDOFUS>>, 5> dofus_sets;
		for(auto &[sc, s] : dofus_sets) sc = numeric_limits<double>::min();
		for(int pm = 0; pm < 3; ++pm)
			sort(dofus[pm].rbegin(), dofus[pm].rend());
		double scdx = 0;
		array<int, NDOFUS> ds;
		for(int x = 0;; ++x) {
			double scdy = scdx;
			for(int y = 0;; ++y) {
				const int xy = x+y;
				const int pm = 1-x+y;
				const int z = NDOFUS-x-y;
				double sc = scdy;
				for(int i = 0; i < z; ++i) {
					sc += dofus[1][i].first;
					ds[xy+i] = dofus[1][i].second;
				}
				if(sc > dofus_sets[pm].first) {
					dofus_sets[pm].first = sc;
					dofus_sets[pm].second = ds;
				}
				if(y >= (int)dofus[2].size()) break;
				scdy += dofus[2][y].first;
				ds[xy] = dofus[2][y].second;
			}
			if(x >= (int)dofus[0].size()) break;
			scdx += dofus[0][x].first;
			ds[x] = dofus[0][x].second;
		}
		return dofus_sets;
	};
	const auto dofus_sets_1 = get_dofus_sets();
	const int BPID = s2i["Bonus_de_panoplies"];
	for(int pm = 0; pm < 3; ++pm) {
		for(int i = 0; i < (int) dofus[pm].size(); ++i) {
			const int c = items[dofus[pm][i].second].condition;
			if(c == -1) continue;;
			if(conditions[c].type < 2) continue;
			if(conditions[c].effect_id != BPID) continue;
			assert(conditions[c].value == 3);
			dofus[pm][i--] = dofus[pm].back();
			dofus[pm].pop_back();
		}
	}
	const auto dofus_sets_2 = get_dofus_sets();

	int best_bonus = -1, best_pm = -1;
	double best_score = numeric_limits<double>::min();
	for(int b = 0; b < 4; ++b) {
		const auto &ds = b == 3 ? dofus_sets_2 : dofus_sets_1;
		for(int pm = 0; pm < 5; ++pm) {
			const double sc = A.back()[b][PM-pm].first + ds[pm].first;
			if(sc > best_score) {
				best_score = sc;
				best_bonus = b;
				best_pm = pm;
				copy_n(A.back()[b][PM-pm].second.begin(), SLOTS_WD, sol.begin());
				copy_n(ds[pm].second.begin(), NDOFUS, sol.begin()+SLOTS_WD);
			}
		}
	}
	
	cerr << "Score: " << best_score << ' ' << best_pm << ' ' << PM-best_pm << endl;
	
	return sol;
}

int main() {
	ifstream data("equipments.txt");
	if(!data) {
		cerr << "can't read database...\n";
		return 1;
	}

	int N_EFFECTS_IN;
	data >> N_EFFECTS_IN;
	assert(N_EFFECTS_IN == EFFECTS);
	for(string &name : effects_name) data >> name;

	int N_ITEMS;
	data >> N_ITEMS;
	items.resize(N_ITEMS);
	for(Item &it : items) {
		int type, n_effects;
		data >> it.name >> it.level >> type >> n_effects;
		it.type = (Item::TYPE) type;
		it.effects.resize(n_effects);
		for(auto &[id, min, max] : it.effects) data >> id >> min >> max;
		it.condition = read_condition(data);
	}

	int N_PANOS;
	data >> N_PANOS;
	panos.resize(N_PANOS);
	N_PANOS = 0;
	for(Pano &p : panos) {
		int n_items, n_effects;
		data >> p.name >> n_items;
		p.items.resize(n_items);
		++ N_PANOS;
		for(int &i : p.items) {
			data >> i;
			items[i].pano = N_PANOS;
		}
		data >> n_effects;
		p.effects.resize(n_effects);
		for(vector<Effect> &es : p.effects) {
			data >> n_effects;
			es.resize(n_effects);
			for(auto &[id, min, max] : es) data >> id >> min >> max;
		}
	}

	for(int i = 0; i < EFFECTS; ++i)
		s2i[effects_name[i]] = i;

	array<double, EFFECTS> W;
	for(const auto &[s, v] : MULTI_DO_FIXE_12_5_165)
		W[s2i[s]] = v;

	array<int, SLOTS> sol = optimize(W);
	for(int i : sol) {
		for(char &c : items[i].name) if(c == '_') c = ' ';
		cout << items[i].name << '\n';
	}

	array<int, EFFECTS> stats;
	stats.fill(0);
	stats[s2i["PA"]] = LVL < 100 ? 6 : 7;
	stats[s2i["PM"]] = 3;
	unordered_map<int, int> pan_set;
	for(int i : sol) {
		for(const auto &[id, min, max] : items[i].effects) stats[id] += max;
		if(items[i].pano) {
			++ pan_set[items[i].pano];
		}
	}
	for(const auto &[p, c] : pan_set) if(c > 1) {
		cout << "Bonus " << panos[p-1].name << " " << c << "\n";
		stats[s2i["Bonus_de_panoplies"]] += c-1;
		for(const Effect &e : panos[p-1].effects[c-1])
			stats[e.id] += e.max;
	}
	for(const string &s : {"PA", "PM"})
		cout << s << ": " << stats[s2i[s]] << '\n';
	
	for(int i : sol) {
		const int c = items[i].condition;
		if(c == -1) continue;
		cout << items[i].name << endl;
		const bool v = check_cond(stats, c);
		cout << "  => " << v << endl;
	}

	return 0;
}