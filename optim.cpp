#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include <array>
#include <limits>
#include <algorithm>
#include <unordered_map>

#include "weights.h"

using namespace std;

const int LVL = 165;
const int PA = 12;
const int PM = 5;

constexpr int EFFECTS = 52;
constexpr int SLOTS = 16;
constexpr int NDOFUS = 6;
constexpr int SLOTS_WD = SLOTS - NDOFUS;

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
	return it.level <= LVL && check_cond_PA_PM(it.condition);
}

array<int, SLOTS> optimize(const array<double, EFFECTS> &W) {
	array<int, SLOTS> sol;

	array<pair<int, double>, SLOTS_WD> best_solo;
	vector<pair<double, int>> dofus;
	best_solo.fill(make_pair(-1, numeric_limits<double>::min()));
	for(int i = 0; i < (int) items.size(); ++i) {
		const Item &it = items[i];
		if(it.pano) continue;
		if(!check_item(it)) continue;
		const double sc = score(W, it);
		switch(it.type) {
		case Item::Anneau:
			if(best_solo[Item::Dofus].first < sc) {
				best_solo[Item::Dofus] = make_pair(sc, i);
				if(best_solo[Item::Anneau].first < sc)
					swap(best_solo[Item::Anneau], best_solo[Item::Dofus]);
			}
			break;
		case Item::Dofus:
			dofus.emplace_back(sc, i);
			break;
		default:
			if(best_solo[it.type].first < sc)
				best_solo[it.type] = make_pair(sc, i);
		}
	}

	array<
		// array<array<
		array<
			pair<double, array<int, SLOTS_WD>>,
			// 4>, PM-2>,
			4>,
		(1<<SLOTS_WD)
	> A;
	for(int m = 0; m < (int) A.size(); ++m) {
		for(int b = 1; b < 4; ++b) A[m][b].first = numeric_limits<double>::min();
		if((m>>Item::Dofus)&1 && !((m>>Item::Anneau)&1)) {
			A[m][0].first = numeric_limits<double>::min();
			continue;
		}
		auto &[sc, a] = A[m][0];
		sc = 0;
		a.fill(-1);
		for(int u = 0; u < SLOTS_WD; ++u) if((m>>u)&1) {
			sc += best_solo[u].first;
			a[u] = best_solo[u].second;
		}
	}

	for(const Pano &p : panos) {
		vector<pair<double, int>> its;
		vector<double> bonus;
		for(int i : p.items) if(check_item(items[i]))
			its.emplace_back(score(W, items[i]), i);
		for(const vector<Effect> &es : p.effects) {
			double sc = 0;
			for(const Effect &e : es) sc += W[e.id] * e.max;
			bonus.push_back(sc);
		}

		struct Subset {
			vector<pair<int, int>> updt;
			double score = 0.;
			int m = 0;
			int bonus;
		};
		vector<Subset> subsets;
		const int nits = its.size();
		const int MP = 1<<nits;
		for(int mp = 1; mp < MP; ++mp) {
			Subset &ss = subsets.emplace_back();
			for(int u = 0; u < nits; ++u) if((mp>>u)&1) {
				const auto &[sc, i] = its[u];
				int t = items[i].type;
				if(ss.m&(1<<t)) {
					if(t == Item::Arme) goto bad_sub;
					if(t == Item::Anneau) {
						t = Item::Dofus;
						if(ss.m&(1<<t)) goto bad_sub;
					} else assert(false);
				}
				ss.score += sc;
				ss.m |= 1<<t;
				ss.updt.emplace_back(t, i);
			}
			ss.score += bonus[ss.bonus = min(ss.updt.size()-1, bonus.size()-1)];
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
			if(A[m][0].first == numeric_limits<double>::min()) continue;
			for(const auto &[updt, score, m2, bonus2] : subsets) {
				if(m&m2) continue;
				for(int bonus = 0; bonus < 4; ++bonus) {
					if(A[m][bonus].first == numeric_limits<double>::min()) continue;
					const double sc2 = A[m][bonus].first + score;
					const int m3 = m|m2;
					const int bonus3 = min(bonus + bonus2, 3);
					if(A[m3][bonus3].first >= sc2) continue;
					A[m3][bonus3].first = sc2;
					A[m3][bonus3].second = A[m][bonus].second;
					for(const auto &[t, i] : updt) A[m3][bonus3].second[t] = i;
				}
			}
		}
	}

	sort(dofus.rbegin(), dofus.rend());
	double dofus_score_1 = 0., dofus_score_2 = 0.;
	array<int, NDOFUS> dofus_1, dofus_2;
	dofus_1.fill(-1);
	dofus_2.fill(-1);
	for(int i = 0; i < NDOFUS; ++i) {
		dofus_score_1 += dofus[i].first;
		dofus_1[i] = dofus[i].second;
	}
	for(int i = 0, j = 0; i < NDOFUS; ++i, ++j) {
		while(j < (int) dofus.size()) {
			const int c = items[dofus[j].second].condition;
			if(c == -1) break;
			if(conditions[c].type < 2) break;
			if(conditions[c].effect_id != s2i["Bonus_de_panoplies"]) break;
			assert(conditions[c].value == 3);
			++j;
		}
		if(j >= (int) dofus.size()) break;
		dofus_score_2 += dofus[j].first;
		dofus_2[i] = dofus[j].second;
	}

	int best_bonus = max_element(A.back().begin(), A.back().begin()+3) - A.back().begin();
	double best_score = A.back()[best_bonus].first + dofus_score_1;
	if(best_score < A.back()[3].first + dofus_score_2) {
		best_bonus = 3;
		best_score = A.back()[3].first + dofus_score_2;
		copy_n(dofus_2.begin(), NDOFUS, sol.begin()+SLOTS_WD);
	} else {
		copy_n(dofus_1.begin(), NDOFUS, sol.begin()+SLOTS_WD);
	}
	for(int i = 0; i < SLOTS_WD; ++i)
		sol[i] = A.back()[best_bonus].second[i];
	
	cerr << "Score: " << best_score << endl;
	
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
			if(pan_set.count(items[i].pano)) ++ stats[s2i["Bonus_de_panoplies"]];
			++ pan_set[items[i].pano];
		}
	}
	for(const auto &[p, c] : pan_set) if(c > 1)
		cout << "Bonus " << panos[p-1].name << " " << c << "\n";
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