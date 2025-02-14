#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include <array>
#include <limits>
#include <algorithm>
#include <unordered_map>

using namespace std;

constexpr int EFFECTS = 52;
constexpr int SLOTS = 16;
constexpr int NDOFUS = 6;
constexpr int SLOTS_WD = SLOTS - NDOFUS;

struct Effect {
	int id, min, max;
};
array<string, EFFECTS> effects_name;

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

array<int, SLOTS> optimize(const array<double, EFFECTS> &W) {
	array<int, SLOTS> sol;

	array<pair<int, double>, SLOTS_WD> best_solo;
	vector<pair<double, int>> dofus;
	best_solo.fill(make_pair(-1, numeric_limits<double>::min()));
	for(int i = 0; i < (int) items.size(); ++i) {
		const Item &it = items[i];
		if(it.pano) continue;
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

	array<pair<double, array<int, SLOTS_WD>>, 1<<SLOTS_WD> A, B;
	for(int m = 0; m < (int) A.size(); ++m) {
		if((m>>Item::Dofus)&1 && !((m>>Item::Anneau)&1)) continue;
		auto &[sc, a] = A[m];
		sc = 0;
		a.fill(-1);
		for(int u = 0; u < SLOTS_WD; ++u) if((m>>u)&1) {
			sc += best_solo[u].first;
			a[u] = best_solo[u].second;
		}
	}

	cerr << A.back().first << endl;

	for(const Pano &p : panos) {
		B = A;
		vector<pair<double, int>> its;
		vector<double> bonus;
		for(int i : p.items) its.emplace_back(score(W, items[i]), items[i].type);
		for(const vector<Effect> &es : p.effects) {
			double sc = 0;
			for(const Effect &e : es) sc += W[e.id] * e.max;
			bonus.push_back(sc);
		}
		vector<pair<int, int>> cur;

		const int ps = its.size();
		const int MP = 1<<ps;
		for(int mp = 1; mp < MP; ++mp) {
			double scp = 0;
			int m2 = 0;
			cur.clear();
			for(int u = 0; u < ps; ++u) if((mp>>u)&1) {
				int t = its[u].second;
				if(m2&(1<<t)) {
					if(t == Item::Arme) goto end_of_pano_sub;
					if(t == Item::Anneau) {
						t = Item::Dofus;
						if(m2&(1<<t)) goto end_of_pano_sub;
					} else assert(false);
				}
				scp += its[u].first;
				m2 |= 1<<t;
				cur.emplace_back(t, p.items[u]);
			}
			scp += bonus[min(cur.size()-1, p.effects.size()-1)];
			
			updt_A:
			for(int m = 0; m < (int) B.size(); ++m) {
				if((m>>Item::Dofus)&1 && !((m>>Item::Anneau)&1)) continue;
				if(m&m2) continue;
				const double sc2 = B[m].first + scp;
				const int m3 = m|m2;
				if(A[m3].first >= sc2) continue;
				A[m3].first = sc2;
				A[m3].second = B[m].second;
				for(const auto &[t, i] : cur) A[m3].second[t] = i;
			}

			if((m2>>Item::Anneau)&1 && !((m2>>Item::Dofus)&1)) {
				m2 ^= (1<<Item::Anneau) | (1<<Item::Dofus);
				for(auto &[t, i] : cur) if(t == Item::Anneau) {
					t = Item::Dofus;
					break;
				}
				goto updt_A;
			}

			// TODO: second slot of ring
			end_of_pano_sub:
			continue;
		}
	}

	cerr << A.back().first << endl;

	sort(dofus.rbegin(), dofus.rend());
	for(int i = 0; i < NDOFUS; ++i)
		sol[SLOTS_WD + i] = dofus[i].second;
	for(int i = 0; i < SLOTS_WD; ++i)
		sol[i] = A.back().second[i];
	
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
	cout << EFFECTS << endl;
	
	int N_ITEMS;
	data >> N_ITEMS;
	items.resize(N_ITEMS);
	array<int, Item::TYPE::SIZE> tc;
	tc.fill(0);
	cout << N_ITEMS << endl;
	for(Item &it : items) {
		int type, n_effects;
		data >> it.name >> it.level >> type >> n_effects;
		it.type = (Item::TYPE) type;
		++ tc[it.type];
		it.effects.resize(n_effects);
		for(auto &[id, min, max] : it.effects) data >> id >> min >> max;
		it.condition = read_condition(data);
	}
	for(int c : tc) cout << c << " ";
	cout << endl;

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
	cout << N_PANOS << endl;

	unordered_map<string, int> s2i;
	for(int i = 0; i < EFFECTS; ++i)
		s2i[effects_name[i]] = i;

	array<double, EFFECTS> W;
	W.fill(0.0);
	W[s2i["Force"]] = 3.;
	W[s2i["Dommage_Terre"]] = 9.;
	W[s2i["Vitalité"]] = 1.;
	W[s2i["\%_Critique"]] = 6.;
	W[s2i["Tacle"]] = 0.2;
	W[s2i["PM"]] = 1.;
	W[s2i["PA"]] = 180.;
	W[s2i["Dommage"]] = 10.;
	W[s2i["Puissance"]] = 4.;
	W[s2i["Dommage_Critiques"]] = 4.;
	W[s2i["Dommage_Poussée"]] = 1.;
	W[s2i["Retrait_PM"]] = 1.;
	W[s2i["Retrait_PA"]] = 1.;
	W[s2i["%_Résistance_Terre"]] = 3.;
	W[s2i["%_Résistance_Neutre"]] = 3.;
	W[s2i["%_Résistance_Eau"]] = 3.;
	W[s2i["%_Résistance_Air"]] = 3.;
	W[s2i["%_Résistance_Feu"]] = 3.;
	W[s2i["Résistance_Terre"]] = 2.;
	W[s2i["Résistance_Neutre"]] = 2.;
	W[s2i["Résistance_Eau"]] = 2.;
	W[s2i["Résistance_Air"]] = 2.;
	W[s2i["Résistance_Feu"]] = 2.;

	array<int, SLOTS> sol = optimize(W);
	for(int i : sol) cout << items[i].name << ' ' << items[i].pano << ' ' << score(W, items[i]) << endl;

	return 0;
}