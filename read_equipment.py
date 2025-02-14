import json
import sys

# TODO:
# consider active of weapons
# consider effects over spells
# consider "quelqu'un vous suit" effect
# consider "Ajouter un sort temporaire" effect

SPECIAL_EFFECTS = [
	'Lié au personnage',
	'Arme de chasse',
	'Échangeable :',
	'-special spell-',
	'Craft coopératif impossible',
	'Reçu le :',
	'Titre :',
	'Change les paroles',
	'Change l\'apparence'
]

OBJECT_TYPES = [
	"Amulette",
	"Anneau",
	"Arme",
	"Bottes",
	"Bouclier",
	"Cape",
	"Ceinture",
	"Chapeau",
	"Dofus",
	"Familier"
]

OBJECT_MAP = {
	"Montilier": "Familier",
	"Prysmaradite": "Dofus",
	"Trophée": "Dofus"
}

MJ_ITEMS = {
	"Annobusé de Maître Jarbo",
	"Surpuissant Chacha de Combat (MJ)",
	"Petit Chacha de Combat (MJ)",
	"Amulette de Jiva",
	"Bague De Sendar"
}

file = open("all_equipment_fr.json", "r")
j = json.load(file)
j = j["items"]

items = []
EtypeIDs = {}
Etypes = []
itemsID = {}

def cond(it, c, root=True):
	global EtypeIDs, Etypes
	if c["is_operand"]:
		assert "children" not in c
		assert "relation" not in c
		c = c["condition"]
		assert c["operator"] in ["<", ">"]
		t = c["element"]
		if t["name"] in ["Kamas", "Niveau d'alignement"]:
			assert root
			it.append("none")
			return
		if t["id"] not in EtypeIDs:
			EtypeIDs[t["id"]] = len(Etypes)
			Etypes.append(t["name"])
		it.append(c["operator"])
		it.append(EtypeIDs[t["id"]])
		it.append(c["int_value"])
	else:
		assert "condition" not in c
		cs = c["children"]
		assert len(cs) == 2
		assert c["relation"] in ["and", "or"]
		it.append(c["relation"])
		cond(it, cs[0], False)
		cond(it, cs[1], False)

def effect(it, es, x=None):
	len_ind = len(it)
	it.append(len(es))
	for e in es:
		t = e["type"]
		if t["name"] in SPECIAL_EFFECTS:
			assert e["ignore_int_min"]
			assert e["ignore_int_max"]
			it[len_ind] -= 1
			continue
		if t["name"] in ["Attitude", "Nombre de victimes :", "Quelqu'un vous suit !", "Ajouter un sort temporaire"]:
			assert e["ignore_int_max"]
			it[len_ind] -= 1
			continue
		if t["name"] == "/":
			it[len_ind] -= 1
			continue
		if t["is_active"]:
			assert x["is_weapon"]
			it[len_ind] -= 1
			continue
		if t["name"][0] == ":":
			it[len_ind] -= 1
			continue
		assert not t["is_meta"]
		assert not e["ignore_int_min"]
		assert not e["ignore_int_max"] or e["int_maximum"] == 0
		if e["ignore_int_max"]: e["int_maximum"] = e["int_minimum"]
		if t["id"] not in EtypeIDs:
			EtypeIDs[t["id"]] = len(Etypes)
			Etypes.append(t["name"])
		it.append(EtypeIDs[t["id"]])
		it.append(e["int_minimum"])
		it.append(e["int_maximum"])
	if not it[len_ind]:
		if x is not None and "parent_set" in x:
			assert x["parent_set"]["name"] in [
				"Panoplie Ankarton", # Ne fonctionne qu'avec les bonus de pano (obligé de porter que des items de la pano)
				"Panoplie Téméraire", # Boost les sorts des iops
				"Panoplie de Gouttière", # Boost les sorts des ecaflips
				"Panoplie Criminelle", # Boost les sorts des scrams
				"Panoplie Intemporelle", # Boost les sorts des xélors
				"Panoplie Altruiste", # Boost les sorts des eniripsas
				"Panoplie Sauvage", # Boost les sorts des sadidas
				"Panoplie Exsangue", # Boost les sorts des sacrieurs
				"Panoplie de l'Innombrable", # Boost les sorts des osamodas
				"Panoplie Ethylique", # Boost les sorts des pandawas
				"Panoplie du Prince des voleurs", # Boost les sorts des cras
				"Panoplie Indestructible", # Boost les sorts des fécas
				"Panoplie Vénérable", # Boost les sorts des enutrofs
				"Panoplie Explosive", # Boost les sorts des roublards
				"Panoplie Lunatique", # Boost les sorts des zobals
				"Panoplie Submersible", # Boost les sorts des steamers
				"Panoplie Transcendante", # Boost les sorts des eliotropes
				"Panoplie Quadramentale", # Boost les sorts des huppermages
				"Panoplie Enragée", # Boost les sorts des ouginaks
				"Panoplie de l'Héritage", # Boost les sorts des forgelances
			]
		else:
			return False
	return True


for x in j:
	if x["name"] in MJ_ITEMS:
		continue

	it = [
		x["name"].replace(" ", "_"),
		x["level"],
		# x["pods"]
	]
	
	# type
	t = "Arme" if x["is_weapon"] else x["type"]["name"]
	if t.startswith("Certificat") or t == "Compagnon" or t.endswith("Percepteur"):
		continue
	if t in OBJECT_MAP: t = OBJECT_MAP[t]
	t_id = OBJECT_TYPES.index(t)
	it.append(OBJECT_TYPES.index(t))
	assert x["is_weapon"] == ("max_cast_per_turn" in x)
	assert x["is_weapon"] == ("ap_cost" in x)
	assert x["is_weapon"] == ("range" in x)
	assert x["is_weapon"] == ("critical_hit_probability" in x)
	assert x["is_weapon"] == ("critical_hit_bonus" in x)
	# if x["is_weapon"]:
		# it.append(x["max_cast_per_turn"])
		# it.append(x["ap_cost"])
		# it.append(x["range"]["min"])
		# it.append(x["range"]["max"])
		# it.append(x["critical_hit_probability"])
		# it.append(x["critical_hit_bonus"])

	if "effects" in x:
		if not effect(it, x["effects"], x):
			continue
	else:
		assert "parent_set" not in x
		continue

	if "conditions" in x:
		cond(it, x["conditions"])
	else:
		it.append("none")

	itemsID[x["ankama_id"]] = len(items)
	items.append(it)

print(len(Etypes))
for name in Etypes:
	print(name.replace(" ", "_"))
print(len(items))
for it in items:
	print(*(x for x in it))

file.close()
file = open("all_sets_fr.json", "r")
j = json.load(file)
j = j["sets"]

sets = []

for x in j:
	if x["contains_cosmetics_only"]:
		continue
	assert not x["contains_cosmetics"]

	s = [
		x["name"].replace(" ", "_"),
		x["items"]
	]

	for i in x["equipment_ids"]:
		i = itemsID[i]
		assert items[i][1] <= x["level"]
		assert items[i][2] != 8, (items[i][2], OBJECT_TYPES[items[i][2]], s[0])
		s.append(i)

	ess = x["effects"]
	s.append(len(ess))
	assert s[-1] <= 8
	
	for i in range(1, s[-1]+1):
		es = ess[str(i)]
		if es is None:
			s.append(0)
		else:
			effect(s, es)
	
	sets.append(s)

print(len(sets))
for s in sets:
	print(*(x for x in s))
