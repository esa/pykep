from ._cassini_mga import cassini1, cassini1_a, cassini1_n
from ._cassini_mga_1dsm import cassini2
from ._eve_mga_1dsm import eve_mga1dsm, eve_mga1dsm_a, eve_mga1dsm_n    
from ._rosetta_mga_1dsm import rosetta
from ._juice_mga_1dsm import juice, juice_mo
from ._messenger import messenger
from ._em_Nimp import em3imp, em5imp, em7imp


# Load TOPS benchmark --------------------------------------------------------
import json
from importlib import resources
def _load_tops_json(name: str):
    with resources.as_file(resources.files(__package__) / "tops" / name) as _path:
        with _path.open("r", encoding="utf-8") as _f:
            return json.load(_f)


zoh_cr3bp = _load_tops_json("_tops_cr3bp.json")
zoh_twobody = _load_tops_json("_tops_twobody.json")
zoh_ss = _load_tops_json("_tops_ss.json")
zoh_mee = _load_tops_json("_tops_mee.json")

del _load_tops_json, resources, json