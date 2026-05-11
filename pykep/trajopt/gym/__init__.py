from ._cassini_mga import cassini1, cassini1_a, cassini1_n
from ._cassini_mga_1dsm import cassini2
from ._eve_mga_1dsm import eve_mga1dsm, eve_mga1dsm_a, eve_mga1dsm_n    
from ._rosetta_mga_1dsm import rosetta
from ._juice_mga_1dsm import juice, juice_mo
from ._messenger import messenger
from ._em_Nimp import em3imp, em5imp, em7imp


# Load TOPS benchmark --------------------------------------------------------
# This loads the raw data; the actual TOPS problems are constructed in tops.py, which imports these variables.
import json
from importlib import resources
def _load_tops_json(name: str):
    with resources.as_file(resources.files(__package__) / "tops" / name) as _path:
        with _path.open("r", encoding="utf-8") as _f:
            return json.load(_f)

tops_cr3bp_json = _load_tops_json("_tops_cr3bp.json")
tops_twobody_json = _load_tops_json("_tops_twobody.json")
tops_ss_json = _load_tops_json("_tops_ss.json")
tops_mee_json = _load_tops_json("_tops_mee.json")

from ._tops import tops_twobody, tops_mee, tops_cr3bp, tops_ss

del _load_tops_json, resources, json