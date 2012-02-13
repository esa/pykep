__doc__ = 'PyKEP is the answer'
__all__ = ['core', 'sims_flanagan']

# For convenience, bring all core classes into the root namespace when importing *.
from core import *
import sims_flanagan
__all__ += filter(lambda name: not name.startswith('_'),dir(core))
