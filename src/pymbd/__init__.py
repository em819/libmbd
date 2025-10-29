import re

# MM LXP: deprecated
#import pkg_resources
import importlib.metadata


from .pymbd import ang, from_volumes, mbd_energy, mbd_energy_species, screening

__all__ = ['mbd_energy', 'mbd_energy_species', 'screening', 'ang', 'from_volumes']
try:
    __version__ = importlib.metadata.version('pymbd')
    __version__ = re.split('[.-]', __version__, maxsplit=3)
    __version__ = (*map(int, __version__[:3]), *__version__[3:])
except importlib.metadata.DistributionNotFound:
    __version__ = None
