__version__ = '0.7.7'

try:
    __ISOCHRONES_SETUP__
except NameError:
    __ISOCHRONES_SETUP__ = False

if not __ISOCHRONES_SETUP__:
    __all__ = ['dartmouth','basti','padova']
    from .isochrone import Isochrone
    from .starmodel import StarModel, BinaryStarModel, TripleStarModel
     
