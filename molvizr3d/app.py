# import packages/modules
import os

# internal
from .config import packageName
from .config import packageShortName
from .config import __version__
from .config import __description__
from .docs import MolParser
from .docs import Compound


def main():
    # short description and version
    _des = f"{packageShortName} {__version__}: {__description__}"
    print(_des)


def td(f):
    '''
    3d visualizer of a molecule

    Parameters
    ----------
    f : str
        molecule file format e.g sdf, mol, json

    Returns
    -------
    None
    '''
    # check file exists
    if os.path.exists(f):
        # parse file
        MolParserC = MolParser(f)
        compound_info = MolParserC.read_file()
        # compound
        compound = Compound(compound_info)
        return compound
    else:
        raise Exception("file path is not valid.")


if __name__ == "__main__":
    main()
