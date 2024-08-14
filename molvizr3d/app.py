# import packages/modules
import os
import pubchemquery as pcq
import pandas as pd
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


def td(file, display_legend=True):
    '''
    3d visualizer of a compound

    Parameters
    ----------
    file : str
        molecule file format (sdf)
    display_legend : bool
        display legend (default True)

    Returns
    -------
    None
        display 3d 
    '''
    # check file exists
    if os.path.exists(file):
        # parse file
        MolParserC = MolParser(file)
        compound_info = MolParserC.read_file()
        # compound
        compound = Compound(compound_info)
        # display 3d
        compound.view3d(display_legend=display_legend)
    else:
        raise Exception("file path is not valid.")


def td_by_inchi(inchi, display_legend=True):
    '''
    3d visualizer of a compound using its InChI identifier

    Parameters
    ----------
    inchi : str
        inchi code
    display_legend : bool
        display legend (default True)

    Returns
    -------
    None
        display 3d
    '''
    # check inchi
    if inchi is not None:
        # get cid
        cid = pcq.get_cid_by_inchi(inchi)
        # check
        if cid is not None:
            # get sdf
            sdf = pcq.get_structure_by_cid(cid)
            # check
            if sdf is not None:
                # display 3d
                # parse file
                MolParserC = MolParser(None)
                compound_info = MolParserC.read_file(
                    sourceContent={
                        'content': sdf,
                        'format': 'sdf'
                    })
                # compound
                compound = Compound(compound_info)
                # display 3d
                compound.view3d(display_legend=display_legend)
            else:
                raise Exception("sdf is not found!")
        else:
            raise Exception("cid is not found!")
    else:
        raise Exception("inchi is not valid.")


if __name__ == "__main__":
    main()
