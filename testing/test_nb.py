import hydrophobicity_explorer as hx
from hydrophobicity_explorer import nb
import pytest, os, json

def test_nb_write():
    nb.main()
    with open(nb.name, "r") as S:
        read_nb = json.load(S)
    
    os.remove(nb.name)
    assert read_nb == nb.notebook