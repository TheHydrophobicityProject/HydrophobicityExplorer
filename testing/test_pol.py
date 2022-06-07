import pytest
from ..LogP.MakePolymer import main

def test_pol():
    assert main() == 0
#doesn't really work now.