from kc.data_structures.logicalterms import *
from kc.data_structures.substitutions import *
from kc.data_structures.literals import *

def test_predicate_equality():
    smokes1 = Predicate('smokes', 4)
    smokes2 = Predicate('smokes', 4)
    smokes_smaller_arity = Predicate('smokes', 3)
    friends = Predicate('friends', 4)

    assert(smokes1 == smokes2)
    assert(smokes1 != smokes_smaller_arity)
    assert(smokes1 != friends)




