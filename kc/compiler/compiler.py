"""This file contains a Compiler class that implements the main
Compile algorithm. This algorithm calls various additional algorithms (rules)
which operate on a CNF, as defined in the data_structures package."""

from kc.data_structures import *

from kc.compiler import KCRule
from kc.compiler import LeafConstruction, UnitPropagation, VacuousConjunction
from kc.compiler import Independence, ShannonDecomposition, ShatteredCompilation
from kc.compiler import IndependentSingleGroundings, IndependentPairedGroundings, AtomCounting
from kc.compiler import Ground


from typing import Dict, Optional, Tuple, Any, Type

# TODO: delta -> theory

class Compiler:
    """A knowledge compilation compiler that takes CNFs and produces 
    FO-da-DNNF (or FO-sda-DNNF) circuits, following VdB's PhD."""

    def __init__(self) -> None:
        """We initialise the compiler with an empty cache of previously-seen theories
        and a list of the compilation rules that can be applied"""
        self._cache: Dict['CNF', 'NNFNode'] = dict()
        self.rules: Tuple[Type['KCRule'], ...] = (LeafConstruction,
                                       UnitPropagation,
                                       VacuousConjunction,
                                       Independence,
                                       ShannonDecomposition,
                                       ShatteredCompilation,
                                       IndependentSingleGroundings,
                                       IndependentPairedGroundings,
                                       AtomCounting,
                                       Ground)

    def compile(self, delta: 'CNF') -> 'NNFNode':  # TODO: implement NNFNode
        """This function follows closely the algorithm described in the PhD and 
        the one used in Forclift"""
        if self.cache_contains(delta):
            return self.cache_get(delta)

        nnf: Optional['NNFNode'] = None

        applicable_rule: Optional[Type['KCRule']]
        stored_data: Optional[Any]
        applicable_rule, stored_data = self.find_rule(delta)

        if applicable_rule is None:
            raise ValueError("Compilation failed - no rule found for {delta}")
            nnf = None
        else:
            nnf = self.apply_rule(delta, applicable_rule, stored_data) 
        self.cache_set(delta, nnf)
        return nnf

    def cache_contains(self, cnf: 'CNF') -> bool:
        """Check if the cache contains a given CNF"""
        return cnf in self._cache.keys()

    # TODO: get_cache instead of cache_get
    def cache_get(self, cnf: 'CNF') -> 'NNFNode':
        """Return the circuit (represented by a single Node) stored in the cache for a
        given CNF.
        NOTE: Assumes that we've already checked for the cnf being present,
        otherwise throws an error if it's not there"""
        return self._cache[cnf]

    def cache_set(self, cnf: Optional['CNF'], nnf: 'NNFNode') -> None:
        """Set the cache value for a given cnf if that cnf compiled AND we haven't
        seen it before (because they should always compile to the same thing)"""
        if not (cnf is None or self.cache_contains(cnf)):
            self._cache[cnf] = nnf

    def find_rule(self, delta: 'CNF') -> Tuple[Optional[Type['KCRule']], Optional[Any]]:
        """Check each compilation rule to see if its preconditions are met
        for the given cnf.
        Return the rule if one is found, plus optional info already computed.
        Otherwise return None and None
        Loops over:
                LeafConstruction,
                UnitPropagation,
                VacuousConjunction,
                Independence,
                ShannonDecomposition,
                ShatteredCompilation,
                IndependentSingleGroundings,
                IndependentPairedGroundings,
                AtomCounting,
                Ground
       """
        for rule in self.rules: 
            applicable: bool
            stored_data: Optional[Any]
            applicable, stored_data = rule.is_applicable(delta)
            if applicable:
                return rule, stored_data
        return None, None

    def apply_rule(self, delta: 'CNF', rule: Type['KCRule'], stored_data: Optional[Any]) -> 'NNFNode':
        """Apply a given compilation rule to a cnf and return the constructed NNF
        NOTE: we also accept precomputed data from find_rule and pass it on if it's not None"""
        # we pass this compiler in so it can be called recursively
        # TODO self first
        nnf = rule.apply(delta, stored_data, self)
        return nnf
