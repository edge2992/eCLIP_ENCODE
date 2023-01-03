from src.util.metrics.condition import ConditionAnd, ConditionEq
from src.eclip.encodecondition import GeneNumberGteCondition, GeneNumberLtCondition

HEPG2_COMPARESET = {
    "Protein_GTE_1000": ConditionAnd(
        [ConditionEq("Biosample name", "HepG2"), GeneNumberGteCondition(1000)]
    ),
    "Protein_LT_1000": ConditionAnd(
        [ConditionEq("Biosample name", "HepG2"), GeneNumberLtCondition(1000)]
    ),
}

K562_COMPARESET = {
    "Protein_GTE_1000": ConditionAnd(
        [ConditionEq("Biosample name", "K562"), GeneNumberGteCondition(1000)]
    ),
    "Protein_LT_1000": ConditionAnd(
        [ConditionEq("Biosample name", "K562"), GeneNumberLtCondition(1000)]
    ),
}

COMPARESET = {"HepG2": HEPG2_COMPARESET, "K562": K562_COMPARESET}
