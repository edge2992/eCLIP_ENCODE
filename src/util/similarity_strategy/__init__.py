from .foldseek import FoldSeekTMScore
from .interaction import Cosine, Dice, Jaccard, Lift, Simpson
from .interface import (
    InteractionSimilarityStrategy,
    ProteinSimilarityStrategy,
    SimilarityStrategy,
    Default,
)
from .peak import PeakStrategy
from .protein import TAPE, BlastP
from .stringdb import DirectStringScore, STRINGDB_SCORE_METRICS
from .uniprot_keyword import KeywordCosine, KeywordAA

INTERACTION_SIMILARITY_STRATEGIES = [
    Jaccard(),
    Simpson(),
    PeakStrategy(metrics="jaccard"),
    PeakStrategy(metrics="union-intersection"),
    PeakStrategy(metrics="intersection"),
    PeakStrategy(metrics="n_intersections"),
]

PROTEIN_SIMILARITY_SYMMETRIC_STRATEGIES = [
    BlastP(symmetric_method="avg"),
    TAPE(symmetric_method="avg"),
    FoldSeekTMScore(symmetric_method="avg"),
    KeywordAA(),
    # DirectStringScore(metrics="score", symmetric_method="avg"),
]

PROTEIN_SIMILARITY_STRATEGIES = [
    BlastP(),
    TAPE(),
    FoldSeekTMScore(),
    KeywordAA(),
    # DirectStringScore(metrics="score"),
]

STRINGDB_SCORE_STRATEGIES = [
    DirectStringScore(metrics=met, fillna=None) for met in STRINGDB_SCORE_METRICS
]

__all__ = [
    "Default",
    "SimilarityStrategy",
    "ProteinSimilarityStrategy",
    "InteractionSimilarityStrategy",
    "BlastP",
    "TAPE",
    "FoldSeekTMScore",
    "PeakStrategy",
    "KeywordCosine",
    "KeywordAA",
    "DirectStringScore",
    "Jaccard",
    "Simpson",
    "Dice",
    "Lift",
    "Cosine",
    "INTERACTION_SIMILARITY_STRATEGIES",
    "PROTEIN_SIMILARITY_STRATEGIES",
    "PROTEIN_SIMILARITY_SYMMETRIC_STRATEGIES",
    "STRINGDB_SCORE_STRATEGIES",
]
