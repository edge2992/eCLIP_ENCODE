from dotenv import load_dotenv
import pandas as pd
from typing import Dict

from src.eclip.sampleset import SampleSetECLIP
from src.util.metrics import Metrics
from src.util.similarity_strategy import (
    INTERACTION_SIMILARITY_STRATEGIES,
    PROTEIN_SIMILARITY_SYMMETRIC_STRATEGIES,
)
from thesis.utils.compareset import COMPARESET

load_dotenv()

report_data: Dict[str, Dict[str, SampleSetECLIP]] = {
    biosample: {
        k: SampleSetECLIP(condition) for k, condition in COMPARESET[biosample].items()
    }
    for biosample in COMPARESET
}


STRATEGIES = PROTEIN_SIMILARITY_SYMMETRIC_STRATEGIES + INTERACTION_SIMILARITY_STRATEGIES

data: Dict[str, Dict[str, pd.DataFrame]] = {
    biosample: {
        key: Metrics(SampleSetECLIP(condition).report)(
            STRATEGIES, add_description=False
        )
        for key, condition in COMPARESET[biosample].items()
    }
    for biosample in COMPARESET.keys()
}  # type: ignore
