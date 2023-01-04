from dotenv import load_dotenv
from typing import Dict

from src.eclip.sampleset import SampleSetECLIP
from thesis.utils.compareset import COMPARESET

load_dotenv()

COMPARE_REPORT_SET: Dict[str, Dict[str, SampleSetECLIP]] = {
    biosample: {
        k: SampleSetECLIP(condition) for k, condition in COMPARESET[biosample].items()
    }
    for biosample in COMPARESET
}
