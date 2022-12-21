# %%
from src.eclip.encodecondition import ProteinKeywordCondition
from src.eclip.sampleset import SampleSetECLIP
from src.eclip.uniprot.keyword import Keyword
from src.util.metrics.condition import ConditionAnd, ConditionEq

keywords = Keyword().keywords()
SKIP_N = 5

# %%
for key in keywords:
    condition = ConditionAnd(
        [ConditionEq("Biosample name", "HepG2"), ProteinKeywordCondition(key)]
    )
    sampleset = SampleSetECLIP(condition)
    if sampleset.report.shape[0] < SKIP_N:
        continue
    print(key, sampleset.report.shape)


# %%
