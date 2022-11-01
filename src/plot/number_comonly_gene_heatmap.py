# 共通して相互作用する遺伝子数をヒートマップで表示する
# replicate1, 2を使用
# 遺伝子数が500件以上のアッセイで行う

import os
from typing import Dict
import pandas as pd
import sys
from joblib import Parallel, delayed
from dotenv import load_dotenv
import matplotlib.pyplot as plt

print(os.getcwd())


load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]
sys.path.append(PROJECT_PATH)
from src.util.bedfile import load_report, read_annotated_bed
from src.util.bed_format_strategy import FormatStrategy
from src.util.get_bed_path import get_formatted_file_path
