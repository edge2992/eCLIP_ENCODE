import os
from typing import List, Tuple

import numpy as np
import pandas as pd
from dotenv import load_dotenv

from src.util.similarity_strategy.interface import ProteinSimilarityStrategy

load_dotenv()

PROJECT_PATH = os.environ["PROJECT_PATH"]


class MSA(ProteinSimilarityStrategy):
    def _protein_similarity(self) -> pd.DataFrame:
        similarities, labels = self._load()
        return pd.DataFrame(similarities, index=labels, columns=labels, dtype=float)

    def _load(self) -> Tuple[List[List[str]], List[str]]:
        similarities: List[List[str]] = []
        labels: List[str] = []

        if self.loadfile is None:  # type: ignore
            self.loadfile = os.environ["CLUSTALO_DISTMAT_PATH"]

        if not os.path.exists(self.loadfile):
            self._execute_msa(self.loadfile)

        with open(self.loadfile, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                if line.startswith("#"):
                    continue
                elif line == "\n":
                    continue
                else:
                    similarities.append(line.split()[1:])
                    labels.append(line.split()[0])
        return similarities, labels

    def _execute_msa(self, path: str):
        """Clustal Omegaを実行してファイルをpathに保存する"""
        raise NotImplementedError()

    def __repr__(self) -> str:
        return "MSA"

    @property
    def lower_better(self) -> bool:
        return False


class BlastP(ProteinSimilarityStrategy):
    FASTAFILE = "/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/reviewed.fasta"
    LOADFILE = "/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/blast/reviewed.tsv"
    DBFILE = "/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/blast/db/RbpHomoDB"

    def _protein_similarity(self) -> pd.DataFrame:
        similarities, labels = self._load()
        return pd.DataFrame(similarities, index=labels, columns=labels, dtype=float)

    def _load(self) -> Tuple[np.ndarray, List[str]]:
        similarities: np.ndarray
        labels: List[str] = []

        if self.loadfile is None:  # type: ignore
            self.loadfile = self.LOADFILE

        if not os.path.exists(self.loadfile):
            proc = self._execute_blastp(self.loadfile)
            proc.wait()
            if proc.returncode != 0:
                raise RuntimeError("blastp failed")

        df = pd.read_table(
            self.loadfile,
            header=None,
            names=["qseqid", "sseqid", "pindent", "evalue", "bitscore", "score"],
        ).drop_duplicates(subset=["qseqid", "sseqid"], keep="first")
        labels = list(df["qseqid"].unique())
        dd = {id_: i for i, id_ in enumerate(labels)}
        similarities = np.zeros((len(labels), len(labels)), dtype=float)
        for _, row in df.iterrows():
            similarities[dd[row["qseqid"]], dd[row["sseqid"]]] = row["bitscore"]

        return similarities, labels

    def _execute_blastp(self, path: str):
        """blastpを実行してファイルをpathに保存する"""
        print("blastpを実行しています...")
        import subprocess

        if not os.path.exists(self.DBFILE + ".phr"):
            proc = self._make_dB()
            proc.wait()
            if proc.returncode != 0:
                raise RuntimeError("makeblastdb failed")

        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))

        cmd_blastp = 'blastp -num_threads 6 -query {} -db {} -outfmt "6 qseqid sseqid pident evalue bitscore score" -out {} -evalue 1000'.format(
            self.FASTAFILE, self.DBFILE, path
        )
        print(cmd_blastp)
        return subprocess.Popen(cmd_blastp, shell=True)

    def _make_dB(self):
        """makeblastdbを実行してDBを作成する"""
        import subprocess

        print("makeblastdbを実行しています...")

        cmd_makeblastdb = (
            "makeblastdb -in {} -dbtype prot -title {} -parse_seqids".format(
                self.FASTAFILE, self.DBFILE
            )
        )
        print(cmd_makeblastdb)
        return subprocess.Popen(cmd_makeblastdb, shell=True)

    def __repr__(self) -> str:
        return "Blastp"

    @property
    def lower_better(self) -> bool:
        return False


class TAPE(ProteinSimilarityStrategy):
    def _protein_similarity(self) -> pd.DataFrame:
        from scipy.spatial.distance import cosine, pdist, squareform

        from src.util.tape import load_embedding

        if self.loadfile is None:  # type: ignore
            self.loadfile = "/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/reviewed_tape.npz"

        if not os.path.exists(self.loadfile):
            self._execute_tape(self.loadfile)

        embedding = load_embedding(how="avg", filepath=self.loadfile)

        return pd.DataFrame(
            squareform(pdist(np.array(list(embedding.values()), dtype=float), cosine)),
            index=list(embedding.keys()),
            columns=list(embedding.keys()),
        )

    def _execute_tape(self, path: str):
        """TAPEを実行してファイルをpathに保存する"""
        raise NotImplementedError()

    def __repr__(self) -> str:
        return "TAPE"

    @property
    def lower_better(self) -> bool:
        return True
