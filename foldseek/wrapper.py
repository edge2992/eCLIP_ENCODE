import os
from dotenv import load_dotenv

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]


class FoldSeek:
    def __init__(
        self,
        target_DB_path: str,
        query_DB_path: str,
        output_tsv_path: str,
        db_dir: str = os.path.join(PROJECT_PATH, "foldseek/db"),
    ):
        self.target_DB_path = target_DB_path
        self.query_DB_path = query_DB_path
        self.db_dir = db_dir
        self.output_tsv_path = output_tsv_path

        if not os.path.exists(self.db_dir):
            os.makedirs(self.db_dir)

    def exec(self, e_value: float = 0.001):
        current_dir = os.getcwd()
        os.chdir(self.db_dir)
        self._create_db(self.query_DB_path, "queryDB")
        self._create_db(self.target_DB_path, "targetDB")
        self._search("queryDB", "targetDB", "aln", "tmpFolder", e_value)
        self._aln2tmscore("queryDB", "targetDB", "aln", "aln_tmscore")
        self._createtsv("queryDB", "targetDB", "aln_tmscore", self.output_tsv_path)
        os.chdir(current_dir)

    def _create_db(self, db_name: str, db_path: str):
        """foldseek createdb rdb queryDB"""
        cmd = "foldseek createdb {} {}".format(db_name, db_path)
        return os.system(cmd)

    def _search(
        self,
        query_db: str,
        target_db: str,
        aln: str,
        tmp_folder: str,
        e_value: float = 0.001,
    ):
        """foldseek search queryDB targetDB aln tmpFolder -a"""
        cmd = "foldseek search {} {} {} {} -a -e {}".format(
            query_db, target_db, aln, tmp_folder, e_value
        )
        return os.system(cmd)

    def _aln2tmscore(self, query_db: str, target_db: str, aln: str, aln_tmscore: str):
        """foldseek aln2tmscore queryDB targetDB aln aln_tmscore"""
        cmd = "foldseek aln2tmscore {} {} {} {}".format(
            query_db, target_db, aln, aln_tmscore
        )
        return os.system(cmd)

    def _createtsv(
        self, query_db: str, target_db: str, aln_tmscore: str, aln_tmscore_tsv: str
    ):
        """foldseek createtsv queryDB targetDB aln_tmscore aln_tmscore.tsv"""
        cmd = "foldseek createtsv {} {} {} {}".format(
            query_db, target_db, aln_tmscore, aln_tmscore_tsv
        )
        return os.system(cmd)


if __name__ == "__main__":
    target_DB_path = os.path.join(PROJECT_PATH, "data/afdb/rdb")
    query_DB_path = os.path.join(PROJECT_PATH, "data/afdb/rdb")

    for e_value in [1e-10, 1e-5, 1e-3]:
        output_tsv_path = os.path.join(
            PROJECT_PATH, "data/afdb/aln_tmscore_{}.tsv".format(e_value)
        )
        fs = FoldSeek(target_DB_path, query_DB_path, output_tsv_path)
        fs.exec()
