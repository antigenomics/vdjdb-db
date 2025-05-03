from collections import defaultdict
import re
import pandas as pd


class VdjdbScoreFactory:
    def __init__(self, master_table: pd.DataFrame):
        self.score_map = {}

        for _, row in master_table.iterrows():
            try:
                sign = self.get_signature(row)

                # Compute score
                if pd.notnull(row["meta.structure.id"]):
                    score = 3  # we have structure
                else:
                    freq_str = row["method.frequency"]
                    freq = self.get_frequency(freq_str)
                    count = self.get_number_of_cells(freq_str)

                    assert freq <= 1.0, "Frequency exceeds 1.0"

                    # Sequencing score
                    seq_score = 1
                    single_cell = row["method.singlecell"].strip().lower() if pd.notnull(row["method.singlecell"]) else row["method.singlecell"]

                    if pd.notnull(row["method.singlecell"]) and single_cell != "no":
                        seq_score = 3
                    else:
                        sequencing_method = row["method.sequencing"].strip().lower() if pd.notnull(row["method.sequencing"]) else row["method.sequencing"]
                        if sequencing_method == "sanger":
                            seq_score = 3 if count >= 2 else 2
                        elif sequencing_method == "amplicon-seq":
                            seq_score = 3 if freq >= 0.01 else 1

                    # Moderate confidence regarding specificity
                    identify_method = row["method.identification"].lower().strip() if pd.notnull(row["method.identification"]) else row["method.identification"]
                    spec_score1 = 0

                    if self.culture_based_identification(identify_method):
                        spec_score1 = 1 if freq >= 0.5 else 0
                    else:
                        if self.sort_based_method(identify_method):
                            spec_score1 = 1 if freq >= 0.05 else 0
                        elif self.stimulation_based_method(identify_method):
                            spec_score1 = 1 if freq >= 0.25 else 0

                    # High confidence regarding specificity
                    verify_method = row["method.verification"].lower() if pd.notnull(row["method.verification"]) else row["method.verification"]
                    spec_score2 = 0

                    if self.direct_verification(verify_method):
                        spec_score2 = 3
                    elif self.stimulation_based_method(verify_method):
                        spec_score2 = 2
                        seq_score = 3  # tcr was cloned
                    elif self.sort_based_method(verify_method):
                        spec_score2 = 1
                        seq_score = 3  # tcr was cloned

                    score = min(seq_score, spec_score1 + spec_score2)

                self.score_map[sign] = max(self.score_map.get(sign, 0), score)

            except Exception as e:
                raise RuntimeError(
                    f"Error: {e}\nin {row['reference.id']} {row['cdr3.alpha']} {row['cdr3.beta']}"
                )

    @staticmethod
    def get_number_of_cells(freq):
        if  pd.isnull(freq):
            return 0
        return int(re.split(r'/+', freq)[0]) if "/" in freq else 0

    @staticmethod
    def get_frequency(freq):
        if pd.isnull(freq):
            return 0.0
        if "/" in freq:
            x = [float(x) for x in re.split(r'/+', freq)]
            return x[0] / x[1]
        elif freq.endswith("%") and len(freq) > 1:
            return float(freq[:-1]) / 100.0
        return float(freq)

    @staticmethod
    def sort_based_method(method):
        return any(x in method for x in ["sort", "beads", "separation", "stain"]) if pd.notnull(method) else False

    @staticmethod
    def stimulation_based_method(method):
        return "targets" in method if pd.notnull(method) else False

    @staticmethod
    def culture_based_identification(method):
        return any(x in method for x in ["culture", "cloning"]) if pd.notnull(method) else False

    @staticmethod
    def direct_verification(method):
        return "direct" in method if pd.notnull(method) else False

    SIGNATURE_COLS = [
        "cdr3.alpha",
        "v.alpha",
        "j.alpha",
        "cdr3.beta",
        "v.beta",
        "j.beta",
        "species",
        "mhc.a",
        "mhc.b",
        "mhc.class",
        "antigen.epitope"
    ]

    def get_score(self, row):
        return self.score_map.get(self.get_signature(row), 0)

    def get_signature(self, row):
        return "\t".join(str(row[col]) for col in self.SIGNATURE_COLS)