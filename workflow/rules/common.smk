from pathlib import Path
import pandas as pd


# Define input for the rules
data = []
with open(config["sumstats_path"], "r") as fp:
    lines = fp.readlines()

for line in lines:
    p = Path(line.strip())
    seqid = ".".join(p.stem.split(".")[:3]) #seqid = p.stem.split("_")[3]
    data.append((seqid, str(p)))

analytes = (
    pd.DataFrame.from_records(data, columns=["seqid", "sumstat_path"])
    .set_index("seqid", drop=False)
    .sort_index()
)


def get_sumstats(wildcards):
    return analytes.loc[wildcards.seqid, "sumstat_path"]


# define the functions generating files' path
def ws_path(file_path):
    return str(Path(config.get("workspace_path"), file_path))

