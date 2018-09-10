import fileinput
import gzip
from glob import glob
from os.path import join, exists
from os import mkdir
import pandas as pd

outpath = "raw"

samples = pd.read_csv("barcodes.csv")


def analyze_group(df):
    """Summarize a single folder."""
    for i in df.Barcode:
        current = df[df.Barcode == i].iloc[0]
        out = join(outpath, current.loc["ISB ID"] + ".fastq.gz")
        print("Writing to %s..." % out)
        path = join(current["Group"],
                    "barcode" + str(current.Barcode).zfill(2),
                    "uploaded")
        files = glob(join(path, "*.fastq"))
        with gzip.open(out, "wt") as gzfile:
            for f in files:
                with open(f) as input:
                    gzfile.write(input.read())


if not exists(outpath):
    mkdir(outpath)

for grp in samples.Group.unique():
    analyze_group(samples[samples.Group == grp])
