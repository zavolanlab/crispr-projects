"""
Public resources used by lincdisease: URLs, cached download and readers. Files are downloaded once into the
data directory (--data-dir; never into the repository) and re-used.
"""
import gzip
import os
import shutil
import urllib.request

import pandas as pd

GIAB = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3/GRCh38@all/OtherDifficult/"
URLS = {
    "segdups": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz",
    # GRCh38 regions falsely duplicated relative to T2T-CHM13 (Aganezov et al. 2022), both copies
    "false_dup_correct": GIAB + "GRCh38_false_duplications_correct_copy.bed.gz",
    "false_dup_incorrect": GIAB + "GRCh38_false_duplications_incorrect_copy.bed.gz",
    # snoRNABase / miRBase small RNAs (includes snoRNAs missing from GENCODE, e.g. those hosted by GAS5)
    "wgrna": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgRna.txt.gz",
    "gencode_gtf": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz",
    "hgnc": "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
    "gwas_catalog": "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
                    "gwas-catalog-associations_ontology-annotated-full.zip",
    # TCGA Pan-Cancer Clinical Data Resource (Liu et al. 2018), survival endpoints
    "tcga_cdr": "https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81",
    # all GTEx samples of recount3 with their tissue (project)
    "gtex_projects": "https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/gtex/metadata/"
                     "gtex.recount_project.MD.gz",
}

# local file names for URLs without one
FILENAMES = {"tcga_cdr": "TCGA-CDR-SupplementalTableS1.xlsx"}

# columns of UCSC genomicSuperDups.txt (bin first)
SEGDUP_COLUMNS = ["bin", "chrom", "chromStart", "chromEnd", "name", "score", "strand", "otherChrom",
                  "otherStart", "otherEnd", "otherSize", "uid", "posBasesHit", "testResult", "verdict",
                  "chits", "ccov", "alignfile", "alignL", "indelN", "indelS", "alignB", "matchB",
                  "mismatchB", "transitionsB", "transversionsB", "fracMatch", "fracMatchIndel", "jcK", "k2K"]


def fetch(key, data_dir, url=None):
    """Download URLS[key] (or `url`) into data_dir unless present; return the local path."""
    url = url or URLS[key]
    os.makedirs(data_dir, exist_ok=True)
    path = os.path.join(data_dir, FILENAMES.get(key) or os.path.basename(url.split("?")[0]))
    if not os.path.exists(path):
        tmp = path + ".part"
        with urllib.request.urlopen(url, timeout=300) as r, open(tmp, "wb") as fh:
            shutil.copyfileobj(r, fh)
        os.replace(tmp, path)
    return path


def read_segdups(data_dir):
    """UCSC genomicSuperDups as chrom/start0/end1/fracMatch."""
    df = pd.read_csv(fetch("segdups", data_dir), sep="\t", header=None, names=SEGDUP_COLUMNS,
                     usecols=["chrom", "chromStart", "chromEnd", "fracMatch"])
    return df.rename(columns={"chromStart": "start0", "chromEnd": "end1"})


def read_bed3(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        rows = [ln.split("\t")[:3] for ln in fh if ln.strip() and not ln.startswith(("#", "track", "browser"))]
    df = pd.DataFrame(rows, columns=["chrom", "start0", "end1"])
    return df.astype({"start0": int, "end1": int})


def read_false_dups(data_dir):
    """Union of the 'correct copy' and 'incorrect copy' BEDs of GRCh38 false duplications."""
    return pd.concat([read_bed3(fetch(k, data_dir)) for k in ("false_dup_correct", "false_dup_incorrect")],
                     ignore_index=True)


def read_wgrna(data_dir):
    """UCSC wgRna (bin, chrom, start, end, name, score, strand, thickStart, thickEnd, type)."""
    df = pd.read_csv(fetch("wgrna", data_dir), sep="\t", header=None, usecols=[1, 2, 3, 4, 6, 9],
                     names=["chrom", "start0", "end1", "name", "strand", "type"])
    return df
