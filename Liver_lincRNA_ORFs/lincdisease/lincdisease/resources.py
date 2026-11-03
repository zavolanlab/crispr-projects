"""
Public resources used by lincdisease: URLs, cached download and readers. Files are downloaded once into the
data directory (--data-dir; never into the repository) and re-used. Every download is recorded in
<directory>/downloads.tsv (date, file, URL, size, md5), as some URLs (GWAS Catalog `latest`) change over time.
"""
import datetime
import gzip
import hashlib
import os
import shutil
import urllib.request

import pandas as pd

GIAB = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3/GRCh38@all/OtherDifficult/"
XENA_GDC = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/"
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
    # TCGA-LIHC masked copy-number segments (GDC, GRCh38) via the UCSC Xena GDC hub
    "tcga_lihc_cn": XENA_GDC + "TCGA-LIHC.masked_cnv_DNAcopy.tsv.gz",
    # TCGA PanCanAtlas ABSOLUTE purity / ploidy (Taylor et al. 2018)
    "tcga_purity": "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5",
    # all GTEx samples of recount3 with their tissue (project)
    "hg19_to_hg38_chain": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
    "gtex_projects": "https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/gtex/metadata/"
                     "gtex.recount_project.MD.gz",
}

# local file names for URLs without one
FILENAMES = {"tcga_cdr": "TCGA-CDR-SupplementalTableS1.xlsx",
             "tcga_purity": "TCGA_mastercalls.abs_tables_JSedit.fixed.txt",
             "eclip_search": "files.json", "rbp_localisation": "vannostrand2020_S3.xlsx"}

# files placed in the data directory by hand (no stable direct URL, or extracted from a larger archive):
# path relative to the data directory -> where to get it
GENCODE42 = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/"
UNIPROT_HUMAN = ("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/"
                 "Eukaryota/UP000005640/")
PRENSNER = "Prensner et al. 2021 (Nat Biotechnol), PMC supplement NIHMS1698134"
MANUAL = {
    # genome for the STAR index of the Govaere 2020 alignments (workflow/rules/govaere.smk)
    "homology/GRCh38.primary_assembly.genome.fa.gz": GENCODE42 + "GRCh38.primary_assembly.genome.fa.gz",
    "homology/gencode.v42.pc_translations.fa.gz": GENCODE42 + "gencode.v42.pc_translations.fa.gz",
    "homology/UP000005640_9606.fasta.gz": UNIPROT_HUMAN + "UP000005640_9606.fasta.gz",
    "homology/UP000005640_9606_additional.fasta.gz": UNIPROT_HUMAN + "UP000005640_9606_additional.fasta.gz",
    "riboseq_orfs/phase1/Ribo-seq_ORFs.peptideins.bed": "https://www.gencodegenes.org/pages/riboseq_orfs/ (phase 1)",
    "riboseq_orfs/phase2/Ribo-seq_ORFs.comprehensive.bed": "https://www.gencodegenes.org/pages/riboseq_orfs/ (phase 2)",
    "Pfam-A.hmm": "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz (gunzip)",
    "crispr_screens/NIHMS1698134-supplement-Supplementary_Tables_1-17.xlsx": PRENSNER,
    "crispr_screens/NIHMS1698134-supplement-Supplementary_Tables_21-32.xlsx": PRENSNER,
    **{f"crispr_screens/aay0262_table_{n}.xlsx": "Chen et al. 2020 (Science, doi 10.1126/science.aay0262), "
                                                 "supplementary tables" for n in ("s1", "s3", "s4")},
    "gtex_eqtl/GTEx_Analysis_v8_eQTL/Liver.v8.egenes.txt.gz":
        "curl -s <eqtl.GTEX_TAR> | tar -x -C <data_dir>/gtex_eqtl --wildcards '*Liver.v8.egenes*' '*Liver.v8.signif*'",
    "gtex_eqtl/GTEx_Analysis_v8_eQTL/Liver.v8.signif_variant_gene_pairs.txt.gz": "as Liver.v8.egenes.txt.gz",
    "mhcflurry/2.3.0/models_class1_presentation":
        "MHCFLURRY_DATA_DIR=<data_dir>/mhcflurry <tool_env>/bin/mhcflurry-downloads fetch models_class1_presentation",
}


def missing_manual(data_dir):
    """(path, source) of the hand-placed files that are not in data_dir."""
    return [(os.path.join(data_dir, p), src) for p, src in MANUAL.items()
            if not os.path.exists(os.path.join(data_dir, p))]


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
        _record(path, url)
    return path


def _record(path, url):
    md5 = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            md5.update(block)
    log = os.path.join(os.path.dirname(path), "downloads.tsv")
    new = not os.path.exists(log)
    with open(log, "a") as fh:
        if new:
            fh.write("date\tfile\turl\tbytes\tmd5\n")
        fh.write(f"{datetime.date.today()}\t{os.path.basename(path)}\t{url}\t{os.path.getsize(path)}\t"
                 f"{md5.hexdigest()}\n")


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
