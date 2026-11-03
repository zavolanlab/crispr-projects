// Build: npm install docx@9 && node methods.js lincdisease_Methods.docx
// Methods document for lincdisease (disease relevance and mechanism annotation of lincRNA ORFs).
const fs = require("fs");
const d = require("docx");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType, Table, TableRow, TableCell, WidthType,
  ShadingType, BorderStyle, LevelFormat, Footer, PageNumber, TableOfContents,
  Math: OMath, MathRun, MathFraction, MathSubScript, MathSuperScript, MathSubSuperScript, MathSum, MathRadical,
  MathRoundBrackets, MathCurlyBrackets, MathSquareBrackets,
} = d;

const out = process.argv[2] || "lincdisease_Methods.docx";
const FONT = "Calibri";
const PAGE_W = 11906, MARGIN = 1300, TEXT_W = PAGE_W - 2 * MARGIN;   // A4

// ---------------------------------------------------------------- text helpers
// p("plain *italic* **bold** `code`")
function runs(s, base = {}) {
  const out = [];
  const re = /(\*\*[^*]+\*\*|\*[^*]+\*|`[^`]+`|_\{[^}]+\})/g;
  let last = 0, m;
  while ((m = re.exec(s))) {
    if (m.index > last) out.push(new TextRun({ text: s.slice(last, m.index), ...base }));
    const t = m[0];
    if (t.startsWith("**")) out.push(new TextRun({ text: t.slice(2, -2), bold: true, ...base }));
    else if (t.startsWith("_{")) out.push(new TextRun({ text: t.slice(2, -1), subScript: true, ...base }));
    else if (t.startsWith("`")) out.push(new TextRun({ text: t.slice(1, -1), font: "Consolas", size: 19, ...base }));
    else out.push(new TextRun({ text: t.slice(1, -1), italics: true, ...base }));
    last = m.index + t.length;
  }
  if (last < s.length) out.push(new TextRun({ text: s.slice(last), ...base }));
  return out;
}
const p = (s, opt = {}) => new Paragraph({ children: runs(s), spacing: { after: 120, line: 276 }, ...opt });
const h1 = (s) => new Paragraph({ heading: HeadingLevel.HEADING_1, children: [new TextRun(s)] });
const h2 = (s) => new Paragraph({ heading: HeadingLevel.HEADING_2, children: [new TextRun(s)] });
const h3 = (s) => new Paragraph({ heading: HeadingLevel.HEADING_3, children: [new TextRun(s)] });
const bullet = (s) => new Paragraph({ numbering: { reference: "bullets", level: 0 }, children: runs(s),
                                      spacing: { after: 60, line: 276 } });

// ---------------------------------------------------------------- math helpers
const R = (t) => new MathRun(t);
const F = (n, den) => new MathFraction({ numerator: [].concat(n), denominator: [].concat(den) });
const Sub = (b, s) => new MathSubScript({ children: [].concat(b), subScript: [].concat(s) });
const Sup = (b, s) => new MathSuperScript({ children: [].concat(b), superScript: [].concat(s) });
const SubSup = (b, s, u) => new MathSubSuperScript({ children: [].concat(b), subScript: [].concat(s),
                                                    superScript: [].concat(u) });
const Sum = (lo, hi, body) => new MathSum({ children: [].concat(body), subScript: lo ? [].concat(lo) : undefined,
                                            superScript: hi ? [].concat(hi) : undefined });
const Sqrt = (x) => new MathRadical({ children: [].concat(x) });
const Par = (x) => new MathRoundBrackets({ children: [].concat(x) });
const Sq = (x) => new MathSquareBrackets({ children: [].concat(x) });
const Cur = (x) => new MathCurlyBrackets({ children: [].concat(x) });
let eqNo = 0;
function eq(...parts) {
  eqNo += 1;
  return new Paragraph({
    alignment: AlignmentType.CENTER, spacing: { before: 80, after: 140 },
    children: [new OMath({ children: parts }), new TextRun({ text: `\t(${eqNo})` })],
    tabStops: [{ type: d.TabStopType.RIGHT, position: TEXT_W }],
  });
}

// ---------------------------------------------------------------- tables
const border = { style: BorderStyle.SINGLE, size: 4, color: "BFC5CC" };
const borders = { top: border, bottom: border, left: border, right: border };
function table(header, rows, widths) {
  const total = widths.reduce((a, b) => a + b, 0);
  const w = widths.map((x) => Math.round((x / total) * TEXT_W));
  w[w.length - 1] += TEXT_W - w.reduce((a, b) => a + b, 0);
  const cell = (s, i, head) => new TableCell({
    borders, width: { size: w[i], type: WidthType.DXA },
    shading: head ? { fill: "E8EDF2", type: ShadingType.CLEAR, color: "auto" } : undefined,
    margins: { top: 60, bottom: 60, left: 100, right: 100 },
    children: [new Paragraph({ children: runs(String(s), head ? { bold: true, size: 19 } : { size: 19 }) })],
  });
  return new Table({
    width: { size: TEXT_W, type: WidthType.DXA }, columnWidths: w,
    rows: [new TableRow({ tableHeader: true, children: header.map((s, i) => cell(s, i, true)) }),
           ...rows.map((r) => new TableRow({ children: r.map((s, i) => cell(s, i, false)) }))],
  });
}
const gap = () => new Paragraph({ children: [], spacing: { after: 60 } });

// ---------------------------------------------------------------- content
const C = [];
C.push(new Paragraph({ heading: HeadingLevel.TITLE, children: [new TextRun("lincdisease: Methods")] }));
C.push(p("Disease relevance and functional mechanism of translated lincRNA ORFs in the liver. This document derives " +
  "every metric reported in the lincdisease tables (the D1-D6 disease components, the M1-M6 mechanism classes, " +
  "and the summary scores) as computed by the `lincdisease` package of the `crispr-projects` repository. The " +
  "specification kept with the code (`LOGIC.md`) is the reference; this document gives the formulas."));
C.push(new TableOfContents("Contents", { hyperlink: true, headingStyleRange: "1-2" }));

// 1
C.push(h1("1. Inputs and background"));
C.push(p("**Input.** ORFs called by the `lincorfz` pipeline on HepG2 and Huh7 Ribo-seq, RNA-seq and polysome " +
  "profiling (run v2). The candidate set is the curated table of 46 ORFs on 33 lincRNA loci. Every measure is " +
  "computed for all 11,910 ORFs on 2,057 expressed lincRNA loci of the same run, which form the **background**: " +
  "each statistic is converted into a percentile among the background, so 0.5 is the value of a typical " +
  "expressed lincRNA and 1 the most extreme."));
C.push(p("**Two levels.** Disease evidence attaches to the lincRNA **locus** (expression data cannot " +
  "separate the ORFs of one transcript); mechanism evidence is computed per **ORF**, with locus-level features " +
  "copied to every ORF of the locus. Locus and ORF scores are reported side by side and are never multiplied."));
C.push(p("**Controls.** 25 annotated genes (for example HULC, H19, COL1A1, PNPLA3, MALAT1, GAS5) are processed by " +
  "the same code as extra features (`ctrl:<gene>`, measured on the union of their GENCODE exons) and scored " +
  "against the lincRNA background without being part of it (Section 9)."));

// 2
C.push(h1("2. Scoring framework"));
C.push(h2("2.1 Percentiles against the background"));
C.push(p("For a statistic with value *x* and the non-missing background values *B* = {*b*_{1}, …, *b*_{N}}, the " +
  "mid-rank empirical percentile is"));
C.push(eq(Sub(R("p"), R("B")), R("(x) = "), F([R("#{b ∈ B : b < x} + "), F(R("1"), R("2")), R(" #{b ∈ B : b = x}")], R("N"))));
C.push(p("so ties count one half and *p* lies in [0, 1]. Missing values stay missing. **Signed statistics** " +
  "(D1 HCC, D1 chronic, D2, D3) are ranked on their magnitude, *p* = *p*_{|B|}(|*x*|), and the sign " +
  "*s* = sgn(*x*) is kept separately to give a direction (Section 2.3): a locus strongly down in tumours is as " +
  "disease-relevant as one strongly up."));
C.push(p("**Centring.** Across the background most loci are somewhat higher in tumours than in their adjacent " +
  "normals, somewhat lower in chronic disease and somewhat higher with worse survival. A locus that only follows " +
  "such a shift is not specifically disease-associated, so D1 HCC, D1 chronic and D2 are ranked, and give their " +
  "direction, as *x* − median(*B*), the background median subtracted. The reported statistics (`stat_` columns) " +
  "are uncentred."));
C.push(p("Where the background is ORFs (mechanism classes, actionability), background ORF *i* carries the weight " +
  "*ω*_{i} = 1 / (number of ORFs of its locus), and the counts in Equation 1 become sums of weights over the total " +
  "weight. Every locus thus counts once: locus-level features are copied to all ORFs of a locus, and loci with many " +
  "ORFs would otherwise dominate the background."));

C.push(h2("2.2 Composite scores and shrinkage for missing data"));
C.push(p("A composite of components *k* with fixed non-negative weights *w*_{k} and percentiles *p*_{k} is the " +
  "weighted mean in which a missing component contributes the neutral value 0.5:"));
C.push(eq(R("S = "), F(Sum(R("k"), null, [Sub(R("w"), R("k")), Sub(R("q"), R("k"))]), Sum(R("k"), null, Sub(R("w"), R("k")))),
  R(",    "), Sub(R("q"), R("k")), R(" = "), Cur([Sub(R("p"), R("k")), R(" if observed;  0.5 if missing")])));
C.push(p("Its **coverage** is the fraction of the weight with data, *C* = Σ_{obs} *w*_{k} / Σ *w*_{k}. Writing " +
  "*S*_{obs} for the weighted mean over the observed components only, the composite is"));
C.push(eq(R("S = C · "), Sub(R("S"), R("obs")), R(" + (1 − C) · 0.5")));
C.push(p("that is, the observed evidence shrunk toward 0.5 in proportion to the missing weight. A locus without " +
  "data therefore never scores low by default, and a score near 0.5 with low coverage means *unknown*, not " +
  "*unremarkable*; coverage is reported with every score."));

C.push(h2("2.3 Directions"));
C.push(p("A direction is read from the signs of the **informative** signed components, those whose own " +
  "percentile is ≥ 0.75 (centred as in Section 2.1 where that applies): all positive → up (oncogene-like), all negative → down (suppressor-like), both → inconsistent, " +
  "none informative → none. Signs are defined so that > 0 means higher expression with disease (D1), worse " +
  "survival with higher expression (D2), or expression that follows a frequent copy gain (D3)."));

C.push(h2("2.4 Multiple testing"));
C.push(p("p-values are adjusted by Benjamini-Hochberg over all tested features (loci and controls) of a test; " +
  "FDR values are reported, while scores use the effect sizes."));

// 3
C.push(h1("3. Quantification in public RNA-seq"));
C.push(p("Many loci are de novo assemblies without a GENCODE gene, so gene-level public tables cannot be used. " +
  "Coverage is read from recount3 base-level coverage bigWigs over the **quantifiable positions** *Q*_{f} of " +
  "each feature *f*: the union of the exons of its isoforms minus the exons of every other GENCODE gene on " +
  "either strand (public RNA-seq is unstranded). Features with |*Q*_{f}| < 100 are flagged " +
  "`tcga_quantifiable = False` and receive no expression-based component. For sample *j* with total base " +
  "coverage *T*_{j} (the summed coverage of all bases, recount3 `bc_auc.all_reads_all_bases`):"));
C.push(eq(Sub(R("ncov"), R("fj")), R(" = "), F(Sum([R("i ∈ "), Sub(R("Q"), R("f"))], null, Sub(R("cov"), R("ij"))), [R("|"), Sub(R("Q"), R("f")), R("| · "), Sub(R("T"), R("j"))]), R(" · "), Sup(R("10"), R("9")),
  R(",    "), Sub(R("e"), R("fj")), R(" = "), Sub(R("log"), R("2")), Par([F(Sub(R("ncov"), R("fj")), Sub(R("s"), R("j"))), R(" + 0.05")])));
C.push(p("**Size factors.** Dividing by the total coverage alone would let the few hepatocyte transcripts that " +
  "dominate liver RNA (ALB, APOA1, …) set the scale; their share falls in tumours and in injured liver, which " +
  "lifts every other transcript. *ncov* is therefore divided by median-of-ratios size factors *s*_{j}, computed " +
  "from the recount3 gene sums *y*_{gj} (scaled by *T*_{j}) of protein-coding genes detected in ≥ 90% of the samples, " +
  "the 1% most abundant left out:"));
C.push(eq(Sub(R("s"), R("j")), R(" = exp "), Sq([Sub(R("median"), R("g")), Par([R("ln "), Sub(R("y"), R("gj")), R(" − "),
  F(R("1"), R("n")), Sum(R("j′"), null, [R("ln "), Sub(R("y"), R("gj′"))])])])));
C.push(p("rescaled to geometric mean 1 over the samples that are compared with each other: TCGA-LIHC and all GTEx " +
  "samples together (tumours, adjacent normals and the GTEx tissues are compared in D1, D6 and tumour " +
  "specificity), and each HCC replication, chronic liver disease and NMD cohort on its own. The recount3 gene " +
  "sums use GENCODE v26 identifiers and the regions GENCODE v42; genes are matched on unversioned Ensembl ids."));
C.push(p("Samples: TCGA-LIHC (424 samples: 371 primary tumours, one per patient, and 50 adjacent normals with a " +
  "matched tumour), GTEx (550 samples: all 251 liver samples and 10 random samples of each of the other 30 " +
  "tissues), and three recount3 SRA chronic liver disease cohorts: Hoang et al. 2019 (GSE130970, 78 NAFLD " +
  "biopsies with fibrosis stage 0-4 and NAFLD activity score), Suppli et al. 2019 (GSE126848, 57: normal-weight, " +
  "obese, NAFL, NASH) and Boldanova et al. 2017 (GSE84346, 28 untreated biopsies, HCV-infected vs non-HCV), " +
  "with phenotypes from the recount3 SRA sample attributes. A fourth cohort, Govaere et al. 2020 (GSE135251, " +
  "206 NAFLD biopsies with fibrosis stage and NAFLD activity score, 10 controls; 215 used, as the read 2 file of " +
  "one NASH F2 run is corrupt at ENA), is not in recount3: its reads " +
  "were aligned with STAR 2.7.11b to GRCh38 / GENCODE v42, the loci counted in the alignments over the same " +
  "regions (samtools bedcov), and gene counts (reverse-stranded library) and technical covariates taken from " +
  "STAR; phenotypes come from the GEO sample characteristics. Loci flagged `tcga_quantifiable = False` are not " +
  "tested in any expression-based statistic and are not part of its BH family."));

// 4
C.push(h1("4. Disease components D1-D6 (per locus)"));
C.push(h2("D1 HCC: tumour versus normal expression"));
C.push(p("A rank-based effect size in [−1, 1] (> 0: higher in tumours). For the *n* = 50 patients with tumour and " +
  "adjacent normal, the differences *d*_{i} = *e*_{tumour} − *e*_{normal}, zeros included (Pratt), are ranked by " +
  "|*d*_{i}|, and the matched-pairs rank-biserial correlation is"));
C.push(eq(Sub(R("D1"), R("HCC")), R(" = "), Sub(R("r"), R("paired")), R(" = "), F([Sub(R("R"), R("+")), R(" − "), Sub(R("R"), R("−"))], [Sub(R("R"), R("+")), R(" + "), Sub(R("R"), R("−")), R(" + "), Sub(R("R"), R("0"))])));
C.push(p("with *R*_{+}, *R*_{−} and *R*_{0} the rank sums of positive, negative and zero differences (Wilcoxon " +
  "signed-rank p-value with Pratt's zero handling). Pairs without a difference, common for lowly expressed " +
  "lincRNAs, thus pull *r* toward 0, and a locus undetected in most patients cannot reach an extreme *r* from its " +
  "few detected pairs. Loci detected (*ncov* > 0) in < 20% of the paired samples are not tested. " +
  "**Tumours vs GTEx liver** (reported, not scored): the 371 tumours against the 251 " +
  "GTEx livers, with the Mann-Whitney statistic *U* of tumours over normals,"));
C.push(eq(Sub(R("r"), R("GTEx")), R(" = "), F(R("2U"), [Sub(R("n"), R("T")), Sub(R("n"), R("G"))]), R(" − 1,    AUC = "), F([Sub(R("r"), R("GTEx")), R(" + 1")], R("2"))));
C.push(p("It is not scored because it mixes two studies: the TCGA adjacent normals and the GTEx livers, both " +
  "non-tumour liver, already differ with |*r*| > 0.3 for most loci, and only the paired design is free of this " +
  "batch. The mean log2 fold changes and the BH-adjusted p-values of both tests are reported (`log2FC_paired`, " +
  "`fdr_paired`, `auc_vs_gtex`, `fdr_vs_gtex`)."));
C.push(p("**Replication** (reported, not scored). The matched-pairs statistic *r*_{paired} is computed in four " +
  "independent recount3 cohorts of HCC tumours and adjacent non-tumour liver of the same patients: GSE77314 " +
  "(49 pairs), GSE77509 (20 HBV-HCC with portal vein tumour thrombus; the runs of a sample averaged, thrombi not " +
  "used), GSE94660 (21 HBV-HCC) and GSE56545 (14), quantified as in Section 3 with size factors per cohort. A feature detected in < 20% of a " +
  "cohort's samples is not tested there. Over the cohorts *c* with a statistic and *m*_{c} pairs, the combined " +
  "*r* is the pair-weighted mean Σ *m*_{c} *r*_{c} / Σ *m*_{c} and the combined p-value is Stouffer's (signed " +
  "*z* from the two-sided p, weights √*m*_{c}, as for D1 chronic below), BH-adjusted (`hcc_rep_fdr`). A locus is " +
  "`hcc_replicated` if the combined FDR < 0.05 and the combined *r* and D1 HCC lie on the same side of their " +
  "background medians (Section 2.1; missing if no cohort tested it)."));

C.push(h2("D1 chronic: expression in chronic liver disease"));
C.push(p("Expression is scaled by size factors per cohort (Section 3; Govaere: STAR gene counts). " +
  "**Per-cohort statistic.** A partial rank correlation of expression with the phenotype *y* (Hoang, Govaere: " +
  "fibrosis stage; Suppli: NAFL + NASH = 1 vs normal-weight + obese = 0; Boldanova: HCV = 1 vs non-HCV = 0), " +
  "adjusted for *k* = 3 technical covariates (log10 input reads, % uniquely mapped reads, exonic fraction = " +
  "sum of the gene sums scaled by total coverage / 10^9 in recount3, reads assigned to genes / uniquely mapped reads in Govaere). With *X* = [1, ranks of the covariates] and the hat matrix *H* = *XX*⁺, the ranks " +
  "of expression and phenotype are residualized (residuals *u*, *v*) and correlated:"));
C.push(eq(R("u = (I − H) rank(e),    v = (I − H) rank(y),    r = "), F(R("u · v"), [R("‖u‖ ‖v‖")])));
C.push(eq(R("t = r "), Sqrt(F(R("n − 2 − k"), [R("1 − "), Sup(R("r"), R("2"))])), R(",    p = 2 "), Sub(R("P"), R("t(n−2−k)")), R("(T > |t|)")));
C.push(p("A feature detected (ncov > 0) in < 20% of a cohort's samples is not tested in that cohort. The NAFLD " +
  "activity score association (Hoang, Govaere) is reported."));
C.push(p("**Combination.** Per feature, over the scored cohorts *c* (Hoang, Govaere, Boldanova) that tested " +
  "it, with *n*_{c} the samples tested, the Stouffer *z* is the statistic, so that each cohort counts by its size:"));
C.push(eq(Sub(R("z"), R("c")), R(" = sgn("), Sub(R("r"), R("c")), R(") "), Sup(R("Φ"), R("−1")), Par([R("1 − "), F(Sub(R("p"), R("c")), R("2"))]),
  R(",    "), Sub(R("D1"), R("chronic")), R(" = Z = "), F(Sum(R("c"), null, [Sqrt(Sub(R("n"), R("c"))), Sub(R("z"), R("c"))]), Sqrt(Sum(R("c"), null, Sub(R("n"), R("c")))))));
C.push(p("The n-weighted mean of the cohort correlations, Σ *n*_{c} *r*_{c} / Σ *n*_{c}, is reported as " +
  "`chronic_r`."));
C.push(p("The p-value of *Z* is 2(1 − Φ(|*Z*|)), BH-adjusted (`chronic_fdr`). Suppli is reported but not " +
  "scored: its cases differ technically from its controls (fewer uniquely mapped reads, lower exonic fraction; " +
  "the 15 samples with the lowest exonic fraction are all cases), and after rescaling and adjustment 40% of loci " +
  "keep p < 0.01 with a median r of +0.14 (Hoang: 9%, 0.00)."));
C.push(p("**Parent genes.** Processed pseudogene loci can pick up reads of their parent genes (e.g. KRT8, KRT18, " +
  "YWHAZ), which rise in injured liver. Let *g* be the homology parent of the locus (best DIAMOND hit, Section 6) " +
  "and ρ the median over cohorts of the Spearman correlation between locus and parent expression. The locus is " +
  "`parent_driven` if ρ ≥ 0.5, the parent itself is associated with chronic disease (FDR < 0.05) and the signs " +
  "agree (parent FDR over the parent genes tested); its cohort statistics, scored and cell-intrinsic, are then " +
  "recomputed with the parent's expression as a (*k* + 1)th covariate, and the FDR is recomputed over all features."));
C.push(p("**Cell-intrinsic statistic (reported, not scored).** `D1_chronic_cell` repeats the computation with the " +
  "liver cell-type marker scores as extra covariates (mean z-score of log expression of MacParland et al. 2018 " +
  "markers for hepatocytes, stellate cells, cholangiocytes, macrophages, lymphocytes and endothelial cells; " +
  "COL1A1, KRT8 and KRT18 left out), i.e. the part of the association not explained by a change in cell mix."));
C.push(p("**D1 component.** D1 = 0.7 · *p*(D1_{HCC}) + 0.3 · *p*(D1_{chronic}), as a composite (Equation 2)."));

C.push(h2("D2 prognosis"));
C.push(p("For the TCGA-LIHC primary tumours with TCGA-CDR endpoints, age and AJCC stage (*n* = 341 for overall " +
  "survival, OS; 342 for the progression-free interval, PFI), Cox proportional hazards models are fitted per " +
  "feature on standardized expression *x* = (*e* − ē)/sd(*e*), adjusted for age, late stage (III-IV vs I-II) and " +
  "the log size factor of the tumour, ln *s* (Section 3; less differentiated tumours lose hepatocyte transcripts, " +
  "and this composition shift itself predicts survival):"));
C.push(eq(R("h(t) = "), Sub(R("h"), R("0")), R("(t) exp("), Sub(R("β"), R("x")), R("x + "), Sub(R("β"), R("age")), R("age + "),
  Sub(R("β"), R("stage")), R("stage + "), Sub(R("β"), R("s")), R("ln s),    "), Sub(R("z"), R("E")), R(" = "), F(Sub(R("b"), R("x")), [R("SE("), Sub(R("b"), R("x")), R(")")]),
  R(",    D2 = "), F([Sub(R("z"), R("OS")), R(" + "), Sub(R("z"), R("PFI"))], R("2"))));
C.push(p("D2 > 0 means higher expression with worse outcome. Hazard ratios per standard deviation, exp(*b*_{x}), " +
  "*b*_{x} being the estimate of β_{x}, and BH-adjusted p-values are reported. Features expressed in < 20% of tumours are not tested. " +
  "OS and PFI are correlated endpoints, so D2 is close to a single endpoint."));

C.push(h2("D3 somatic copy number"));
C.push(p("TCGA-LIHC primary tumours (GDC GRCh38 masked copy-number segments from the UCSC Xena GDC hub): segment " +
  "value (log2 ratio) at the midpoint of the locus; ρ_{CN} = Spearman correlation of CN with expression, partial " +
  "on tumour purity (ABSOLUTE calls; rank residuals after linear regression on the purity rank), because " +
  "non-tumour cells shift expression and CN amplitude together; gain and loss frequencies *f*_{gain}, *f*_{loss} " +
  "= fractions of tumours with CN > +0.3 and < −0.3:"));
C.push(p("HCC copy-number changes span whole chromosome arms (1q and 8q gains), where most expressed genes " +
  "follow the copy number to some degree, so ρ_{CN} alone mostly says where a locus lies. The dosage response is " +
  "therefore measured relative to the neighbourhood: the same purity-adjusted correlation is computed for the " +
  "protein-coding genes (recount3 gene sums of the same tumours, scaled by the size factors, detected in ≥ 20% of " +
  "them), and ρ_{nb} is its median over the 20 protein-coding genes nearest to the locus on the same chromosome:"));
C.push(eq(R("D3 = max("), Sub(R("ρ"), R("CN")), R(" − "), Sub(R("ρ"), R("nb")), R(", 0) · ("), Sub(R("f"), R("gain")), R(" − "), Sub(R("f"), R("loss")), R(")")));
C.push(p("set to 0 where the CN-expression correlation is not significant (BH-FDR ≥ 0.05). D3 is thus large only " +
  "when expression follows a frequent copy change more closely than the genes around it (D3 > 0: gain, < 0: " +
  "loss). The unadjusted correlation, ρ_{nb} (`cn_rho_neighbours`) and the excess (`cn_rho_excess`) are reported."));

C.push(h2("D4 germline genetics"));
C.push(p("GWAS Catalog associations with whitelisted liver traits (EFO labels of hepatocellular carcinoma, " +
  "MASLD/MASH, fibrosis, cirrhosis, alcoholic liver disease, viral hepatitis, and the liver enzymes ALT, AST, " +
  "GGT) within the locus span ± 50 kb (1-based positions), a proxy for an LD window:"));
C.push(eq(R("D4 = "), Sub(R("max"), R("a ∈ A(f)")), R(" "), Par([R("−"), Sub(R("log"), R("10")), R(" "), Sub(R("p"), R("a"))]), R(",    D4 = 0 if A(f) = ∅")));
C.push(p("For the summary scores D4 is split into D4_{HCC} (the hepatocellular carcinoma trait) and D4_{chronic} " +
  "(the other liver disease traits and the enzymes). GTEx v8 liver eQTLs and their link to a GWAS hit are " +
  "reported, not scored, as only annotated loci can have GTEx eQTLs."));

C.push(h2("D5 prior knowledge"));
C.push(p("Europe PMC title/abstract counts for the exact symbol or up to five HGNC previous symbols or aliases, " +
  "AND an RNA term, AND a liver term. Usable aliases have ≥ 5 characters with a digit and are not another " +
  "gene's approved symbol, not listed by HGNC for more than one gene, not the stem of another gene's symbol " +
  "(EEF1B, the protein name of EEF1B2) and not a generic token (a Greek letter or 1-3 letters with a number). " +
  "Europe PMC splits phrases at hyphens, so a name *G*-AS, -IT or -DT is dropped when it has ≥ 0.5 times the " +
  "liver papers of *G* itself (*G* with ≥ 20); a locus without a usable name has missing counts. Exact-phrase " +
  "matching keeps, for example, LINC01089 papers from being credited to LINC01090; loci named only by an " +
  "Ensembl, StringTie or NONCODE identifier count 0."));
C.push(eq(R("D5 = "), Sub(R("log"), R("2")), R("(1 + "), Sub(R("n"), R("liver")), R(" − "), Sub(R("n"), R("retracted")), R("),    "),
  Sub(R("D5"), R("HCC")), R(" = "), Sub(R("log"), R("2")), R("(1 + "), Sub(R("n"), R("liver ∧ cancer")), R("),    "),
  Sub(R("D5"), R("chronic")), R(" = "), Sub(R("log"), R("2")), R("(1 + "), Sub(R("n"), R("liver ∧ chronic")), R(")")));
C.push(p("Cancer terms: cancer, carcinoma, tumour, HCC, hepatoma, hepatocellular, oncogen*, metasta*. Chronic " +
  "disease terms: NAFLD, MASLD, NASH, MASH, steatosis, steatohepatitis, fatty liver, fibrosis, cirrhosis, " +
  "hepatitis, liver injury/failure, cholestasis. A paper can count for both."));

C.push(h2("D6 liver context"));
C.push(p("Equal-weight composite (Equation 2) of the percentiles of three statistics: median expression in GTEx " +
  "liver, *m*_{liver}; liver specificity, *m*_{liver} − max over the other GTEx tissues *t* of *m*_{t}; and " +
  "max(log2 TPM in HepG2, Huh7) from `lincorfz`."));

// 5
C.push(h1("5. Summary scores"));
C.push(p("All summary scores are composites (Equation 2) of component percentiles with weights fixed in advance; " +
  "each is reported with its coverage."));
C.push(table(["Score", "Components (weight)", "Meaning"], [
  ["`disease_score`", "D1 0.25, D2 0.25, D3 0.20, D4 0.10, D5 0.10, D6 0.10", "overall liver disease relevance; mostly HCC (chronic disease expression carries 0.25 × 0.3 = 7.5%)"],
  ["`disease_score_noD5`", "as `disease_score` without D5", "relevance without prior literature (unbiased by how well studied a locus is); orders the candidates"],
  ["`hcc_score`", "D1 HCC 0.25, D2 0.25, D3 0.20, D4 HCC 0.10, D5 HCC 0.10, D6 0.10", "hepatocellular carcinoma"],
  ["`chronic_score`", "D1 chronic 0.50, D4 chronic 0.20, D5 chronic 0.15, D6 0.15", "chronic liver disease (MASLD/MASH, fibrosis, viral hepatitis)"],
  ["`actionability_score`", "tumour specificity, predicted HLA-I presentation, lincORF-specific HLA-I peptides (1/3 each)", "per ORF: suitability as a tumour antigen"],
], [1.4, 2.6, 2.4]));
C.push(gap());
C.push(p("Directions: `disease_direction` from D1 HCC, D1 chronic, D2 and D3; `hcc_direction` from D1 HCC, D2 and " +
  "D3; `chronic_direction` (up/down with disease) from D1 chronic (Section 2.3)."));
C.push(p("**Actionability** (per ORF). Tumour specificity = median *e* in TCGA-LIHC primary tumours − max over " +
  "GTEx tissues except testis of the median *e*. Predicted HLA class I presentation: number of distinct 8-11-mer " +
  "peptides of the ORF with an MHCflurry 2 presentation percentile ≤ 2 for any of the HepG2 (A*02:01, A*24:02, " +
  "B*35:14, B*51:01, C*04:01, C*16:02), Huh7 (A*11:01, B*54:01, C*01:02; Cellosaurus) or 16 common alleles, " +
  "leaving out self-peptides, those that also occur in an annotated protein (GENCODE v42 translations, UniProt " +
  "human with isoforms), to which T cells are tolerant. Identified peptides: lincORF-specific HLA-I " +
  "immunopeptides (lincorfz immunopeptidomics). Percentiles are taken among all ORFs, locus-weighted (Section 2.1)."));
C.push(p("**Fitted weights (reported, not used).** Non-negative logistic regression of literature labels " +
  "(loci with ≥ 5 liver papers vs named loci without one; loci named only by an Ensembl, StringTie or NONCODE " +
  "identifier are not labelled) on the component percentiles, with an L2 penalty chosen by 5-fold " +
  "cross-validated AUC, gives `disease_score_fitted`. On run v2 (164 positives, 998 negatives) the fitted weights " +
  "are D3 0.46, D1 0.36 and D4 0.18 (D2 and D6 0), with a cross-validated AUC of 0.73 against 0.64 for the default " +
  "weights. Because the labels mark lncRNAs studied in liver rather than disease relevance, the default weights are " +
  "kept and the fit is reported as a sensitivity analysis."));

// 6
C.push(h1("6. Quality checks and copy attribution (per ORF)"));
C.push(p("**Flags.** Overlap of the ORF's coding exons with segmental duplications (UCSC genomicSuperDups, identity `fracMatch`) or with " +
  "GRCh38 regions falsely duplicated relative to T2T-CHM13; peptide homology (DIAMOND blastp --more-sensitive " +
  "against GENCODE v42 translations and the UniProt human reference proteome; hits to genes overlapping the ORF " +
  "discarded; homolog = best hit with E ≤ 1e-5 and query coverage ≥ 30%); the lincorfz flags antisense to " +
  "coding exons, small-RNA overlap, multimapper-dominated and possible coding-gene extension; " +
  "`tcga_quantifiable` (Section 3)."));
C.push(p("**Copy-attribution confidence** *a* ∈ [0, 1], for ORFs in duplicated regions or with a protein homolog " +
  "(*a* = 1 otherwise), is the mean of the available parts (0.5 if none):"));
C.push(bullet("expression: fraction of the ORF's P-sites from uniquely mapped reads;"));
C.push(bullet("peptides: 1 if an identified peptide (proteomics or immunopeptidomics) maps to this ORF only; missing otherwise, as absence of peptides is no evidence against the copy;"));
C.push(bullet("duplication age: clip((0.98 − identity)/0.08, 0, 1), identity = `fracMatch`, else protein identity to the homolog; old duplications make the conservation alignment more likely to hold orthologs of this copy."));
C.push(p("The first two parts are evidence that the copy itself is translated; the age only says how far the " +
  "conservation alignment belongs to the copy."));
C.push(p("**Confidence tier**: `likely_artefact` = duplicated or homologous with < 50% unique P-sites (or " +
  "multimapper-dominated), and no copy-specific peptide; `expressed_paralog` = " +
  "duplicated or homologous otherwise; `caution` = any lincorfz flag or not quantifiable; `clean` = none. Flags " +
  "and tiers do not remove ORFs and do not change scores."));

// 7
C.push(h1("7. Mechanism classes M1-M6 (per ORF)"));
C.push(p("Each feature is converted to a percentile among all 11,910 ORFs, locus-weighted (Section 2.1; locus " +
  "features are copied to every ORF of the locus). The **class score** *S*_{k} is the equal-weight composite " +
  "(Equation 2) of its feature percentiles, so missing features count 0.5. Raw class scores have different " +
  "spreads (few or binary features give wide ones), so classes are compared on **rank-based standardized " +
  "scores**, *z*_{k} = Φ^{−1}(*p*_{k}(*S*_{k})), with *p*_{k} the locus-weighted percentile of the class score " +
  "among the ORFs with data for class *k* (coverage > 0), clipped to [1/2*N*, 1 − 1/2*N*]. The top of every class " +
  "maps to the same quantiles, so a class with one rare binary feature cannot win by its spread alone; an ORF " +
  "without data for a class has no *z* for it. The **dominant class** is the class with the largest *z* if " +
  "*z* > 0 and it exceeds the second largest by ≥ 0.5, else none. Classes are not mutually exclusive: the scores " +
  "are absolute evidence, the *z* relative to the other ORFs."));
C.push(p("The M1 features that a paralogous copy can supply as well as this one (translation strength, P-site " +
  "share, sequence features, PhyloCSF, PhyloP contrast) are weighted by the copy-attribution confidence *a*, so " +
  "that they count only as far as they can be attributed to this copy:"));
C.push(eq(Sub(R("p"), R("w")), R(" = 0.5 + (p − 0.5) · a")));
C.push(p("Features that rise with ORF length (translation strength, P-site share, sequence features and " +
  "PhyloCSF; locus-weighted Spearman 0.19-0.56 with length) are ranked among the background ORFs of the same " +
  "length bin (20 quantile bins of the background ORF lengths), so that M1 does not reward length as such; length " +
  "itself is not a feature. Features that rise or fall with HepG2 expression (PhyloCSF, NMD sensitivity, low " +
  "polysome association; |Spearman| 0.34-0.44 with log2 TPM HepG2) are ranked among the ORFs of the same " +
  "expression bin (5 quantile bins; no HepG2 expression is a bin of its own), PhyloCSF in joint length × " +
  "expression bins."));
C.push(table(["Class", "Features (each one percentile)"], [
  ["M1 Encoded peptide", "translation strength (pooled log10 Bayes factor, lincorfz); fraction of the transcript's P-sites in the ORF; identified peptides (proteomics peptides absent from canonical proteins + lincORF-specific immunopeptides); independent translation (same stop as a GENCODE Ribo-seq ORF, +1 with phase 1 peptide support); PhyloCSF; ORF − non-ORF exon PhyloP100way mean over quantifiable positions (≥ 30 positions each); sequence features = 0.5 · [TM segment, signal-like N-terminus or Pfam-A domain] + 0.5 · (1 − disordered fraction, metapredict); CRISPR screen value (Section 8); translation, P-site share, sequence features, PhyloCSF and PhyloP contrast *a*-weighted"],
  ["M2 Translation-coupled RNA regulation", "NMD substrate (stop > 50 nt upstream of the last exon junction); NMD sensitivity (below)"],
  ["M3 Cis-regulation of a neighbour", "divergent (TSSs ≤ 5 kb, head-to-head) or antisense position to a protein-coding gene; co-expression excess with the nearest protein-coding gene (≤ 100 kb): Spearman ρ over TCGA-LIHC + GTEx liver minus the median ρ of nearest protein-coding pairs in the same distance bin (0-5, 5-20, 20-50, 50-100 kb); enhancer-like TSS chromatin, log2((H3K4me1 + 1)/(H3K4me3 + 1)) at active TSSs (H3K27ac ≥ 2-fold; inactive lowest); low polysome association (−MRL percentile)"],
  ["M4 Trans-acting RNA, nuclear", "abundance (log2 TPM HepG2); eCLIP enrichment of nuclear RBPs"],
  ["M5 Trans-acting RNA, cytoplasmic", "cytoplasmic enrichment (−nuclear); eCLIP enrichment of cytoplasmic RBPs; no data (coverage 0) for loci not expressed in HepG2 and without fractionation data"],
  ["M6 Small-RNA host", "number of snoRNA, scaRNA or miRNA genes within the locus span on the same strand (GENCODE + UCSC wgRna)"],
], [1.4, 5]));
C.push(gap());
C.push(p("**ENCODE HepG2 features** use stranded data, measured on the exons of the feature minus same-strand " +
  "exons of other GENCODE genes. Nuclear enrichment from CeFra-seq fractionation (mean signal per position, two " +
  "replicates averaged): log2[(nucleus + 0.01)/(cytosol + 0.01)], missing if both < 0.05; its negative, " +
  "cytoplasmic enrichment, is a feature of M5 only, as cis-acting (M3) and nuclear trans-acting (M4) RNAs are " +
  "both nuclear and the feature cannot tell them apart. eCLIP: reproducible " +
  "(IDR) peaks of 105 RBPs on the locus strand with log2 enrichment ≥ 3 and −log10 p ≥ 3; RBPs grouped by " +
  "immunofluorescence (Van Nostrand et al. 2020) into nuclear-only (17) and cytoplasmic-only (15). For group *G* " +
  "with share *f*_{G} of the 105 RBPs, *n*_{G} bound RBPs of the group and *n* bound RBPs in total:"));
C.push(eq(Sub(R("E"), R("G")), R(" = "), Sub(R("log"), R("2")), Sq(F(F([Sub(R("n"), R("G")), R(" + "), Sub(R("f"), R("G"))], R("n + 1")), Sub(R("f"), R("G"))))));
C.push(p("the over-representation of the group among the bound RBPs, shrunk toward 0 by one pseudo-RBP split by " +
  "the group shares: *E*_{G} = 0 without binding, and the further from 0 the more RBPs support it. Loci expressed " +
  "in HepG2 (TPM ≥ 1) without a bound RBP thus have *E*_{G} = 0 (no evidence either way); loci not expressed in " +
  "HepG2 are missing."));
C.push(p("**NMD sensitivity** is measured in the HeLa knock-downs of Colombo et al. 2017 (recount3 SRP083135): " +
  "6 scrambled-shRNA controls, UPF1 knock-down and SMG6 + SMG7 double knock-down (3 replicates each). Loci are " +
  "quantified as in Section 3, and *ncov* is divided by median-of-ratios size factors of the protein-coding gene " +
  "sums (Equation 8). With *ē*_{c} the mean of log2(*ncov*/*s* + 0.05) over the samples of condition *c*:"));
C.push(eq(R("NMD = "), F([Par([Sub(R("ē"), R("UPF1")), R(" − "), Sub(R("ē"), R("ctrl"))]), R(" + "),
  Par([Sub(R("ē"), R("SMG6+7")), R(" − "), Sub(R("ē"), R("ctrl"))])], R("2"))));
C.push(p("(> 0: up when NMD is inhibited), tested as this contrast in a linear model with the three conditions " +
  "(t test with the residual variance pooled within conditions, 12 − 3 df, so that a difference between the two " +
  "knock-downs does not inflate the variance; BH-FDR); missing for loci detected in < 50% of the samples. The " +
  "feature is NMD shrunk toward 0 by its standard error, NMD · τ² / (τ² + se²) with τ² = var(NMD) − median se² " +
  "over the tested loci (empirical Bayes, normal prior), so that noisy, lowly covered loci do not reach the top " +
  "by chance. Known NMD targets rise as expected " +
  "(ATF4 +1.4, GADD45B +2.0, GAS5 +1.3, ZFAS1 +2.0 log2), and loci with a predicted NMD-substrate ORF rise more " +
  "than the others (median +0.41 vs −0.60, Mann-Whitney p = 2 × 10⁻⁸⁶)."));

// 8
C.push(h1("8. CRISPR screen evidence"));
C.push(p("ORFs are matched to published knock-out screens by their stop codon. **Prensner et al. 2021** (553 " +
  "non-canonical ORFs, 8 cell lines including HepG2, secondary tiling screen of the hits): 2 = viability " +
  "phenotype specific to the ORF; 1 = phenotype not resolved by the secondary screen; 0 = no phenotype or one " +
  "attributed to the transcript, the region or an adjacent coding gene. **Chen et al. 2020** (~2,300 ORFs, " +
  "K562 and iPSC, coordinates lifted from hg19): 1 = hit, 0 = screened without a hit. Their hit calls are not " +
  "published, and their null (non-targeting sgRNAs) makes nearly every ORF significant because Cas9 cutting " +
  "alone costs a phenotype of about −0.015; hits are re-called as BH-FDR < 0.05 in a one-sided Mann-Whitney test " +
  "of the ORF's sgRNA phenotypes against the non-targeting sgRNAs shifted by −0.015, with a top-3 sgRNA phenotype " +
  "≤ −0.02 (146 hits in both cell lines vs 401 reported: stricter)."));
C.push(p("`crispr_screen_value` = Prensner's secondary-screen verdict where it exists (the only test of ORF " +
  "specificity), otherwise the larger of the two values; missing if neither screen tested the ORF. ORFs with a " +
  "value ≥ 1 form the screen-supported list."));

// 9
C.push(h1("9. Positive controls"));
C.push(p("Each score has positive controls, chosen in advance and processed like the candidates. Controls that " +
  "fail are reported, not fixed by tuning weights. The mechanism controls are feature checks: they show that " +
  "the features of a class are defined as intended, not that the mechanism is detected, as several pass on " +
  "annotation features alone (divergent position, small-RNA host status) and M4 on the very high abundance of " +
  "MALAT1 and NEAT1. The checks are defined once (`controls.PANEL`); a check with a missing value is " +
  "\"not evaluated\", not \"fail\". Result on run v2 (25 pass, 8 fail, 2 not evaluated):"));
C.push(table(["Score", "Control", "Expectation", "Result"], [
  ["D1 HCC", "HULC (up), LINC00261 (down)", "log2FC_paired in that direction, FDR < 0.05", "LINC00261 pass; HULC fails (log2FC 0.31, FDR 0.11)"],
  ["D1 HCC replication", "HULC (up), LINC00261 (down)", "combined *r* in that direction, FDR < 0.05", "LINC00261 pass; HULC fails (down, r −0.46)"],
  ["D1 chronic", "COL1A1", "Hoang *r* > 0, p < 0.05", "pass (Hoang r 0.40)"],
  ["D2", "BIRC5, CYTOR", "HR of overall survival > 1", "pass (HR 1.35 and 1.27)"],
  ["D3", "PVT1", "ρ_{CN} > 0 (FDR < 0.05) and D3 > 0", "fails: ρ 0.49 equals that of the neighbouring genes, D3 = 0"],
  ["D4", "PNPLA3, TM6SF2", "genome-wide significant liver-trait GWAS hit (p < 5e-8)", "pass"],
  ["D5", "HULC; LINC01090", "D5 percentile ≥ 0.9; ≤ 1 liver paper (no LINC01089 papers credited)", "pass"],
  ["D6", "HNF1A-AS1", "liver-specific (D6 ≥ 0.75)", "fails (D6 0.71)"],
  ["disease / fitted disease / HCC score", "HULC, H19", "> 90th percentile of the background", "H19 pass; HULC fails (0.45, 0.44, 0.41)"],
  ["chronic score", "COL1A1, PNPLA3, TM6SF2", "> 90th percentile", "pass"],
  ["actionability", "GPC3 (FVGEFFTDV, HLA-A*02:01)", "tumour specificity > 0; presentation percentile ≤ 2", "pass (percentile 0.05)"],
  ["M1", "MTLN, SPAAR", "top class (highest *z*)", "not evaluated (not lincorfz loci)"],
  ["M2", "GAS5, ZFAS1", "M2 score ≥ 0.6", "pass"],
  ["M3", "HNF1A-AS1, LINC00261", "top class", "pass"],
  ["M4", "MALAT1, NEAT1", "top class", "pass"],
  ["M5", "NORAD", "top class", "pass"],
  ["M6", "SNHG1, MIR22HG", "top class", "SNHG1 pass; MIR22HG fails (top class M4)"],
], [1.2, 1.7, 1.9, 2.2]));
C.push(gap());

// 10
C.push(h1("10. Data and software"));
C.push(table(["Resource", "Use"], [
  ["GENCODE v42 (GRCh38)", "annotation, protein-coding translations"],
  ["recount3 (TCGA-LIHC, GTEx, SRA SRP197353, SRP186450, SRP078450, SRP068976, SRP069212, SRP099053, SRP040998)", "base-level coverage, gene sums, QC metadata"],
  ["TCGA-CDR (Liu et al. 2018); ABSOLUTE purity (PanCanAtlas)", "survival endpoints; tumour purity"],
  ["UCSC Xena GDC hub", "TCGA-LIHC masked copy number"],
  ["GWAS Catalog; GTEx v8 liver eQTLs", "D4"],
  ["Europe PMC REST API; HGNC", "D5"],
  ["ENCODE HepG2 (CeFra-seq, eCLIP, histone ChIP-seq); Van Nostrand et al. 2020", "M3-M5 features"],
  ["Colombo et al. 2017 HeLa NMD knock-downs (GSE86148, recount3 SRP083135)", "M2 NMD sensitivity"],
  ["UCSC genomicSuperDups, wgRna, PhyloP100way, hg19→hg38 chain; GIAB GRCh38 false duplications", "QC, M1, M6"],
  ["UniProt human reference proteome; Pfam-A; GENCODE Ribo-seq ORFs", "homology, M1"],
  ["Prensner et al. 2021; Chen et al. 2020 supplementary tables", "CRISPR screens"],
  ["Python 3.11 (numpy, pandas, scipy, statsmodels, lifelines, pyBigWig); DIAMOND; HMMER; metapredict; MHCflurry 2.1", "computation"],
  ["Snakemake 9.27.1 (cluster-generic executor)", "workflow"],
], [3, 3]));

// ---------------------------------------------------------------- document
const doc = new Document({
  creator: "lincdisease", title: "lincdisease: Methods",
  styles: {
    default: { document: { run: { font: FONT, size: 21 } } },
    paragraphStyles: [
      { id: "Title", name: "Title", basedOn: "Normal", run: { size: 40, bold: true, color: "1F3A5F" },
        paragraph: { spacing: { after: 200 } } },
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 30, bold: true, color: "1F3A5F" }, paragraph: { spacing: { before: 320, after: 140 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 24, bold: true, color: "2E5A88" }, paragraph: { spacing: { before: 220, after: 100 }, outlineLevel: 1 } },
    ],
  },
  numbering: { config: [{ reference: "bullets", levels: [{ level: 0, format: LevelFormat.BULLET, text: "•",
    alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 540, hanging: 270 } } } }] }] },
  features: { updateFields: true },
  sections: [{
    properties: { page: { size: { width: PAGE_W, height: 16838 },
                          margin: { top: 1300, bottom: 1300, left: MARGIN, right: MARGIN } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER,
      children: [new TextRun({ children: [PageNumber.CURRENT], size: 18, color: "808080" })] })] }) },
    children: C,
  }],
});
Packer.toBuffer(doc).then((buf) => { fs.writeFileSync(out, buf); console.log(`wrote ${out}`); });
