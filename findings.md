# Findings — Unknown Bacterial Isolate (ONT WGS)

**Input:** `unknown_isolate.fastq.gz` — raw Oxford Nanopore reads, single bacterial isolate, submitted
by Prof. Kılıç (Molecular Biology), recovered from a patient whose infection was not responding to
standard antibiotics.

**Bottom line up front:** this is *Klebsiella pneumoniae*, sequence type **ST258** — the single most
significant globally-disseminated KPC-producing high-risk clone — carrying a plasmid-borne **KPC-3
carbapenemase**, three further resistance plasmids, chromosomal porin-loss mutations that compound the
carbapenem resistance, and target-site mutations conferring fluoroquinolone resistance on top of an
acquired resistance gene for the same drug class. This is about as multidrug-resistant as a clinical
*K. pneumoniae* isolate gets, though notably it does **not** carry the hypervirulence markers seen in
the more invasive emerging lineages.

---

## Note to Prof. Kılıç (plain-language summary)

The isolate is **_Klebsiella pneumoniae_**, confirmed three independent ways from the genome we
assembled from your sequencing reads: a 16S ribosomal RNA match (99.7% identity to the *K. pneumoniae*
type strain), and two separate strain-typing tools that both landed on the same specific strain type,
**ST258**. If that number means anything to your clinical microbiology colleagues, it should — ST258
is the single most notorious carbapenem-resistant *K. pneumoniae* lineage worldwide; it's been behind
a large share of hospital carbapenemase outbreaks globally since the mid-2000s. Finding it here, on a
strain that also carries a carbapenemase, is not a coincidence — it's the textbook pairing.

The most important single finding is a gene called **KPC-3**, a carbapenemase — an enzyme that
destroys carbapenems, the antibiotics doctors reach for when everything else has failed. This sits on
a small, circular, independently-mobile piece of DNA (a plasmid) separate from the bacterium's main
chromosome, meaning it carries a real risk of spreading to other bacteria in the same patient or on
the same ward, not just staying put in this one isolate.

There's more, though, and it's worth understanding because it changes how confident we should be that
standard "resistant" antibiotics genuinely won't work here:

- **The bacterium has also lost two of its normal pore proteins** (called OmpK35 and OmpK36) that
  antibiotics use to get into the cell. On their own these mutations raise resistance to several
  antibiotic classes; stacked on top of a working carbapenemase, they compound it further.
- **It carries mutations in its own DNA-replication machinery** (genes called *gyrA* and *parC*) that
  independently confer resistance to fluoroquinolones (the ciprofloxacin/levofloxacin family) — on
  top of an *already-acquired* gene doing the same thing. Based on this exact genetic profile, the
  analysis software predicts ciprofloxacin resistance with very high confidence (it's seen this same
  combination before in over 2,100 other genomes, and it was resistant essentially every time).
- Three further resistance plasmids are present, adding resistance to aminoglycosides, sulfonamides,
  chloramphenicol, trimethoprim, and macrolides.

On the virulence side — separate from resistance — I specifically checked for the genetic markers of
the "hypervirulent" *K. pneumoniae* strains that have been an emerging global concern (these combine
resistance *and* unusual invasiveness). **This isolate does not carry any of them.** It scores zero on
that specific virulence panel. So this is a maximally drug-resistant strain, but not, on this evidence,
an unusually invasive one beyond what *K. pneumoniae* normally is.

**My recommendation:**
1. **Phenotypic confirmation is still essential**, especially given how many resistance mechanisms are
   stacking here — a full antibiogram, not just a carbapenem screen, since fluoroquinolone resistance
   is independently predicted with high confidence too.
2. **Colistin remains a plausible treatment option based on the genotype** — no colistin-resistance
   genes or mutations were detected — but this absolutely needs phenotypic confirmation before it's
   relied on, since it's often the last resort once carbapenems are off the table.
3. **Infection control**: this meets the profile of a high-risk CRE (carbapenem-resistant
   Enterobacteriaceae) clone. I'd recommend contact precautions and considering screening of other
   patients on the ward, per your infection control team's protocol.

---

## 1. Data quality assessment

| Metric | Value |
|---|---|
| Total reads | 260,294 |
| Total bases | 576,590,333 (~577 Mb) |
| Read N50 | 15,932 bp |
| Mean read quality (Q) | 20.25 |
| Q20 / Q30 | 91.37% / 84.49% |
| GC content | 55.85% |
| Reads ≥1kb used for assembly | 83,661 reads / 506.8 Mb |

A mean Q of ~20 across called bases is a signature of modern Q20+/R10.4.1-era Nanopore chemistry,
which is why I assembled in `--nano-hq` mode rather than `--nano-raw`. GC of 55.85% is in range for
*Klebsiella pneumoniae* (~57% genomic GC), consistent with the species call below. At ~577 Mb of raw
data against a ~5.8 Mb assembled genome, that's roughly 99x raw coverage — well beyond what an
ONT-only assembly needs.

See `read_qc.png` for the length/quality distribution plots.

---

## 2. Assembly

Flye (`--nano-hq`) assembled the filtered reads into **8 contigs, 5,843,794 bp total, mean coverage
84x**. The chromosome fully circularised at a size right in the expected range for the species.

| Contig | Length (bp) | Coverage | Circular | Call |
|---|---|---|---|---|
| contig_3 | 5,306,069 | 75x | Y | **Chromosome** |
| contig_8 | 214,836 | 117x | Y | Plasmid — IncFIB(K)/IncFII(K) |
| contig_4 | 125,021 | 188x | Y | Plasmid — IncFII(Yp)/repB(R1701)* |
| contig_9 | 79,498 | 122x | Y | **Plasmid — IncI2 — carries KPC-3** |
| contig_7 | 41,154 | 67x | Y | Plasmid — IncR |
| contig_1 | 38,533 | 824x | N | Unresolved — high-copy, no genes matched |
| contig_6 | 35,481 | 112x | N | Unresolved — no genes matched |
| contig_2 | 3,202 | 313x | N | Unresolved — tiny, high-copy, likely repeat unit |

Full table: `contig_summary.tsv`. Assembly log and per-contig stats: `flye_run.log`,
`flye_asm/assembly_info.txt`.

**Chromosome size confirmed clean.** 5.31 Mb is squarely within the typical *K. pneumoniae*
chromosome range (5.2–5.9 Mb) — Kleborate's independent contig-stats module reported the identical
figure. (An earlier pass through this same pipeline had shown a chromosome contig inflated to 7.3 Mb;
that turned out to be specific to that run and did not reproduce here — flagging that it's resolved
rather than silently dropping it, since I raised it as an open concern previously.)

**contig_4 caveat still stands.** Flye still flags it as a repeat region (`mult.=2`), and it still
carries a distinct SHV-12 ESBL allele plus IncFII(Yp)/repB(R1701) replicon markers arguing for a real,
independent plasmid. I'm calling it a plasmid at moderate confidence, same as before — this one wasn't
resolved by the assembly re-run.

---

## 3. Species and strain identification — now fully confirmed

**Species: *Klebsiella pneumoniae*, high confidence, three independent lines of evidence:**

1. **16S rRNA BLAST** against NCBI's curated database: best hit **99.7% identity** to *Klebsiella
   pneumoniae* strain DSM 30104 (the species type strain), across multiple 16S copies on the
   chromosome. (`16S_hits.tsv`)
2. **MLST** (genus-agnostic `mlst` tool, pubMLST *klebsiella* scheme): 7/7 loci exact allele matches
   (gapA-3, infB-3, mdh-1, pgi-1, phoE-1, rpoB-1, tonB-79) → **ST258**. (`mlst_result.tsv`)
3. **Kleborate** (species+typing tool purpose-built for the *K. pneumoniae* species complex): species
   call "Klebsiella pneumoniae" at "strong" match confidence, and independently confirms **ST258**
   from its own MLST module. (`kleborate_out/`)

*(Note on process: my first pass through this pipeline had both the 16S and MLST steps fail due to
two Colab-environment tooling bugs — the fixes are documented in the notebook, and this is the
successful re-run.)*

**Strain type: ST258.** This is not a neutral detail — ST258 is the dominant globally-disseminated
clone behind the KPC carbapenemase pandemic in Enterobacteriaceae, first characterised in the US
mid-2000s and since found on every populated continent. Its pairing with a KPC-3 carbapenemase here is
the archetypal combination for this clone, not a coincidental co-occurrence.

**Capsule/O-antigen typing (Kaptive, via Kleborate)** — supplementary detail, not clinically central
here: K-locus **KL107** (typeable, 99.97% identity, serotype not yet mapped to a K-antigen in
Kaptive's reference table), O-locus **OL13/O13** (flagged "untypeable" — several expected O-antigen
biosynthesis genes, including *wbbB* and *wzt*, were not found intact, so the O-serotype call is
lower-confidence than the K-locus call).

---

## 4. Antimicrobial resistance genes

Screened against CARD, ResFinder, and ARG-ANNOT (≥90% identity, ≥80% reference coverage), then
cross-checked against Kleborate's dedicated *K. pneumoniae* AMR module (CARD v3.2.9 curated calls),
which agreed on every gene and additionally caught two categories of resistance my BLAST-only screen
could not: porin-loss frameshift mutations and fluoroquinolone target-site point mutations (neither is
a horizontally-acquired gene, so a straightforward gene-presence BLAST screen structurally can't see
them — this is exactly why I ran Kleborate as a second, purpose-built pass rather than relying on one
generic method).

### contig_9 (plasmid, IncI2) — the carbapenemase plasmid

| Gene | %id | Drug class |
|---|---|---|
| **KPC-3** | **100.0%** | **Carbapenem, cephalosporin, monobactam, penicillin — last-line concern** |
| OXA-9 | 100.0% | Penicillin |
| TEM-90 | 99.9% | Cephalosporin, monobactam, penicillin |
| AAC(6')-Ib-cr | 99.6% | Aminoglycoside **and** fluoroquinolone (dual-action, acquired) |
| aadA | 99.4% | Aminoglycoside |

### contig_8 (plasmid, IncFIB(K)/IncFII(K)) — class 1 integron MDR plasmid

sul1 (100%, sulfonamide), dfrA12 (100%, trimethoprim), qacEdelta1/qacE (100%, disinfectant/antiseptic
resistance), mphA (100%, macrolide), catA1 (99.8%, phenicol/chloramphenicol), aadA2 (99.9%,
aminoglycoside), APH(3')-Ia (99.4%, aminoglycoside).

### contig_7 (plasmid, IncR)

APH(4)-Ia (100%, hygromycin), sul3 (100%, sulfonamide), aadA2 (100%, aminoglycoside), AAC(3)-IVa
(100%, aminoglycoside), cmlA1 (99.9%, chloramphenicol).

### contig_3 (chromosome) — intrinsic resistome, plus two acquired-resistance mechanisms

- **Intrinsic:** oqxA/oqxB efflux pump (100%, fluoroquinolone/tigecycline/nitrofuran-relevant),
  **SHV-11** (Kleborate's curated call — corrects my earlier BLAST-only estimate of "SHV-182"; SHV-11
  is the standard *K. pneumoniae* chromosomal narrow-spectrum beta-lactamase), LEN family, OmpA and
  OmpK37 porins, KpnEF/MdtQ efflux components.
- **SHV chromosomal substitutions worth flagging:** Kleborate detected three amino-acid substitutions
  on the chromosomal SHV allele — p.Leu35Gln, **p.Gly238Ser, p.Glu240Lys**. The latter two (G238S,
  E240K) are classic ESBL-conferring substitutions when found in SHV — so even the "intrinsic" copy
  may have extended-spectrum activity, on top of the separately-acquired SHV-12 ESBL on contig_4.
- **Porin loss — new finding, not visible to a gene-presence-only BLAST screen:** frameshift mutations
  in **OmpK35** (p.Glu42fs) and **OmpK36** (p.Val59fs). Loss of these two major porins is a
  well-established carbapenem-resistance-*compounding* mechanism — on a strain that already has a
  functional carbapenemase, porin loss pushes MICs higher still and is part of why phenotypic
  susceptibility testing matters more than a single-gene genomic call here.
- **Fluoroquinolone target-site mutations — also a new finding:** **GyrA p.Ser83Ile, GyrA p.Asp87Asn,
  ParC p.Ser80Ile** (classic QRDR mutations). Combined with the acquired AAC(6')-Ib-cr gene above,
  Kleborate's curated genotype-to-phenotype database predicts **ciprofloxacin resistance with high
  confidence** (profile matched 2,147/2,167 = 99.21% of genomes with this exact combination being
  resistant; predicted MIC 2 mg/L, within the 2–4 mg/L resistant range).

### contig_4 (plasmid, moderate confidence — see assembly caveat above)

**SHV-12** (100%, a genuine ESBL-conferring SHV variant, distinct from the chromosomal SHV-11) plus a
second copy of LEN/AAK/OHIO/OKP-C family genes, consistent with the repeat-region caveat noted above.

Full filtered hit table: `resistance_summary.tsv`. Kleborate's structured AMR call table (includes the
mutation-based findings above): `kleborate_out/klebsiella_pneumo_complex_output.txt`. Raw (unfiltered)
BLAST hits: `blast/*_hits.tsv`.

**Kleborate's summary AMR score: 2/3** (carbapenemase-positive, no colistin resistance detected genomic
ally — score 3 would additionally require colistin resistance). **10 resistance classes, 19 resistance
genes/mutations total** across the isolate.

**Clinically notable, in order of importance:** (1) KPC-3 carbapenemase — plasmid-borne, mobilisable;
(2) porin loss compounding that carbapenem resistance; (3) dual-mechanism fluoroquinolone resistance
(acquired gene + target-site mutations) with a specific high-confidence resistant-phenotype prediction;
(4) no colistin resistance detected, so it remains a genotype-plausible option pending confirmation.

---

## 5. Chromosome vs. plasmid — spread risk, gene by gene

- **KPC-3, OXA-9, TEM-90, AAC(6')-Ib-cr** → contig_9, **circular, IncI2 replicon confirmed** → high
  confidence, genuinely mobilisable plasmid. This is the transmission-risk finding.
- **sul1, dfrA12, mphA, catA1, aadA2, qacEdelta1** → contig_8, circular, IncFIB(K)/IncFII(K) replicons
  confirmed → high confidence, independent plasmid.
- **APH(4)-Ia, sul3, aadA2, AAC(3)-IVa, cmlA1** → contig_7, circular, IncR replicon confirmed → high
  confidence, independent plasmid.
- **SHV-11 (+ ESBL-associated substitutions), LEN, oqxAB, OmpK35/36 porin loss, GyrA/ParC mutations**
  → all on contig_3, the chromosome → fixed in this strain's own genome; won't hop to unrelated
  bacteria the way a plasmid-borne gene can, though fully transmissible if this exact strain (ST258)
  spreads person-to-person, which — given ST258's track record — is a real possibility worth taking
  seriously for infection control.
- **SHV-12** → contig_4 → moderate-confidence plasmid, per the caveat above.

Plasmid replicon typing (PlasmidFinder): `plasmid_replicons.tsv`.

---

## 6. Virulence factors

**Kleborate virulence score: 0/5** — the isolate does not carry any of the acquired virulence loci
Kleborate specifically screens for: no yersiniabactin (*ybt*), no colibactin (*clb*), no complete
aerobactin operon (*iucABCD* — the standalone *iutA* receptor from the VFDB screen below doesn't count
as a functional acquired system without the biosynthesis genes), no salmochelin (*iro*), and no
*rmpADC*/*rmpA2* hypermucoviscosity regulators. This is the strongest single piece of evidence that
this is a "classical" resistant *K. pneumoniae* rather than a convergent hypervirulent-and-resistant
strain — worth taking at face value given the size and conservation of the gene clusters a screen at
this depth would have caught if present.

Separately, 49 distinct VFDB gene hits on the chromosome cover the standard *K. pneumoniae* housekeeping
virulence toolkit (present in essentially all clinical isolates, not indicative of anything unusual):

- **Type 1 fimbriae** (fimA–fimK, full operon) — adhesion.
- **Type 3 fimbriae** (mrkA–mrkJ, full operon) — biofilm formation.
- **Enterobactin siderophore system** (entA/B/C/E/F/S, fepA/B/C/D/G, fes) — core iron acquisition.
- **Capsule regulation** (galF, gndA, ugd, rcsA/rcsB) — standard capsule biosynthesis machinery.
- **Type VI secretion system** (full T6SS cluster) — inter-bacterial competition.
- **iutA** (aerobactin receptor only, no biosynthesis genes — see virulence-score note above).

Full table: `virulence_summary.tsv`.

---

## 7. Overall confidence and limitations

- **High confidence, multiply confirmed:** species (*K. pneumoniae*, 3 independent methods); strain
  type (ST258, 2 independent methods); plasmid-borne KPC-3 carbapenemase on a circular,
  replicon-confirmed IncI2 plasmid; three further resistance plasmids with confirmed replicons;
  chromosomal porin-loss and fluoroquinolone target-site mutations (Kleborate, curated calls); absence
  of hypervirulence markers (Kleborate virulence score 0/5).
- **Moderate confidence:** contig_4's status as a genuinely independent plasmid rather than a
  duplicated chromosome segment — Flye's repeat-region flag and the distinct SHV-12 allele pull in
  opposite directions on this one, and it wasn't resolved by the successful re-run.
- **What this analysis does not establish:** phenotypic susceptibility. Genotype-to-phenotype
  concordance for well-characterised mechanisms like KPC-3 and the GyrA/ParC quinolone mutations is
  generally high — Kleborate's own curated database gives 99.21% concordance for the ciprofloxacin
  profile specifically — but a lab susceptibility panel is what should actually drive treatment
  decisions, especially given how many resistance mechanisms are stacking on this isolate. This screen
  also can't rule out resistance mechanisms not yet catalogued in CARD/ResFinder/ARG-ANNOT/Kleborate's
  reference set.

---

## 8. Tools and database versions (reproducibility)

| Tool | Version / notes |
|---|---|
| seqkit | v2.8.2 (GitHub release binary) |
| minimap2 | via apt, Ubuntu 22.04 |
| samtools | via apt, Ubuntu 22.04 |
| Flye | installed via `pip install git+https://github.com/fenderglass/Flye.git` (no apt package on Ubuntu 22.04) |
| BLAST+ | via apt, Ubuntu 22.04 |
| mlst | github.com/tseemann/mlst, pubMLST *klebsiella* scheme |
| Kleborate | v3.2.4, using CARD v3.2.9 for AMR calls, Kaptive for K/O-locus typing |

| Database | Source | Note |
|---|---|---|
| CARD | abricate mirror (github.com/tseemann/abricate) | thresholds: ≥90% id, ≥80% coverage |
| ResFinder | abricate mirror | same thresholds |
| ARG-ANNOT | abricate mirror | same thresholds |
| PlasmidFinder | abricate mirror | ≥90% id, ≥60% coverage |
| VFDB | abricate mirror | ≥90% id, ≥80% coverage |
| 16S_ribosomal_RNA | NCBI FTP direct download | 99.7% id to *K. pneumoniae* type strain |
| CARD (via Kleborate) | v3.2.9, bundled with Kleborate v3.2.4 | independent cross-check of BLAST-based CARD calls above |

*Resistance-gene calls are only as good as the database snapshot they're screened against — these
reads screened against a newer CARD/ResFinder release could plausibly return a different allele call
on a borderline hit. The KPC-3 call here was confirmed independently by two separate database
snapshots (the abricate CARD mirror and Kleborate's bundled CARD v3.2.9), which is about as solid as
this kind of allele call gets without going to phenotypic testing.*

---

## Appendix: raw output files referenced above

```
seqkit_stats_raw.txt / seqkit_stats_filtered.txt   — QC summary tables
read_qc.png                                        — length/quality distribution plots
reads_vs_card.paf                                  — quick read-based AMR triage (not used for final calls)
flye_run.log / flye_asm/assembly_info.txt          — assembly log and per-contig stats
contig_summary.tsv                                 — assembly + AMR/VF/replicon join (chromosome/plasmid table)
16S_hits.tsv                                       — species ID (16S BLAST)
mlst_result.tsv                                    — strain typing (ST258)
kleborate_out/                                     — species/strain/AMR/virulence/K-O-locus typing (Kleborate)
resistance_summary.tsv / blast/*_hits.tsv          — AMR gene screening (filtered / raw)
plasmid_replicons.tsv                              — PlasmidFinder replicon typing
virulence_summary.tsv                              — VFDB screening
```
