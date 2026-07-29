# Unknown Isolate — Species ID & AMR Profiling

Case analysis for the 15th Generation Long-Term Internship Program (Massive Bioinformatics),
built around a single input file: `unknown_isolate.fastq.gz` (raw Oxford Nanopore reads).

**Status: complete, fully confirmed.** *Klebsiella pneumoniae*, sequence type **ST258** (the
dominant globally-disseminated KPC clone), carrying a plasmid-borne **KPC-3 carbapenemase** plus
three further resistance plasmids, chromosomal porin-loss mutations, and dual-mechanism
fluoroquinolone resistance. No hypervirulence markers detected. See `findings.md` for the full
write-up, including the note to Prof. Kılıç.

## What's in here

```
.
├── README.md          <- you are here
├── findings.md         <- the final report: species, AMR, plasmids, virulence, clinician note
├── results/            <- curated final outputs (see note below)
└── code/
    ├── unknown_isolate_analysis.ipynb <- pipeline, runs clean top-to-bottom, all bugs fixed
    └── run_pipeline.sh                <- same pipeline as a plain bash script, non-notebook route
```

`results/` is a **curated** subset of what the notebook actually produced — the final assembly,
every summary table `findings.md` cites, the logs, and the species-ID outputs (16S, MLST,
Kleborate). It's deliberately *not* everything Colab zipped up: I left out Flye's disposable
intermediate working directories (draft/consensus/polished assemblies at each internal stage,
alignment BAMs) and the filtered read set, which together were ~460MB of regenerable data with
no verification value. What's here (~15MB) is enough to check every claim in `findings.md`
against the actual output without re-running the ~2hr assembly yourself.

## What broke during the real run, and what's fixed now

Getting `mlst` and 16S species ID working took several rounds of debugging against Colab's actual
environment. All of it is now folded into the notebook's install cell and the relevant analysis
cells, so a fresh run shouldn't hit any of this — but documenting it here since it's a realistic
picture of what "runs the first time" actually took:

1. **16S species-ID download** — `update_blastdb.pl` isn't on PATH in Colab's `ncbi-blast+` apt
   package. Fixed: pull the `16S_ribosomal_RNA` BLAST DB tarball directly from NCBI's FTP instead.
2. **`mlst` PATH** — a `%%bash` cell's `export PATH=...` doesn't persist to later `%%bash` cells in
   Colab (each is its own subprocess). Fixed: symlink the `mlst` binaries into `/usr/local/bin` at
   install time.
3. **`mlst` Perl dependencies** — `mlst` is a Perl tool, and Colab's base Perl is missing several
   modules it needs: `JSON`, `List::MoreUtils`, `Moo`, `Type::Tiny`. Fixed via a mix of apt packages
   (`libjson-perl`, `liblist-moreutils-perl`) and `cpanm` (`Moo`, `Type::Tiny` — not reliably
   available as correctly-named apt packages on Ubuntu 22.04).
4. **`any2fasta` missing** — a small standalone script `mlst` shells out to internally, not bundled
   with the `mlst` git repo and not in apt. Fixed: fetched directly from its own GitHub repo.
5. **`mlst` scheme auto-detection picked the wrong species** — with no `--scheme` flag, `mlst`
   scores every scheme by BLAST hits and picks the best match. On this genome it tied *exactly*
   between `ecoli_achtman_4` and `klebsiella` (conserved housekeeping genes like *gyrB*/*mdh*/*recA*
   score similarly across related Enterobacteriaceae schemes) and silently picked the wrong one
   (E. coli) on the first pass. Since species was already confirmed *K. pneumoniae* by that point
   (16S, and later Kleborate), the fix was to force `--scheme klebsiella` explicitly rather than
   trust the tie-break. **This one is worth flagging generally**: `mlst`'s auto-detect is not
   reliable proof of species on its own when close relatives are in play — always cross-check
   against an independent species call (16S, or a purpose-built tool like Kleborate) before trusting
   an auto-detected ST.
6. **Kleborate confirmed everything independently** and additionally caught two resistance
   mechanisms a gene-presence BLAST screen structurally can't see: porin-loss frameshift mutations
   (OmpK35, OmpK36) and fluoroquinolone target-site mutations (GyrA, ParC) — both are point mutations
   in native genes, not acquired genes, so they don't show up in a CARD/ResFinder BLAST screen no
   matter how the thresholds are tuned. This is the actual reason the pipeline runs both a generic
   BLAST screen *and* a species-specific tool rather than relying on either alone.

## How to run it (Colab)

1. Upload `unknown_isolate.fastq.gz` to your Google Drive first, via the normal Drive web UI (not
   Colab's in-notebook upload widget — at ~460 MB it's unreliable and silently truncates too often,
   and `/content` gets wiped on every session restart anyway, so you'd have to re-upload each time).
2. Open a new Colab notebook, `Runtime -> Upload notebook`, pick `unknown_isolate_analysis.ipynb`.
3. `Runtime -> Run all`. The first code cell mounts your Drive and will prompt you to authorize
   access — accept it. The next cells locate the file and set `FASTQ` to its path (edit if your
   file isn't at the top level of My Drive).
4. Expect roughly:
   - Environment setup: ~3–4 min (Flye built from source via pip; `mlst`'s Perl dependencies and
     `any2fasta` installed alongside it — this is the step that took several debugging rounds to get
     right, now consolidated into one cell)
   - QC + read filtering: <1 min
   - Assembly (Flye, `--nano-hq`): the bulk of the runtime — on the actual run this produced a clean
     8-contig, 5.84 Mb assembly (5.31 Mb chromosome + 4 plasmids + 3 unresolved small contigs) in a
     bit under 2 hours on Colab's standard 2-vCPU runtime.
   - Species ID (16S + MLST + Kleborate) + AMR/virulence/plasmid database downloads and BLAST
     screening: a few minutes total.
5. The Kleborate cell runs unconditionally now (`IS_KLEBSIELLA = True`) since species is confirmed —
   if you're running this pipeline on a different sample where the organism is genuinely unknown,
   check the 16S output first and set that flag accordingly before running the Kleborate cell.
6. The last two cells zip everything in `/content/results` and copy it to your Drive as
   `results_bundle.zip`.

**On long Flye runs:** Colab free tier disconnects after ~90 min idle and has a hard session cap
around 12 hrs. A ~2 hr run is within budget, but keep the tab open and check in every 20–30 min
rather than closing the laptop — a lid-close or long idle period can kill the connection and lose
progress mid-assembly. If it does die partway through consensus/polishing, re-run the same Flye
cell with `--resume` appended rather than starting over.

## How to run it (plain bash, no notebook)

```bash
# install what you need first (Ubuntu/Debian shown; adjust for your OS)
sudo apt-get install -y minimap2 samtools ncbi-blast+ mash
# flye has no PyPI wheel, and is only in apt on Ubuntu 24.04+ (noble) — on
# anything older (incl. Ubuntu 22.04/jammy, which Colab itself runs) install
# it straight from source instead:
pip install git+https://github.com/fenderglass/Flye.git
pip install biopython pandas kleborate

# mlst + its Perl dependencies + any2fasta (see "what broke" above for why each of these is here)
git clone https://github.com/tseemann/mlst /opt/mlst && ln -sf /opt/mlst/bin/* /usr/local/bin/
sudo apt-get install -y libjson-perl liblist-moreutils-perl cpanminus
cpanm --notest Moo Type::Tiny
curl -sL https://raw.githubusercontent.com/tseemann/any2fasta/master/any2fasta -o /usr/local/bin/any2fasta
chmod +x /usr/local/bin/any2fasta

curl -sL https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_amd64.tar.gz \
  | tar xz && sudo mv seqkit /usr/local/bin/

./code/run_pipeline.sh /path/to/unknown_isolate.fastq.gz ./results
```

The script covers raw output generation (assembly, BLAST hit tables, 16S/MLST calls — with
`--scheme klebsiella` forced, per the fix above). The identity/coverage filtering and
chromosome-vs-plasmid join logic lives in the notebook's pandas cells, with the reasoning for each
threshold documented right next to the code.

## Why this approach

The notebook's markdown cells carry my reasoning for every methodological choice next to the code
that implements it — why `--nano-hq` over `--nano-raw`, why filter reads before assembly, why three
AMR databases instead of one, why the chromosome-vs-plasmid call rests on two independent lines of
evidence, why `mlst`'s scheme auto-detect isn't trustworthy on its own, and why Kleborate ran as a
second, species-specific pass rather than relying on the generic BLAST screen alone (it caught
resistance mutations the BLAST screen structurally cannot see).

## A note on database currency

AMR/virulence databases (CARD, ResFinder, ARG-ANNOT, PlasmidFinder, VFDB) are pulled from the
`tseemann/abricate` GitHub mirror at run time; the pipeline logs both the fetch timestamp and the
last commit hash touching `db/` in that repo (`db_provenance.txt`). Kleborate additionally cross-
checked the AMR calls against its own bundled CARD v3.2.9 snapshot — the KPC-3 call is confirmed by
two independent database snapshots, not just one. Re-running this notebook later can legitimately
produce different allele-level gene calls if any of those databases have been updated upstream since
— that's expected, not a bug, and it's exactly why the snapshot is logged.
