#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# unknown_isolate.fastq.gz — species ID + AMR/virulence/plasmid profiling
#
# This is the command-line twin of unknown_isolate_analysis.ipynb. I keep both
# because the notebook is what I actually ran and annotated my reasoning in,
# but a plain script is easier to drop into a Snakemake/Nextflow wrapper later
# or run unattended on a bigger machine. Same commands, same order, same
# thresholds — if the two ever disagree, the notebook is the one I trust,
# this one just hasn't been re-synced.
#
# Note on the notebook version: if you're running in Colab instead of here,
# use Google Drive (not the in-browser upload widget) to get the ~460 MB fastq
# onto the instance — the upload widget silently truncates large files often
# enough that it's not worth the risk, and /content gets wiped on every
# session restart anyway. Not relevant to this script since it just takes a
# local path.
#
# Usage:
#   ./run_pipeline.sh /path/to/unknown_isolate.fastq.gz /path/to/outdir
# ---------------------------------------------------------------------------
set -euo pipefail

FASTQ="${1:?Usage: $0 <fastq.gz> <outdir>}"
OUT="${2:?Usage: $0 <fastq.gz> <outdir>}"
THREADS="${THREADS:-2}"

mkdir -p "$OUT/dbs" "$OUT/blast"
cd "$OUT"

echo "[1/7] Tool check (install via apt/pip/GitHub-release first if missing)"
# NOTE: flye has no PyPI wheel and is only in apt on Ubuntu 24.04+ (noble).
# On older Ubuntu (e.g. 22.04/jammy, which Colab uses), install it with:
#   pip install git+https://github.com/fenderglass/Flye.git
# NOTE: mlst additionally needs libjson-perl, liblist-moreutils-perl (apt), Moo and
# Type::Tiny (cpanm — not reliably packaged correctly for apt on 22.04), and any2fasta
# (a standalone script mlst shells out to, not bundled or in apt):
#   sudo apt-get install -y libjson-perl liblist-moreutils-perl cpanminus
#   cpanm --notest Moo Type::Tiny
#   curl -sL https://raw.githubusercontent.com/tseemann/any2fasta/master/any2fasta -o /usr/local/bin/any2fasta
#   chmod +x /usr/local/bin/any2fasta
for t in seqkit minimap2 samtools blastn makeblastdb mash flye; do
  command -v "$t" >/dev/null || { echo "MISSING: $t — see README.md for install commands"; exit 1; }
done
{
  echo "seqkit:    $(seqkit version)"
  echo "minimap2:  $(minimap2 --version)"
  echo "samtools:  $(samtools --version | head -1)"
  echo "blast+:    $(blastn -version | head -1)"
  echo "mash:      $(mash --version)"
  echo "flye:      $(flye --version)"
} | tee tool_versions.txt

echo "[2/7] Raw-read QC"
seqkit stats -a "$FASTQ" | tee seqkit_stats_raw.txt
# Drop very short reads before assembly — cuts graph noise and runtime; at typical
# ONT isolate coverage (>>50x) this costs essentially nothing in contiguity.
seqkit seq -m 1000 "$FASTQ" -o filtered_1kb.fastq.gz
seqkit stats -a filtered_1kb.fastq.gz | tee seqkit_stats_filtered.txt

echo "[3/7] Quick read-based AMR triage (fast, low-confidence — sanity check only)"
curl -sL https://raw.githubusercontent.com/tseemann/abricate/master/db/card/sequences -o dbs/card.fasta
minimap2 -x map-ont -t "$THREADS" --secondary=no -c dbs/card.fasta filtered_1kb.fastq.gz \
  > reads_vs_card.paf 2> reads_vs_card.log
echo "Distinct CARD genes with a raw-read hit:"
awk '{print $6}' reads_vs_card.paf | sort -u | wc -l

echo "[4/7] Assembly (Flye, --nano-hq — see README for why, and for the memory caveat)"
flye --nano-hq filtered_1kb.fastq.gz --out-dir flye_asm --threads "$THREADS" \
  2>&1 | tee flye_run.log
echo "If this dies at 'samtools sort: couldn't allocate memory', re-run with --resume appended."

echo "[5/7] Species ID (16S BLAST) + genus-agnostic MLST"
mkdir -p dbs/16S_db && cd dbs/16S_db
curl -sL https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz -o 16S.tar.gz
tar xzf 16S.tar.gz
cd - > /dev/null
blastn -query flye_asm/assembly.fasta -db dbs/16S_db/16S_ribosomal_RNA \
  -out 16S_hits.tsv \
  -outfmt "6 qseqid sseqid pident length mismatch evalue bitscore stitle" \
  -max_target_seqs 5 -evalue 1e-50
sort -k1,1 -k7,7nr 16S_hits.tsv | sort -u -k1,1 --merge | tee 16S_best_per_contig.tsv

if command -v mlst >/dev/null; then
  # NOTE: plain `mlst` (no --scheme) auto-picks a scheme by BLAST score. On a real K. pneumoniae
  # run this tied exactly against ecoli_achtman_4 (conserved housekeeping genes score similarly
  # across related Enterobacteriaceae schemes) and silently picked the wrong one. If the species
  # is already known/suspected from the 16S step above, force the scheme explicitly rather than
  # trust the auto-detect tie-break, e.g.:
  #   mlst --scheme klebsiella flye_asm/assembly.fasta | tee mlst_result.tsv
  # Left as auto-detect by default here since this script doesn't know the species in advance —
  # but treat an auto-detected ST as provisional until cross-checked (16S / Kleborate / manual).
  mlst flye_asm/assembly.fasta | tee mlst_result.tsv
else
  echo "mlst not installed — 'git clone https://github.com/tseemann/mlst /opt/mlst' and add to PATH"
fi
# If the 16S call points to Klebsiella pneumoniae species complex, also run:
#   kleborate -a flye_asm/assembly.fasta -o kleborate_out -p kpsc
# (pip install kleborate first). Not run unconditionally — it's genus-specific.

echo "[6/7] AMR / plasmid-replicon / virulence screening"
for db in card resfinder plasmidfinder vfdb argannot; do
  curl -sL "https://raw.githubusercontent.com/tseemann/abricate/master/db/${db}/sequences" -o "dbs/${db}.fasta"
  makeblastdb -in "dbs/${db}.fasta" -dbtype nucl -out "dbs/${db}_db" > /dev/null
  blastn -query flye_asm/assembly.fasta -db "dbs/${db}_db" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle" \
    -evalue 1e-20 -out "blast/${db}_hits.tsv"
done
curl -s "https://api.github.com/repos/tseemann/abricate/commits?path=db&per_page=1" \
  | python3 -c "import sys,json; d=json.load(sys.stdin); print('abricate db/ last commit:', d[0]['sha'][:10], d[0]['commit']['committer']['date'])" \
  | tee db_provenance.txt

echo "[7/7] Done. Now run the pandas cells in the notebook (or write a small python"
echo "      script yourself) to filter blast/*_hits.tsv at >=90% identity / >=80%"
echo "      coverage, and to join hits back onto flye_asm/assembly_info.txt for the"
echo "      chromosome-vs-plasmid call. That logic is the last few cells of"
echo "      unknown_isolate_analysis.ipynb — I kept it in Python there because a"
echo "      pandas groupby is a lot more readable than the awk equivalent."

echo "All raw outputs are in $OUT — zip it as-is for the code/ submission folder."
