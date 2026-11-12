version 1.0

## GTDBTK_ANI_REP
##
## Runs `gtdbtk ani_rep` once over a whole genome panel: every input assembly is
## compared with skani against all GTDB species representatives, and each is
## assigned to the species cluster of its closest representative.
##
## This answers one question directly -- do genomes that separate into different
## ANI bands against a single reference also fall into different GTDB species
## clusters? -- without paying for `classify_wf`, whose pplacer step needs about
## 140 GB of RAM for bacteria.
##
## The reference package is the expensive part (about 100 GB unpacked). Copy it
## into your workspace bucket ONCE and pass it as `gtdbtk_data_tar`:
##
##   gsutil cp <tarball from
##     https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/ > \
##     gs://<your-workspace-bucket>/refdata/
##
## Localising from the bucket is far faster than re-downloading from the
## internet on every attempt, and it is what makes a re-run cheap. If you leave
## `gtdbtk_data_tar` empty the task falls back to downloading `gtdbtk_data_url`.
##
## Docker: the default is the Seqera Wave build of GTDB-Tk 2.7.2 that nf-core
## uses, because it is NOT on Docker Hub. Docker Hub rate-limits anonymous pulls
## from cloud IP ranges, and a throttled pull kills the job before the container
## starts -- no stderr, no stdout, no rc, which is exactly what the first run of
## this workflow produced with ecogenomic/gtdbtk. This image is also conda-built,
## so it has no ENTRYPOINT to fight Cromwell's script. 2.7.0 is the minimum
## version supporting GTDB release R232.

workflow gtdbtk_ani_rep {
  input {
    Array[File] assemblies
    String      output_prefix        = "panel"

    File?       gtdbtk_data_tar
    String      gtdbtk_data_url      = "https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz"

    Float       min_af               = 0.5
    String      docker               = "community.wave.seqera.io/library/gtdbtk:2.7.2--64b0fd171db01270"
    Int         cpu                  = 16
    Int         memory_gb            = 64
    Int         disk_gb              = 320
    # 0, not 1. The reference download runs for about 100 minutes before any
    # useful work happens, and a preemption in that window throws all of it
    # away. This task was preempted twice at roughly 1h50m before this default
    # was changed. A preemptible VM is false economy for a long unattended
    # download; set this to 1 only if the package is already staged in a bucket
    # and passed as gtdbtk_data_tar.
    Int         preemptible          = 0
  }

  call ani_rep {
    input:
      assemblies      = assemblies,
      output_prefix   = output_prefix,
      gtdbtk_data_tar = gtdbtk_data_tar,
      gtdbtk_data_url = gtdbtk_data_url,
      min_af          = min_af,
      docker          = docker,
      cpu             = cpu,
      memory_gb       = memory_gb,
      disk_gb         = disk_gb,
      preemptible     = preemptible
  }

  output {
    File  ani_closest      = ani_rep.ani_closest
    File  ani_summary      = ani_rep.ani_summary
    File  species_counts   = ani_rep.species_counts
    File  per_genome_table = ani_rep.per_genome_table
    File  run_log          = ani_rep.run_log
    String gtdbtk_version  = ani_rep.gtdbtk_version
    String gtdb_release    = ani_rep.gtdb_release
  }
}

task ani_rep {
  input {
    Array[File] assemblies
    String      output_prefix
    File?       gtdbtk_data_tar
    String      gtdbtk_data_url
    Float       min_af
    String      docker
    Int         cpu
    Int         memory_gb
    Int         disk_gb
    Int         preemptible
  }

  command <<<
    set -euo pipefail

    echo "### gtdbtk ani_rep -- started $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    gtdbtk --version 2>&1 | tee version.txt || echo "unknown" > version.txt
    head -n 1 version.txt | tr -d '\n' > gtdbtk_version.txt

    # The conda-built GTDB-Tk images ship no curl and no wget, which is what
    # killed the first run of this workflow (exit 127). Every external tool this
    # task depends on is therefore checked up front, so a missing one is named
    # here rather than discovered forty lines later.
    echo "### tool availability"
    for tool in python3 tar gzip cp basename grep head tr; do
      printf '  %-10s %s\n' "${tool}" "$(command -v "${tool}" || echo MISSING)"
    done

    # ---------------------------------------------------------------- reference
    mkdir -p refdata
    TARBALL="~{default='' gtdbtk_data_tar}"
    if [ -n "${TARBALL}" ] && [ -s "${TARBALL}" ]; then
      echo "### unpacking reference package localised from the bucket"
    else
      echo "### no gtdbtk_data_tar given; downloading the reference package"
      TARBALL="gtdbtk_data.tar.gz"
      # The URL goes through a file rather than into the Python source, so the
      # heredoc stays quoted and WDL does not interpolate inside it.
      printf '%s' "~{gtdbtk_data_url}" > data_url.txt
      cat > fetch.py <<'FETCHEOF'
import shutil
import sys
import time
import urllib.request

url = open("data_url.txt").read().strip()
dest = "gtdbtk_data.tar.gz"

for attempt in range(1, 4):
    try:
        print(f"### download attempt {attempt}: {url}", flush=True)
        start = time.time()
        with urllib.request.urlopen(url, timeout=120) as response:
            total = response.headers.get("Content-Length")
            total = int(total) if total else None
            with open(dest, "wb") as out:
                shutil.copyfileobj(response, out, length=16 * 1024 * 1024)
        elapsed = time.time() - start
        import os
        size = os.path.getsize(dest)
        rate = size / elapsed / 1e6 if elapsed else 0
        print(f"### downloaded {size} bytes in {elapsed:.0f}s ({rate:.1f} MB/s)",
              flush=True)
        if total is not None and size != total:
            raise IOError(f"size mismatch: got {size}, expected {total}")
        break
    except Exception as exc:
        print(f"### attempt {attempt} failed: {exc}", flush=True)
        if attempt == 3:
            sys.exit("reference package could not be downloaded")
        time.sleep(15)
FETCHEOF
      python3 fetch.py
    fi

    # The package unpacks to a single release directory; --strip-components=1
    # puts its contents directly under refdata/ so GTDBTK_DATA_PATH is stable
    # whatever the directory inside the tarball happens to be called.
    if command -v tar >/dev/null 2>&1; then
      tar -xzf "${TARBALL}" -C refdata --strip-components=1
    else
      echo "### tar not present; extracting with Python"
      printf '%s' "${TARBALL}" > tarball_path.txt
      cat > extract.py <<'EXTRACTEOF'
import os
import tarfile

src = open("tarball_path.txt").read().strip()
with tarfile.open(src, "r:gz") as tf:
    for member in tf:
        parts = member.name.split("/", 1)
        if len(parts) == 1:
            continue
        member.name = parts[1]
        tf.extract(member, "refdata")
print("### extraction complete")
EXTRACTEOF
      python3 extract.py
    fi
    rm -f gtdbtk_data.tar.gz

    export GTDBTK_DATA_PATH="$(pwd)/refdata"
    echo "### GTDBTK_DATA_PATH=${GTDBTK_DATA_PATH}"
    ls "${GTDBTK_DATA_PATH}" > refdata_listing.txt
    head -20 refdata_listing.txt

    # The release is recorded so it can be cited in the Methods rather than
    # recalled later.
    if [ -s "${GTDBTK_DATA_PATH}/metadata/metadata.txt" ]; then
      grep -i -m1 'VERSION_DATA' "${GTDBTK_DATA_PATH}/metadata/metadata.txt" \
        | tr -d '\n' > gtdb_release.txt || true
    fi
    if [ ! -s gtdb_release.txt ]; then
      basename "~{gtdbtk_data_url}" | tr -d '\n' > gtdb_release.txt
    fi

    # ---------------------------------------------------------------- genomes
    # Every assembly is staged under one extension, because ani_rep selects
    # input files by extension and the panel may mix .fasta, .fna and .gz.
    # NOT write_lines(): on an Array[File] that writes the gs:// paths as they
    # stand before localisation, so every cp fails with "No such file or
    # directory" on a path that looks correct. Interpolating the array into the
    # command gives the localised paths Cromwell actually put on this disk.
    mkdir -p genomes
    n=0
    for src in ~{sep=' ' assemblies}; do
      [ -z "${src}" ] && continue
      [ -s "${src}" ] || { echo "missing or empty input: ${src}" >&2; exit 1; }
      base="$(basename "${src}")"
      base="${base%.gz}"
      base="${base%.fasta}"
      base="${base%.fna}"
      base="${base%.fa}"
      case "${src}" in
        *.gz) python3 -c "import gzip,shutil,sys; shutil.copyfileobj(gzip.open(sys.argv[1],'rb'), open(sys.argv[2],'wb'))" "${src}" "genomes/${base}.fna" ;;
        *)    cp "${src}" "genomes/${base}.fna" ;;
      esac
      n=$((n + 1))
    done
    echo "### staged ${n} assemblies"
    [ "${n}" -gt 0 ] || { echo "no assemblies staged" >&2; exit 1; }

    # ---------------------------------------------------------------- run
    gtdbtk ani_rep \
      --genome_dir genomes \
      --out_dir out \
      --extension fna \
      --prefix "~{output_prefix}" \
      --min_af ~{min_af} \
      --cpus ~{cpu} \
      2>&1 | tee -a gtdbtk_run.log

    # The prefixed names are what the documentation gives, but after a two-hour
    # download this must not fall over on a filename. Fall back to a glob.
    echo "### ani_rep output directory"
    ls -l out/ || true
    if [ -s "out/~{output_prefix}.ani_closest.tsv" ]; then
      cp "out/~{output_prefix}.ani_closest.tsv" ani_closest.tsv
    else
      found="$(ls out/*ani_closest*.tsv 2>/dev/null | head -1 || true)"
      [ -n "${found}" ] || { echo "no ani_closest output produced" >&2; exit 1; }
      cp "${found}" ani_closest.tsv
    fi
    if [ -s "out/~{output_prefix}.ani_summary.tsv" ]; then
      cp "out/~{output_prefix}.ani_summary.tsv" ani_summary.tsv
    else
      found="$(ls out/*ani_summary*.tsv 2>/dev/null | head -1 || true)"
      if [ -n "${found}" ]; then cp "${found}" ani_summary.tsv; else
        echo "no ani_summary produced" > ani_summary.tsv; fi
    fi
    cat out/gtdbtk.log >> gtdbtk_run.log 2>/dev/null || true

    # ---------------------------------------------------------------- summarise
    cat > summarise.py <<'PYEOF'
import collections
import csv
import sys

# ani_closest.tsv holds one row per input genome: its closest GTDB species
# representative, the ANI to it and the alignment fraction.
rows = []
with open("ani_closest.tsv", newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        rows.append({k.strip(): (v or "").strip() for k, v in row.items()})

if not rows:
    sys.exit("ani_closest.tsv has no rows")

cols = rows[0].keys()


def pick(*names):
    # exact match first, then a suffix match: GTDB-Tk 2.7 names the ANI and
    # alignment-fraction columns skani_ani and skani_af, and an exact-only
    # lookup silently produced empty columns in the first run of this workflow.
    for n in names:
        for c in cols:
            if c.lower() == n:
                return c
    for n in names:
        for c in cols:
            if c.lower().endswith("_" + n):
                return c
    return None


c_genome = pick("user_genome", "genome", "query")
c_ref = pick("reference_genome", "reference", "closest_genome")
c_tax = pick("taxonomy", "classification", "reference_taxonomy")
c_ani = pick("skani_ani", "ani", "ani_percent")
c_af = pick("skani_af", "af", "alignment_fraction")

with open("per_genome_table.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["genome", "closest_gtdb_representative", "species", "ani", "af"])
    for r in rows:
        tax = r.get(c_tax, "") if c_tax else ""
        species = ""
        for part in tax.split(";"):
            if part.startswith("s__"):
                species = part[3:].strip()
        w.writerow([
            r.get(c_genome, ""),
            r.get(c_ref, "") if c_ref else "",
            species or tax,
            r.get(c_ani, "") if c_ani else "",
            r.get(c_af, "") if c_af else "",
        ])

counts = collections.Counter()
with open("per_genome_table.tsv", newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        counts[r["species"] or "(unassigned)"] += 1

with open("species_counts.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["species", "genomes"])
    for name, n in counts.most_common():
        w.writerow([name, n])

print("### closest GTDB species representative, by count")
for name, n in counts.most_common():
    print(f"  {n:4d}  {name}")
if len(counts) > 1:
    print("### more than one species cluster is represented in this panel")
else:
    print("### every genome falls in a single GTDB species cluster")
PYEOF

    python3 summarise.py 2>&1 | tee -a gtdbtk_run.log

    echo "### finished $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  >>>

  output {
    File   ani_closest      = "ani_closest.tsv"
    File   ani_summary      = "ani_summary.tsv"
    File   per_genome_table = "per_genome_table.tsv"
    File   species_counts   = "species_counts.tsv"
    File   run_log          = "gtdbtk_run.log"
    String gtdbtk_version   = read_string("gtdbtk_version.txt")
    String gtdb_release     = read_string("gtdb_release.txt")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " SSD"
    preemptible: preemptible
    maxRetries: 1
  }
}
