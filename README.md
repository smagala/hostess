# hostess

Dehost short reads in a single alignment pass, and keep the removed reads as a
coordinate-sorted BAM so you can look at them.

A rework of [hostile](https://github.com/bede/hostile). The dehosting method is
deliberately unchanged — bowtie2 against a host index, keep the pairs that align
to nothing — but the plumbing around it differs.

## Why

hostile discards the host reads. You get counts of what was removed and nothing
else, so "did the dehosting over-remove?" cannot be answered from its output.
Getting the removed reads means a second full alignment pass with `--invert`,
which on a large sample is another hour and another index staging.

hostess writes the alignment once and derives everything from it.

## What it produces

| file | purpose |
| --- | --- |
| `<prefix>.clean_1.fastq.gz`, `.clean_2.fastq.gz` | dehosted reads |
| `<prefix>.host.bam` + `.bai` | removed reads, coordinate-sorted and indexed |
| `<prefix>.host.flagstat.txt`, `.host.stats.txt` | BAM statistics |
| `<prefix>.hostess.json` | run report |

Load the BAM in IGV against the host reference to see exactly what was removed
and where it aligned.

## MultiQC

Everything is parsed by stock MultiQC. No plugin, no custom module, no config.

`<prefix>.hostess.json` matches the shape MultiQC's **Hostile** module expects:
a JSON list whose entries carry `fastq1_in_name`, `index`, `reads_in`,
`reads_out`, `reads_removed`, `reads_removed_proportion` and `version`. Its
search pattern is any `*.json` containing `"reads_removed_proportion"` in the
first 100 lines.

`flagstat` and `stats` are written verbatim, which is what MultiQC's **samtools**
module recognises — it keys off content signatures, not filenames.

`tests/test_report.py` encodes that contract, so a change that would break an
existing report fails there first.

## Single pass, not `tee`

hostile duplicates the alignment stream with `tee` so several consumers can each
filter it. hostess splits it with `samtools view -U`, which writes records
matching the filter to one destination and everything else to another:

```
bowtie2 ─┬─ selected: pairs where both mates are unmapped ─→ samtools fastq ─→ clean FASTQs
         └─ unselected: everything that hit the host       ─→ BAM ─→ sort ─→ BAM + BAI
```

Read counts come from the flagstat that is generated for MultiQC anyway, rather
than from `tee >(samtools view -c)` stages inserted to count.

## Speed

Measured, on 533k read pairs, varying only `samtools fastq --threads`:

| threads | time |
| --- | --- |
| 1 | 2456 ms |
| 2 | 1286 ms |
| 3 | 936 ms |
| 4 | 772 ms |
| 6 | 696 ms |

Compression starved at one thread costs more than anything else in the pipeline,
so it gets a real share of the budget rather than the leftovers.

**Lowering the gzip level is the wrong lever.** `-c 1` at one thread ran in
679 ms against `-c 6` at 2456 ms, but `-c 6` at four threads reached 772 ms with
output 15 percent smaller. Threads buy the same speed without the size penalty.

**Writing uncompressed and gzipping afterwards is worse**, at 1704 ms and the
largest output of any variant. The extra I/O pass costs more than it saves.

`--reorder` is off by default. bowtie2's `--reorder` buffers finished alignments
so output order matches input order; measured peak RSS rose from under 12 GB at
one alignment thread to 43.6 GB at four. Nothing here depends on read order —
the host reads are coordinate-sorted anyway. Pass `--reorder` if you need
byte-reproducible FASTQ ordering.

### It is not faster than hostile

Head to head, same data, same thread count, both without `--reorder`:

| dataset | hostess | hostile |
| --- | --- | --- |
| 533k pairs, 0 percent host | 1618 ms | 1475 ms |
| 40k pairs, 100 percent host | 915 ms | 717 ms |

hostess is roughly 25 percent slower because it is doing strictly more work:
sorting, indexing and computing statistics over the host BAM. That is the
feature, not overhead to be optimised away — the alternative is a second
alignment pass with `--invert`, which costs far more.

Both datasets are small enough that fixed costs dominate. The ratio on
multi-hour samples has not been measured, and the two things that would help
there — one pass instead of two, and `--reorder` off — both favour hostess.

## Usage

```bash
hostess \
  --fastq1 sample_R1.fastq.gz \
  --fastq2 sample_R2.fastq.gz \
  --index /path/to/bowtie2/index_prefix \
  --output results/ \
  --threads 8
```

`--dry-run` prints the pipeline without running it. Single-end input is
supported by omitting `--fastq2`.

Requires `bowtie2` and `samtools` on `PATH`.

## Status

Early. The pipeline, report and thread allocation are covered by tests and have
been run end to end, but this has not been used in production. Long-read
support via minimap2 is not implemented.

## Licence

Apache-2.0. hostile is MIT, © Bede Constantinides.
