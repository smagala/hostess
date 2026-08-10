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

By default, the same things hostile produces:

| file | purpose |
| --- | --- |
| `<prefix>.clean_1.fastq.gz`, `.clean_2.fastq.gz` | dehosted reads |
| `<prefix>.hostess.json` | run report |

With `--bam`, the removed reads are kept as well:

| file | purpose |
| --- | --- |
| `<prefix>.host.bam` + `.bai` | removed reads, coordinate-sorted and indexed |
| `<prefix>.host.flagstat.txt`, `.host.stats.txt` | BAM statistics |

Load that BAM in IGV against the host reference to see exactly what was removed
and where it aligned.

## Aligners

`--aligner auto` (the default) resolves the way hostile does: bowtie2 given two
FASTQs, minimap2 given one. That infers read length from file count, so
**single-end short reads resolve to minimap2 and should be given
`--aligner bowtie2` explicitly** — hostile's own help says the same.

| input | aligner | preset |
| --- | --- | --- |
| paired short reads | bowtie2 | `-k 1 --mm` |
| long reads | minimap2 | `-ax map-ont --secondary no` |
| paired short reads, forced | minimap2 | `-ax sr --secondary no` |

## Relationship to hostile's command line

`hostess clean` takes hostile's arguments where they mean the same thing:
`--fastq1`, `--fastq2`, `--aligner`, `--index`, `-o/--output`, `-t/--threads`,
`--aligner-args`, `--force`, `--reorder`.

Two deliberately differ:

- `--airplane` is accepted and ignored. hostile uses it to stop the tool
  downloading a reference index; hostess never downloads one, so it is always
  effectively on. Accepted so existing command lines keep working.
- `--invert` is rejected, with a message pointing at `--bam`. hostile inverts
  because it throws the host reads away, so the only way to see them is to
  throw the other set away instead. `--bam` keeps both.

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

### Compared with hostile

Producing only clean reads, hostile is slightly quicker, because hostess is also
sorting, indexing and computing statistics over the host BAM:

| dataset | hostess | hostile |
| --- | --- | --- |
| 533k pairs, 0 percent host | 1618 ms | 1475 ms |
| 40k pairs, 100 percent host | 915 ms | 717 ms |

That comparison is not like for like, though, because hostile has not produced
the removed reads. Getting them out of hostile means running `clean` again with
`--invert`, which yields FASTQs rather than a BAM, and then aligning those a
third time to get something sortable and viewable.

Producing the **same artefacts** — clean FASTQs, a coordinate-sorted indexed BAM
of the removed reads, flagstat and stats — the ordering reverses:

| dataset | hostess (1 alignment) | hostile route (3 alignments) | |
| --- | --- | --- | --- |
| 40k pairs, 100 percent host | **909 ms** | 2158 ms | 2.4x |
| 533k pairs, 0 percent host | **1616 ms** | 2278 ms | 1.4x |

Verified equivalent: both routes yield 80000 reads in the host BAM and the same
clean read count.

The margin widens with host content, which is the direction real samples go. It
should also widen with sample size, since hostess pays the index load and the
read decompression once instead of three times — though that has not been
measured on a multi-hour sample.

### Compression is already BGZF

`samtools fastq` writes BGZF through htslib, not plain gzip, confirmed by the
`BC` subfield in the header and by `bgzip -t` validating the block structure.
So the output is already block-compressed: gzip-compatible for anything that
reads it, while remaining indexable and parallel-decompressible. There is
nothing to gain by piping through `bgzip` separately, and doing so would add a
process and a copy.

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

## Dependencies

`bowtie2`, `minimap2` and `samtools` on `PATH`, and Python 3.10+. Nothing
outside the standard library.

All three binaries are already present in the published hostile container
(`bowtie2`, `minimap2 2.30`, `samtools 1.22.1`, Python 3.13), so hostess can be
added to it as an alternative entry point without changing its dependencies.

## Status

Early, and not yet used in production. Covered by 100 tests and verified end to
end in all four combinations of paired/long-read and with/without `--bam`,
including that the read counts agree between the two counting paths.

## Licence

Apache-2.0. hostile is MIT, © Bede Constantinides.
