# Aligning
_Updated: May 11, 2017_

Before any alignment, a graph needs to be generated with [vargas define](define.md). To use default parameters and produce a SAM file,

```
vargas align -g <graph_def> -U <reads> -S <aligns_out.sam>
```

where `<reads>` can be a FASTA/Q or SAM file.

## Scoring options

- `--ete` Use end to end alignment. This is faster than full local.
- `--ma <INT>` Match bonus.
- `--mp <INT>` Mismatch penalty (provide a positive number).
- `--rdg <INT,INT>` Read gap penalty <Open,Extend>.
- `--rfg <INT,INT>` Reference gap penalty <Open,Extend>.
- `--np <INT>` Penalty for non-A/C/G/T bases.

## Alignment targets

By default, reads will align to the base graph. If the input is in SAM format, specific read groups can be aligned to specific graphs with the `-a` option. The argument can either (a) be a graph to align to, or (b) A list of target mappings with the format:

```
vargas align -g <graph_def> -r <reads> -t <out.sam> -a "RG:<tag>:<val>,<subgraph>; ..."
```

where `<tag>` can be `ID` or some auxiliary tag in the SAM Read Group line. RG's with `<val>` in the `<tag>` field will be aligned to `<subgraph>`. Reads that are not associated with a read group are assigned to `VAUGRP`.

Using a SAM input where an alignment is already defined will enable the reporting of the `cf` and `ts` flags.

## Assess

If a SAM read file is provided, `-s` can attempt to match a previous scoring function. Currently Bowtie2, HISAT2, and BWA MEM are supported.


## Output

Alignments are written as SAM files. A `CIGAR` and `POS` are reported for the maximum-scoring
alignment. The SIMD scorer records both the max-scoring position and the exact graph node the max
cell was in, and the traceback is specialized by graph structure (selected automatically, or forced
with `--traceback`):

- **Linear reference / trie / forest** (every node has in-degree ≤ 1): the recorded end node's
  unique ancestor chain *is* the reference, so the alignment is recovered by walking parents and
  running one linear traceback — no search. This is exact.
- **General (re-convergent) graph**: the traceback is recomputed over a small window subgraph
  upstream of the recorded end node by enumerating that window's candidate paths (see
  `--max-trace-paths`) and keeping the best-scoring one.

Pass `--notraceback` to skip CIGAR computation. Because the traceback needs the max-scoring
position, a graph alignment must be run with `--maxonly` (not `--msonly`, which reports scores
only). `--traceback auto|linear|trie|graph` (default `auto`) forces a specialization; `linear` and
`trie` error out if the graph does not actually have that structure. Fields relevant to Vargas are
listed below.

- `POS` Start position of the maximum-scoring alignment.
- `RNAME` Maximum scoring sequence.
- `FLAG` Reverse complement flag bit if aligned to opposite strand.
- `CIGAR` Alignment of the maximum-scoring alignment (unless `--notraceback`).
- `AS` Maximum score.
- `vp` Nodes/alleles traversed by a graph alignment; non-reference alleles are marked with `*`.

### Graph traceback limitations

Linear-reference, trie, and forest alignments are exact. For a **general (re-convergent) graph** the
traceback enumerates candidate paths through a local window and re-scores each with the same affine
DP; in some cases the re-scored best path does not reach the SIMD-reported max score, in which case a
`[WARNING]` is printed and the best path found is still emitted as a valid (if slightly sub-optimal)
CIGAR. Known cases:

- **Reads shorter than the batch's longest read.** The SIMD scorer pads short reads and may score
  the padding, so its `AS` can exceed the best score achievable by the true read. (Applies to all
  structures.)
- **Ambiguous (`N`) bases.** The SIMD scorer and the traceback DP treat `N` differently.
- **Re-convergent nodes reached by paths of different lengths.** When a node has multiple
  predecessors whose paths differ in length, its single linear coordinate is ambiguous and the
  upstream enumeration may miss the true path.

Recovering these exactly requires a full partial-order (graph) DP with per-node traceback, which
is future work.
- `mp` Position of the maximum scoring cell.
- `ss` Second best score.
- `st` Strand of second best score.
- `sp` Position of the second highest scoring cell.
- `su` Sequence name of the second best score.
- `mc` Number of max-score occurrences.
- `sc` Number of second-best score occurrences.
- `gd` Read group tag. Subgraph aligned to.

`vargas convert` can be used to extract these fields into a CSV file.

## Coordinates

When a graph diverges, parallel nodes can have different lengths. To project the position onto the reference, nodes are anchored to the end of the node. For example:

```
    012 1234 5678
        ||||
        ACAC
       /    \
    AAC--GT--AAAA
         ||
    012  34  5678
```
Note how `C` and `T`, the ends of the nodes, have the same position.
