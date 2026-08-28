alamem is an approximate local alignment tool written specifically to replicate many of the functions of BLAT, but allowing for streaming the database rather than the query to save memory, and also to be faster by avoiding base-level extension and using ideas from probabilistic sketching and pseudo-matching statistics instead.

The basic algorithm is just seed and chain without full extension. Major points are as follows:
 - alamem indexes the query and streams the database across the query index, so that we can search very large databases (like the GTDB), whose indices won't fit into memory.
 - alamem uses a default seed and stride size of 11, matching that of BLAT, though this is configurable, and then extends anchors using exact match to get MEMs.
 - MEMs are chained together, with a maximum gap length of 100 (configurable).
 - Unlike algorithms like BLAT, alamem doesn't do any alignment after chaining. Instead, we use a combination of average MEM length, aligned gap lengths, and unaligned gap lengths to approximate ANI (see paper).
 - Only hits with a lower-bound chain ANI of at least 90% (configurable) are kept.
 - Only hits of length at least 40 (configurable) are returned.

## Install and run instructions
After cloning the repo, compile via 

```
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

Then run
```
target/release/alamem reference_list.txt query.fna[.gz] out.txt
```

reference_list.txt should be a newline-delimited list of all genomes (gzipped is fine)

out.txt will be where results are written.

Note that alamem is almost symmetric with respect to reference and query in terms of output, but you should have the query be the shorter sequence, so that indexing is fast. alamem is fast enough that we can just stream the much larger reference database.

