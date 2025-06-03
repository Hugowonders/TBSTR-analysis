#!/bin/bash
### Run GAT for genomic association and enrichment 
### GAT v1.3.4
wd=$your_working_directory

### Segment file
segFile=$wd/your_segment_file.txt

### Annotation file
annotFile=$wd/your_annotatin_file.txt

### Workspace file
### Gaps in GRCh38 were removed
### Refer to https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=gap
ws=$wd/your_workspace_file.txt

### Run GAT with 2000 permutations
gat-run.py --segments=$segFile --annotations=$annotFile --workspace=$ws \
	--num-samples=2000 --log=$wd/gat.log > $wd/gat_output_file.txt
