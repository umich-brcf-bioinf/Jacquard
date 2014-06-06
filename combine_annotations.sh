#!/bin/bash

EXPECTED_ARGS=2
if [ $# -ne $EXPECTED_ARGS ]; then
    echo "Usage: `basename $0` {source_dir} {destination_dir}"
    echo "Example `basename $0` foo/epee_output foo/combined_output"
    exit 1
fi

SOURCE_DIR=$1
DEST_DIR=$2

PRESERVED_COLUMNS=2,3,4,5,6,7,10,16,17,18,19,20,21,22,23,24,25,27,28,29,30,31,32,33,34,35,36,37,38,39,40
find $SOURCE_DIR -name '*_All_Genes.txt' | xargs -I {} grep -v '^#' {} >> $DEST_DIR/all_annotation.tmp
find $SOURCE_DIR -name '*_All_Genes.txt' | head -n1 | xargs -I {} sh -c "grep -P '^#.*CHROM' -m 1 {} | cut -f$PRESERVED_COLUMNS  > $DEST_DIR/combined_annotation.txt"
cat $DEST_DIR/all_annotation.tmp | cut -f$PRESERVED_COLUMNS | sort | uniq >> $DEST_DIR/combined_annotation.txt
