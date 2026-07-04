#!/bin/bash

BLOCK_FILE="blocks.txt"

while IFS= read -r block; do
    # Skip empty lines and comments
    [[ -z "$block" || "$block" =~ ^# ]] && continue

    echo "Creating rule for $block"

    rucio add-rule "$block" 1 T2_CH_CERN \
        --lifetime 259200 \
		--activity "User AutoApprove" --ask-approval \
        --comment "For urgent EGM NGT studies"
done < "$BLOCK_FILE"
