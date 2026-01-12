#!/bin/bash

# A script to quickly convert the depth per base files to the bed format

INPUT=$1

OUTPUT=$2

awk '{print $1"\t"$2"\t"$2"\t"$3}' $INPUT > $OUTPUT
