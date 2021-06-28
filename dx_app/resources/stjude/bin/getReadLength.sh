#!/bin/bash

bamfile=$1

samtools view $bamfile | head -n 1000 | cut -f 10 | awk '{ print length($0) }' | sort -n | tail -n 1
