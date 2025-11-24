# Hair_Fragmentomics
Scripts for parsing termini occurences around CpG sites from pileups. 

# CpG Pileup Termini Parser

## Overview
This repository contains a command-line Python script for parsing **5′ and 3′ DNA termini** on the forward and reverse strands from a **samtools pileup** file. The script identifies termini based on the `^` (start) and `$` (end) characters that appear in the pileup readbase string.

The parser is optimized for downstream analysis of DNA breakpoints surrounding **CpG sites**.

---

## Features
- Extracts 5′ and 3′ termini from samtools pileup files.
- Identifies **12-bp windows** surrounding each CpG site  
  (5 positions upstream of the C and 5 positions downstream of the G).
- Outputs termini counts at each position around the CpG to a tsv file.

---

## Requirements
polars
csv  # standard library
