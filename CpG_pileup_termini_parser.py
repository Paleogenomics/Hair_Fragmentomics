import csv
import polars as pl

class CommandLine():
    def __init__(self, inOpts=None):
        import argparse
        self.parser = argparse.ArgumentParser(
            description="Get Fragmentation contexts for a list of CpG's from a mpileup file",
            add_help=True,  # default is True
            prefix_chars='-',
            usage='python CpG_pileup_termini_parser.py -i my_pileup.pileup -o my_output.tsv'
            )

        self.parser.add_argument('-i', '--input', nargs='?', type=str, default="", action='store',
                                 help='Pileup File')
        self.parser.add_argument('-o', '--outfile', nargs='?', type=str, default="output.csv", action='store',
                                 help='Name for output file')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

def writeTSV(CpG_fragKons, outdir):
    with open(outdir, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(CpG_fragKons)

def parsePileup(mpileup):
    myDF = pl.read_csv(mpileup, separator='\t', has_header=False)  # Read in pileup
    myDF = myDF.rename({"column_1": "CHR", "column_2": "Position", "column_3": "Base", 
                        "column_4": "Coverage", "column_5": "Pileup", "column_6": "Quality"})
    myDF = myDF.with_columns([
        pl.col("Position").cast(pl.Int64),
        pl.col("CHR").cast(pl.Utf8)
    ])

    print("Parsing pileup file...")

    myDF = myDF.with_columns([
        pl.col("CHR").cast(pl.Utf8),
        pl.col("Position").cast(pl.Int64),
        pl.col("Coverage").cast(pl.Float64),

        # Termini matches
        pl.col("Pileup").str.count_matches(r"\^.\.").alias("ForStrandStarts"),
        pl.col("Pileup").str.count_matches(r"\.\$").alias("ForStrandEnds"),
        pl.col("Pileup").str.count_matches(r"\^.,").alias("RevStrandStarts"),
        pl.col("Pileup").str.count_matches(r",\$").alias("RevStrandEnds"),

        # Coverage Normalized
        (pl.col("Pileup").str.count_matches(r"\^.\." ) / pl.col("Coverage")).alias("ForStrandStartN"),
        (pl.col("Pileup").str.count_matches(r"\.\$") / pl.col("Coverage")).alias("ForStrandEndN"),
        (pl.col("Pileup").str.count_matches(r"\^.,") / pl.col("Coverage")).alias("RevStrandStartN"),
        (pl.col("Pileup").str.count_matches(r",\$") / pl.col("Coverage")).alias("RevStrandEndN"),
    ])

    print(myDF['ForStrandStarts'].max())
    print(myDF['ForStrandEnds'].max())
    print(myDF['RevStrandStarts'].max())
    print(myDF['RevStrandEnds'].max())
    print(myDF['Coverage'].max())

    # Return relevant columns (saves memory by not storing pileup and quality scores)
    return myDF.select([
        "CHR", "Position", "Base", "Coverage", "ForStrandStarts", "ForStrandEnds", 
        "RevStrandStarts", "RevStrandEnds", "ForStrandStartN", "ForStrandEndN", 
        "RevStrandStartN", "RevStrandEndN"
    ])
             
def generateFragkons(parsedDF):
    # Generating column labels for positions and initialize the data structure that will store termini information for each CpG
    CpG_fragKons = []  # Storing CpG frag contexts as list of lists
    positionMap = [1, 2, 3, 4, 5, "C", "G", 8, 9, 10, 11, 12]
    forStartHeaders = ["forStart" + str(x) for x in positionMap]
    forEndHeaders = ["forEnd" + str(x) for x in positionMap]
    revStartHeaders = ["revStart" + str(x) for x in positionMap]
    revEndHeaders = ["revEnd" + str(x) for x in positionMap]
    forStartNHeaders = ["forStartN" + str(x) for x in positionMap]
    forEndNHeaders = ["forEndN" + str(x) for x in positionMap]
    revStartNHeaders = ["revStartN" + str(x) for x in positionMap]
    revEndNHeaders = ["revEndN" + str(x) for x in positionMap]
    coverages = ["Coverage" + str(x) for x in positionMap]

    headerList = ["CpGID", "SeqContext",
                  *forStartHeaders, *forEndHeaders, *revStartHeaders, *revEndHeaders,
                  *forStartNHeaders, *forEndNHeaders, *revStartNHeaders, *revEndNHeaders,
                  *coverages]

    # Obtain row indices of CpGs by detecting CG sequences
    parsedDF = parsedDF.with_columns(
        pl.col("Base").shift(-1).alias("NextBase")
    )
    indexDF = parsedDF.with_row_index()
    CpG_rowIndices = indexDF.filter((pl.col("Base") == "C") & (pl.col("NextBase") == "G"))['index'].to_list()

    print(f'Iterating through {len(CpG_rowIndices)} CpG sites')

    CpG_fragKons = [headerList] 
    for i, coord in enumerate(CpG_rowIndices):

        if i % 250000 == 0:
            print(f'Processing CpG {i}/{len(CpG_rowIndices)} at position {coord}')

        startPos = coord - 5
        endPos = coord + 7
        window_len = endPos - startPos
        if startPos < 0:
            continue

        # Slice the polars DataFrame
        CpG_DF = parsedDF.slice(startPos, window_len)
        if CpG_DF.height != window_len:
            continue

        # Check for adjacent bases
        pos_series = CpG_DF["Position"].to_list()
        if pos_series[-1] - pos_series[0] != 11:
            continue

        base_series = CpG_DF["Base"].to_list()
        upperSeqContext = ''.join(base_series).upper()
        chr_series = CpG_DF["CHR"].to_list()
        cpgid = str(chr_series[coord - startPos]) + "_" + str(pos_series[coord - startPos])
    
        CpG_fragKons.append([cpgid, upperSeqContext,
            *CpG_DF["ForStrandStarts"].to_list(), *CpG_DF["ForStrandEnds"].to_list(),
            *CpG_DF["RevStrandStarts"].to_list(), *CpG_DF["RevStrandEnds"].to_list(),
            *CpG_DF["ForStrandStartN"].to_list(), *CpG_DF["ForStrandEndN"].to_list(),
            *CpG_DF["RevStrandStartN"].to_list(), *CpG_DF["RevStrandEndN"].to_list(),
            *CpG_DF["Coverage"].to_list()])

    print("Total Number of CpGs Identified:", len(CpG_rowIndices))
    print("Total Number of CpGs Successfully Processed:", len(CpG_fragKons))
    print("Total Number of CpGs Skipped:", len(CpG_rowIndices) - len(CpG_fragKons))
    return CpG_fragKons

def main(options=None):
    cl = CommandLine(options)  # setup the command line
    pileup = cl.args.input
    outDir = cl.args.outfile

    CpG_frags = parsePileup(pileup)
    fragKons = generateFragkons(CpG_frags)
    writeTSV(fragKons, outDir)

if __name__ == '__main__':
    main()
