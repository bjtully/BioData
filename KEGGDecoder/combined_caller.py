import os
from plumbum import local
from .arg_parse import ArgParse

"""
Run all three python scripts in BioData workflow to generate combined results
Display options provided in KEGG_decoder will be available for final step only 

"""

decoder = local[os.path.join(__file__, "KEGG_decoder.py")]
expander = local[os.path.join(__file__, "KEGG_expander.py")]


def run():
    ap = ArgParse(
        (
            (("-i", "--input"),
             {"help": "Input KOALA output (from blastKOALA, ghostKOALA, or kofamscan)", "required": True}),
            (("-o", "--output"),
             {"help": "Prefix for all output files", "required": True}),
            (("-v", "--vizoption"),
             {"help": "Options: static, interactive, tanglegram", "default": "interactive"}),
            (("-n", "--newick"),
             {"help": "Required input tree for tanglegram visualization", "default": "None"}),

        ),
        description="Decoder and expander workflow to parse through a KEGG-Koala outputs to determine the "
                    "completeness of various metabolic pathways."
    )
