#!/usr/bin/env python3
import os
import glob
import argparse
import pandas as pd
from pathlib import Path
from plumbum import local
from argparse import RawTextHelpFormatter
from plumbum.commands.processes import ProcessExecutionError

"""
Run all three python scripts in BioData workflow to generate combined results
Display options provided in KEGG_decoder will be available for final step only 

"""

decoder = local["KEGG-decoder"]
expander = local["KEGG-expander"]


class ArgParse:

    def __init__(self, arguments_list, description, *args, **kwargs):
        self.arguments_list = arguments_list
        self.args = []
        # Instantiate ArgumentParser
        self.parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description=description,
                                              *args, **kwargs)
        # Add all arguments stored in self.arguments_list
        self._parse_arguments()
        # Parse arguments
        try:
            self.args = self.parser.parse_args()
        except:
            exit(1)

    def _parse_arguments(self):
        """ Protected method for adding all arguments stored in self.arguments_list
            Checks value of "require" and sets accordingly

        """
        for args in self.arguments_list:
            self.parser.add_argument(*args[0], **args[1])

    @staticmethod
    def description_builder(header_line, help_dict, flag_dict):
        """ Static method provides summary of programs/requirements

        :param header_line:
        :param help_dict:
        :param flag_dict:
        :return:
        """
        assert set(help_dict.keys()) == set(flag_dict.keys()), "Program names do not match in key/help dictionaries"
        to_return = header_line + "\n\nAvailable Programs:\n\n"
        programs = sorted(flag_dict.keys())
        for program in programs:
            to_return += program + ": " + help_dict[program] + "\n\t" + \
                         "\t(Flags: {})".format(" --" + " --".join(flag_dict[program])) + "\n"
        to_return += "\n"
        return to_return


def default_viz(genome_df, outfile_name):
    import seaborn as sns
    import matplotlib.pyplot as plt
    sns.set(font_scale=1.2)
    sns.set_style({"savefig.dpi": 200})
    ax = sns.heatmap(genome_df, cmap=plt.cm.YlOrRd, linewidths=2,
                     linecolor='k', square=True, xticklabels=True,
                     yticklabels=True, cbar_kws={"shrink": 0.1})
    ax.xaxis.tick_top()
    # ax.set_yticklabels(ax.get_yticklabels(), rotation=90)
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    # get figure (usually obtained via "fig,ax=plt.subplots()" with matplotlib)
    fig = ax.get_figure()
    # specify dimensions and save
    # xLen = len(genome_df.columns.values.tolist())*20
    # yLen = len(genome_df.index.tolist())*20
    fig.set_size_inches(100, 100)
    fig.savefig(outfile_name, bbox_inches='tight', pad_inches=0.1)


def prefix(_path):
    """ Get prefix of file

    :param _path: Path, possibly relative
    :return:
    """
    return os.path.basename(os.path.splitext(Path(_path).resolve())[0])


def try_which():
    """ Locate hmmsearch on path, if possible

    :return:
    """
    try:
        return str(local["which"]["hmmsearch"]()).rstrip("\r\n")
    except ProcessExecutionError:
        return "None"


def print_run(cmd):
    """
    :param cmd: plumbum local object
    :return:
    """
    print(cmd)
    cmd()


def write_final_output_file(koala_list, hmm_list, output_file, ap):
    koala = pd.read_csv(open(koala_list, "r"), index_col=0, sep="\t")
    hmm = pd.read_csv(open(hmm_list, "r"), index_col=0, sep="\t")
    output_df = koala.merge(hmm, left_index=True, right_index=True)

    # Reorganize column orientation to put like pathways together
    cols = output_df.columns.tolist()
    retinal_index = cols.index('Retinal biosynthesis')
    cols.insert(retinal_index + 1, cols.pop(int(cols.index('beta-carotene 15,15-monooxygenase'))))
    cols.insert(retinal_index + 2, cols.pop(int(cols.index('rhodopsin'))))
    trans_urea = cols.index('transporter: urea')
    cols.insert(trans_urea + 1, cols.pop(int(cols.index('transporter: ammonia'))))
    nifH_index = cols.index('nitrogen fixation')
    cols.insert(nifH_index + 1, cols.pop(int(cols.index('Vanadium-only nitrogenase'))))
    cols.insert(nifH_index + 2, cols.pop(int(cols.index('Iron-only nitrogenase'))))
    dmsplyase_index = cols.index('DMSP demethylation')
    cols.insert(dmsplyase_index, cols.pop(int(cols.index('DMSP lyase (dddLQPDKW)'))))
    cols.insert(dmsplyase_index + 1, cols.pop(int(cols.index('DMSP synthase (dsyB)'))))
    sulfitereductase_index = cols.index('dissimilatory sulfite < > sulfide')
    cols.insert(sulfitereductase_index + 1, cols.pop(int(cols.index('DsrD dissimilatory sulfite reductase'))))
    output_df = output_df[cols]
    output_df.index = [idx.replace(".protein", "") + ".fna" for idx in output_df.index]
    output_df.index.name = "ID"
    print(output_df)
    output_df.to_csv(output_file, sep="\t")

    file_in = open(output_file, "r")
    genome = pd.read_csv(file_in, index_col=0, sep='\t')
    rearrange = False
    if ap.args.myorder != 'None' and os.path.exists(ap.args.myorder):
        rearrange = True
        leaf_order = []
        for line in open(str(ap.args.myorder), "r"):
            line = line.rstrip("\r\n")
            leaf_order.append(line)
        genome = genome.reindex(leaf_order)

    if ap.args.vizoption == 'static':
        from .KEGG_clustering import hClust_euclidean
        if len(genome.index) >= 2 and not rearrange:
            genome = hClust_euclidean(genome)
        default_viz(genome, prefix(ap.args.input) + ".svg")
    if ap.args.vizoption == 'interactive':
        from .Plotly_viz import plotly_viz
        plotly_viz(genome, prefix(ap.args.input) + ".html")
    if ap.args.vizoption == 'tanglegram':
        from .MakeTanglegram import make_tanglegram
        if len(genome.index) >= 3:
            make_tanglegram(genome, ap.args.newick, prefix(ap.args.input) + ".tanglegram.svg",
                            int(prefix(ap.args.tangleopt)))
        else:
            raise ValueError("Tanglegram mode requires three or more genomes")


def run():
    """ Primary calling script for flit

    :return:
    """
    possible_hmmsearch_path = try_which()
    ap = ArgParse(
        (
            (("-i", "--input"),
             {"help": "Input KOALA output (from blastKOALA, ghostKOALA, or kofamscan)", "required": True}),
            (("-p", "--protein_fasta"),
             {"help": "Protein FASTA file with ids matching KOALA output", "required": True}),
            (("-o", "--output"),
             {"help": "Prefix for all output files", "required": True}),
            (("-m", "--myorder"),
             {"help": "Orders output as specified by file", "default": "None"}),
            (("-v", "--vizoption"),
             {"help": "Options: static, interactive, tanglegram", "default": "interactive"}),
            (("-n", "--newick"),
             {"help": "Required input tree for tanglegram visualization", "default": "None"}),
            (("-t", "--tangleopt"),
             {"help": "Number of tree iterations for minimizing tangles in tanglegram", "default": "1000"}),
            (("-s", "--hmmsearch_path"),
             {"help": "Path to hmmsearch%s" %
                      ("; default: %s" % possible_hmmsearch_path if possible_hmmsearch_path else "None"),
              "default": possible_hmmsearch_path}),
        ),
        description="Decoder and expander workflow to parse through a KEGG-Koala outputs to determine the "
                    "completeness of various metabolic pathways."
    )
    # Require hmmsearch
    assert ap.args.hmmsearch_path != "None" and os.path.exists(ap.args.hmmsearch_path), "Provide valid path to " \
                                                                                        "hmmsearch "
    # Run KEGG-decoder
    decoder_outfile = prefix(ap.args.output) + ".decoder.tsv"
    print_run(
        decoder["--input", ap.args.input,
                "--output", decoder_outfile]
    )
    # Run hmmsearch
    hmmsearch_results = prefix(ap.args.input) + ".expander.tbl"
    print_run(
        local[ap.args.hmmsearch_path]
        ["--tblout", hmmsearch_results, "-T", "75",
         glob.glob(os.path.join(os.path.dirname(os.path.dirname(__file__)),
                                "/lib/*/site-packages/KEGGDecoder/HMM_Models/expander_dbv0.7.hmm"))[0],
         ap.args.protein_fasta]
    )
    # Run expander
    expander_outfile = prefix(ap.args.input) + ".expander.tsv"
    print_run(
        expander[
            decoder_outfile, expander_outfile
        ]
    )
    # Run decoder-expander and save final image
    write_final_output_file(decoder_outfile, expander_outfile, prefix(ap.args.output) + ".decode-expand.tsv", ap)


if __name__ == "__main__":
    run()
