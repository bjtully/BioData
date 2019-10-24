#!/usr/bin/env python3

import os
import random
from string import ascii_lowercase, ascii_uppercase, digits


class FastaHash:

    @staticmethod
    def create_fasta_hash(fasta_file, hash_len=25):
        """

        :param fasta_file:
        :param hash_len:
        :return:
        """
        hash_data = {}
        _fasta_file = open(fasta_file, "r")
        hash_file = open(fasta_file + ".hash", "w")
        hash_len = hash_len
        for _line in _fasta_file:
            if _line.startswith(">"):
                line = _line.rstrip("\r\n").split(" ")
                line[0] = line[0].replace(">", "")
                data = hash_data.get(line[0], None)
                if data is None:
                    hash_data[line[0]] = ''.join(random.choice(ascii_uppercase + digits + ascii_lowercase)
                                                 for _ in range(hash_len))
                    hash_file.write("%s\t%s\n" % (line[0], hash_data[line[0]]))
        hash_file.close()
        return fasta_file + ".hash"

    @staticmethod
    def create_tsv_hash(tsv_file, hash_len=25):
        hash_data = {}
        _tsv_file = open(tsv_file, "r")
        hash_file = open(tsv_file + ".hash", "w")
        hash_len = hash_len
        for _line in _tsv_file:
            line = _line.rstrip("\r\n").split("\t")
            data = hash_data.get(line[0], None)
            if data is None:
                hash_data[line[0]] = ''.join(random.choice(ascii_uppercase + digits + ascii_lowercase)
                                             for _ in range(hash_len))
                hash_file.write("%s\t%s\n" % (line[0], hash_data[line[0]]))
        hash_file.close()
        return tsv_file + ".hash"

    @staticmethod
    def rewrite_fasta(fasta_file, hash_file, out_fasta):
        """

        :param hash_file:
        :return:
        """
        _fasta_file = open(fasta_file, "r")
        out_fasta = open(out_fasta, "w")
        dict_data = FastaHash._generate_hash_dict(hash_file)
        for _line in _fasta_file:
            if _line.startswith(">"):
                line = _line.split(" ")
                line[0] = line[0].replace(">", "")
                new_val = dict_data.get(line[0], None)
                if new_val is not None:
                    out_fasta.write(_line.replace(line[0], dict_data[line[0]]))
                else:
                    out_fasta.write(_line)
            else:
                out_fasta.write(_line)
        out_fasta.close()

    @staticmethod
    def rewrite_tsv(tsv_file, hash_file, out_tsv):
        """

        :param tsv_file:
        :param hash_file:
        :param out_tsv:
        :return:
        """
        data_dict = FastaHash._generate_hash_dict(hash_file)
        W = open(out_tsv, "w")
        for _line in open(tsv_file, "r"):
            line = _line.split("\t")
            new_val = data_dict.get(line[0], None)
            if new_val is not None:
                W.write(_line.replace(line[0], new_val))
            else:
                W.write(_line)

    @staticmethod
    def _generate_hash_dict(hash_file):
        """

        :param hash_file:
        :return:
        """
        dict_data = {}
        assert os.path.exists(hash_file), "Hash file does not exist"
        hash_file = open(hash_file, "r")
        for _line in hash_file:
            line = _line.rstrip("\r\n").split("\t")
            dict_data[line[0]] = line[1]
            dict_data[line[1]] = line[0]
        return dict_data


if __name__ == '__main__':
    import sys

    # Create a hash file
    hash_file_fasta = FastaHash.create_fasta_hash("TOBG-CPC-31.fna")
    hash_file_tsv = FastaHash.create_tsv_hash("combined.ko")
    # Test with fasta file
    FastaHash.rewrite_fasta("TOBG-CPC-31.fna", hash_file_fasta, "TOBG-CPC-31.fna.rewrite")
    FastaHash.rewrite_fasta("TOBG-CPC-31.fna.rewrite", hash_file_fasta, "TOBG-CPC-31.fna.original")

    FastaHash.rewrite_tsv("combined.ko", hash_file_tsv, "combined.ko.rewrite")
    FastaHash.rewrite_tsv("combined.ko.rewrite", hash_file_tsv, "combined.ko.original")
    # Test with tsv file
