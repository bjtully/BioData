import subprocess


class SeqMap:
    """ Class creates has hash function for mapping user input file to
    non-redundant hash values

    """
    def __init__(self, in_file, out_ext=".hash"):
        """ Creates dict with mapping, sets length of each hash

        :param in_file: File to read
        :param hash_len: Length of each hash to generate
        :param out_ext:
        """
        self.data = {}
        self.hash_len = int(subprocess.getoutput("wc -l < %s" % in_file))
        self.out_ext = out_ext

    def add(self, genome_id):
        """ Adds additional genome id to the dict

        :param genome_id:
        :return:
        """
        value = self.data.get(genome_id, None)
        if value is None:
            self.data[genome_id] = SeqMap._create_hash(genome_id, self.hash_len)

    def write(self, outfile):
        """ Writes all items in stored dict to file

        :param outfile:
        :return:
        """
        with open(outfile, "w") as W:
            for k, v in self.data.items():
                W.write("%s\t%s\n" % (k, v))

    @staticmethod
    def _create_hash(genome_id, length):
        """

        :param genome_id:
        :param length:
        :return:
        """
        return genome_id

    @staticmethod
    def rewrite_fasta(fasta_file, out_file, hash_file):
        out = SeqMap(fasta_file)
        new_fasta = open(out_file + ".tmp")
        with open(fasta_file, "r") as R:
            for line in R:
                if line.startswith(">"):
                    line = line.strip("\r\n").split(" ")
                    out.add(line[0].replace(">", ""))

        out.write(hash_file)

