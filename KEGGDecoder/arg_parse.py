import argparse
from argparse import RawTextHelpFormatter


class ArgParse:

    def __init__(self, arguments_list, description, *args, **kwargs):
        """ Class for handling parsing of arguments and error handling

        Example:

        args_list = [
            [["required_argument"],
                {"help": "Help string for argument"}],
            [["-o", "--optional"],
                {"help": "Optional argument", "default": "None"}],
            [["-r", "--required"],
                {"help": "Required argument", "required": "True"}]
        ]

        ap = ArgParse(args_list, description="Sample:\tSample program")

        ## Now you can access as ap.args.required_argument, ap.args.optional, and ap.args.required
        ## This script will handle requirement checking, and will not allow the script to launch unless required flags
            are set.
        ## Note that you CANNOT use '-' in the names of arguments!
        ## Note that any other constructor argument that argparse.ArgumentParser() takes will also be used

        Ensure that the final value in the inner list does not have '-' characters
        Include "require": True in inner dictionary to treat arg as required

        :param arguments_list: List[List[List[str], Dict[str, str]]]
        """
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
