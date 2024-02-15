import logging
import argparse
import os
import NPS_generation
<<<<<<< HEAD
from NPS_generation.commands import (preprocess, create_training_sets, inner_train_models_RNN,
                                     inner_sample_molecules_RNN,inner_tabulate_molecules,
                                     inner_collect_tabulated_molecules, inner_process_tabulated_molecules,
                                     inner_write_structural_prior_CV, inner_write_formula_prior_CV)
=======
from NPS_generation.commands import preprocess, create_training_sets, inner_train_models_RNN
>>>>>>> master


logger = logging.getLogger("NPS_generation")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=NPS_generation.__version__)

<<<<<<< HEAD
    modules = (preprocess, create_training_sets, inner_train_models_RNN, inner_sample_molecules_RNN,
               inner_tabulate_molecules, inner_collect_tabulated_molecules, inner_process_tabulated_molecules,
               inner_write_structural_prior_CV, inner_write_formula_prior_CV)
=======
    modules = (preprocess, create_training_sets, inner_train_models_RNN)
>>>>>>> master

    subparsers = parser.add_subparsers(title="Choose a command")
    subparsers.required = True

    def get_str_name(module):
        return os.path.splitext(os.path.basename(module.__file__))[0]

    for module in modules:
        this_parser = subparsers.add_parser(
            get_str_name(module), description=module.__doc__
        )
        this_parser.add_argument(
            "-v", "--verbose", action="store_true", help="Increase verbosity"
        )
        module.add_args(this_parser)
        this_parser.set_defaults(func=module.main)

    args = parser.parse_args()
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    args.func(args)


if __name__ == "__main__":
    main()