import os
import hashlib
import tempfile
from collections import namedtuple
from NPS_generation.clean_SMILES import main as clean_SMILES_main
from NPS_generation.augment_SMILES import main as augment_SMILES_main
from NPS_generation.train_model import main as train_model_main
import NPS_generation.data as data_folder

test_dir = os.path.join(os.path.dirname(__file__), "test_data")
data_dir = data_folder.__path__[0]


def test_clean_SMILES():
    input_file = os.path.join(data_dir, "chembl_28_2000.smi")

    with tempfile.TemporaryDirectory() as temp_dir:
        output_file = os.path.join(temp_dir, "output1.smi")
        clean_SMILES_main(input_file, output_file)

        baseline_file = os.path.join(test_dir, "output_step1.smi")
        baseline_checksum = hashlib.md5(''.join(sorted(open(baseline_file, "r").readlines())).encode('utf8')).hexdigest()
        output_checksum = hashlib.md5(''.join(sorted(open(output_file, "r").readlines())).encode('utf8')).hexdigest()

        assert baseline_checksum == output_checksum


def test_augment_SMILES():
    input_file = os.path.join(test_dir, "output_step1.smi")
    with tempfile.TemporaryDirectory() as temp_dir:
        output_file = os.path.join(temp_dir, "output2.smi")
        augment_SMILES_main(input_file, output_file, enum_factor=1)



# def test_train_model():
#     input_file = os.path.join(test_dir, "test_data", "output2.smi")
#     output_dir = os.path.join(test_dir, "test_data")
#
#     # Argv = namedtuple('Argv', ["smiles_file", "selfies", "output_dir", "rnn_type", "embedding_size", "hidden_size", "n_layers",
#     # "dropout", "bidirectional", "nonlinearity", "tie_weights", "learning_rate", "learning_rate_decay",
#     # "learning_rate_decay_steps", "gradient_clip", "seed", "batch_size", "max_epochs", "patience",
#     # "sample_idx", "sample_every_epochs", "sample_every_steps", "log_every_epochs", "log_every_steps",
#     # "sample_size", "pretrain_model", "vocab_file", "stop_if_exists"])
#     #
#     # argv = Argv("smiles_file": input_file, "selfies": False, "output_dir": output_dir, "rnn_type": "RNN",
#     #         "embedding_size": 128,"hidden_size": 16, "n_layers": 1, "dropout": 0.0, "bidirectional": False,
#     #         "nonlinearity": "tanh", "tie_weights": False,"learning_rate": 0.001, "learning_rate_decay": None,
#     #         "learning_rate_decay_steps": None, "gradient_clip": None,"seed": 0, "batch_size": 128, "max_epochs": 1000,
#     #         "patience": 100, "sample_idx": 0, "sample_every_epochs": 200,"sample_every_steps": 200,
#     #         "log_every_epochs": 1, "log_every_steps": None, "sample_size": 100000, "pretrain_model": None,
#     #         "vocab_file": None, "stop_if_exists": False)
#
#     argv = {
#         "smiles_file": input_file, "selfies": False, "output_dir": output_dir, "rnn_type": "RNN",
#         "embedding_size": 128, "hidden_size": 16, "n_layers": 1, "dropout": 0.0, "bidirectional": False,
#         "nonlinearity": "tanh", "tie_weights": False, "learning_rate": 0.001, "learning_rate_decay": None,
#         "learning_rate_decay_steps": None, "gradient_clip": None, "seed": 0, "batch_size": 128, "max_epochs": 1000,
#         "patience": 100, "sample_idx": 0, "sample_every_epochs": 200, "sample_every_steps": 200,
#         "log_every_epochs": 1, "log_every_steps": None, "sample_size": 100000, "pretrain_model": None,
#         "vocab_file": None, "stop_if_exists": False
#     }
#
#     ArgsTuple = namedtuple('ArgsTuple', argv.keys())
#
#     args_namedtuple = ArgsTuple(**argv)
#
#     train_model_main(args_namedtuple)
#     assert True
