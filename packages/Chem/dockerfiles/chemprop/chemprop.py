#!/usr/bin/env python

from utils import *
from engine import *
import rdkit
from rdkit import Chem

class ChemProp(Engine):
    type = 'Chemprop'
    description = '<a target="_blank" href="https://github.com/swansonk14/chemprop">Chemprop: Molecular Property Prediction</a>'

    def __init__(self):
        super(ChemProp, self).__init__()

    def train_impl(self, id: str, table: pd.DataFrame, predict: str, parameter_values: dict):
        tmp_dir = Engine.get_temporary_directory(id)
        table_path = os.path.join(tmp_dir, 'table.csv')
        table.columns = [n if n == predict else 'smiles' for n in table.columns]
        # Convert molblocks to SMILES if needed
        table['smiles'] = table['smiles'].apply(self.convert_to_smiles)
        ChemProp._save_table(table, table_path)
        params = [
            'chemprop',
            'train',
            '--data-path', table_path,
            '--output-dir', tmp_dir,
        ]
        params.extend(self.parameter_values_to_shell_params_string(parameter_values))
        index_of_type = params.index('--dataset-type')
        params[index_of_type] = '--task-type'
        log = call_process(params)
        model_blob = ''
        try:
            model_blob = open(ChemProp._get_model_blob_path(tmp_dir), 'rb').read()
        except:
            pass
        return model_blob, log

    def predict_impl(self, id: str, model_blob, table: pd.DataFrame, estimate_performance: bool = False):
        model_path = ChemProp._get_model_blob_path(Engine.get_temporary_directory(id, False))
        save_blob = not os.path.exists(model_path)
        tmp_dir = Engine.get_temporary_directory(id, save_blob)
        if save_blob:
            model_dir = os.path.dirname(model_path)
            shutil.rmtree(model_dir, ignore_errors=True)
            os.makedirs(model_dir)
            open(model_path, 'wb').write(model_blob)
        table_path = os.path.join(tmp_dir, 'table.csv')
        if not estimate_performance:
            table.columns = ['smiles']
            # Convert molblocks to SMILES if needed
            table['smiles'] = table['smiles'].apply(self.convert_to_smiles)
            ChemProp._save_table(table, table_path)
        predictions_path = os.path.join(tmp_dir, 'table_preds_0.csv')
        params = [
            'chemprop',
            'predict',
            '--test-path', table_path,
            '--model-path', model_path,
        ]
        call_process(params)
        prediction = pd.read_csv(predictions_path, na_values=['Invalid SMILES'])
        if not estimate_performance:
            prediction = prediction.filter([c for c in list(prediction) if c not in list(table)])
        return prediction

    @staticmethod
    def _get_model_blob_path(tmp_dir: str):
        return os.path.join(tmp_dir, 'model_0', 'best.pt')

    @staticmethod
    def _save_table(table: pd.DataFrame, table_path: str):
        table = table[['smiles'] + [col for col in table.columns if col != 'smiles']]
        table.to_csv(table_path, index=False)

    @staticmethod
    def convert_to_smiles(molblock):
        try:
            mol = Chem.MolFromMolBlock(molblock)
            return Chem.MolToSmiles(mol) if mol else molblock
        except:
            return molblock

    options = {
        'features_enabled': True,
        'features_sem_types': ['Molecule']
    }

    parameters = {

        'dataset_type': {
            'type': Types.STRING,
            'description': 'Type of dataset, e.g. classification or regression.'
                           'This determines the loss function used during training.',
            'choices': ['regression', 'classification'],#, 'classification', 'multiclass'],
            'category': 'general',
            'default_value': 'regression',
            'required': True
        },

         'metric': {
            'type': Types.STRING,
            'choices': ['auc', 'prc-auc', 'rmse', 'mae', 'mse', 'r2', 'accuracy', 'cross_entropy'],
            'description': 'Metric to use during evaluation.'
                           'Note: Does NOT affect loss function used during training'
                           '(loss is determined by the `dataset_type` argument).'
                           'Note: Defaults to "auc" for classification and "rmse" for regression.',
            'category': 'general',
            'default_value': None
        },

        'multiclass_num_classes': {
            'type': Types.INT,
            'description': 'Number of classes when running multiclass classification',
            'category': 'general',
            'default_value': 3
        },

        'num_folds': {
            'type': Types.INT,
            'description': 'Number of folds when performing cross validation',
            'category': 'general',
            'default_value': 1
        },

        'data_seed': {
            'type': Types.INT,
            'description': 'Random seed to use when splitting data into train/val/test sets.'
                           'When `num_folds` > 1, the first fold uses this seed and all'
                           'subsequent folds add 1 to the seed.',
            'category': 'general',
            'default_value': 0
        },

        'split_sizes': {
            'type': Types.LIST,
            'description': 'Split proportions for train/validation/test sets',
            'category': 'general',
            'default_value': [0.8, 0.1, 0.1]
        },

        'split_type': {
            'type': Types.STRING,
            'choices': ['random', 'scaffold_balanced', 'predetermined', 'crossval', 'index_predetermined'],
            'description': 'Method of splitting the data into train/val/test',
            'category': 'general',
            'default_value': 'random'
        },

        'activation': {
            'type': Types.STRING,
            'description': 'Activation function',
            'choices': ['ReLU', 'LeakyReLU', 'PReLU', 'tanh', 'SELU', 'ELU'],
            'category': 'model',
            'default_value': 'ReLU',
        },

        'atom_messages': {
            'type': Types.BOOL,
            'description': 'Use messages on atoms instead of messages on bonds',
            'category': 'model',
            'default_value': False
        },

        'message_bias': {
            'type': Types.BOOL,
            'description': 'Whether to add bias to linear layers',
            'category': 'model',
            'default_value': False
        },

        'ensemble_size': {
            'type': Types.INT,
            'description': 'Number of models in ensemble',
            'category': 'model',
            'default_value': 1
        },

        'message_hidden_dim': {
            'type': Types.INT,
            'description': 'Dimensionality of hidden layers in MPN',
            'category': 'model',
            'default_value': 300
        },

        'depth': {
            'type': Types.INT,
            'description': 'Number of message passing steps',
            'category': 'model',
            'default_value': 3
        },

        'dropout': {
            'type': Types.FLOAT,
            'description': 'Dropout probability',
            'category': 'model',
            'default_value': 0.0
        },

        'undirected': {
            'type': Types.BOOL,
            'description': 'Undirected edges (always sum the two relevant bond vectors)',
            'category': 'model',
            'default_value': False
        },

        'ffn_hidden_dim': {
            'type': Types.INT,
            'description': 'Hidden dim for higher-capacity FFN (defaults to hidden_size)',
            'category': 'model',
            'default_value': 300
        },

        'ffn_num_layers': {
            'type': Types.INT,
            'description': 'Number of layers in FFN after MPN encoding',
            'category': 'model',
            'default_value': 2
        },

        'epochs': {
            'type': Types.INT,
            'description': 'Number of epochs to run',
            'category': 'training',
            'default_value': 50
        },

        'batch_size': {
            'type': Types.INT,
            'description': 'Batch size',
            'category': 'training',
            'default_value': 64
        },

        'warmup_epochs': {
            'type': Types.FLOAT,
            'description': 'Number of epochs during which learning rate increases linearly from'
                           'init_lr to max_lr. Afterwards, learning rate decreases exponentially'
                           'from max_lr to final_lr.',
            'category': 'training',
            'default_value': 2.0
        },

        'init_lr': {
            'type': Types.FLOAT,
            'description': 'Initial learning rate',
            'category': 'training',
            'default_value': 0.0001
        },

        'max_lr': {
            'type': Types.FLOAT,
            'description': 'Maximum learning rate',
            'category': 'training',
            'default_value': 0.001
        },

        'final_lr': {
            'type': Types.FLOAT,
            'description': 'Final learning rate',
            'category': 'training',
            'default_value': 0.0001
        },

        'no_descriptor_scaling': {
            'type': Types.BOOL,
            'description': 'Turn off scaling of features',
            'category': 'training',
            'default_value': False
        }
    }
