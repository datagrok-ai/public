#!/usr/bin/env python

from utils import *
from engine import *

class ChemProp(Engine):
    type = 'Chemprop'
    description = '<a target="_blank" href="https://github.com/swansonk14/chemprop">Chemprop: Molecular Property Prediction</a>'

    def __init__(self):
        super(ChemProp, self).__init__()

    def train_impl(self, id: str, table: pd.DataFrame, predict: str, parameter_values: dict):
        tmp_dir = Engine.get_temporary_directory(id)
        table_path = os.path.join(tmp_dir, 'table.csv')
        table.columns = [n if n == predict else 'smiles' for n in table.columns]
        ChemProp._save_table(table, table_path)
        params = [
            'chemprop_train',
            '--data_path', table_path,
            '--save_dir', tmp_dir]
        params.extend(self.parameter_values_to_shell_params_string(parameter_values))
        log = call_process(params)
        model_blob = ''
        try:
            model_blob = open(ChemProp._get_model_blob_path(tmp_dir), 'rb').read()
        except:
            pass
        return model_blob, log

    def predict_impl(self, id: str, model_blob, table: pd.DataFrame):
        model_path = ChemProp._get_model_blob_path(Engine.get_temporary_directory(id, False))
        save_blob = not os.path.exists(model_path)
        tmp_dir = Engine.get_temporary_directory(id, save_blob)
        if save_blob:
            model_dir = os.path.dirname(model_path)
            shutil.rmtree(model_dir, ignore_errors=True)
            os.makedirs(model_dir)
            open(model_path, 'wb').write(model_blob)
        table_path = os.path.join(tmp_dir, 'table.csv')
        save_table = not os.path.exists(table_path)
        if save_table:
            table.columns = ['smiles']
            ChemProp._save_table(table, table_path)
        predictions_path = os.path.join(tmp_dir, 'predictions.csv')
        params = [
            'chemprop_predict',
            '--test_path', table_path,
            '--checkpoint_path', model_path,
            '--preds_path', predictions_path
        ]
        call_process(params)
        prediction = pd.read_csv(predictions_path, na_values=['Invalid SMILES'])
        if save_table:
            prediction = prediction.filter([c for c in list(prediction) if c not in list(table)])
        shutil.rmtree(tmp_dir, ignore_errors=True)
        return prediction

    @staticmethod
    def _get_model_blob_path(tmp_dir: str):
        return os.path.join(tmp_dir, 'fold_0', 'model_0', 'model.pt')

    @staticmethod
    def _save_table(table: pd.DataFrame, table_path: str):
        table = table[['smiles'] + [col for col in table.columns if col != 'smiles']]
        table.to_csv(table_path, index=False)

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

        'log_frequency': {
            'type': Types.INT,
            'description': 'The number of batches between each logging of the training loss',
            'category': 'general',
            'default_value': 10
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

        'no_cache': {
            'type': Types.BOOL,
            'description': 'Turn off caching mol2graph computation',
            'category': 'general',
            'default_value': False
        },

        'num_folds': {
            'type': Types.INT,
            'description': 'Number of folds when performing cross validation',
            'category': 'general',
            'default_value': 1
        },

        'seed': {
            'type': Types.INT,
            'description': 'Random seed to use when splitting data into train/val/test sets.'
                           'When `num_folds` > 1, the first fold uses this seed and all'
                           'subsequent folds add 1 to the seed.',
            'category': 'general',
            'default_value': 0
        },

        'show_individual_scores': {
            'type': Types.BOOL,
            'description': 'Show all scores for individual targets, not just average, at the end',
            'category': 'general',
            'default_value': False
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

        'test': {
            'type': Types.BOOL,
            'description': 'Whether to skip training and only test the model',
            'category': 'general',
            'default_value': False
        },

        'use_compound_names': {
            'type': Types.BOOL,
            'description': 'Use when test data file contains compound names in addition to SMILES strings',
            'category': 'general',
            'default_value': False
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

        'bias': {
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

        'hidden_size': {
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

        'ffn_hidden_size': {
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
            'default_value': 30
        },

        'batch_size': {
            'type': Types.INT,
            'description': 'Batch size',
            'category': 'training',
            'default_value': 50
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

        'no_features_scaling': {
            'type': Types.BOOL,
            'description': 'Turn off scaling of features',
            'category': 'training',
            'default_value': False
        },

        'max_data_size': {
            'type': Types.INT,
            'description': 'Maximum number of data points to load',
            'category': 'predict',
            'default_value': None
        }
    }
