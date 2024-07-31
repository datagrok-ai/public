#!/usr/bin/env python

import os
import math
import shutil
import grok
import numpy as np
import pandas as pd
from sklearn import metrics
from utils import Types

class Engine(object):
    type = None
    description = None
    options = {}
    parameters = {}

    def __init__(self):
        pass

    def train(self, id: str, table_server_url: str, table_token: str, predict: str, parameter_values: dict):
        """
        Trains model.

        :param id: New model ID.
        :param table_server_url: Table server URL.
        :param table_token: Table token.
        :param predict: Predict column name.
        :param parameter_values: Dictionary with model specific parameter values.
        :return: Model blob, log and performance.
        """
        grok.grok_init(table_server_url)
        table = grok.grok_read(table_token)
        model_blob, log = self.train_impl(id, table, predict, parameter_values)
        # TODO Use validation dataset for performance estimation
        performance = self.estimate_performance_impl(id, model_blob, table, predict,
                                                     True) if model_blob == '' else self.estimate_performance_impl(id,
                                                                                                                   model_blob,
                                                                                                                   table,
                                                                                                                   predict,
                                                                                                                   False)
        return model_blob, log, performance

    def predict(self, id: str, model_blob, table_server_url: str, table_token: str) -> pd.DataFrame:
        """
        Predicts outcome.

        :param id: Model ID.
        :param model_blob: Model blob.
        :param table_server_url: Table server URL.
        :param table_token: Table token.
        :return: Outcome, dataframe.
        """
        grok.grok_init(table_server_url)
        table = grok.grok_read(table_token)
        return self.predict_impl(id, model_blob, table)

    def estimate_performance(self, id: str, model_blob, table_server_url: str, table_token: str, predict: str) -> dict:
        """
        Predicts outcome.

        :param id: Model ID.
        :param model_blob: Model blob.
        :param table_server_url: Table server URL.
        :param table_token: Table token.
        :param predict: Predict column name.
        :return: Dictionary with model performance.
        """
        grok.grok_init(table_server_url)
        table = grok.grok_read(table_token)
        performance = self.estimate_performance_impl(id, model_blob, table, predict)
        return performance

    def train_impl(self, id: str, table: pd.DataFrame, predict: str, parameter_values: dict):
        """
        Trains model, engine specific implementation.

        :param id: New model ID.
        :param table: Table.
        :param predict: Predict column name.
        :param parameter_values: Dictionary with model specific parameter values.
        :return: Model blob and train log.
        """
        return None, None

    def predict_impl(self, id: str, model_blob, table: pd.DataFrame) -> pd.DataFrame:
        """
        Predicts outcome, engine specific implementation.

        :param id: Model ID.
        :param model_blob: Model blob.
        :param table: Table.
        :return: Outcome, dataframe.
        """
        return

    def estimate_performance_impl(self, id: str, model_blob, table: pd.DataFrame, predict: str,
                                  cancelled: bool) -> dict:
        """
        Predicts outcome, engine specific implementation.

        :param id: Model ID.
        :param model_blob: Model blob.
        :param table: Table.
        :param predict: Predict column name.
        :return: Dictionary with model performance.
        """
        metrics_names = ['mse', 'rmse', 'corr', 'nobs', 'r2', 'logloss', 'auc']
        if cancelled:
            return dict(zip(metrics_names, ['null'] * 7))
        predicted = self.predict_impl(id, model_blob, table, estimate_performance=True)
        y_true = np.array(table[predict])
        y_pred = np.array(predicted[predict])

        index = (y_true != np.array(None)) & (y_pred != np.array(None)) & \
                (~np.isnan(y_true)) & (~np.isnan(y_pred))
        y_true = y_true[index]
        y_pred = y_pred[index]

        try:
            fpr, tpr = metrics.roc_curve(y_true, y_pred, pos_label=2)
            auc = metrics.auc(fpr, tpr)
        except:
            auc = None

        try:
            log_loss = metrics.log_loss(y_true, y_pred)
        except:
            log_loss = None

        r2 = metrics.r2_score(y_true, y_pred)

        return {
            'mse': metrics.mean_absolute_error(y_true, y_pred),
            'rmse': metrics.mean_squared_error(y_true, y_pred),
            'corr': math.sqrt(abs(r2)),
            'nobs': len(y_true),
            'r2': r2,
            'logloss': log_loss,
            'auc': auc,
        }

    def parameter_values_to_shell_params_string(self, parameter_values: dict) -> list:
        """
        Converts parameters values into list of shell parameters strings.

        :param parameter_values: Parameter values.
        :return: List of shell parameters strings.
        """
        params = []
        for param in parameter_values.keys():
            value = parameter_values[param]
            if value != self.parameters[param]['default_value'] or \
                    ('required' in self.parameters[param] and self.parameters[param]['required']):
                if self.parameters[param]['type'] == Types.LIST:
                    param = param.replace('_', '-')
                    params.append('--' + param)
                    for number in value:
                        params.append(str(number))
                    continue
                str_value = str(value)
                param = param.replace('_', '-')
                params.append('--' + param)
                params.append(str_value)
        return params

    @staticmethod
    def get_temporary_directory(id: str, create_new=True):
        """
        Gets temporary directory for train or predict.

        :param id: Model ID.
        :param create_new: Creates new folder and delete previous, if exists.
        :return: Path to temporary directory.
        """
        tmp_dir = os.path.join('tmp', 'models', id)
        if create_new:
            shutil.rmtree(tmp_dir, ignore_errors=True)
            os.makedirs(tmp_dir)
        return tmp_dir
