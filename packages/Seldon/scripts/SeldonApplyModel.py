#name: SeldonApplyModelPy
#language: python
#input: string seldonHost
#input: string seldonOIDCServer
#input: string seldonClientID
#input: string seldonUser
#input: string seldonPassword
#input: string seldonNamespace
#input: string seldonDeploymentName
#input: dataframe seldonInput
#output: dataframe seldonResult
#output: string seldonErrorState

# Algorithm

# 1. By the namespace and deployment names, fetch the corresponding model's metadata
# 2. Disregard of columns order, match (case-insensitive, treat dash and space equally)
#    columns by name to the columns from the model metadata. For now, ignore the type (we
#    will have to pass the type info to the function separately)
# 3. Iterpret the result, either as a successful set out outputs, or an error state
# 4. Outputs be a dataframe with names, formed based on the Seldon outputs types
# 5. On error, prepare a human-readible interpretation of this error; if no error, seldonErrorState == ""

seldonErrorState = ""

import seldon_deploy_sdk
from seldon_deploy_sdk.auth import OIDCAuthenticator
from seldon_deploy_sdk import Configuration, ApiClient, SeldonDeploymentsApi
from seldon_deploy_sdk.rest import ApiException
from urllib3.exceptions import InsecureRequestWarning
from seldon_core.seldon_client import SeldonClient

import numpy as np
import warnings
import contextlib
import requests
import time
import json
import pandas as pd

old_merge_environment_settings = requests.Session.merge_environment_settings

@contextlib.contextmanager
def no_ssl_verification():
  opened_adapters = set()
  def merge_environment_settings(self, url, proxies, stream, verify, cert):
    # Verification happens only once per connection so we need to close
    # all the opened adapters once we're done. Otherwise, the effects of
    # verify=False persist beyond the end of this context manager.
    opened_adapters.add(self.get_adapter(url))
    settings = old_merge_environment_settings(self, url, proxies, stream, verify, cert)
    settings['verify'] = False
    return settings
  requests.Session.merge_environment_settings = merge_environment_settings
  try:
    with warnings.catch_warnings():
      warnings.simplefilter('ignore', InsecureRequestWarning)
      yield
  finally:
    requests.Session.merge_environment_settings = old_merge_environment_settings
    for adapter in opened_adapters:
      try:
        adapter.close()
      except:
        pass

with no_ssl_verification():

  config = Configuration()
  config.host = seldonHost
  config.oidc_server = seldonOIDCServer
  config.oidc_client_id = seldonClientID
  config.verify_ssl = False
  config.username = seldonUser
  config.password = seldonPassword

  def auth():
    auth = OIDCAuthenticator(config)
    config.access_token = auth.authenticate()
    api_client = ApiClient(config)
    return api_client

  
  # 1. Get the model's metadata object
  
  deployment_api = seldon_deploy_sdk.ModelMetadataServiceApi(auth())
  responseDeploymentToModel = deployment_api.model_metadata_service_list_runtime_metadata_for_model(deployment_namespace = seldonNamespace, deployment_name = seldonDeploymentName)
  responseModelMetadata = deployment_api.model_metadata_service_list_model_metadata(uri=responseDeploymentToModel.runtime_metadata[0].model_uri)
  metadataObject = responseModelMetadata.models[0].prediction_schema
  
  # 2. Form input and output columns names as sequences
  
  inputColumnNames = []
  inputColumnTypes = []
  inputColumnMeta = []
  requestsObject = metadataObject.requests
  for i in range(0, len(requestsObject)):
    inputColumnNames += [requestsObject[i].name]
    columnType = requestsObject[i].type
    inputColumnTypes += [columnType]
    inputColumnMeta += [requestsObject[i].category_map]
  
  responsesObject = metadataObject.responses
  outputColumnNames = []
  for i in range(0, len(responsesObject)):
    responseObject = responsesObject[i]
    if responseObject.type == "PROBA":
      prefix = responseObject.name
      schema = responseObject.schema
      for j in range(len(schema)):
        outputColumnNames += [prefix + ":" + schema[j].name]

  # 3. Predict
  deployment_api = SeldonDeploymentsApi(auth())
  actualInput = seldonInput[inputColumnNames].values.tolist()
  for j in range(0, len(inputColumnTypes)):
    if inputColumnTypes[j] == "CATEGORICAL":
      categoriesObject = inputColumnMeta[j]
      categoryNamesToInput = {}
      for key, value in categoriesObject.items():
        categoryNamesToInput[value] = key
      for i in range(0, len(actualInput)):
        value = actualInput[i][inputColumnNames[j]]
        actualInput[i][inputColumnNames[j]] = categoryNamesToInput[value]
  obj = { "data" : { "ndarray": actualInput } }
  api_instance = seldon_deploy_sdk.PredictApi(seldon_deploy_sdk.ApiClient(config))
  api_response = api_instance.predict_seldon_deployment(name=seldonDeploymentName, namespace=seldonNamespace, prediction=obj)
  seldonResult = pd.DataFrame(api_response['data']['ndarray'], columns = outputColumnNames)