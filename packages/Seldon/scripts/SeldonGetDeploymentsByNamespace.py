#name: SeldonGetDeploymentsByNamespacePy
#language: python
#input: string seldonHost
#input: string seldonOIDCServer
#input: string seldonClientID
#input: string seldonUser
#input: string seldonPassword
#input: string seldonNamespace
#output: dataframe seldonDeployments

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
      
  # Obtain a list of deployments in the namespace
  
  deployment_api = SeldonDeploymentsApi(auth())
  api_response = deployment_api.list_seldon_deployments(seldonNamespace)
  itemsList = api_response.items
  deploymentsNames = []
  for i in range(0, len(itemsList)):
    deploymentsNames += [itemsList[i].spec.name]
  seldonDeployments = pd.DataFrame(deploymentsNames, columns=['seldonDeployments'])