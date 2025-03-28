#name: esmfold
#description: Predicts the 3D structure of a protein from its amino acid sequence
#language: python
#meta.cache: all
#meta.cache.invalidateOn: 0 * * * *
#input: string sequence
#input: string api_key
#output: string response

import requests

invoke_url = "https://health.api.nvidia.com/v1/biology/nvidia/esmfold"

headers = {
    "Authorization": f"Bearer {api_key}",
    "Accept": "application/json",
}

payload = {
  "sequence": sequence
}

# re-use connections
session = requests.Session()

response = session.post(invoke_url, headers=headers, json=payload)

response.raise_for_status()
response_body = response.json()
response = response_body['pdbs'][0]