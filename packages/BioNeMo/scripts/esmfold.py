#name: esmfoldPython
#description: Predicts the 3D structure of a protein from its amino acid sequence
#language: python
#meta.cache: all
#meta.cache.invalidateOn: 0 * * * *
#input: string sequence
#input: string api_key
#output: string result

import requests
import json

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

try:
  response = session.post(invoke_url, headers=headers, json=payload)
  response.raise_for_status()
  response_body = response.json()
  response = response_body['pdbs'][0]
  result = {"success": True, "pdb": response}
except requests.exceptions.RequestException as e:
  result = {"success": False, "error": str(e)}

result = json.dumps(result)