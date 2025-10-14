#name: MolMIMGenerate
#description: MolMIM performs controlled generation, finding molecules with the right properties
#language: python
#input: string algorithm
#input: int num_molecules
#input: string property_name
#input: bool minimize
#input: double min_similarity
#input: int particles
#input: int iterations
#input: string smi
#input: string api_key
#output: string response

import requests

def get_query_url_and_headers(api_key):
  if api_key.strip():
    return (
      "https://health.api.nvidia.com/v1/biology/nvidia/molmim/generate",
      {"Authorization": f"Bearer {api_key}", "Accept": "application/json"},
    )
  return (
    "http://0.0.0.0:8000/generate",
    {"Content-Type": "application/json"},
  )

query_url, headers = get_query_url_and_headers(api_key)

data = {
  "algorithm": algorithm,
  "num_molecules": num_molecules,
  "property_name": property_name,
  "minimize": minimize,
  "min_similarity": min_similarity,
  "particles": particles,
  "iterations": iterations,
  "smi": smi,
}

response = requests.post(query_url, headers=headers, json=data).json()