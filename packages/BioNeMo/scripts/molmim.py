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
#output: string response
import requests

invoke_url = "https://health.api.nvidia.com/v1/biology/nvidia/molmim/generate"

headers = {
    "Authorization": "Bearer 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",
    "Accept": "application/json",
}

payload = {
  "algorithm": algorithm,
  "num_molecules": num_molecules,
  "property_name": property_name,
  "minimize": minimize,
  "min_similarity": min_similarity,
  "particles": particles,
  "iterations": iterations,
  "smi": smi
}

session = requests.Session()
response = session.post(invoke_url, headers=headers, json=payload)
response.raise_for_status()
response_body = response.json()
response = response_body