#name: MolMIMGenerate
#description: MolMIM performs controlled generation, finding molecules with the right properties
#language: python
#input: string algorithm = "CMA-ES"
#input: int num_molecules = 30 {min: 1; max: 100}
#input: string property_name = "QED" {choices: ["QED", "plogP"]}
#input: bool minimize = False
#input: double min_similarity = 0.3 {min: 0; max: 1}
#input: int particles = 30 {min: 2; max: 1000}
#input: int iterations = 10 {min: 1; max: 1000}
#input: string smi = "[H][C@@]12Cc3c[nH]c4cccc(C1=C[C@H](NC(=O)N(CC)CC)CN2C)c34" {semType: Molecule}
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