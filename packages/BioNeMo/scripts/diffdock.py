#name: diffdockPython
#friendlyName: DiffDock
#description: Predicts the 3D structure of how a small molecule ligand docks into a protein
#language: python
#input: string protein [Protein receptor structure (PDB text)]
#input: string ligand [Ligand structure (SDF/molblock or SMILES)]
#input: int num_poses [Number of docking poses to generate]
#input: string api_key [NVIDIA BioNeMo API key (blank uses the local endpoint)]
#output: blob result

import requests

def get_query_url_and_headers(api_key):
  if api_key.strip():
    return (
      "https://health.api.nvidia.com/v1/biology/mit/diffdock",
      {"Authorization": f"Bearer {api_key}", "Accept": "application/json"},
    )
  return (
    "http://0.0.0.0:8000/molecular-docking/diffdock/generate",
    {"Content-Type": "application/json"},
  )

query_url, headers = get_query_url_and_headers(api_key)

data = {
  "ligand": ligand,
  "ligand_file_type": "sdf",
  "protein": protein,
  "num_poses": num_poses,
  "time_divisions": 20,
  "steps": 18,
  "save_trajectory": False,
  "is_staged": False
}

response = requests.post(query_url, headers=headers, json=data)

result = response.text.encode("utf-8")