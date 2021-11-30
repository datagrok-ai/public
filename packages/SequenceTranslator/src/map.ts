export const map: {[synthesizer: string]: {[technology: string]: {[code: string]: {"name": string, "weight": number, "normalized": string, "SMILES"?: string}}}} = {
  "BioSpring Codes": {
    "For ASO Gapmers": {
      "5": {
        "name": "2'MOE-5Me-rU",
        "weight": 378.27,
        "normalized": "rU"
      },
      "6": {
        "name": "2'MOE-rA",
        "weight": 387.29,
        "normalized": "rA"
      },
      "7": {
        "name": "2'MOE-5Me-rC",
        "weight": 377.29,
        "normalized": "rC"
      },
      "8": {
        "name": "2'MOE-rG",
        "weight": 403.28,
        "normalized": "rG"
      },
      "9": {
        "name": "5-Methyl-dC",
        "weight": 303.21,
        "normalized": "dC"
      },
      "*": {
        "name": "ps linkage",
        "weight": 16.07,
        "normalized": ""
      },
      "A": {
        "name": "Adenine",
        "weight": 313.21,
        "normalized": "dA"
      },
      "C": {
        "name": "Cytosine",
        "weight": 289.18,
        "normalized": "dC"
      },
      "G": {
        "name": "Guanine",
        "weight": 329.21,
        "normalized": "dG"
      },
      "T": {
        "name": "Tyrosine",
        "weight": 304.2,
        "normalized": "dT"
      }
    },
    "For 2'-OMe and 2'-F modified siRNA": {
      "1": {
        "name": "2'-fluoro-U",
        "weight": 308.16,
        "normalized": "rU"
      },
      "2": {
        "name": "2'-fluoro-A",
        "weight": 331.2,
        "normalized": "rA"
      },
      "3": {
        "name": "2'-fluoro-C",
        "weight": 307.18,
        "normalized": "rC"
      },
      "4": {
        "name": "2'-fluoro-G",
        "weight": 347.19,
        "normalized": "rG"
      },
      "5": {
        "name": "2'OMe-rU",
        "weight": 320.2,
        "normalized": "rU"
      },
      "6": {
        "name": "2'OMe-rA",
        "weight": 343.24,
        "normalized": "rA"
      },
      "7": {
        "name": "2'OMe-rC",
        "weight": 319.21,
        "normalized": "rC"
      },
      "8": {
        "name": "2'OMe-rG",
        "weight": 359.24,
        "normalized": "rG"
      },
      "*": {
        "name": "ps linkage",
        "weight": 16.07,
        "normalized": ""
      }
    }
  },
  "Axolabs Codes": {
    "For 2'-OMe and 2'-F modified siRNA": {
      "Uf": {
        "name": "2'-fluoro-U",
        "weight": 308.16,
        "normalized": "rU"
      },
      "Af": {
        "name": "2'-fluoro-A",
        "weight": 331.2,
        "normalized": "rA"
      },
      "Cf": {
        "name": "2'-fluoro-C",
        "weight": 307.18,
        "normalized": "rC"
      },
      "Gf": {
        "name": "2'-fluoro-G",
        "weight": 347.19,
        "normalized": "rG"
      },
      "u": {
        "name": "2'OMe-rU",
        "weight": 320.2,
        "normalized": "rU"
      },
      "a": {
        "name": "2'OMe-rA",
        "weight": 343.24,
        "normalized": "rA"
      },
      "c": {
        "name": "2'OMe-rC",
        "weight": 319.21,
        "normalized": "rC"
      },
      "g": {
        "name": "2'OMe-rG",
        "weight": 359.,
        "normalized": "rG"
      },
      "s": {
        "name": "ps linkage",
        "weight": 16.07,
        "normalized": ""
      }
    }
  },
  "Janssen GCRS Codes": {
    "For ASO Gapmers": {
      "moeT": {
        "name": "2'MOE-5Me-rU",
        "weight": 378.27,
        "normalized": "rU",
        "SMILES": "OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))[C@H](OCCOC)[C@@H]1O"
      },
      "moeA": {
        "name": "2'MOE-rA",
        "weight": 387.29,
        "normalized": "rA",
        "SMILES": "OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OCCOC)[C@@H]1O"
      },
      "moe5mC": {
        "name": "2'MOE-5Me-rC",
        "weight": 377.29,
        "normalized": "rC",
        "SMILES": "OC[C@H]1O[C@@H](N2C=C(C)C(N)=NC2(=O))[C@H](OCCOC)[C@@H]1O"
      },
      "(5m)moeC": {
        "name": "2'MOE-5Me-rC",
        "weight": 377.29,
        "normalized": "rC",
        "SMILES": "OC[C@H]1O[C@@H](N2C=C(C)C(N)=NC2(=O))[C@H](OCCOC)[C@@H]1O"
      },
      "moeG": {
        "name": "2'MOE-rG",
        "weight": 403.28,
        "normalized": "rG",
        "SMILES": "OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OCCOC)[C@@H]1O"
      },
      "5mC": {
        "name": "5-Methyl-dC",
        "weight": 303.28,
        "normalized": "dC",
        "SMILES": "OC[C@H]1O[C@@H](N2C=C(C)C(N)=NC2(=O))C[C@@H]1O"
      },
      "(5m)C": {
        "name": "5-Methyl-dC",
        "weight": 303.28,
        "normalized": "dC",
        "SMILES": "OC[C@H]1O[C@@H](N2C=C(C)C(N)=NC2(=O))C[C@@H]1O"
      },
      "ps": {
        "name": "ps linkage",
        "weight": 16.07,
        "normalized": "",
        "SMILES": "OP(=O)(O)S"
      },
      "A": {
        "name": "Adenine",
        "weight": 313.21,
        "normalized": "dA",
        "SMILES": "OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)C[C@@H]1O"
      },
      "dA": {
        "name": "Adenine",
        "weight": 313.21,
        "normalized": "dA",
        "SMILES": "OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)C[C@@H]1O"
      },
      "C": {
        "name": "Cytosine",
        "weight": 289.18,
        "normalized": "dC",
        "SMILES": "OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))C[C@@H]1O"
      },
      "dC": {
        "name": "Cytosine",
        "weight": 289.18,
        "normalized": "dC",
        "SMILES": "OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))C[C@@H]1O"
      },
      "G": {
        "name": "Guanine",
        "weight": 329.21,
        "normalized": "dG",
        "SMILES": "OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)C)[C@@H]1O"
      },
      "dG": {
        "name": "Guanine",
        "weight": 329.21,
        "normalized": "dG",
        "SMILES": "OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)C)[C@@H]1O"
      },
      "T": {
        "name": "Tyrosine",
        "weight": 304.2,
        "normalized": "dT",
        "SMILES": "OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))C[C@@H]1O"
      },
      "dT": {
        "name": "Tyrosine",
        "weight": 304.2,
        "normalized": "dT",
        "SMILES": "OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))C[C@@H]1O"
      },
      "rA": {
        "name": "Adenine",
        "weight": 329.21,
        "normalized": "rA",
        "SMILES": "OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](O)[C@@H]1O"
      },
      "rC": {
        "name": "Cytosine",
        "weight": 305.18,
        "normalized": "rC",
        "SMILES": "OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](O)[C@@H]1O"
      },
      "rG": {
        "name": "Guanine",
        "weight": 345.21,
        "normalized": "rG",
        "SMILES": "OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](O)[C@@H]1O"
      },
      "rU": {
        "name": "Uracil",
        "weight": 306.17,
        "normalized": "rU",
        "SMILES": "OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](O)[C@@H]1O"
      }
    },
    "For 2'-OMe and 2'-F modified siRNA": {
      "fU": {
        "name": "2'-fluoro-U",
        "weight": 308.16,
        "normalized": "rU",
        "SMILES": "OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1O"
      },
      "fA": {
        "name": "2'-fluoro-A",
        "weight": 331.2,
        "normalized": "rA",
        "SMILES": "OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1O"
      },
      "fC": {
        "name": "2'-fluoro-C",
        "weight": 307.18,
        "normalized": "rC",
        "SMILES": "OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1O"
      },
      "fG": {
        "name": "2'-fluoro-G",
        "weight": 347.19,
        "normalized": "rG",
        "SMILES": "OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1O"
      },
      "mU": {
        "name": "2'OMe-rU",
        "weight": 320.2,
        "normalized": "rU",
        "SMILES": "OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1O"
      },
      "mA": {
        "name": "2'OMe-rA",
        "weight": 343.24,
        "normalized": "rA",
        "SMILES": "OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1O"
      },
      "mC": {
        "name": "2'OMe-rC",
        "weight": 319.21,
        "normalized": "rC",
        "SMILES": "OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1O"
      },
      "mG": {
        "name": "2'OMe-rG",
        "weight": 359.24,
        "normalized": "rG",
        "SMILES": "OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1O"
      }
    }
  }
};