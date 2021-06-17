const rnaColor = "rgb(255,230,153)";
export const axolabsMap: {[index: string]: {symbols: [string, string, string, string], color: string}} = {
  "RNA": {
    symbols: ["A", "C", "G", "U"],
    color: rnaColor
  },
  "DNA": {
    symbols: ["dA", "dC", "dG", "dT"],
    color: "rgb(197,224,180)"
  },
  "2'-Fluoro": {
    symbols: ["Af", "Cf", "Gf", "Uf"],
    color: "rgb(68,114,196)"
  },
  "2'-O-Methyl": {
    symbols: ["a", "c", "g", "u"],
    color: "rgb(166,166,166)"
  },
  "2'-O-MOE": {
    symbols: ["Am", "Cm", "Gm", "Tm"],
    color: "rgb(112,48,160)"
  },
  "GNA": { // Glycol nucleic acid
    symbols: ["(GNA-A)", "(GNA-C)", "(GNA-G)", "(GNA-T)"],
    color: "rgb(255,192,0)"
  },
  "LNA": {
    symbols: ["Ab", "Cb", "Gb", "Tb"],
    color: "rgb(54,229,238)"
  },
  "UNA": { // Unlocked nucleotides
    symbols: ["Ao", "Co", "Go", "Uo"],
    color: "rgb(255,192,0)"
  },
  "A": {
    symbols: ["a", "a", "a", "a"],
    color: rnaColor
  },
  "C": {
    symbols: ["c", "c", "c", "c"],
    color: rnaColor
  },
  "G": {
    symbols: ["g", "g", "g", "g"],
    color: rnaColor
  },
  "U": {
    symbols: ["u", "u", "u", "u"],
    color: rnaColor
  },
  "X-New": {
    symbols: ["X", "X", "X", "X"],
    color: "rgb(108,0,0)"
  },
  "Y-New": {
    symbols: ["Y", "Y", "Y", "Y"],
    color: "rgb(210,146,146)"
  },
  "Z-New": {
    symbols: ["Z", "Z", "Z", "Z"],
    color: "rgb(84,130,53)"
  },
  "InvAbasic": {
    symbols: ["(invabasic)", "(invabasic)", "(invabasic)", "(invabasic)"],
    color: "rgb(240,62,202)"
  },
  "5'-vinylps": {
    symbols: ["(vinu)", "(vinu)", "(vinu)", "(vinu)"],
    color: "rgb(0,0,139)"
  }
};