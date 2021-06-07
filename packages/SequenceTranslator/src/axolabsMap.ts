export const axolabsMap: {[index: string]: {symbols: [string, string, string, string], color: string}} = {
  "RNA": {
    symbols: ["A", "C", "G", "U"],
    color: "purple"
  },
  "DNA": {
    symbols: ["dA", "dC", "dG", "dT"],
    color: "magenta"
  },
  "2'-Fluoro": {
    symbols: ["Af", "Cf", "Gf", "Uf"],
    color: "blue"
  },
  "2'-O-Methyl": {
    symbols: ["a", "c", "g", "u"],
    color: "silver"
  },
  "2'-O-MOE": {
    symbols: ["Am", "Cm", "Gm", "Tm"],
    color: "cyan"
  },
  "GNA": { // Glycol nucleic acid
    symbols: ["(GNA-A)", "(GNA-C)", "(GNA-G)", "(GNA-T)"],
    color: "yellow"
  },
  "LNA": {
    symbols: ["Ab", "Cb", "Gb", "Tb"],
    color: "aquamarine"
  },
  "UNA": { // Unlocked nucleotides
    symbols: ["Ao", "Co", "Go", "Uo"],
    color: "green"
  },
  "A": {
    symbols: ["A", "A", "A", "A"],
    color: "darkviolet"
  },
  "U": {
    symbols: ["U", "U", "U", "U"],
    color: "indigo"
  }
};