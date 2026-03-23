const CHEMICAL_MIXTURE_SEM_TYPE = 'ChemicalMixture';
const MIX_FILE_VERSION = 'mixfileVersion';

class ChemPackageDetectors extends DG.Package {
  static likelyNames = [
    'structure', 'mol', 'molecule', 'smiles', 'rdkit',
    'canonical_smiles', 'core', 'scaffold'];

  // typical r-group names, such as "R1" or "R100-family"
  static likelyRegexps = [/[r|R][0-9]+(\W|$)/];

  static likelyChemicalName(s) {
    return ChemPackageDetectors.likelyNames.some((likelyName) => s.includes(likelyName)) ||
      ChemPackageDetectors.likelyRegexps.some((regexp) => regexp.test(s));
  }

  static validSmilesChars = new Uint32Array([0, 805305144, 1073741823, 1073741822, 0, 0, 0, 0]);
  static validFirstSmilesChars = new Uint32Array([0, 256, 134857308, 573448, 0, 0, 0, 0]);
  static numberOfMolsToCheck = 100;
  static numberOfMolsToCheckRdKit = 10;

  static validChar(pos) {
    return (ChemPackageDetectors.validSmilesChars[Math.floor(pos / 0x20)] & (1 << (pos & 0x1f))) !== 0;
  }

  static validFirstChar(pos) {
    return (ChemPackageDetectors.validFirstSmilesChars[Math.floor(pos / 0x20)] & (1 << (pos & 0x1f))) !== 0;
  }

  static likelyValidSmiles(s) {
    if (!ChemPackageDetectors.validFirstChar(s[0].charCodeAt(0)))
      return false;
    for (let j = 1; j < s.length; j++) {
      const ch = s[j];
      if (ch === ' ' && j + 1 < s.length && s[j + 1] === '|')
        continue;
      if (!ChemPackageDetectors.validChar(ch.charCodeAt(0)))
        return false;
    }
    return true;
  }

  //name: detectMolecules
  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  //meta.skipTest: GROK-17630
  detectMolecules(col) {
    if (DG.Detector.sampleCategories(col, (s) => s.includes('M  END'), 1)) {
      col.meta.units = DG.UNITS.Molecule.MOLBLOCK;
      col.meta.cellRenderer = 'Molecule';
      return DG.SEMTYPE.MOLECULE;
    }

    const lowerCaseName = col.name.toLowerCase();
    const likelyMolName = ChemPackageDetectors.likelyChemicalName(lowerCaseName);
    const minUnique = likelyMolName ? 1 : 3;
    let longest = '';
    try {
      longest = col.aggregate('longest') ?? '';
    } catch (x) {}
    if (!likelyMolName && longest.length < 5)
      return null;

    // temporary workaround to make it understand r-groups like "CC[*:2]"
    if (likelyMolName) {
      grok.functions.call('Chem:detectSmiles', {col: col, min: minUnique}).then(() => {});
      return null;
    }

    if (DG.Detector.sampleCategories(col, ChemPackageDetectors.likelyValidSmiles, minUnique, 10, 0.8))
      grok.functions.call('Chem:detectSmiles', {col: col, min: minUnique}).then(() => {});
  }

  static likelyReactionNames = ['reaction', 'rxn'];

  // Check if a string looks like a reaction SMILES/SMARTS.
  // Supports three formats:
  //   "reactants>>products"           — standard reaction SMILES
  //   "reactants>agents>products"      — SMIRKS with agent/reagent field
  //   "A>>B>>C>>D"                     — multi-step reaction chain
  static isLikelyReaction(s) {
    if (s.includes('M  END')) return false;

    // Standard/multi-step: contains ">>"
    if (s.includes('>>')) {
      const idx = s.indexOf('>>');
      const left = s.substring(0, idx);
      const rightFull = s.substring(idx + 2);
      // For multi-step, validate just the first product group
      const right = rightFull.split('>>')[0];
      if (left.length === 0 || right.length === 0) return false;
      // Validate the last reactant and first product molecule
      const lastReactant = left.split('.').pop();
      const firstProduct = right.split('.')[0];
      return lastReactant.length > 0 && firstProduct.length > 0 &&
        ChemPackageDetectors.likelyValidSmiles(lastReactant) &&
        ChemPackageDetectors.likelyValidSmiles(firstProduct);
    }

    // SMIRKS format: "reactants>agents>products" (exactly two single ">")
    const parts = s.split('>');
    if (parts.length === 3 && parts[0].length > 0 && parts[2].length > 0) {
      const lastReactant = parts[0].split('.').pop();
      const firstProduct = parts[2].split('.')[0];
      return lastReactant.length > 0 && firstProduct.length > 0 &&
        ChemPackageDetectors.likelyValidSmiles(lastReactant) &&
        ChemPackageDetectors.likelyValidSmiles(firstProduct);
    }

    return false;
  }

  //name: detectReactions
  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectReactions(col) {
    const lowerName = col.name.toLowerCase();
    const likelyName = ChemPackageDetectors.likelyReactionNames.some((n) => lowerName.includes(n));
    const minUnique = likelyName ? 1 : 2;

    if (DG.Detector.sampleCategories(col, ChemPackageDetectors.isLikelyReaction, minUnique, 10, 0.8)) {
      col.semType = DG.SEMTYPE.CHEMICAL_REACTION;
      col.meta.cellRenderer = 'ChemicalReaction';
      return DG.SEMTYPE.CHEMICAL_REACTION;
    }
    return null;
  }

  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  //meta.skipTest: GROK-17630
  detectMixture(col) {
    if (DG.Detector.sampleCategories(col, (s) => {
      return s.includes(MIX_FILE_VERSION);
    }, 1)) {
      col.semType = CHEMICAL_MIXTURE_SEM_TYPE;
      return col.semType;
    }
    return null;
  }
}
