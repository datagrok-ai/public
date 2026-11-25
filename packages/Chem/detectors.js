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

  static validSmilesChars = new Uint32Array([0, 805306168, 1073741823, 1073741822, 0, 0, 0, 0]);
  static validFirstSmilesChars = new Uint32Array([0, 1280, 134857308, 573448, 0, 0, 0, 0]);
  static numberOfMolsToCheck = 100;
  static numberOfMolsToCheckRdKit = 10;

  static validChar(pos) {
    return (ChemPackageDetectors.validSmilesChars[Math.floor(pos / 0x20)] & (1 << (pos & 0x1f))) !== 0;
  }

  static validFirstChar(pos) {
    return (ChemPackageDetectors.validFirstSmilesChars[Math.floor(pos / 0x20)] & (1 << (pos & 0x1f))) !== 0;
  }

  // characters SMARTS holds past an atom's first character and SMILES does not: `!` and `&` in atom and bond
  // expressions ([C;!R], [N&H2], C!@C), `~` as the any-bond
  static smartsOnlyChars = '!&~';

  // with smarts set, smartsOnlyChars are accepted too (a reaction template's sides)
  static likelyValidSmiles(s, smarts = false) {
    if (!ChemPackageDetectors.validFirstChar(s[0].charCodeAt(0)))
      return false;
    for (let j = 1; j < s.length; j++) {
      const ch = s[j];
      if (ch === ' ' && j + 1 < s.length && s[j + 1] === '|')
        continue;
      if (!ChemPackageDetectors.validChar(ch.charCodeAt(0)) &&
        !(smarts && ChemPackageDetectors.smartsOnlyChars.includes(ch)))
        return false;
    }
    return true;
  }

  //name: detectMolecules
  //meta.role: semTypeDetector
  //input: column col {type: string}
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

  /** A CDXML page (as cdxml-reaction-cells.ts isCdxmlText reads one). */
  static isCdxml(s) {
    return typeof s === 'string' && /^\s*(<\?xml[^>]*>\s*)?(<!DOCTYPE[^>]*>\s*)?<CDXML[\s>]/.test(s.slice(0, 4096));
  }

  // the share of sampled values that must look like reactions: 80 %, passed as 0.79 because js-api's
  // sampleCategories (js-api/src/data.ts:473-481) allows 10 * (1 - 0.8) = 1.9999... failures in ten samples, so 0.8
  // rejects a column at its second failure (90 %); 0.79 allows two
  static reactionSampleRatio = 0.79;

  // a reaction column's units: the notation its values are written in
  static reactionUnits = {smiles: 'smiles', smarts: 'smarts', rxn: 'rxn', route: 'route'};
  // how many of a column's distinct values the units are read from (evenly spaced)
  static reactionUnitsSample = 100;

  // the notation of one reaction value: an RXN block; a route (steps joined by --**--, or a chain A>>B>>C); reaction
  // SMARTS when SMARTS syntax shows outside a CX extension (the test of _isSmarts, src/utils/mol-creation_rdkit.ts);
  // else reaction SMILES
  static reactionNotation(s) {
    const units = ChemPackageDetectors.reactionUnits;
    if (s.trimStart().startsWith('$RXN'))
      return units.rxn;
    if (s.includes('--**--') || s.split('>>').length > 2)
      return units.route;
    if (/\[.?#\d|\$|&|;|,|!.?]|~/.test(s.replace(/\s+\|[^|]*\|\s*$/, '')))
      return units.smarts;
    return units.smiles;
  }

  // a column's units: the first of route, rxn, smarts and smiles that one of its sampled values is written in
  static reactionColumnUnits(col) {
    const categories = col.categories;
    const n = Math.min(ChemPackageDetectors.reactionUnitsSample, categories.length);
    const found = new Set();
    for (let i = 0; i < n; i++) {
      const s = categories[Math.floor(categories.length * i / n)];
      if (s)
        found.add(ChemPackageDetectors.reactionNotation(s));
    }
    const units = ChemPackageDetectors.reactionUnits;
    return [units.route, units.rxn, units.smarts, units.smiles].find((u) => found.has(u)) ?? units.smiles;
  }

  // Check if a string looks like a reaction SMILES/SMARTS or RXN block.
  // Supports four formats:
  //   "reactants>>products"           — standard reaction SMILES
  //   "reactants>agents>products"      — SMIRKS with agent/reagent field
  //   "A>>B>>C>>D"                     — multi-step reaction chain
  //   "A.B>C>D--**--D.E>>F"            — a route, its steps joined by --**-- (the enumerator's form): its first
  //                                      step decides
  //   "$RXN..."                        — MDL RXN format (molblock-based)
  static isLikelyReaction(value) {
    // RXN format: structured block starting with $RXN
    if (value.trimStart().startsWith('$RXN')) return true;
    const s = value.split('--**--')[0];

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
        ChemPackageDetectors.likelyValidSmiles(lastReactant, true) &&
        ChemPackageDetectors.likelyValidSmiles(firstProduct, true);
    }

    // SMIRKS format: "reactants>agents>products" (exactly two single ">")
    const parts = s.split('>');
    if (parts.length === 3 && parts[0].length > 0 && parts[2].length > 0) {
      const lastReactant = parts[0].split('.').pop();
      const firstProduct = parts[2].split('.')[0];
      return lastReactant.length > 0 && firstProduct.length > 0 &&
        ChemPackageDetectors.likelyValidSmiles(lastReactant, true) &&
        ChemPackageDetectors.likelyValidSmiles(firstProduct, true);
    }

    return false;
  }

  //name: detectReactions
  //meta.role: semTypeDetector
  //input: column col {type: string}
  //output: string semType
  detectReactions(col) {
    // a column of CDXML pages (ChemDraw schemes): Chem reads them as reaction drawings in Crux mode
    // (cdxml-reaction-cells.ts); the engine setting is Chem's to check, so it is asked as detectSmiles is
    if (DG.Detector.sampleCategories(col, ChemPackageDetectors.isCdxml, 1, 10, 0.8)) {
      grok.functions.call('Chem:detectCdxmlReactions', {col: col}).then(() => {});
      return null;
    }
    const lowerName = col.name.toLowerCase();
    const likelyName = ChemPackageDetectors.likelyReactionNames.some((n) => lowerName.includes(n));
    const minUnique = likelyName ? 1 : 2;

    if (DG.Detector.sampleCategories(col, ChemPackageDetectors.isLikelyReaction, minUnique, 10,
      ChemPackageDetectors.reactionSampleRatio)) {
      col.semType = DG.SEMTYPE.CHEMICAL_REACTION;
      col.meta.cellRenderer = 'ChemicalReaction';
      col.meta.units = ChemPackageDetectors.reactionColumnUnits(col);
      return DG.SEMTYPE.CHEMICAL_REACTION;
    }
    return null;
  }

  //meta.role: semTypeDetector
  //input: column col {type: string}
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
