/**
 * Foot-gun rules: what a value written to a cell, or held in one, would break downstream (in reactions, the editors,
 * the molecule sketcher), and how loud to say it, so that it is flagged before it does. Every rule's level and words
 * are here, in RULES, so that the dialog footer's middle slot (footer-slot.ts), the context panel's Checks section
 * (widgets/checks-panel.ts) and the tests all say the same. Three levels: a block when Chem itself cannot read what
 * would be written (it would be drawn crossed out or as text, and every consumer skips it), which disables OK; a warn
 * when it reads but something downstream breaks or changes meaning, OK enabled; a note for a harmless oddity. The rules
 * read with the module they are given (Crux's: checks.ts hands it in Crux mode, and to the reaction editor in both
 * modes) and import nothing of the platform, so that they run in node too. Each rule names its id in the foot-gun
 * catalogue (MOL-01, TPL-07, ...), which the package tests and the BDD scenarios name too.
 */
import type {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';

/** What the rules do about an error of a Crux call: `noteError` hears every one (a
 * trap the module could not recover from turns Crux off: crux-module.ts noteCruxError), and no call is made while
 * `failed()` says Crux is off. checks.ts sets them; in node they do nothing. */
const cruxGuards = {noteError: (_e: unknown): void => {}, failed: (): boolean => false};

/** Sets what the rules do about Crux's errors (checks.ts: crux-module.ts's noteCruxError and cruxLoadFailed). */
export function setCruxGuards(guards: {noteError: (e: unknown) => void, failed: () => boolean}): void {
  Object.assign(cruxGuards, guards);
}

/** An error of a Crux call heard (a trap turns Crux off), for a catch that goes on without the call's result. */
function heard(e: unknown): void {
  cruxGuards.noteError(e);
}

/** How loud a flag is. */
export type Level = 'block' | 'warn' | 'note';

/** Below this confidence an atom mapping is flagged. */
export const MAPPING_CONFIDENCE_LIMIT = 0.5;
/** The deepest nesting of recursive SMARTS (`$(` in `$(`) Crux's reader takes (found by reading deeper and deeper
 * nestings). */
export const CRUX_RECURSION_LIMIT = 8;
/** The atoms a V2000 molblock holds at most (its counts line has three digits). */
export const V2000_MAX_ATOMS = 999;

/** Every rule: its level and its words, `d` filling in the specifics. The words say what is wrong and what it breaks,
 * in the user's terms; the tests look for their key words. */
export const RULES = {
  // ---- molecules written by a sketcher, or held in a molecule cell
  /** MOL-01, MOL-02: a molecule Chem can't read back as one (`d`: every problem, "an impossible valence on C (atom 1)",
   * "an aromatic ring with no Kekulé structure"). */
  moleculeUnreadable: {level: 'block', words: (d: string) =>
    `Chem can't read this molecule back: ${d}. It would be read as a query at best, and every calculation, search ` +
    'and canonical SMILES skips it'},
  /** MOL-03: a query drawn into a SMILES column: no SMILES holds it, and OK would write an empty cell. */
  queryInSmilesColumn: {level: 'block', words: (d: string) =>
    `A query (${d}) has no SMILES: this column holds SMILES, and OK would empty the cell`},
  /** MOL-17: a ChemDraw reaction pasted into a molecule editor: its arrow left out, every component one molecule. */
  arrowLeftOut: {level: 'warn', words: (_d: string) =>
    'The pasted reaction\'s arrow was left out: the cell gets one molecule of every component, reactants and ' +
    'products together'},
  /** A Crux getter that the engine could not answer (it failed on that read twice, on a fresh engine too: crux-sketch
   * API-069; or it never loaded): OK writes nothing, so the cell keeps its value (the never-write-empty guard, the
   * engine-poison foot-gun); OK is not blocked and nothing asks for a reload: the engine recovers by itself. */
  engineStopped: {level: 'warn', words: (_d: string) =>
    'The drawing could not be read back just now: OK keeps the cell as it was'},
  /** A Crux getter that threw because the drawing cannot be written in the column's notation: OK would write ''. */
  notWritable: {level: 'block', words: (d: string) =>
    `The drawing cannot be written as ${d}: OK would empty the cell, so it is off and the cell keeps its value`},
  /** RXN-15: a value the reaction editor reads as a molecule (its engine detects no reaction in it). */
  editorReadsMolecule: {level: 'block', words: (d: string) =>
    `The editor reads this value as a molecule, not a reaction (${d}): OK would write a molecule into a reaction ` +
    'column'},
  /** RXN-06: a reaction the editor's sketcher does not open (a convergent route): nothing to edit, OK writes nothing.
   * */
  editorRefused: {level: 'block', words: (d: string) =>
    `The editor does not open this reaction (${d}): there is nothing to edit, and the cell stays as it is`},
  /** RXN-16: no arrow left in the reaction editor. */
  reactionNoArrow: {level: 'block', words: (_d: string) =>
    'No reaction arrow is left: OK would write a molecule into a reaction column'},
  /** TPL-12: query features an edit of a template lost (the editor writes an edited component anew). */
  queryLost: {level: 'warn', words: (d: string) =>
    `The edit loses query features of the template: ${d}; the template would match more than it did`},
  /** RXN-17: an edited RXN block written in another notation (crux-sketch has no RXN writer for an edit). */
  rxnNotationLost: {level: 'warn', words: (d: string) =>
    `The RXN block is written as ${d}: its coordinates, its layout and the column's RXN notation are lost`},
  /** Atoms drawn as labels (a reagent name such as HATU) that reaction SMILES or SMARTS write as plain atoms, "*" for
   * a name, which no one would otherwise notice (`d`: "An atom|reaction SMILES"; reactions/reaction-losses.ts). */
  labelAsStar: {level: 'warn', words: (d: string) => {
    const [atoms, notation] = d.split('|');
    return `${atoms} drawn as a label (a name such as HATU) ${atoms === 'An atom' ? 'is' : 'are'} written in ` +
      `${notation} as a plain atom, "*" for a name: calculations, search and mapping see no structure there. Draw ` +
      'the structure instead';
  }},
  /** The drawing's objects (texts, compound numbers, arrow styles) that the notation OK writes does not keep (`d`:
   * "reaction SMILES: the texts (conditions and captions)"; reactions/reaction-losses.ts). */
  drawingNotKept: {level: 'warn', words: (d: string) =>
    `Not kept in ${d}: the reaction is written without them`},
  /** MOL-04: a query in a molecule column (a molblock column takes one). */
  queryInMoleculeColumn: {level: 'warn', words: (d: string) =>
    `A query (${d}) in a molecule column: calculators, canonical SMILES and duplicate checks ignore its query ` +
    'features'},
  /** MOL-05: atom map numbers in a molecule column. */
  mapsInMoleculeColumn: {level: 'warn', words: (d: string) =>
    `Atom map numbers (${d}) in a molecule column: canonical SMILES keep them, so the molecule is no duplicate of ` +
    'its ' +
    'unmapped self'},
  /** MOL-06: R groups the edit adds to a plain molecule. */
  rGroupsAdded: {level: 'warn', words: (d: string) =>
    `R groups (${d}) in a plain molecule: properties and searches count each as an atom`},
  /** MOL-06: attachment points the edit adds. */
  attachmentPointsAdded: {level: 'warn', words: (d: string) =>
    `An attachment point (${d}) in a plain molecule: properties and searches count it as an atom`},
  /** MOL-06: R groups a cell holds (an R-group core column holds them on purpose). */
  rGroupsInCell: {level: 'note', words: (d: string) =>
    `R groups (${d}): properties and searches count each as an atom`},
  /** MOL-06: attachment points a cell holds. */
  attachmentPointsInCell: {level: 'note', words: (d: string) =>
    `Attachment points (${d}): properties and searches count each as an atom`},
  /** MOL-08: enhanced stereo kept in a cell. */
  enhancedStereo: {level: 'note', words: (d: string) =>
    `Enhanced stereo (${d}) is kept in the cell, but canonical SMILES and duplicate checks ignore it`},
  /** MOL-09: S-groups kept in a cell. */
  sgroupsKept: {level: 'note', words: (_d: string) =>
    'S-groups (a polymer bracket, data) are kept in the cell, but canonical SMILES and duplicate checks ignore them'},
  /** MOL-10: a 3D molblock edited flat. */
  lostZ: {level: 'warn', words: (_d: string) =>
    'The 3D coordinates are lost: the edited molecule is written flat (every z = 0)'},
  /** MOL-11: more atoms than a V2000 molblock holds. */
  overV2000: {level: 'warn', words: (d: string) =>
    `${d} atoms: a V2000 molfile holds at most ${V2000_MAX_ATOMS}, so the cell is written as V3000`},
  /** MOL-12: hydrogens drawn as atoms. */
  explicitHydrogens: {level: 'note', words: (d: string) =>
    `${d} hydrogen${d === '1' ? '' : 's'} drawn as atoms ([H]): the canonical SMILES drops them, so the molecule ` +
    'stays a ' +
    'duplicate of its plain form'},
  /** MOL-13: several molecules in one cell. */
  fragments: {level: 'note', words: (d: string) =>
    `${d} fragments in one cell: calculators take them as one molecule`},

  // ---- reactions held in a reaction cell, or written by a reaction editor
  /** TPL-07: a `${…}` placeholder. */
  placeholder: {level: 'block', words: (d: string) =>
    `${d} is a placeholder: Chem reads no reaction until it is replaced`},
  /** TPL-15: recursive SMARTS nested deeper than Crux reads. */
  recursionTooDeep: {level: 'block', words: (d: string) =>
    `Chem can't read this reaction: a recursive SMARTS nested ${d} deep (${CRUX_RECURSION_LIMIT} levels at most)`},
  /** TPL-17, RXN-01, RXN-02: a text Crux reads no reaction from (its reader's own words in `d`): Chem's cells and
   * panel show it as text. */
  reactionUnreadable: {level: 'block', words: (d: string) =>
    `Chem reads no reaction from this text: ${d}`},
  /** TPL-06: a route in a template column. */
  routeAsTemplate: {level: 'warn', words: (d: string) =>
    `A route of ${d} steps written as a template: enumeration and Transform take one reaction a template`},
  /** TPL-13, RXN-12: a template in a column Chem's detector took for something else (one route among templates makes
   * a route column; templates in SMILES syntax a SMILES column): checked as a template all the same. */
  readAsTemplate: {level: 'note', words: (d: string) =>
    `This column was taken for ${d} (its units), but this value is a template: it is checked as one`},
  /** RXN-05: a chain with agents over one of its arrows. */
  chainWithAgents: {level: 'note', words: (_d: string) =>
    'A chain with agents over an arrow (A>X>B>>C): Chem reads its steps, but a tool that takes one reaction ' +
    'refuses it'},
  /** TPL-14, RXN-04: a reaction with no product. */
  noProduct: {level: 'warn', words: (_d: string) => 'No products: running it makes nothing'},
  /** TPL-14, RXN-04: a reaction with no reactant. */
  noReactant: {level: 'warn', words: (_d: string) => 'No reactants: running or mapping it has nothing to start from'},
  /** RXN-03: a component of a reaction (not a template) Chem can't read as a molecule. */
  componentUnreadable: {level: 'warn', words: (d: string) =>
    `${d}: calculators, search and atom mapping skip it`},
  /** TPL-03: a map number twice among the reactants. */
  mapTwiceReactants: {level: 'warn', words: (d: string) =>
    `Map ${d} used twice in the reactants: every run of this template fails`},
  /** TPL-04: a map number twice in the products. */
  mapTwiceProducts: {level: 'warn', words: (d: string) =>
    `Map ${d} used twice in the products: the run joins the atoms into one`},
  /** TPL-01: a map number on the product side only. */
  mapProductOnly: {level: 'warn', words: (d: string) =>
    `Map ${d} on the product side only: such an atom is made new, and the reactant's is deleted`},
  /** TPL-02: map numbers on the reactant side only (how a leaving group is written). */
  mapReactantOnly: {level: 'note', words: (d: string) =>
    `Map ${d} on the reactant side only: such atoms are removed by the reaction`},
  /** TPL-10: a reactant with no mapped atom (an agent written as a reactant). */
  agentAsReactant: {level: 'warn', words: (d: string) =>
    `This template takes ${d} with no mapped atom: it looks like an agent; write it over the arrow`},
  /** TRF-01: a template whose reactants are not the ones a Transform dialog gives (two columns, or one). */
  reactantCount: {level: 'warn', words: (d: string) => `This reaction takes ${d}: every run would fail`},
  /** TPL-11: agents over the arrow in a template. */
  agentsOverArrow: {level: 'warn', words: (d: string) =>
    `Agents over the arrow (${d}): enumeration skips this template`},
  /** TPL-09: query features on product atoms. */
  productQuery: {level: 'warn', words: (d: string) =>
    `Query features in the product are ignored when the template runs: ${d}`},
  /** TPL-05: a product atom with no map number. */
  productAtomUnmapped: {level: 'note', words: (d: string) =>
    `${d} no map number: made new by the reaction`},
  /** TRF-05: a template's blocking SMARTS no engine reads: the enumerator drops it, and blocks nothing. */
  blockingUnreadable: {level: 'block', words: (d: string) =>
    `Blocking SMARTS ${d} cannot be read: the enumerator would drop it and block nothing`},
  /** TRF-05's kin: a blocking SMARTS holding ';' (its low-precedence AND), which the enumerator cuts at ';' into pieces
   * no engine reads: it drops them, and blocks nothing. */
  blockingCut: {level: 'block', words: (d: string) =>
    `Blocking SMARTS ${d} is cut at ";" by the enumerator, which splits a blocking cell at ";" and "|": no engine ` +
    'reads its pieces, so it would block nothing. Write "&" for AND'},
  /** TRF-06: an exclusion SMARTS no engine reads: the enumerator drops it, and excludes nothing. */
  exclusionUnreadable: {level: 'block', words: (d: string) =>
    `Exclusion SMARTS ${d} cannot be read: the enumerator would drop it and exclude nothing`},
  /** TRF-07: a blocking SMARTS that matches none of the building blocks: it blocks nothing. */
  blockingMatchesNothing: {level: 'warn', words: (d: string) =>
    `Blocking SMARTS ${d} matches none of the building blocks: it blocks nothing`},
  /** MAP-02: Map atoms on a template: RXNMapper reads reaction SMILES, and a template's map numbers are its author's.
   * */
  templateMapping: {level: 'warn', words: (_d: string) =>
    'This is a template: its map numbers are its author\'s, and every run follows them; Map atoms maps reactions, ' +
    'not templates'},
  /** MAP-01: a mapping the mapper is not sure of. */
  lowConfidenceMapping: {level: 'warn', words: (d: string) =>
    `Check this mapping: the mapper is only ${Math.round(Number(d) * 100)}% sure ` +
    `(under ${Math.round(MAPPING_CONFIDENCE_LIMIT * 100)}%)`},
  /** MAP-01: products with heavy atoms no reactant or reagent drawn holds (an N-acetylation with no
   * acetyl source): the mapper maps them onto something else, sure or not. `d`: the atoms, "C 2, O 1", with "|repeated"
   * after them when they are one reactant's atoms once more or more times over (crux-reaction-drawing.ts
   * productAtomsMissing; without "|repeated", a self-condensation drawn once was told to add a reagent). */
  unbalancedMapping: {level: 'warn', words: (d: string) => {
    const [atoms, repeated] = d.split('|');
    return `Check this mapping: the reactants drawn don't hold every product atom (${atoms} more in the products): ` +
      (repeated ? 'is a reactant used more than once?' : 'add the reagent that brings them');
  }},
} as const;

export type RuleId = keyof typeof RULES;

/** One thing a value would break, worded for the person about to write it. */
export interface Flag {
  /** The rule that raised it (a key of RULES). */
  rule: RuleId;
  level: Level;
  /** The words: what is wrong, and what it breaks. */
  text: string;
  /** The specifics the words were filled in with. */
  detail?: string;
}

/** The flag `rule` raises, its words filled in with `detail`. */
export function flag(rule: RuleId, detail = ''): Flag {
  const r = RULES[rule];
  return {rule, level: r.level, text: r.words(detail), detail};
}

/** The order flags are shown in: a block first, a note last. */
export const LEVEL_ORDER: {[level in Level]: number} = {block: 0, warn: 1, note: 2};

/** The worst level among `flags`, or null for none. */
export function worstLevel(flags: readonly Flag[]): Level | null {
  return flags.reduce<Level | null>((w, f) => w === null || LEVEL_ORDER[f.level] < LEVEL_ORDER[w] ? f.level : w,
    null);
}

/** `flags` the worst first, their words joined with `separator`. */
export function flagsText(flags: readonly Flag[], separator: string): string {
  return [...flags].sort((a, b) => LEVEL_ORDER[a.level] - LEVEL_ORDER[b.level]).map((f) => f.text).join(separator);
}

// ---------------------------------------------------------------------------------------------------------------------
// Reading SMILES and SMARTS text

/** Whether a string is a molblock. */
export function isMolBlockText(s: string): boolean {
  return s.includes('M  END');
}

/** A SMILES/SMARTS string without its CX extension (` |…|`). */
export function withoutCx(s: string): string {
  return s.replace(/\s+\|[^|]*\|\s*$/, '');
}

/** The top-level bracket atoms of a SMILES or SMARTS (the text between `[` and `]`), a recursive SMARTS inside one kept
 * whole. */
export function bracketAtoms(s: string): string[] {
  const atoms: string[] = [];
  let depth = 0;
  let start = -1;
  for (let i = 0; i < s.length; i++) {
    if (s[i] === '[') {
      if (depth++ === 0)
        start = i + 1;
    } else if (s[i] === ']' && depth > 0 && --depth === 0)
      atoms.push(s.slice(start, i));
  }
  return atoms;
}

/** The atom map number of a bracket atom's text (`C:3` → 3), or 0. */
export function mapNumber(atom: string): number {
  const m = /:(\d+)$/.exec(atom);
  return m ? +m[1] : 0;
}

/** The atoms written outside brackets (the organic subset, `*`, SMARTS' `A` and `a`): they carry no map number. */
export function bareAtomCount(s: string): number {
  let n = 0;
  let depth = 0;
  for (let i = 0; i < s.length; i++) {
    const c = s[i];
    if (c === '[')
      depth++;
    else if (c === ']')
      depth--;
    else if (depth === 0) {
      if ((c === 'C' && s[i + 1] === 'l') || (c === 'B' && s[i + 1] === 'r'))
        i++;
      if ('BCNOPSFIbcnopsAa*'.includes(c))
        n++;
    }
  }
  return n;
}

/** The element symbols a SMILES bracket atom takes, and the aromatic ones (as mol-creation_rdkit.ts's SMILES_ELEMENT).
 * */
const ELEMENT = 'H|He|Li|Be|B|C|N|O|F|Ne|Na|Mg|Al|Si|P|S|Cl|Ar|K|Ca|Sc|Ti|V|Cr|Mn|Fe|Co|Ni|Cu|Zn|Ga|Ge|As|Se|Br|Kr|' +
  'Rb|Sr|Y|Zr|Nb|Mo|Tc|Ru|Rh|Pd|Ag|Cd|In|Sn|Sb|Te|I|Xe|Cs|Ba|La|Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu|Hf|Ta|W|Re|' +
  'Os|Ir|Pt|Au|Hg|Tl|Pb|Bi|Po|At|Rn|Fr|Ra|Ac|Th|Pa|U|Np|Pu|Am|Cm|Bk|Cf|Es|Fm|Md|No|Lr|Rf|Db|Sg|Bh|Hs|Mt|Ds|Rg|Cn|Nh|' +
  'Fl|Mc|Lv|Ts|Og|as|se|si|te|[bcnops*]';
/** What a bracket atom holds after its element in SMILES syntax: chirality, H count, charge, map number. */
const ATOM_PROPERTIES = '(?:@{1,2}|@TH[12]|@AL[12]|@SP[1-3]|@TB\\d{1,2}|@OH\\d{1,2})?(?:H\\d?)?' +
  '(?:[+-]\\d{0,2}|\\+\\+|--)?(?::\\d+)?';
/** A bracket atom in SMILES syntax alone: isotope, element, chirality, H count, charge, map number. */
const SMILES_ATOM = new RegExp(`^\\d{0,3}(?:${ELEMENT})${ATOM_PROPERTIES}$`);
/** A molecule's own atom as SMARTS writes it (get_smarts of a molecule, not of a query): the element as its atomic
 * number, then the atom's isotope, chirality, H count, charge and map number, joined with no `&` or `;` (`[#6@H]`,
 * `[#7+]`, `[2#1]`). Those are the molecule's features, not query features: the same atom as a query is written
 * `[#6&H1]`, `[#7&+]`. */
const NUMBERED_ATOM = new RegExp(`^\\d{0,3}#\\d{1,3}${ATOM_PROPERTIES}$`);

/** What makes a bracket atom a query, as the words name it; null for an atom SMILES writes. `#n`, alone or with the
 * atom's own properties ({@link NUMBERED_ATOM}), counts here (no SMILES writes it); the molecule checks (a molblock
 * read as SMARTS) and the product check leave it out. */
export function atomQuery(atom: string): string | null {
  if (SMILES_ATOM.test(atom))
    return null;
  const a = atom.replace(/:\d+$/, '');
  if (a.includes('$('))
    return 'a recursive SMARTS';
  if (a.includes(','))
    return 'an atom list';
  if (/[;&!]/.test(a))
    return 'an atom query';
  if (NUMBERED_ATOM.test(a))
    return 'an atomic-number atom';
  return 'an atom query';
}

/** The query features of a SMILES/SMARTS text, as the words name them (none for plain SMILES). */
export function queryFeatures(s: string): string[] {
  const body = withoutCx(s);
  const found = new Set<string>();
  for (const atom of bracketAtoms(body)) {
    const q = atomQuery(atom);
    if (q)
      found.add(q);
  }
  const outside = body.replace(/\[[^\]]*\]/g, '');
  if (outside.includes('~'))
    found.add('an any-bond');
  if (/[,;!&@]/.test(outside))
    found.add('a bond query');
  if (/(^|[^lr])[Aa](?![sglmurtcde])/.test(outside.replace(/Cl|Br/g, '')))
    found.add('an any-atom');
  return [...found];
}

/** Whether a reaction text is a template by its own words (MAP-02): query features no reaction SMILES writes
 * (`[N;H2:3]`, `[OH,O-]`); an RXN block or a template in SMILES syntax tells nothing (its column's units do). */
export function isTemplateText(text: string): boolean {
  return !text.trimStart().startsWith('$RXN') && !isMolBlockText(text) &&
    queryFeatures(text.replace(/>/g, '.')).length > 0;
}

/** The deepest nesting of recursive SMARTS in a text. */
export function recursionDepth(s: string): number {
  let deepest = 0;
  const open: boolean[] = []; // per open parenthesis: whether it opened a `$(`
  let depth = 0;
  for (let i = 0; i < s.length; i++) {
    if (s[i] === '(') {
      const recursive = s[i - 1] === '$';
      open.push(recursive);
      if (recursive)
        deepest = Math.max(deepest, ++depth);
    } else if (s[i] === ')' && open.pop())
      depth--;
  }
  return deepest;
}

// ---------------------------------------------------------------------------------------------------------------------
// Molecules (MOL-*)

/** What a molecule column holds, for the checks: its units (DG.chem.Notation's values: smiles, molblock, v3Kmolblock,
 * smarts). */
export interface MoleculeTarget {
  units: string;
}

/** Each element's highest neutral valence that Chem's reader takes, for naming the atom of a valence it rejects. */
const MAX_VALENCE: {[z: number]: number} = {5: 3, 6: 4, 7: 3, 8: 2, 9: 1, 17: 1, 35: 1, 53: 1};
const SYMBOL: {[z: number]: string} = {5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 17: 'Cl', 35: 'Br', 53: 'I'};

/** The atom of an unsanitized molecule whose valence is too high, as "C (atom 1)" (numbered from `offset` + 1: a
 * fragment's atoms after the fragments before it), or a general phrase. */
function valenceAtom(mol: RDMol, offset = 0): string {
  try {
    const json = JSON.parse(mol.get_json());
    const defaults = json.defaults ?? {};
    const m = json.molecules[0];
    const bonds = new Array(m.atoms.length).fill(0);
    for (const b of m.bonds ?? []) {
      const order = b.bo ?? defaults.bond?.bo ?? 1;
      for (const a of b.atoms)
        bonds[a] += order;
    }
    for (let i = 0; i < m.atoms.length; i++) {
      const atom = {...defaults.atom, ...m.atoms[i]};
      const max = MAX_VALENCE[atom.z];
      if (max !== undefined && bonds[i] + (atom.impHs ?? 0) - Math.abs(atom.chg ?? 0) > max)
        return `${SYMBOL[atom.z]} (atom ${offset + i + 1})`;
    }
  } catch {
    // the general phrase
  }
  return 'an atom';
}

/** Whether `module` reads `s` with `details`; the molecule handed to `use` and deleted after. */
function reads<T>(module: RDModule, s: string, details: object, use?: (mol: RDMol) => T): T | boolean {
  if (cruxGuards.failed())
    return false;
  let mol: RDMol | null = null;
  try {
    mol = module.get_mol(s, JSON.stringify(details)) ?? null;
    return mol ? (use ? use(mol) : true) : false;
  } catch (e) {
    heard(e);
    return false;
  } finally {
    mol?.delete();
  }
}

/** The block of a molecule Chem cannot read back as one, as its readers read it (getMolSafe: sanitized, then without
 * kekulizing; a query only is no molecule), or null. */
export function unreadableMolecule(module: RDModule, s: string): Flag | null {
  if (reads(module, s, {}))
    return null;
  if (reads(module, s, {kekulize: false}))
    return flag('moleculeUnreadable', 'an aromatic ring with no Kekulé structure');
  const atom = reads(module, s, {sanitize: false}, (mol) => valenceAtom(mol));
  if (atom !== false)
    return flag('moleculeUnreadable', `an impossible valence on ${atom}`);
  return flag('moleculeUnreadable', 'text that doesn\'t parse');
}

/** Whether a molblock holds query atoms or bonds: an atom list, a query bond, a query property. */
export function molblockQuery(module: RDModule, molblock: string): string | null {
  const smarts = reads(module, molblock, {sanitize: false}, (mol) => mol.get_smarts());
  if (typeof smarts === 'string') {
    const f = queryFeatures(smarts).filter((x) => x !== 'an atomic-number atom');
    if (f.length)
      return f.join(', ');
    // a molecule read from it shows none, though the molblock holds some (a ring query, as a package test found):
    // the query read names them
    if (!molblockQueryMarks(molblock))
      return null;
  }
  if (cruxGuards.failed())
    return null;
  let q: RDMol | null = null;
  try {
    q = module.get_qmol(molblock);
    if (!q)
      return null;
    const f = queryFeatures(q.get_smarts()).filter((x) => x !== 'an atomic-number atom');
    return f.length ? f.join(', ') : 'query atoms or bonds';
  } catch (e) {
    heard(e);
    return null;
  } finally {
    q?.delete();
  }
}

/**
 * Whether a molblock carries query features that a molecule read from it drops: V2000's query bond types 5 to 8 (single
 * or double, single or aromatic, double or aromatic, any), its atom-list, ring-bond-count, substitution and
 * unsaturation lines, and a query atom's SMARTS kept as data (an `M  SED` holding a bracket atom, as Crux and RDKit
 * write `[#6;R]`); V3000's bond types 5 to 8, its RBCNT, SUBST and UNSAT fields and its atom lists.
 */
function molblockQueryMarks(molblock: string): boolean {
  if (molblock.includes('V3000')) {
    return /^M {2}V30 \d+ [5-8] \d+ \d+/m.test(molblock) || /(?:RBCNT|SUBST|UNSAT)=/.test(molblock) ||
      /^M {2}V30 \d+ (?:NOT)?\[/m.test(molblock);
  }
  const lines = molblock.split('\n');
  const atoms = molblockAtoms(molblock);
  const bonds = +(lines[3] ?? '').slice(3, 6) || 0;
  if (lines.slice(4 + atoms, 4 + atoms + bonds).some((l) => +l.slice(6, 9) >= 5 && +l.slice(6, 9) <= 8))
    return true;
  return lines.some((l) => /^M {2}(?:ALS|RBC|SUB|UNS) /.test(l) || /^M {2}SED +\d+ +\[.*\]\s*$/.test(l));
}

/** A molecule's atom count as its text says it (a V2000 counts line, a V3000 COUNTS line), 0 for a SMILES. */
function molblockAtoms(molblock: string): number {
  const v3 = /M {2}V30 COUNTS (\d+)/.exec(molblock);
  return v3 ? +v3[1] : +(molblock.split('\n')[3] ?? '').slice(0, 3) || 0;
}

/** Whether a molblock's atoms have a z coordinate other than 0 (a conformer's, a 3D structure's). */
export function molblockHasZ(molblock: string): boolean {
  if (molblock.includes('V3000'))
    return [...molblock.matchAll(/^M {2}V30 \d+ \S+ \S+ \S+ (\S+)/gm)].some((m) => Math.abs(+m[1]) > 1e-4);
  const lines = molblock.split('\n');
  return lines.slice(4, 4 + molblockAtoms(molblock)).some((l) => Math.abs(+l.slice(20, 30)) > 1e-4);
}

/** A molecule's text as SMILES for the text checks: a SMILES as it is, a molblock as Crux writes it (unsanitized). */
function asSmiles(module: RDModule, value: string): string {
  if (!isMolBlockText(value))
    return value;
  const cx = reads(module, value, {sanitize: false}, (mol) => mol.get_cxsmiles());
  return typeof cx === 'string' ? cx : '';
}

/** The atom map numbers of a SMILES (its bracket atoms', R groups' `[*:n]` aside). */
function atomMaps(smiles: string): number[] {
  return bracketAtoms(withoutCx(smiles)).filter((a) => !a.startsWith('*')).map(mapNumber).filter((n) => n > 0);
}

/** The R groups (`[*:n]` as Rn, CX `_Rn` labels) and attachment points (CX `_APn`) of a SMILES. */
function rGroups(smiles: string): {rs: string[], aps: string[]} {
  const rs = new Set<string>();
  for (const a of bracketAtoms(withoutCx(smiles))) {
    const m = /^\*:(\d+)$/.exec(a);
    if (m)
      rs.add(`R${m[1]}`);
  }
  for (const m of smiles.matchAll(/_R(\d+)/g))
    rs.add(`R${m[1]}`);
  const aps = new Set([...smiles.matchAll(/_AP(\d+)/g)].map((m) => `_AP${m[1]}`));
  return {rs: [...rs], aps: [...aps]};
}

/** The enhanced stereo groups a CXSMILES or molblock holds, in words, or null. */
function enhancedStereo(value: string): string | null {
  if (isMolBlockText(value))
    return /MDLV30\/STE(RAC|REL)/.test(value) ? 'AND/OR groups' : null;
  const cx = /\|([^|]*)\|\s*$/.exec(value)?.[1] ?? '';
  const groups = [...cx.matchAll(/(?:^|,)([&o])(\d+):/g)].map((m) => `${m[1] === '&' ? 'and' : 'or'}${m[2]}`);
  return groups.length ? groups.join(', ') : null;
}

/** Whether a CXSMILES or molblock holds S-groups (a polymer bracket, a data S-group: CX `Sg:`/`SgD:`, molblock `STY`).
 * */
function sgroups(value: string): boolean {
  if (isMolBlockText(value))
    return /^M {2}STY/m.test(value) || /^M {2}V30 BEGIN SGROUP/m.test(value);
  const cx = /\|([^|]*)\|\s*$/.exec(value)?.[1] ?? '';
  return /(?:^|,)SgD?:/.test(cx);
}

/** What Chem can't read back in each fragment of a SMILES, as one phrase (every problem, so that the footer's tooltip
 * names all), or null when every fragment reads. */
function fragmentProblems(module: RDModule, smiles: string): string | null {
  if (reads(module, smiles, {}))
    return null;
  const body = withoutCx(smiles);
  // a ring closure across a dot joins the fragments: then the whole is read as one
  const fragments = /\.\d|\.%/.test(body) ? [body] : body.split('.').filter((f) => f);
  const problems: string[] = [];
  const valences: string[] = [];
  let offset = 0;
  for (const f of fragments) {
    if (!reads(module, f, {})) {
      if (reads(module, f, {kekulize: false}))
        problems.push('an aromatic ring with no Kekulé structure');
      else {
        const atom = reads(module, f, {sanitize: false}, (mol) => valenceAtom(mol, offset));
        if (atom !== false)
          valences.push(atom as string);
        else
          problems.push('text that doesn\'t parse');
      }
    }
    offset += reads(module, f, {sanitize: false}, (mol) => mol.get_num_atoms()) as number || 0;
  }
  if (valences.length)
    problems.unshift(`an impossible valence on ${valences.join(' and on ')}`);
  return problems.length ? [...new Set(problems)].join('; ') : 'text that doesn\'t parse';
}

/**
 * The flags of a molecule sketcher's OK: writing `written` into a column of `target`'s units over the cell's
 * `previous` value. `smilesQuery` is non-null when the sketcher holds a query its SMILES can't write (Crux's SMILES
 * of a query is empty, so OK would empty a SMILES column's cell: MOL-03): the query's features in words, or '' to read
 * them from `drawing`, the sketcher's molblock. What the user adds that the cell did not have is a warn; what it had
 * already keeps the level the Checks section gives it.
 */
export function sketcherWriteFlags(module: RDModule, written: string, target: MoleculeTarget,
  smilesQuery: string | null = null, drawing = '', previous = ''): Flag[] {
  if (target.units === 'smarts')
    return [];
  if (target.units === 'smiles' && smilesQuery !== null) {
    const features = smilesQuery || (drawing ? molblockQuery(module, drawing) : null) || 'query features';
    return [flag('queryInSmilesColumn', features)];
  }
  // another sketcher's empty SMILES of a query drawing (it says nothing of its own): the same block
  const drawnQuery = target.units === 'smiles' && !written?.trim() && drawing ? molblockQuery(module, drawing) : null;
  if (drawnQuery)
    return [flag('queryInSmilesColumn', drawnQuery)];
  if (!written?.trim())
    return [];
  const before = previous ? moleculeFlags(module, previous, target) : [];
  const flags = moleculeFlags(module, written, target, true).map((f) =>
    before.some((b) => b.rule === CELL_RULE[f.rule]) ? flag(CELL_RULE[f.rule]!, detailOf(f)) : f);
  if (isMolBlockText(written) && isMolBlockText(previous) && molblockHasZ(previous) && !molblockHasZ(written))
    flags.push(flag('lostZ'));
  if (target.units === 'molblock' && isMolBlockText(written) && molblockAtoms(written) > V2000_MAX_ATOMS)
    flags.push(flag('overV2000', `${molblockAtoms(written)}`));
  return flags;
}

/** The detail a flag was raised with (kept beside its words). */
function detailOf(f: Flag): string {
  return f.detail ?? '';
}

/**
 * The flags of a molecule `value` in a column of `target`'s units: a block where Chem can't read it back as a molecule;
 * a query in a molecule column (a warn); atom maps (a warn); R groups and attachment points (a note in a cell, a warn
 * when `edit`: the sketcher adds them); notes for enhanced stereo, explicit hydrogens and several fragments.
 */
export function moleculeFlags(module: RDModule, value: string, target: MoleculeTarget, edit = false): Flag[] {
  if (!value?.trim() || target.units === 'smarts')
    return [];
  const query = isMolBlockText(value) ? molblockQuery(module, value) :
    queryFeatures(value).filter((f) => f !== 'an any-atom').join(', ') || null;
  if (query)
    return [flag('queryInMoleculeColumn', query)];
  const flags: Flag[] = [];
  const unreadable = isMolBlockText(value) ? unreadableMolecule(module, value) :
    (() => {
      const problems = fragmentProblems(module, value);
      return problems === null ? null : flag('moleculeUnreadable', problems);
    })();
  if (unreadable)
    flags.push(unreadable);
  const smiles = asSmiles(module, value);
  const maps = [...new Set(atomMaps(smiles))];
  if (maps.length)
    flags.push(flag('mapsInMoleculeColumn', maps.join(', ')));
  const {rs, aps} = rGroups(smiles);
  if (rs.length)
    flags.push(flag(edit ? 'rGroupsAdded' : 'rGroupsInCell', rs.join(', ')));
  if (aps.length)
    flags.push(flag(edit ? 'attachmentPointsAdded' : 'attachmentPointsInCell', aps.join(', ')));
  const stereo = enhancedStereo(value);
  if (stereo)
    flags.push(flag('enhancedStereo', stereo));
  if (sgroups(value))
    flags.push(flag('sgroupsKept'));
  if (!unreadable) {
    const hs = bracketAtoms(withoutCx(smiles)).filter((a) => /^\d*H(:\d+)?$/.test(a)).length;
    if (hs)
      flags.push(flag('explicitHydrogens', `${hs}`));
    const fragments = withoutCx(smiles).split('.').filter((f) => f).length;
    if (fragments > 1)
      flags.push(flag('fragments', `${fragments}`));
  }
  return flags;
}

/** The rule a cell already holding what an edit adds raises instead: the edit keeps what the cell had. */
const CELL_RULE: {[rule in RuleId]?: RuleId} = {
  rGroupsAdded: 'rGroupsInCell',
  attachmentPointsAdded: 'attachmentPointsInCell',
};

// ---------------------------------------------------------------------------------------------------------------------
// Reactions and templates (TPL-*, TRF-01)

/** What a reaction column holds, for the checks: whether its values are templates (a reaction column of units
 * `smarts`, as Chem's detector writes them, detectors.js reactionNotation). */
export interface ReactionTarget {
  template: boolean;
  /** The column's units (`smiles`, `route`, `rxn`, `smarts`), when the value is a column's. */
  units?: string;
}

/** A reaction column's units as the words say them. */
const UNITS_WORDS: {[units: string]: string} = {smiles: 'reaction SMILES', route: 'routes', rxn: 'RXN blocks',
  smarts: 'reaction SMARTS'};

/** Whether a molecule text reads with an atom short of its bonds (a radical): `[C:1]` written without its H count. */
function hasRadical(module: RDModule, s: string): boolean {
  return reads(module, s, {}, (mol) => {
    const json = JSON.parse(mol.get_json());
    const defaults = json.defaults?.atom ?? {};
    return (json.molecules?.[0]?.atoms ?? []).some((a: any) => ({...defaults, ...a}).nRad > 0);
  }) === true;
}

/**
 * Whether a reaction text is a template by what it holds, whatever its column's units (TPL-13, RXN-12): query
 * features, or map numbers with a component that is no molecule as written (an atom written without its H count,
 * `[C:1]`, `[NH2:3]` alone, reads as a radical; an aromatic atom outside a ring reads as nothing). A mapped reaction
 * (a mapper's: `[CH3:1][C:2](=[O:3])[OH:4]`) is none.
 */
export function readsAsTemplate(module: RDModule, text: string): boolean {
  if (isTemplateText(text))
    return true;
  const body = withoutCx(text);
  if (body.trimStart().startsWith('$RXN') || isMolBlockText(body) || !bracketAtoms(body).some((a) => mapNumber(a) > 0))
    return false;
  const components = body.split(/>|--\*\*--/).flatMap((side) => side.split('.')).map((c) => c.trim())
    .filter((c) => c);
  return components.some((c) => unreadableMolecule(module, c) !== null || hasRadical(module, c));
}

/** Whether a chain has a step with agents over its arrow (`A>X>B>>C`, RXN-05); a step list's steps (`A>X>B--**--…`)
 * are no chain. */
export function chainHasAgents(text: string): boolean {
  return withoutCx(text).split('--**--').some((b) => {
    const parts = b.split('>');
    return parts.length > 3 && parts.some((p, i) => i % 2 === 1 && p.trim() !== '');
  });
}

/** The steps of a route (`--**--` between steps, or a chain `A>>B>>C`); 1 for one reaction or an RXN block. */
export function routeStepCount(text: string): number {
  if (text.trimStart().startsWith('$RXN'))
    return 1;
  return text.split('--**--').reduce((n, b) => n + Math.max(1, b.split('>>').length - 1), 0);
}

/** A one-step reaction's reactant, agent and product components (`A.B>C>D`, its CX extension dropped); null when it is
 * not three sides. */
export function reactionSides(text: string): {reactants: string[], agents: string[], products: string[]} | null {
  const parts = withoutCx(text.trim()).split('>');
  if (parts.length !== 3)
    return null;
  const split = (s: string) => s.split('.').filter((c) => c);
  return {reactants: split(parts[0]), agents: split(parts[1]), products: split(parts[2])};
}

/** The map numbers of components, in order, repeats kept. */
function mapsOf(components: string[]): number[] {
  return components.flatMap((c) => bracketAtoms(c).map(mapNumber)).filter((n) => n > 0);
}

/** Map numbers as the words say them, with their verb: "number 5 is", "numbers 1, 2, 3 are". */
function mapNumbers(numbers: number[]): string {
  return numbers.length === 1 ? `number ${numbers[0]} is` : `numbers ${numbers.join(', ')} are`;
}

/** The numbers that occur more than once. */
function repeated(numbers: number[]): number[] {
  return [...new Set(numbers.filter((n, i) => numbers.indexOf(n) !== i))].sort((a, b) => a - b);
}

/** Why Crux reads no reaction from `text` (its scheme's own words, the reader Chem's cells and panel use), or null when
 * it reads one. `module` is Crux's (a ChemReactionScheme class and get_rxn). */
export function unreadableReaction(module: any, text: string): Flag | null {
  const placeholders = [...new Set([...text.matchAll(/\$\{[^}]*\}/g)].map((m) => m[0]))];
  if (placeholders.length)
    return flag('placeholder', placeholders.join(', '));
  const depth = recursionDepth(text);
  if (cruxGuards.failed())
    return null;
  let scheme: any = null;
  try {
    scheme = new module.ChemReactionScheme(text, '');
  } catch (e) {
    heard(e);
    if (depth > CRUX_RECURSION_LIMIT)
      return flag('recursionTooDeep', `${depth}`);
    const why = String((e as any)?.message ?? e).replace(/^ChemReactionScheme:\s*/, '');
    return flag('reactionUnreadable', why === 'not a reaction' ? 'the text is not a reaction it can parse' : why);
  } finally {
    scheme?.free();
  }
  return null;
}

/** The template checks of one-step reaction `smarts`: its sides, its map numbers, its agents, its product
 * atoms' queries and map numbers. */
export function templateFlags(smarts: string): Flag[] {
  const sides = reactionSides(smarts);
  if (!sides)
    return [];
  const {reactants, agents, products} = sides;
  const flags: Flag[] = [];
  const left = mapsOf(reactants);
  const right = mapsOf(products);
  if (repeated(left).length)
    flags.push(flag('mapTwiceReactants', mapNumbers(repeated(left))));
  if (repeated(right).length)
    flags.push(flag('mapTwiceProducts', mapNumbers(repeated(right))));
  // one side's numbers missing from a side that is there (an empty side is its own warn)
  const productOnly = [...new Set(right.filter((n) => !left.includes(n)))].sort((a, b) => a - b);
  if (productOnly.length && reactants.length)
    flags.push(flag('mapProductOnly', mapNumbers(productOnly)));
  const reactantOnly = [...new Set(left.filter((n) => !right.includes(n)))].sort((a, b) => a - b);
  if (reactantOnly.length && products.length)
    flags.push(flag('mapReactantOnly', mapNumbers(reactantOnly)));
  const unmappedReactants = reactants.filter((c) => !mapsOf([c]).length);
  if (left.length && unmappedReactants.length) {
    flags.push(flag('agentAsReactant', `${reactants.length} reactants, ${unmappedReactants.length === 1 ? 'one' :
      'some'} of them ${unmappedReactants.join(' and ')}`));
  }
  if (agents.length)
    flags.push(flag('agentsOverArrow', agents.join(', ')));
  const queries = [...new Set(products.flatMap((c) => bracketAtoms(c))
    .filter((a) => {
      const q = atomQuery(a);
      return q !== null && q !== 'an atomic-number atom' && !/^\*(:\d+)?$/.test(a);
    }))].map((a) => `[${a}]`);
  if (queries.length)
    flags.push(flag('productQuery', queries.join(', ')));
  if (right.length && reactants.length) {
    const unmapped = products.reduce((n, c) => n + bareAtomCount(c) +
      bracketAtoms(c).filter((a) => !mapNumber(a)).length,
    0);
    if (unmapped)
      flags.push(flag('productAtomUnmapped', unmapped === 1 ? '1 product atom has' : `${unmapped} product atoms have`));
  }
  return flags;
}

/**
 * The flags of reaction `text` in a column of `target` (a template column, or a reaction editor's template, runs the
 * template checks): a block where Crux reads no reaction; a route in a template column; then the template's own.
 */
export function reactionFlags(module: any, text: string, target: ReactionTarget): Flag[] {
  if (!text?.trim())
    return [];
  const unreadable = unreadableReaction(module, text);
  if (unreadable)
    return [unreadable];
  const flags: Flag[] = [];
  // a template in a column taken for something else is checked as one, and noted (TPL-13, RXN-12)
  if (!target.template && readsAsTemplate(module, text)) {
    if (target.units && target.units !== 'smarts')
      flags.push(flag('readAsTemplate', UNITS_WORDS[target.units] ?? target.units));
    target = {...target, template: true};
  }
  // a chain with agents over an arrow (RXN-05): a note, whatever it is read as
  if (!target.template && chainHasAgents(text))
    return [...flags, flag('chainWithAgents')];
  const steps = routeStepCount(text);
  if (steps > 1)
    return target.template ? [...flags, flag('routeAsTemplate', `${steps}`)] : flags;
  const sides = text.trimStart().startsWith('$RXN') ? null : reactionSides(text);
  if (sides && !sides.reactants.length)
    flags.push(flag('noReactant'));
  if (sides && !sides.products.length)
    flags.push(flag('noProduct'));
  if (!target.template) {
    // a reaction's components read as molecules (a SMARTS component is a template's: none of this applies)
    const named: [string[], string][] = [[sides?.reactants ?? [], 'Reactant'], [sides?.agents ?? [], 'Agent'],
      [sides?.products ?? [], 'Product']];
    for (const [components, role] of named) {
      components.forEach((c, i) => {
        if (queryFeatures(c).length)
          return;
        const problem = unreadableMolecule(module, c);
        if (problem) {
          flags.push(flag('componentUnreadable',
            `${role} ${i + 1} (${c}) can't be read as a molecule: ${problem.detail}`));
        }
      });
    }
    return flags;
  }
  let smarts = text;
  if (text.trimStart().startsWith('$RXN')) {
    let rxn: any = null;
    try {
      rxn = cruxGuards.failed() ? null : module.get_rxn(text);
      smarts = rxn?.get_smarts() ?? text;
    } catch (e) {
      heard(e); // the text as it is
    } finally {
      rxn?.delete();
    }
  }
  return [...flags, ...templateFlags(smarts)];
}

/** The flag of a reaction run by a Transform dialog that gives it `given` reactants (Two-Component: 2, Transformation:
 * 1) when its template takes another number (TRF-01). */
export function reactantCountFlags(smarts: string, given: number): Flag[] {
  const n = reactionSides(smarts)?.reactants.length ?? given;
  return n === given ? [] : [flag('reactantCount', `${n} reactant${n === 1 ? '' : 's'}; the dialog gives ${given}`)];
}

/** The flag of Map atoms on a template (MAP-02). */
export function templateMappingFlags(): Flag[] {
  return [flag('templateMapping')];
}

/** The flag of an atom mapping at `confidence`: a warn under the limit, else none. */
export function mappingFlags(confidence: number): Flag[] {
  return confidence < MAPPING_CONFIDENCE_LIMIT ? [flag('lowConfidenceMapping', confidence.toFixed(2))] : [];
}

// ---------------------------------------------------------------------------------------------------------------------
// The reaction cell editor (RXN-14..17, TPL-12)

/** Why crux-sketch's reaction document reads a reaction text as a molecule (its engine detects no reaction in it, while
 * Chem's cells read one: sweep r-space, t-dative), or null: a space beside the arrow, a dative bond `->`. */
export function readAsMoleculeWhy(text: string): string | null {
  if (text.trimStart().startsWith('$RXN') || !text.includes('>'))
    return null;
  const body = withoutCx(text);
  if (/->/.test(body))
    return 'a dative bond "->" in it, which the editor takes for no reaction arrow';
  if (/\s>|>\s/.test(body.trim()))
    return 'a space beside the arrow';
  return null;
}

/** The query primitives of a template's bracket atoms (`NX3;H2,H1;!$(N=*)` → NX3, H2,H1, !$(N=*)), each once. */
function queryPrimitives(text: string): string[] {
  const out = new Set<string>();
  for (const atom of bracketAtoms(withoutCx(text))) {
    if (atomQuery(atom) === null)
      continue;
    let depth = 0;
    let start = 0;
    const a = atom.replace(/:\d+$/, '');
    for (let i = 0; i <= a.length; i++) {
      const c = a[i];
      if (c === '(' || c === '[')
        depth++;
      else if (c === ')' || c === ']')
        depth--;
      else if ((c === ';' || c === '&' || i === a.length) && depth === 0) {
        const p = a.slice(start, i);
        if (p && atomQuery(p) !== null)
          out.add(p);
        start = i + 1;
      }
    }
  }
  return [...out];
}

/**
 * The flags of the reaction cell editor's OK (RXN-14..17, TPL-12; the template checks of 13 too): `original` the
 * cell's text as the editor opened, `text` what OK would write now, `edited` whether a person changed it. Unedited: a
 * block where the editor reads the value as a molecule (RXN-15), nothing else (OK writes nothing). Edited: a block when
 * no arrow is left (RXN-16), else the value's own checks, the query features the edit lost (TPL-12) and an RXN block
 * written in another notation (RXN-17).
 */
export function reactionEditFlags(module: RDModule, original: string, text: string, target: ReactionTarget,
  edited: boolean): Flag[] {
  if (!edited) {
    const why = readAsMoleculeWhy(original);
    return why ? [flag('editorReadsMolecule', why)] : [];
  }
  const rxnBlock = text.trimStart().startsWith('$RXN');
  if (!rxnBlock && !text.includes('>'))
    return [flag('reactionNoArrow')];
  const flags = reactionFlags(module, text, target);
  const lost = queryPrimitives(original).filter((p) => !text.includes(p));
  if (lost.length)
    flags.push(flag('queryLost', lost.join(', ')));
  if (original.trimStart().startsWith('$RXN') && !rxnBlock)
    flags.push(flag('rxnNotationLost', queryFeatures(text).length ? 'reaction SMARTS' : 'reaction SMILES'));
  return flags;
}
