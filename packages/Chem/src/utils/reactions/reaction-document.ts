/**
 * The Reaction Sketcher's document (utils/reactions/reactions-app.ts): ONE reaction or route, a title and the text
 * Crux reads it from (reaction SMILES or SMARTS with a CX part; a route as a chain, a step list or route JSON; an RXN
 * block; an RD file). This module holds the document, the default example the app opens on the first time (a five-step
 * mapped route with agents: the app never opens empty, but on a drawing that shows what it can do), the nine titled
 * examples of the Examples menu, the copy kept in the browser's localStorage after a person's change, the pictures,
 * and the export to SVG and PNG. Every picture is drawn by Crux's scheme painter (`ChemReactionScheme` in Chem's Crux
 * file), in either engine mode, with Chem's molecule draw options (as reaction cells are).
 */
import {ChemWasm, cruxLoadFailed, noteCruxError} from '../../engine/crux-module';
import {CHEM_CENTRE, chemDrawOptions, tooLargeToDraw} from '../../rendering/crux-reaction-drawing';

/** The browser's localStorage key that holds the document, and the version of what is stored under it: 2, one
 * reaction (`{version, title, text, mapNumbers}`); a page of version 1 is read as none (the default example). */
export const DOCUMENT_KEY = 'Chem.ReactionsApp.document';
const DOCUMENT_VERSION = 2;
/** The page: the most it widens to, and its padding (CSS px). */
export const PAGE_MAX_WIDTH = 1000;
export const PAGE_PADDING = 24;
/** The pictures' bond length (px; the cells draw at Crux's own 22). */
const BOND_LENGTH = 25;
/** A reaction wider than the page is broken into rows by Crux: its box leaves room for this many rows (Crux's most),
 * and the blank around what it drew is cut off. */
const MAX_ROWS = 3;
const ROW_ROOM = 1.5;
/** The room between a picture's ink and its frame (px). */
const PICTURE_MARGIN = 6;
/** A text Crux reads no reaction from: the height of its line in the page and the export (px), font and colour. */
const UNREAD_HEIGHT = 28;
const UNREAD_FONT = '12px monospace';
const UNREAD_COLOUR = '#7f7f7f';
/** The titles in the export: font, colour, the line they take (px). */
const TITLE_FONT_FAMILY = 'Roboto, Arial, sans-serif';
const TITLE_FONT_SIZE = 14;
const TITLE_COLOUR = '#2a2a2a';
export const TITLE_HEIGHT = 24;
/** The gap between two reactions in the export (px). */
const EXPORT_GAP = 18;
/** The PNG export's resolution (dpi; a CSS px is 1/96 in), and the largest side a canvas takes (device px). */
const PNG_DPI = 300;
const CANVAS_MAX_SIDE = 16000;

/** The app's document, and each example: one reaction or route. */
export interface DocumentReaction {
  title: string;
  /** The reaction in any notation Crux reads (module docs). */
  text: string;
  /** Map numbers drawn beside their atoms (Crux's `mapNumbers` option). */
  mapNumbers: boolean;
  /** The drawing as a CDXML page (Crux Sketch's writer), kept beside the text while it holds what no reaction text
   * keeps (reaction-losses.ts: a catalytic cycle, texts, arrow styles, a label such as HATU): the app reopens on it,
   * the text staying as written (from the text alone, a cycle came back as a straight route, two cycles as a refusal).
   * Absent when the text says it all, and dropped when the text is typed or a document loaded. */
  drawing?: string;
}

/**
 * The document the app opens on until a person changes it (an enumerated route, written as a step list of its mapped
 * steps with their agents): five numbered steps, agents over every arrow, every step's reaction centre coloured, a
 * stereocentre; the map numbers hidden (the toggle shows them).
 * Numbered through (crux-core scheme::number_through): one number per atom in every step, the same atom keeping its
 * number from step to step, the steps' own numbering renamed so (the source numbers each step from 1).
 */
const DEFAULT_EXAMPLE: DocumentReaction = {
  title: 'A route of five mapped steps: amide coupling, Boc removal, SNAr, Suzuki coupling, reductive amination',
  text: [
    '[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])[N:8]1[CH2:9][CH2:10][C@H:11]([NH2:12])[CH2:13]1.O[C:15]' +
      '(=[O:14])[c:17]1[cH:18][cH:19][cH:20][n:21][cH:22]1>CN(C)C(On1nnc2cccnc21)=[N+](C)C.F[P-](F)(F)(F)(F)F.C' +
      'CN(C(C)C)C(C)C.CN(C)C=O>[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])[N:8]1[CH2:9][CH2:10][C@H:11]([N' +
      'H:12][C:15](=[O:14])[c:17]2[cH:18][cH:19][cH:20][n:21][cH:22]2)[CH2:13]1',
    '[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])[N:8]1[CH2:9][CH2:10][C@H:11]([NH:12][C:15](=[O:14])[c:1' +
      '7]2[cH:18][cH:19][cH:20][n:21][cH:22]2)[CH2:13]1>O=C(O)C(F)(F)F.ClCCl>[NH:8]1[CH2:9][CH2:10][C@H:11]([NH' +
      ':12][C:15](=[O:14])[c:17]2[cH:18][cH:19][cH:20][n:21][cH:22]2)[CH2:13]1',
    'Cl[c:36]1[n:37][cH:38][c:39]([Br:40])[cH:41][n:42]1.[NH:8]1[CH2:9][CH2:10][C@H:11]([NH:12][C:15](=[O:14]' +
      ')[c:17]2[cH:18][cH:19][cH:20][n:21][cH:22]2)[CH2:13]1>CCN(C(C)C)C(C)C.CS(C)=O>[N:8]1([c:36]2[n:37][cH:38' +
      '][c:39]([Br:40])[cH:41][n:42]2)[CH2:9][CH2:10][C@H:11]([NH:12][C:15](=[O:14])[c:17]2[cH:18][cH:19][cH:20' +
      '][n:21][cH:22]2)[CH2:13]1',
    'OB(O)[c:43]1[cH:44][cH:45][c:46]([CH:47]=[O:48])[n:49][cH:50]1.[N:8]1([c:36]2[n:37][cH:38][c:39]([Br:40]' +
      ')[cH:41][n:42]2)[CH2:9][CH2:10][C@H:11]([NH:12][C:15](=[O:14])[c:17]2[cH:18][cH:19][cH:20][n:21][cH:22]2' +
      ')[CH2:13]1>[Pd].O=C([O-])[O-].[K+].[K+].C1COCCO1.O>[N:8]1([c:36]2[n:37][cH:38][c:39](-[c:43]3[cH:44][cH:' +
      '45][c:46]([CH:47]=[O:48])[n:49][cH:50]3)[cH:41][n:42]2)[CH2:9][CH2:10][C@H:11]([NH:12][C:15](=[O:14])[c:' +
      '17]2[cH:18][cH:19][cH:20][n:21][cH:22]2)[CH2:13]1',
    '[CH3:51][N:52]1[CH2:53][CH2:54][NH:55][CH2:56][CH2:57]1.[N:8]1([c:36]2[n:37][cH:38][c:39](-[c:43]3[cH:44' +
      '][cH:45][c:46]([CH:47]=[O:48])[n:49][cH:50]3)[cH:41][n:42]2)[CH2:9][CH2:10][C@H:11]([NH:12][C:15](=[O:14' +
      '])[c:17]2[cH:18][cH:19][cH:20][n:21][cH:22]2)[CH2:13]1>CC(=O)O[BH-](OC(C)=O)OC(C)=O.[Na+].CC(=O)O.ClCCCl' +
      '>[N:8]1([c:36]2[n:37][cH:38][c:39](-[c:43]3[cH:44][cH:45][c:46]([CH2:47][N:55]4[CH2:54][CH2:53][N:52]([C' +
      'H3:51])[CH2:57][CH2:56]4)[n:49][cH:50]3)[cH:41][n:42]2)[CH2:9][CH2:10][C@H:11]([NH:12][C:15](=[O:14])[c:' +
      '17]2[cH:18][cH:19][cH:20][n:21][cH:22]2)[CH2:13]1',
  ].join('--**--'),
  mapNumbers: false,
};

/**
 * The Examples menu: one feature or two an example (chosen from pictures of a collection of reactions and routes and
 * of Chem's templates; the routes are paracetamol from phenol and imatinib, the templates Chem's enumerator's).
 * Choosing one loads it into the document in place of what it held.
 */
const EXAMPLES: DocumentReaction[] = [
  {
    title: 'A route of three steps, the agents over each arrow: paracetamol from phenol',
    text: 'Oc1ccccc1.O=[N+]([O-])O>O>O=[N+]([O-])c1ccc(O)cc1>[H][H].[Pd].CCO>Nc1ccc(O)cc1' +
      '>CC(=O)OC(C)=O>CC(=O)Nc1ccc(O)cc1',
    mapNumbers: true,
  },
  {
    title: 'Atom-mapped: a Diels-Alder reaction, the bonds made in green, changed in blue',
    text: '[CH2:1]=[C:2]([CH3:7])[CH:3]=[CH2:4].[CH2:5]=[CH:6][C:8](=[O:9])[CH3:10]>>' +
      '[CH2:1]1[C:2]([CH3:7])=[CH:3][CH2:4][CH2:5][CH:6]1[C:8](=[O:9])[CH3:10]',
    mapNumbers: true,
  },
  {
    title: 'A reaction template drawn as a generic scheme: amide coupling, its constraints as footnotes',
    text: '[NX3;!$(N[C]=[S,O,N]);$(N[#6,S]);!$(N[#7]);H1,H2:1].[#6:4][C:5](=[O:6])[OH,O-]>>[N:1][C:5](=[O:6])[#6:4]',
    mapNumbers: true,
  },
  {
    title: 'Atom-mapped, the map numbers hidden: a Suzuki coupling with palladium and potassium carbonate',
    text: 'OB(O)[c:8]1[cH:7][cH:6][cH:5][c:4]([N+:2](=[O:1])[O-:3])[cH:12]1.' +
      'Br[c:25]1[cH:24][cH:23][c:22]([CH2:21][NH:20][C:18]([O:17][C:14]([CH3:13])([CH3:15])[CH3:16])=[O:19])' +
      '[cH:28][cH:27]1>[Pd].O=C([O-])[O-].[K+].[K+]>[O:1]=[N+:2]([O-:3])[c:4]1[cH:5][cH:6][cH:7][c:8]' +
      '(-[c:25]2[cH:24][cH:23][c:22]([CH2:21][NH:20][C:18]([O:17][C:14]([CH3:13])([CH3:15])[CH3:16])=[O:19])' +
      '[cH:28][cH:27]2)[cH:12]1 |f:2.3.4|',
    mapNumbers: false,
  },
  {
    title: 'A template with aryl groups and a halide: Suzuki coupling',
    text: '[c:1]-B(O)O.[Cl,Br,I][c:2]>>[c:1]-[c:2]',
    mapNumbers: true,
  },
  {
    title: 'Stereochemistry: a Mitsunobu esterification inverts the stereocentre',
    text: 'C[C@H](O)CCc1ccccc1.OC(=O)c1ccccc1>c1ccc(cc1)P(c1ccccc1)c1ccccc1.CC(C)OC(=O)N=NC(=O)OC(C)C>' +
      'C[C@@H](OC(=O)c1ccccc1)CCc1ccccc1.O',
    mapNumbers: true,
  },
  {
    title: 'Four components in one step: the Ugi reaction',
    text: 'O=Cc1ccccc1.Nc1ccccc1.CC(O)=O.CC(C)(C)[N+]#[C-]>CO>CC(=O)N(c1ccccc1)C(C(=O)NC(C)(C)C)c1ccccc1.O',
    mapNumbers: true,
  },
  {
    title: 'Isotopes, and a salt as one agent: phenol O-methylated with iodomethane-d3',
    text: '[2H]C([2H])([2H])I.Oc1ccccc1>[K+].[K+].[O-]C([O-])=O>[2H]C([2H])([2H])Oc1ccccc1.I |f:2.3.4|',
    mapNumbers: true,
  },
  {
    title: 'A convergent route of five steps, broken into rows: imatinib',
    text: [
      'CC(=O)c1cccnc1.COC(OC)N(C)C>>CN(C)/C=C/C(=O)c1cccnc1',
      'Cc1ccc([N+](=O)[O-])cc1N.N#CN>O=[N+]([O-])O.CCO>Cc1ccc([N+](=O)[O-])cc1NC(=N)N',
      'CN(C)/C=C/C(=O)c1cccnc1.Cc1ccc([N+](=O)[O-])cc1NC(=N)N>[Na+].[OH-].CC(C)O>' +
        'Cc1ccc([N+](=O)[O-])cc1Nc1nccc(-c2cccnc2)n1',
      'Cc1ccc([N+](=O)[O-])cc1Nc1nccc(-c2cccnc2)n1>[H][H].[Pd].C1CCOC1>Cc1ccc(N)cc1Nc1nccc(-c2cccnc2)n1',
      'CN1CCN(Cc2ccc(C(=O)Cl)cc2)CC1.Cc1ccc(N)cc1Nc1nccc(-c2cccnc2)n1>c1ccncc1>' +
        'Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1',
    ].join('--**--'),
    mapNumbers: true,
  },
];

/** A fresh copy of the default example: the document the app opens on, and Reset to example brings back. */
export function defaultExample(): DocumentReaction {
  return {...DEFAULT_EXAMPLE};
}

/** Fresh copies of the Examples menu's nine examples, in its order. */
export function examples(): DocumentReaction[] {
  return EXAMPLES.map((r) => ({...r}));
}

/** Whether `value` is a stored document: a title, a text and the map numbers' switch. */
function isReaction(value: any): value is DocumentReaction {
  return !!value && typeof value.title === 'string' && typeof value.text === 'string' &&
    typeof value.mapNumbers === 'boolean';
}

/** The document kept in this browser, or null: none kept, unreadable, of another version (a page of several reactions,
 * version 1), no storage, or one with no reaction typed in it: the app never opens empty. */
export function loadDocument(): DocumentReaction | null {
  try {
    const stored = JSON.parse(window.localStorage.getItem(DOCUMENT_KEY) ?? 'null');
    if (stored?.version !== DOCUMENT_VERSION || !isReaction(stored) || !stored.text.trim())
      return null;
    const doc: DocumentReaction = {title: stored.title, text: stored.text, mapNumbers: stored.mapNumbers};
    if (typeof stored.drawing === 'string' && stored.drawing.trim())
      doc.drawing = stored.drawing;
    return doc;
  } catch {
    return null;
  }
}

/** Keeps `doc` in this browser (the app reopens with it); false when the browser keeps nothing (storage blocked or
 * full). A drawing too large for the storage left is dropped, the text kept. */
export function saveDocument(doc: DocumentReaction): boolean {
  const kept = (drawing?: string) => JSON.stringify({version: DOCUMENT_VERSION, title: doc.title, text: doc.text,
    mapNumbers: doc.mapNumbers, ...(drawing ? {drawing} : {})});
  try {
    window.localStorage.setItem(DOCUMENT_KEY, kept(doc.drawing));
    return true;
  } catch {
    try {
      if (!doc.drawing)
        return false;
      window.localStorage.setItem(DOCUMENT_KEY, kept());
      return true;
    } catch {
      return false;
    }
  }
}

/** Forgets the document kept in this browser: the app opens on the example again. */
export function forgetDocument(): void {
  try {
    window.localStorage.removeItem(DOCUMENT_KEY);
  } catch {
    // nothing kept to forget
  }
}

type Scheme = InstanceType<ChemWasm['ChemReactionScheme']>;

/** Where a picture sits in its frame on the page: Crux's box (px; Crux fits the scheme into it, never enlarging it),
 * and the box's top relative to the frame's (the blank around the ink cut off). */
interface Frame {
  boxWidth: number;
  boxHeight: number;
  boxTop: number;
  height: number;
}

/** What the status says when Crux is off after a trap it could not recover from (engine/crux-module.ts
 * `noteCruxError`): its instance may be broken, so nothing more is asked of it, not even a `free()`. */
export const CRUX_OFF = 'Drawing stopped after an error in this session: the reaction is not drawn';
/** What the status says of a reaction too large to draw in time (a 1,000-residue reaction took 7 s, the tab frozen),
 * by the rule cells and the Reaction panel follow (crux-reaction-drawing.ts `tooLargeToDraw`). */
export const TOO_LARGE = 'The reaction is too large to draw';
/** What the status and the preview say of a text with no reaction arrow, a lone molecule (a molecule drawn alone
 * turned the preview into an error dump): one short line saying what to do. */
export const NO_ARROW = 'A reaction needs an arrow: draw one with the arrow tool, or type reactants>>products';

/**
 * Whether reaction texts `a` and `b` are one reaction (Arrange in the app must leave the Text exactly as it was): the
 * same steps as Crux reads them, each step's reactants, agents and products the same molecules
 * in any order within a side (canonical SMILES, map numbers kept). A drawing changed only in its layout (Arrange,
 * Clean Up, a move) gives its reaction anew in Crux Sketch's own writing, and the Text keeps the person's. A text with
 * a CX part, an RXN block or an RD file, or one too large to compare quickly, is compared as written.
 */
export function sameReactionText(crux: ChemWasm, a: string, b: string): boolean {
  if (a === b)
    return true;
  // a CX part that only groups ions into one agent (`|f:2.3|`, as a writer may group HATU's) holds no chemistry
  const body = (t: string) => t.replace(/\s+\|f:[\d.,]+\|\s*$/, '');
  [a, b] = [body(a), body(b)];
  const comparable = (t: string) => t.trim() !== '' && !t.includes('|') && !t.trimStart().startsWith('$') &&
    !tooLargeToDraw(t) && !noReactionArrow(t);
  if (!comparable(a) || !comparable(b) || cruxLoadFailed())
    return false;
  const canonical = (smiles: string): string | null => {
    let mol: ReturnType<ChemWasm['get_mol']> = null;
    try {
      mol = crux.get_mol(smiles);
      return mol ? mol.get_smiles() : null;
    } catch (e) {
      noteCruxError(e);
      return null;
    } finally {
      if (!cruxLoadFailed()) // a trapped instance is asked nothing more, not even a delete
        mol?.delete();
    }
  };
  const stepsOf = (t: string): string[][] | null => {
    let scheme: Scheme | null = null;
    try {
      scheme = new crux.ChemReactionScheme(t, '');
      const out: string[][] = [];
      for (const step of scheme.get_steps()) {
        const sides = step.split('>');
        if (sides.length !== 3)
          return null;
        const read = sides.map((side) => side.split('.').filter((c) => c.trim()).map(canonical));
        if (read.some((side) => side.some((c) => c === null)))
          return null;
        out.push(read.map((side) => (side as string[]).sort().join('.')));
      }
      return out;
    } catch (e) {
      noteCruxError(e);
      return null;
    } finally {
      if (!cruxLoadFailed())
        scheme?.free();
    }
  };
  const [x, y] = [stepsOf(a), stepsOf(b)];
  return x !== null && y !== null && JSON.stringify(x) === JSON.stringify(y);
}

/** `smiles` as canonical SMILES with its map numbers left out (a mapped atom keeps its hydrogens written, [CH3:1]:
 * the number dropped, the atom is as unmapped); null when Crux reads no molecule from it. */
function unmappedSmiles(crux: ChemWasm, smiles: string): string | null {
  let mol: ReturnType<ChemWasm['get_mol']> = null;
  try {
    mol = crux.get_mol(smiles.replace(/(\[[^\]:]+):\d+\]/g, '$1]'));
    return mol?.get_smiles() || null;
  } catch (e) {
    noteCruxError(e);
    return null;
  } finally {
    if (!cruxLoadFailed())
      mol?.delete();
  }
}

/** A side of a step's reaction SMILES as its molecules' texts (a CX part after it dropped). */
const sideParts = (side: string) => side.replace(/\s+\|.*$/, '').split('.').filter((c) => c.trim());

/** The steps of reaction text `text` as Crux reads them, each as `reactants>agents>products`; null when it reads none,
 * the text is too large to read quickly, or Crux is off. */
function readSteps(crux: ChemWasm, text: string): string[][] | null {
  if (!text.trim() || tooLargeToDraw(text) || cruxLoadFailed())
    return null;
  let scheme: Scheme | null = null;
  try {
    scheme = new crux.ChemReactionScheme(text, '');
    return scheme.get_steps().map((step) => step.split('>'));
  } catch (e) {
    noteCruxError(e);
    return null;
  } finally {
    if (!cruxLoadFailed())
      scheme?.free();
  }
}

/** The molecules of reaction text `text` as Crux reads them: every step's reactants and products (a lone molecule's
 * parts), each its canonical SMILES with its map numbers left out; the agents not (a solvent or a reagent two reactions
 * share makes them no kin); null when Crux reads none from it (a template's SMARTS), it is too large to read quickly,
 * or Crux is off. */
export function reactionMolecules(crux: ChemWasm, text: string): Set<string> | null {
  if (!text.trim() || tooLargeToDraw(text) || cruxLoadFailed())
    return null;
  const t = text.trim();
  // a lone molecule: a MOL block read whole, a SMILES by its parts
  const texts = noReactionArrow(text) ? (t.includes('\n') ? [text] : sideParts(t)) :
    readSteps(crux, text)?.flatMap((sides) => [...sideParts(sides[0] ?? ''), ...sideParts(sides[2] ?? '')]) ?? [];
  const out = new Set<string>();
  for (const m of texts) {
    const smiles = unmappedSmiles(crux, m);
    if (smiles)
      out.add(smiles);
  }
  return out.size ? out : null;
}

/** The last product of reaction text `text`: the last molecule of its last step's products, the one a next step
 * starts from (crux-sketch's addMolecule role 'next'), as canonical SMILES with its map numbers left out; null when it
 * reads none. */
export function lastProduct(crux: ChemWasm, text: string): string | null {
  const steps = readSteps(crux, text);
  const products = steps && steps.length ? sideParts(steps[steps.length - 1][2] ?? '') : [];
  return products.length ? unmappedSmiles(crux, products[products.length - 1]) : null;
}

/** Whether reaction texts `a` and `b` have a molecule in common ({@link reactionMolecules}); null when either reads as
 * none. */
export function shareAMolecule(crux: ChemWasm, a: string, b: string): boolean | null {
  const [x, y] = [reactionMolecules(crux, a), reactionMolecules(crux, b)];
  return x === null || y === null ? null : [...x].some((m) => y.has(m));
}

/** Whether `text` holds no reaction arrow: neither `>` nor an RXN or RD file, nor route JSON (a lone molecule, as
 * SMILES or a MOL block). */
export function noReactionArrow(text: string): boolean {
  const t = text.trimStart();
  return t !== '' && !t.includes('>') && !t.startsWith('$RXN') && !t.startsWith('$RDFILE') &&
    !/^(\{|\[\s*\{)/.test(t);
}

/** A document reaction read by Crux once, drawn as often as the page, the SVG or the PNG asks. A Crux call that throws
 * is reported to the module (`noteCruxError`: a trap it could not recover from turns Crux off), and once Crux is off
 * no call reaches it again: the picture is written as text, and what Crux held is let go without a `free()` on a
 * broken instance (which would be a double free). */
export class ReactionPicture {
  readonly reaction: DocumentReaction;
  // The reaction Crux read, or null with {@link error} saying why (it read none, could not draw it, it is too large,
  // or Crux is off).
  private scheme: Scheme | null = null;
  error: string | null = null;
  private readonly frames = new Map<number, Frame>();

  constructor(crux: ChemWasm, reaction: DocumentReaction) {
    this.reaction = reaction;
    if (!reaction.text.trim())
      this.error = 'Type or paste a reaction';
    else if (tooLargeToDraw(reaction.text))
      this.error = TOO_LARGE;
    else if (noReactionArrow(reaction.text))
      this.error = NO_ARROW;
    else if (cruxLoadFailed())
      this.error = CRUX_OFF;
    else {
      try {
        this.scheme = new crux.ChemReactionScheme(reaction.text, '');
      } catch (e) {
        noteCruxError(e);
        this.error = cruxLoadFailed() ? CRUX_OFF : cruxMessage(e);
      }
    }
  }

  // The scheme while Crux may be asked: null once Crux is off (what it held dropped, not freed).
  private live(): Scheme | null {
    if (this.scheme !== null && cruxLoadFailed()) {
      this.scheme = null;
      this.error = CRUX_OFF;
      this.frames.clear();
    }
    return this.scheme;
  }

  // A call on the scheme threw `e`: reported, and the scheme let go (freed only while Crux is still on).
  private lost(e: unknown): void {
    noteCruxError(e);
    const scheme = this.scheme;
    this.scheme = null;
    this.frames.clear();
    if (cruxLoadFailed())
      this.error = CRUX_OFF;
    else {
      this.error = cruxMessage(e);
      scheme?.free();
    }
  }

  /** The steps Crux reads (`get_steps`), none when it reads no reaction. */
  steps(): string[] {
    const scheme = this.live();
    if (scheme === null)
      return [];
    try {
      return scheme.get_steps();
    } catch (e) {
      this.lost(e);
      return [];
    }
  }

  // Crux's drawing options for this picture, besides its box: Chem's look, its reaction centre's widths too.
  private options(): object {
    return {bondLength: BOND_LENGTH, mapNumbers: this.reaction.mapNumbers ? 'show' : 'hide',
      drawOptions: chemDrawOptions(), ...CHEM_CENTRE};
  }

  /** The frame of the picture on a page `width` px wide: Crux's box is the page's width and the scheme's natural
   * height, or room for {@link MAX_ROWS} rows when the scheme is wider than the page; its ink measured on a canvas. */
  frame(width: number): Frame {
    const kept = this.frames.get(width);
    if (kept)
      return kept;
    let frame: Frame = {boxWidth: width, boxHeight: UNREAD_HEIGHT, boxTop: 0, height: UNREAD_HEIGHT};
    const scheme = this.live();
    if (scheme) {
      try {
        const natural = scheme.get_svg(JSON.stringify(this.options()))
          .match(/width='([\d.]+)px' height='([\d.]+)px'/);
        const [w, h] = natural ? [+natural[1], +natural[2]] : [width, UNREAD_HEIGHT];
        const rows = Math.min(MAX_ROWS, Math.ceil(w / width));
        const boxHeight = Math.ceil(rows === 1 ? h : h * rows * ROW_ROOM);
        const canvas = new OffscreenCanvas(width, boxHeight);
        const g = canvas.getContext('2d', {willReadFrequently: true})! as unknown as CanvasRenderingContext2D;
        scheme.draw_to_canvas(canvas, JSON.stringify({width, height: boxHeight, ...this.options()}));
        const [top, bottom] = inkRows(g.getImageData(0, 0, width, boxHeight));
        frame = bottom < top ? {boxWidth: width, boxHeight, boxTop: 0, height: boxHeight} :
          {boxWidth: width, boxHeight, boxTop: PICTURE_MARGIN - top, height: bottom - top + 1 + 2 * PICTURE_MARGIN};
      } catch (e) { // read but not drawn: written as text, as a text Crux reads nothing from
        this.lost(e);
      }
    }
    this.frames.set(width, frame);
    return frame;
  }

  /** Draws the picture into the frame at `x, y` (CSS px) of `g`'s canvas, `width` px wide, at `ratio` device px per
   * CSS px; a text Crux reads no reaction from is written as it is, in grey. */
  draw(g: CanvasRenderingContext2D, x: number, y: number, width: number, ratio: number): void {
    const f = this.frame(width);
    g.save();
    try {
      g.setTransform(ratio, 0, 0, ratio, 0, 0);
      g.beginPath();
      g.rect(x, y, width, f.height);
      g.clip();
      const scheme = this.live();
      if (scheme) {
        try {
          scheme.draw_to_canvas(g.canvas, JSON.stringify({width: f.boxWidth, height: f.boxHeight,
            offsetx: x, offsety: y + f.boxTop, ...this.options()}));
          return;
        } catch (e) {
          this.lost(e);
        }
      }
      g.font = UNREAD_FONT;
      g.fillStyle = UNREAD_COLOUR;
      g.textBaseline = 'middle';
      g.fillText(this.reaction.text.split('\n')[0], x + 4, y + f.height / 2);
    } finally {
      g.restore();
    }
  }

  /** The picture as an SVG element for a page `width` px wide, its frame's top left at `x, y`: Crux's SVG of the box,
   * moved so that the frame's part shows. */
  svg(x: number, y: number, width: number): string {
    const f = this.frame(width);
    const scheme = this.live();
    if (scheme) {
      try {
        const svg = scheme.get_svg(JSON.stringify({width: f.boxWidth, height: f.boxHeight, ...this.options()}));
        return svg.replace(/^<\?xml[^>]*>\s*/, '').replace('<svg', `<svg x='${x}' y='${y + f.boxTop}'`);
      } catch (e) {
        this.lost(e);
      }
    }
    return `<text x='${x + 4}' y='${y + f.height / 2}' dominant-baseline='middle' font-family='monospace' ` +
      `font-size='12' fill='${UNREAD_COLOUR}'>${escapeXml(this.reaction.text.split('\n')[0])}</text>`;
  }

  /** Frees what Crux holds for the picture: never on an instance a trap may have broken (Crux off), where it is only
   * dropped. */
  free(): void {
    const scheme = this.scheme;
    this.scheme = null;
    if (scheme !== null && !cruxLoadFailed())
      scheme.free();
  }
}

/** The first and last rows of `image` with ink (alpha above 0); last < first when it has none. */
export function inkRows(image: ImageData): [number, number] {
  const {width, height, data} = image;
  const inked = (row: number) => {
    for (let i = row * width * 4 + 3, end = (row + 1) * width * 4; i < end; i += 4) {
      if (data[i] > 0)
        return true;
    }
    return false;
  };
  let top = 0;
  while (top < height && !inked(top))
    top++;
  let bottom = height - 1;
  while (bottom > top && !inked(bottom))
    bottom--;
  return top < height ? [top, bottom] : [0, -1];
}

/** What the app says of a value whose `|f:…|` groups Crux does not read (its reader's "the fragment group f:[0, 1, 2]
 * names components out of range or of two roles", "cannot read the fragment group f:…"): in words, as the app shows
 * no raw engine text (the first once showed in its right pane). */
export const FRAGMENT_GROUPS_UNREAD = 'Its salt and mixture groups (the |f:…| part after the reaction) do not match ' +
  'its molecules';

/** What Crux says when it throws, without its binding's prefix; a message in Crux's own notation in words. */
export function cruxMessage(e: unknown): string {
  const message = String(e instanceof Error ? e.message : e).replace(/^ChemReactionScheme: /, '');
  if (/fragment group f:/.test(message)) {
    console.warn(`Chem | Crux did not read the reaction: ${message}`);
    return FRAGMENT_GROUPS_UNREAD;
  }
  return message;
}

function escapeXml(text: string): string {
  return text.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/'/g, '&apos;');
}

/** The export's layout: each picture's top, after its title, on a page `width` px wide; and the page's height. */
function exportLayout(pictures: ReactionPicture[], width: number): {tops: number[], height: number} {
  const tops: number[] = [];
  let y = PAGE_PADDING;
  for (const p of pictures) {
    tops.push(y + TITLE_HEIGHT);
    y += TITLE_HEIGHT + p.frame(width).height + EXPORT_GAP;
  }
  return {tops, height: y - EXPORT_GAP + PAGE_PADDING};
}

/** The document as one SVG, `width` px wide (the pictures' width, the page's padding around): each reaction's title
 * as text and its picture, one under another, on white. */
export function documentSvg(pictures: ReactionPicture[], width: number): string {
  const {tops, height} = exportLayout(pictures, width);
  const total = width + 2 * PAGE_PADDING;
  const parts = pictures.map((p, i) =>
    `<text x='${PAGE_PADDING}' y='${tops[i] - TITLE_HEIGHT / 2}' dominant-baseline='middle' ` +
    `font-family='${TITLE_FONT_FAMILY}' font-size='${TITLE_FONT_SIZE}' font-weight='500' fill='${TITLE_COLOUR}'>` +
    `${escapeXml(p.reaction.title)}</text>\n${p.svg(PAGE_PADDING, tops[i], width)}`);
  return `<?xml version='1.0' encoding='UTF-8'?>\n<svg xmlns='http://www.w3.org/2000/svg' version='1.1' ` +
    `width='${total}px' height='${height}px' viewBox='0 0 ${total} ${height}'>\n` +
    `<rect x='0' y='0' width='${total}' height='${height}' fill='#FFFFFF'/>\n${parts.join('\n')}\n</svg>\n`;
}

/** The document as a PNG at {@link PNG_DPI} (less when the page is too tall for a canvas at that resolution, its
 * pHYs chunk saying which), laid out as {@link documentSvg}. */
export async function documentPng(pictures: ReactionPicture[], width: number): Promise<Blob> {
  const {tops, height} = exportLayout(pictures, width);
  const total = width + 2 * PAGE_PADDING;
  const ratio = Math.min(PNG_DPI / 96, CANVAS_MAX_SIDE / Math.max(total, height));
  const canvas = document.createElement('canvas');
  canvas.width = Math.round(total * ratio);
  canvas.height = Math.round(height * ratio);
  const g = canvas.getContext('2d')!;
  g.fillStyle = '#FFFFFF';
  g.fillRect(0, 0, canvas.width, canvas.height);
  pictures.forEach((p, i) => {
    g.setTransform(ratio, 0, 0, ratio, 0, 0);
    g.font = `500 ${TITLE_FONT_SIZE}px ${TITLE_FONT_FAMILY}`;
    g.fillStyle = TITLE_COLOUR;
    g.textBaseline = 'middle';
    g.fillText(p.reaction.title, PAGE_PADDING, tops[i] - TITLE_HEIGHT / 2);
    p.draw(g, PAGE_PADDING, tops[i], width, ratio);
  });
  const png = await new Promise<Blob>((resolve, reject) =>
    canvas.toBlob((b) => b ? resolve(b) : reject(new Error('the PNG was not written')), 'image/png'));
  return withResolution(png, ratio * 96);
}

/** `png` with a pHYs chunk after its header saying `dpi` (as pixels per metre), so that a document it is put into
 * sizes it at that resolution. */
export async function withResolution(png: Blob, dpi: number): Promise<Blob> {
  const bytes = new Uint8Array(await png.arrayBuffer());
  const ihdrEnd = 8 + 25; // the signature, then IHDR: length, type, 13 bytes of data, CRC
  const chunk = new Uint8Array(21);
  const view = new DataView(chunk.buffer);
  const perMetre = Math.round(dpi / 0.0254);
  view.setUint32(0, 9);
  chunk.set([0x70, 0x48, 0x59, 0x73], 4); // pHYs
  view.setUint32(8, perMetre);
  view.setUint32(12, perMetre);
  chunk[16] = 1; // the unit: the metre
  view.setUint32(17, crc32(chunk.subarray(4, 17)));
  return new Blob([bytes.subarray(0, ihdrEnd), chunk, bytes.subarray(ihdrEnd)], {type: 'image/png'});
}

/** PNG's CRC-32 of `bytes` (ISO 3309, the PNG specification's annex D). */
function crc32(bytes: Uint8Array): number {
  let c = 0xffffffff;
  for (const b of bytes) {
    c ^= b;
    for (let k = 0; k < 8; k++)
      c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
  }
  return (c ^ 0xffffffff) >>> 0;
}
