/**
 * Apps | Chem | Sketcher: the Crux sketcher filling its view in molecule mode, the context panel open beside it and
 * the help panel shut (as the Reaction Sketcher shuts it).
 * - The sketcher is the platform's sketcher host (`DG.chem.Sketcher`) given its own sketcher, Crux, whatever the
 *   session's, and no input row (`{sketcherType: 'Crux', showInput: false}`): the host makes the drawing the current
 *   object at every change (its `syncCurrentObject`), so the context panel shows the drawn molecule's panes, and it
 *   shows a molecule it cannot read as "Malformed molecule" (Chem's `validateMolecule`, as the molecule cell editor).
 *   Crux opens in molecule (normal) mode, or for a link holding a query in query mode (`openCruxInQueryMode`); the
 *   gear's Query mode switches it.
 * - The link: the app's function takes the query parameter `smiles` (core's routing hands an app its URL's parameters
 *   by name), a SMILES, CXSMILES or CXSMARTS, and opens on it; with none, on the example (atorvastatin). Every change
 *   of the drawing writes it back through the view's `path`, which core turns into the address (`history.replaceState`:
 *   the address changes in place, Back leaves the app): the drawing's CXSMILES without coordinates, or for a drawing
 *   with query features no SMILES holds, its CXSMARTS (the SMARTS with what it cannot hold: an enhanced-stereo group,
 *   an attachment point, a link node) or its SMARTS ({@link linkValue}, {@link queryLink}).
 * - The ribbon's Export... opens the Export dialog on the drawing (widgets/molecule-export.ts).
 */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';
import type {SketcherHandle} from '../../../vendor/crux-sketch/crux-sketch';
import {CruxSketcher, openCruxInQueryMode} from '../../crux/crux-sketcher';
import {cruxModuleSync, loadCrux} from '../../engine/crux-module';
import {MoleculeExport, openMoleculeExport} from '../../widgets/molecule-export';
import {sketcherHelpUrl} from './sketcher-help';
import {PackageFunctions} from '../../package';

/** The app's name: its view's, and Apps' (package.ts sketcherApp). */
export const APP_NAME = 'Sketcher';
/** The query parameter the link holds the drawing in (the app function's input of that name: core's routing passes
 * an app's URL parameters by name). */
export const LINK_PARAMETER = 'smiles';
/** The molecule the app opens on without a link: atorvastatin, its two stereocentres as drawn (3R,5R). */
export const EXAMPLE = 'CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O';
/** The view's class (what the tests find the app's view by). */
export const APP_CLASS = 'chem-sketcher-app';
/** The ribbon's button that opens the Export dialog. */
export const EXPORT_BUTTON = 'Export...';
/** The sketcher the app shows. */
export const CRUX = 'Crux';
/** The longest link value written, URL-encoded (characters): a server takes a request line of 8 KB by default
 * (nginx), so a longer address would open as an error. A drawing over it leaves the address without a molecule, and
 * says so once. */
export const MAX_LINK_LENGTH = 6000;
/** What the app says when the drawing grows past {@link MAX_LINK_LENGTH}. */
export function tooLongForLink(length: number): string {
  return `This drawing is too large for the page address (${length.toLocaleString('en-US')} characters, ` +
    `${MAX_LINK_LENGTH.toLocaleString('en-US')} at most), so the address no longer holds it: use Export... to share it`;
}
/** What the app says when the Crux sketcher does not start (the error goes to the console). */
export const SKETCHER_FAILED = 'The Crux sketcher could not start here.';

/** The link value of the drawing `sketch` holds: '' for an empty drawing; its CXSMILES without coordinates (Crux's
 * `cxsmiles` getter, as RDKit.js's `get_cxsmiles` writes it with every field but the coordinates); a drawing with a
 * query feature no SMILES holds (an atom list, a query bond, a custom query: `cxsmiles` is '') its {@link queryLink},
 * by Chem's Crux module `module` (null: the SMARTS alone). Null when Crux cannot write it (the address is then left as
 * it is). */
export function linkValue(sketch: SketcherHandle, module: CxSmartsModule | null = null): string | null {
  try {
    if (sketch.isEmpty)
      return '';
    const cxsmiles = sketch.cxsmiles;
    return cxsmiles !== '' ? cxsmiles : queryLink(sketch, module);
  } catch (e) {
    console.warn('Chem | the drawing was not written for the page address', e);
    return null;
  }
}

/** The calls of Chem's Crux module a query's link takes: chem-wasm's `get_cxsmarts(details)` reads the `CX_*` keys
 * (crux-core e881b8e98; RDKit.js's reads none and writes the coordinates too). */
export interface CxSmartsModule {
  get_qmol(input: string): {get_cxsmarts(details?: string): string, delete(): void} | null;
}

/** The CXSMILES fields a query's link carries: every one but the coordinates and the atom properties. RDKit reads a
 * query molfile's SMARTSQ into the atom property "MRV SMA" and writes it back unescaped: a comma or a colon in it makes
 * the extension unreadable, by RDKit.js and Crux alike, and a text such as `a:6` in it reads as a field of its own (an
 * ABS group), so the properties go but for the attachment points ({@link ATTACHMENT_POINT}). */
export const QUERY_LINK_FIELDS = JSON.stringify(Object.fromEntries(['CX_ATOM_LABELS', 'CX_MOLFILE_VALUES',
  'CX_RADICALS', 'CX_LINKNODES', 'CX_ENHANCEDSTEREO', 'CX_SGROUPS', 'CX_POLYMER', 'CX_BOND_CFG', 'CX_BOND_ATROPISOMER',
  'CX_COORDINATE_BONDS', 'CX_HYDROGEN_BONDS', 'CX_ZERO_BONDS'].map((field) => [field, true])));
/** An attachment point among a query's atom properties as RDKit writes them (`<atom>.molAttchpt.<1, 2 or -1>`). */
const ATTACHMENT_POINT = /(?:^|:)(\d+)\.molAttchpt\.(-?\d+)(?=:|$)/g;

/** A CXSMILES's or CXSMARTS's extension without its bars; '' for a text without one. */
function extensionOf(cx: string): string {
  const bar = cx.indexOf(' |');
  return bar < 0 || !cx.endsWith('|') ? '' : cx.slice(bar + 2, -1);
}

/** A query's link (main's F1, run19): its CXSMARTS, Chem's Crux module's `get_qmol` of Crux's `molV3000` written by
 * `get_cxsmarts` with {@link QUERY_LINK_FIELDS} and the attachment points added back as RDKit writes them, so that the
 * link keeps what SMARTS cannot hold (an AND or OR group, an attachment point, a link node). Where that SMARTS part is
 * not Crux's own `smarts` (which folds drawn hydrogens into their atoms as Datagrok's search reads them, writes an R
 * label as `[*:1]`, perceives aromaticity as it reads its own molfile back), or without the module, the SMARTS alone:
 * the link holds what Crux writes, and reopens as it. */
export function queryLink(sketch: SketcherHandle, module: CxSmartsModule | null): string {
  const smarts = sketch.smarts;
  if (module === null)
    return smarts;
  let query: ReturnType<CxSmartsModule['get_qmol']> = null;
  try {
    query = module.get_qmol(sketch.molV3000);
    if (query === null)
      return smarts;
    const cx = query.get_cxsmarts(QUERY_LINK_FIELDS);
    const bar = cx.indexOf(' |');
    if ((bar < 0 ? cx : cx.slice(0, bar)) !== smarts)
      return smarts;
    const props = extensionOf(query.get_cxsmarts(JSON.stringify({CX_ATOM_PROPS: true}))).replace(/^atomProp:/, '');
    const points = [...props.matchAll(ATTACHMENT_POINT)].map((m) => `${m[1]}.molAttchpt.${m[2]}`);
    if (points.length === 0)
      return cx;
    const fields = extensionOf(cx);
    return `${smarts} |atomProp:${points.join(':')}${fields ? `,${fields}` : ''}|`;
  } catch (e) {
    console.warn('Chem | the CXSMARTS of the query was not written for the page address, its SMARTS is', e);
    return smarts;
  } finally {
    query?.delete();
  }
}

/** The view's path for link value `value` (the app's own address before it, from core): `?smiles=<value, URL-encoded>`,
 * `?smiles=` for an empty drawing; '' (the app's address alone) for null. */
export function linkPath(value: string | null): string {
  return value === null ? '' : `?${LINK_PARAMETER}=${encodeURIComponent(value)}`;
}

/** The app's view, on `smiles` (the link's parameter: SMILES, CXSMILES or CXSMARTS; undefined or null, no parameter:
 * the example; '', an empty drawing). */
export async function sketcherAppView(smiles?: string | null): Promise<DG.View> {
  showSidePanels();
  const view = DG.View.create();
  view.name = APP_NAME;
  view.box = true;
  view.root.classList.add(APP_CLASS);
  view.helpUrl = sketcherHelpUrl();
  new SketcherApp(view, smiles ?? null);
  return view;
}

/** The shell's panels beside the view as the app wants them: the context panel open (it shows the drawing's panes),
 * the help panel shut, as the Reaction Sketcher shuts it (reactions-app.ts `shutSidePanels`); opened again, the help
 * panel shows the view's page (its `helpUrl`). A person opens or shuts them again as ever. */
export function showSidePanels(): void {
  grok.shell.windows.showContextPanel = true;
  grok.shell.windows.showHelp = false;
}

/** The key an app hangs on its view's root by, the tests reach an app the platform opened by: a registered symbol,
 * as Chem's package tests are a bundle of their own (a module's WeakMap there is not the app's). */
const APP_KEY = Symbol.for('datagrok.chem.sketcher-app');
type AppRoot = HTMLElement & {[APP_KEY]?: SketcherApp};

/** The app whose view is `view`, or undefined. */
export function sketcherAppOf(view: DG.ViewBase | null | undefined): SketcherApp | undefined {
  return view ? (view.root as AppRoot)[APP_KEY] : undefined;
}

/** The app on its view: the sketcher host, the link, the ribbon. */
export class SketcherApp {
  readonly host: DG.chem.Sketcher;
  /** Crux's handle once the sketcher is up. */
  sketch: SketcherHandle | null = null;
  // whether the drawing's link is over MAX_LINK_LENGTH (said once each time it gets there)
  private tooLong = false;
  private readonly subs = new Subscription();

  constructor(readonly view: DG.View, smiles: string | null) {
    (view.root as AppRoot)[APP_KEY] = this;
    // Crux whatever the session's sketcher (which the host leaves as it is), and no input row: the sketcher takes the
    // whole view
    this.host = new DG.chem.Sketcher(undefined, (s: string) => PackageFunctions.validateMolecule(s),
      {sketcherType: CRUX, showInput: false});
    Object.assign(this.host.root.style, {width: '100%', height: '100%'});
    // the example, or the link's molecule (a query given as SMARTS is set as SMARTS: the host's setMolecule asks Chem's
    // isSmarts); an empty parameter, an empty drawing
    const start = smiles === null ? EXAMPLE : smiles.trim();
    // molecule (normal) mode (the owner, 2026-11-19); a link holding a query, which the host sets as SMARTS (Chem's
    // isSmarts, as its setMolecule asks), opens in query mode, so that a shared query reopens as one
    if (smiles !== null && start !== '' && PackageFunctions.isSmarts(start))
      openCruxInQueryMode(this.host);
    if (start !== '')
      this.host.setMolecule(start);
    // the address as given while the app opens (the platform writes the view's path into the address as it shows the
    // view): a link Crux cannot read stays in the address as given until the first change
    if (smiles !== null)
      this.view.path = this.written = linkPath(smiles);
    this.given = start;
    this.showCrux();
    view.append(this.host.root);
    view.subs.push(this.subs);
    view.setRibbonPanels([[
      ui.button(EXPORT_BUTTON, () => this.openExport(), 'The molecule as SVG, PNG, CDXML, MOL, SDF, SMILES, ' +
        'CXSMILES, SMARTS, InChI or InChIKey: copied or downloaded'),
    ]]);
    this.subs.add(this.host.onChanged.subscribe(() => this.writeLink()));
    this.subs.add(grok.events.onViewRemoved.subscribe((removed) => {
      if (removed?.root === view.root)
        this.dispose();
    }));
  }

  // Crux once the host has it up: its handle, the panels as the app wants them, the opening molecule's link
  private showCrux(): void {
    this.ready = this.host.sketcherReady().then((ready) => {
      if (!(ready instanceof CruxSketcher) || ready.handle === null) {
        console.error(`Chem | the Sketcher app shows ${ready?.type ?? 'no sketcher'}, not Crux`);
        if (!(ready instanceof CruxSketcher)) // a CruxSketcher that did not start says so itself
          this.view.append(ui.divText(SKETCHER_FAILED));
        return null;
      }
      this.sketch = ready.handle;
      // the panels again, now that the view is up: the shell ignores the settings while it is opening a view
      showSidePanels();
      // the opening molecule's link (a link written otherwise, as OCC for CCO, is rewritten); one Crux did not read
      // leaves the address as given
      if (this.given === '' || !this.sketch.isEmpty)
        this.writeLink();
      return this.sketch;
    });
  }

  // the molecule the app opened on (the link's, or the example; '' for none)
  private given = '';

  /** Crux's handle once the host's sketcher is up and the link written; null when it is not Crux or did not start. */
  ready: Promise<SketcherHandle | null> = Promise.resolve(null);

  /** The drawing's link written through the view's path ({@link linkValue}): over {@link MAX_LINK_LENGTH}, the app's
   * address alone, said once. Nothing while the sketcher is not up, nor when Crux cannot write it. A query's CXSMARTS
   * takes Chem's Crux module, which Crux mode loads with the package and RDKit mode on the first query drawn here:
   * until it is up the link is the SMARTS, then written again. */
  writeLink(): void {
    if (this.sketch === null)
      return;
    const module = cruxModuleSync() as unknown as CxSmartsModule | null;
    if (module === null && !this.cruxAsked && this.sketch.hasQuery) {
      this.cruxAsked = true;
      loadCrux('Sketcher app: the link of a query').then(() => this.writeLink(), () => {});
    }
    const value = linkValue(this.sketch, module);
    if (value === null)
      return;
    let path = linkPath(value);
    const length = encodeURIComponent(value).length;
    if (length > MAX_LINK_LENGTH) {
      path = '';
      if (!this.tooLong)
        grok.shell.warning(tooLongForLink(length));
      this.tooLong = true;
    } else
      this.tooLong = false;
    if (path !== this.written) {
      this.written = path;
      this.view.path = path;
    }
  }

  // the path last written (writeLink); null before the first
  private written: string | null = null;
  // whether writeLink asked for Chem's Crux module (RDKit mode, a query drawn)
  private cruxAsked = false;

  /** The drawing's exports, as the Export dialog makes them (each format in each style); null before Crux is up. Made
   * by the app's own code, which the package tests reach through the app: their bundle's copy of Crux Sketch cannot
   * make the offscreen sketcher the ACS style lays the molecule out on (package.js's has defined its element). */
  exporter(): MoleculeExport | null {
    return this.sketch === null ? null : new MoleculeExport(this.sketch);
  }

  private openExport(): void {
    if (this.sketch === null) {
      grok.shell.warning('The sketcher is not ready yet');
      return;
    }
    openMoleculeExport(this.sketch);
  }

  // The view closed: Crux let go (the host's detach destroys its sketcher), the subscriptions ended. The platform ends
  // no plain view's `subs` when it closes, so the app goes when the platform says the view was removed.
  private dispose(): void {
    this.subs.unsubscribe();
    this.host.detach();
  }
}
