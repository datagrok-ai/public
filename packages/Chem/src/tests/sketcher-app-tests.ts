/**
 * Category `Sketcher: app`: Apps | Chem | Sketcher (utils/sketcher-app/sketcher-app.ts) opened as the platform opens
 * it, from its address (`grok.shell.route`: the address pushed, core's routing runs the app with the URL's parameters):
 * the example or the link's molecule in Crux whatever the session's sketcher (left as it is), in molecule (normal) mode
 * unless the link holds a query, the gear's Query mode toggling it, the
 * sketcher filling the view with no input row (the host's `showInput: false`), the context panel open on the drawing
 * (the current object) and the help panel shut, and the address rewritten by core from the view's path at each change;
 * and the Export dialog (widgets/molecule-export.ts), every format in both styles.
 */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {after, awaitCheck, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import type {SketcherHandle} from '../../vendor/crux-sketch/crux-sketch';
import {APP_CLASS, APP_NAME, EXAMPLE, EXPORT_BUTTON, LINK_PARAMETER, MAX_LINK_LENGTH, SketcherApp, linkPath,
  sketcherAppOf, tooLongForLink} from '../utils/sketcher-app/sketcher-app';
import {ACS_1996, AS_DRAWN, CDXML, COPY, CXSMILES, EXPORT_TITLE, FORMATS, FORMAT_FILES, INCHI, INCHIKEY, MOL_V2000,
  MOL_V3000, NOTE_CLASS, PNG, PNG_DPI, PREVIEW_CLASS, QUERY_HAS_NO_SMILES, SAVE, SDF, SMARTS, SMILES,
  STYLED_FORMATS, SVG, TEXT_CLASS, moleculeExportOf} from '../widgets/molecule-export';
import {SKETCHER_HELP_FILE, sketcherHelpUrl} from '../utils/sketcher-app/sketcher-help';

/** The app's address. */
const ADDRESS = '/apps/Chem/Sketcher';
/** A molecule with an AND group, which CXSMILES keeps and SMILES leaves out. */
const ENHANCED = 'C[C@H](O)[C@@H](C)Cl |&1:1,3|';
/** A query: an atom in a ring bonded to N or O; its SMARTS as Crux (RDKit's MolToSmarts) writes it. */
const QUERY = '[#6;R]-[#7,#8]';
const QUERY_SMARTS = '[#6&R]-[#7,#8]';
/** Queries holding what SMARTS cannot: an AND group on a stereocentre beside an atom list (its link writes the centre's
 * H count, as Crux's `smarts` does), and an attachment point on a query's atom (its link as given). */
const QUERY_AND = '[#6]-[#6@H](-[#8])-[#7,#8] |&1:1|';
const QUERY_AND_LINK = '[#6]-[#6@&H1](-[#8])-[#7,#8] |&1:1|';
const QUERY_ATTACHED = '[#6]-[#6]-[#7,#8] |atomProp:2.molAttchpt.1|';
/** Atorvastatin's standard InChIKey. */
const EXAMPLE_INCHIKEY = 'XUKUURHRXDUEBC-KAYWLYCHSA-N';
/** 3-Aminopropanol drawn as no layout draws it: its five atoms on one vertical line, 1.5 apart. */
const STRAIGHT_UP = ['', '  straight', '', '  5  4  0  0  0  0  0  0  0  0999 V2000',
  ...['N', 'C', 'C', 'C', 'O'].map((e, i) => `${'0.0000'.padStart(10)}${(i * 1.5).toFixed(4).padStart(10)}` +
    `${'0.0000'.padStart(10)} ${e.padEnd(3)} 0  0  0  0  0  0  0  0  0  0  0  0`),
  '  1  2  1  0', '  2  3  1  0', '  3  4  1  0', '  4  5  1  0', 'M  END', ''].join('\n');

/** Test timeouts: the app opens from its address in about 4 s on the stand (the platform's routing, the view, Crux),
 * so a test opening it once gets a minute; three times, or exporting every format, two. */
const OPENS_ONE = {timeout: 60000};
const OPENS_THREE = {timeout: 120000};
const EXPORTS = {timeout: 120000};

let rdkit: any = null;
/** `smiles` as RDKit writes it canonically (no CX part), or '' when RDKit does not read it. */
function canonical(smiles: string): string {
  const mol = rdkit.get_mol(smiles);
  try {
    return mol?.is_valid() ? mol.get_smiles() : '';
  } finally {
    mol?.delete();
  }
}

/** Opens the app from `address` as a link does (the address pushed, core's routing runs the app with its
 * parameters); resolves with the app once its Crux sketcher is up. */
async function openAt(address: string): Promise<SketcherApp> {
  const before = new Set(Array.from(grok.shell.views).map((v) => v.root));
  grok.shell.route(address);
  let app: SketcherApp | undefined;
  await awaitCheck(() => {
    const v = grok.shell.v;
    app = v && v.name === APP_NAME && !before.has(v.root) ? sketcherAppOf(v) : undefined;
    return !!app;
  }, `the app opens from ${address}`, 25000);
  const sketch = await app!.ready;
  expect(sketch !== null, true, 'Crux is up in the app');
  return app!;
}

/** The address bar's path and query. */
const address = () => `${location.pathname}${location.search}`;
/** The address of link value `value`. */
const addressOf = (value: string) => `${ADDRESS}${linkPath(value)}`;

/** Waits until the address is `want`. */
async function addressBecomes(want: string, what: string): Promise<void> {
  await awaitCheck(() => address() === want, `${what}: the address is ${want}, not ${address()}`, 10000);
}

/** The molecule the context panel shows: the current object's value when it is a Molecule semantic value. */
function currentMolecule(): string | null {
  const o = grok.shell.o;
  return o instanceof DG.SemanticValue && o.semType === DG.SEMTYPE.MOLECULE ? String(o.value) : null;
}

function close(app: SketcherApp | null): void {
  app?.view.close();
}

category('Sketcher: app', () => {
  let panel = false;
  let help = false;
  let sketcherType = '';
  before(async () => {
    rdkit = await grok.functions.call('Chem:getRdKitModule');
    panel = grok.shell.windows.showContextPanel;
    help = grok.shell.windows.showHelp;
    sketcherType = DG.chem.currentSketcherType;
  });

  after(async () => {
    grok.shell.windows.showContextPanel = panel;
    grok.shell.windows.showHelp = help;
    DG.chem.currentSketcherType = sketcherType;
  });

  test('registered: Apps | Chem | Sketcher, an app with the link\'s smiles input', async () => {
    const f = DG.Func.find({package: 'Chem', name: 'sketcherApp'})[0];
    expect(!!f, true, 'Chem:sketcherApp');
    expect(f.friendlyName, APP_NAME);
    expect(f.options['role'], 'app');
    expect(f.options['browsePath'], 'Chem', 'in Apps | Chem');
    expect(f.options['icon'], 'images/sketcher.png');
    expect(f.inputs.map((p) => p.name).join(','), LINK_PARAMETER, 'the link\'s parameter is the app\'s input');
  });

  test('opened from its address alone: the example in Crux, molecule (normal) mode, filling the view with no input ' +
    'row and no scrollbar, the context panel open on it and the help panel shut, the address its link', async () => {
    let app: SketcherApp | null = null;
    grok.shell.windows.showContextPanel = false;
    grok.shell.windows.showHelp = true;
    try {
      app = await openAt(ADDRESS);
      const sketch = app.sketch!;
      expect(sketch.cxsmiles, EXAMPLE, 'the example, atorvastatin, its stereo kept');
      expect(sketch.mode, 'molecule', 'molecule (normal) mode: the owner, 2026-11-19 18:0x');
      expect(grok.shell.windows.showContextPanel, true, 'the context panel opened');
      expect(grok.shell.windows.showHelp, false, 'the help panel shut, as the Reaction Sketcher shuts it');
      await awaitCheck(() => currentMolecule() !== null && canonical(currentMolecule()!) === canonical(EXAMPLE),
        'the drawing is the current object', 10000);
      await addressBecomes(addressOf(EXAMPLE), 'opened');
      // the host's input row (its molecule field and ≡ menu) left out by the host (showInput: false), not hidden by
      // the app; Crux takes the view, nothing scrolls
      const root = app.view.root;
      expect(root.classList.contains(APP_CLASS), true);
      expect(root.querySelector('.grok-sketcher-input'), null, 'the host\'s input row is not in the page');
      expect(app.host.molInput.isConnected, false, 'the host\'s molecule field is not in the page');
      const crux = root.querySelector('crux-sketch') as HTMLElement;
      const [v, c] = [root.getBoundingClientRect(), crux.getBoundingClientRect()];
      expect(Math.abs(v.width - c.width) <= 2 && Math.abs(v.height - c.height) <= 2, true,
        `Crux fills the view: ${c.width}x${c.height} in ${v.width}x${v.height}`);
      for (const e of [root, ...Array.from(root.querySelectorAll<HTMLElement>('*')).filter((x) => crux.contains(x) ||
        x.contains(crux))].filter((x) => !crux.contains(x) || x === crux)) {
        expect(e.scrollHeight <= e.clientHeight + 1 && e.scrollWidth <= e.clientWidth + 1, true,
          `no scrollbar on ${e.className || e.tagName}: ${e.scrollWidth}x${e.scrollHeight} in ` +
          `${e.clientWidth}x${e.clientHeight}`);
      }
      expect(app.view.helpUrl, sketcherHelpUrl());
    } finally {
      close(app);
    }
  }, OPENS_ONE);

  // The owner (2026-11-19, decision 11): no workaround: the app once set the session's sketcher to Crux while its host
  // built its sketcher, and gave it back after. The host's own sketcher (js-api ISketcherOptions.sketcherType) leaves
  // the session's alone: read every 5 ms while the app opens, it never changes, and a molecule editor opened beside the
  // app shows the session's.
  test('the session\'s sketcher is untouched by opening the app: OpenChemLib before, during and after, the app ' +
    'showing Crux all the same', async () => {
    let app: SketcherApp | null = null;
    let dialog: DG.Dialog | null = null;
    DG.chem.currentSketcherType = 'OpenChemLib';
    const seen = new Set<string>();
    const watch = setInterval(() => seen.add(DG.chem.currentSketcherType), 5);
    try {
      app = await openAt(ADDRESS);
      expect(app.host.sketcher!.root.classList.contains('crux-sketcher'), true, 'the app shows Crux');
      const beside = new DG.chem.Sketcher();
      dialog = DG.Dialog.create('beside the app').add(beside).show();
      const other = await beside.sketcherReady();
      expect(other.root.classList.contains('crux-sketcher'), false,
        'a molecule editor beside the app shows the session\'s sketcher (OpenChemLib), not Crux');
      dialog.close();
      dialog = null;
      close(app);
      app = null;
      expect([...seen].join(', '), 'OpenChemLib', 'the session\'s sketcher, read every 5 ms while the app opened');
      expect(DG.chem.currentSketcherType, 'OpenChemLib', 'the session\'s sketcher after the app closed');
    } finally {
      clearInterval(watch);
      dialog?.close();
      close(app);
      DG.chem.currentSketcherType = sketcherType;
    }
  }, OPENS_ONE);

  test('opened from a link: a SMILES, a CXSMILES with an AND group and a CXSMARTS query each drawn, the address ' +
    'rewritten as the drawing\'s link', async () => {
    // a molecule's link opens in molecule mode, a query's in query mode (main's call, 2026-11-19 18:0x: a shared query
    // reopens as one)
    const cases: [string, (s: SketcherHandle) => boolean, string, string][] = [
      ['OCC', (s) => s.cxsmiles === 'CCO' && !s.hasQuery, 'CCO', 'molecule'],
      [ENHANCED, (s) => s.cxsmiles === ENHANCED, ENHANCED, 'molecule'],
      [QUERY, (s) => s.hasQuery && s.smarts === QUERY_SMARTS && s.cxsmiles === '', QUERY_SMARTS, 'query'],
    ];
    for (const [given, drawn, link, mode] of cases) {
      let app: SketcherApp | null = null;
      try {
        app = await openAt(`${ADDRESS}?${LINK_PARAMETER}=${encodeURIComponent(given)}`);
        expect(drawn(app.sketch!), true, `${given}: drawn (${app.sketch!.cxsmiles || app.sketch!.smarts})`);
        expect(app.sketch!.mode, mode, `${given}: the mode it opens in`);
        await addressBecomes(addressOf(link), given);
      } finally {
        close(app);
      }
    }
  }, OPENS_THREE);

  // The owner's "query options enabled" (2026-11-19 13:05) holds as the gear's Query mode: the app opens in molecule
  // mode (18:0x), and Query mode is a switch away, changing neither the drawing nor its link
  test('Query mode switches the app\'s sketcher to the query tools and back, the drawing and the address as they were',
    async () => {
      let app: SketcherApp | null = null;
      try {
        app = await openAt(`${ADDRESS}?${LINK_PARAMETER}=CCO`);
        await addressBecomes(addressOf('CCO'), 'opened');
        const sketch = app.sketch!;
        expect(sketch.mode, 'molecule', 'opened in molecule mode');
        sketch.mode = 'query';
        expect(sketch.mode, 'query', 'Query mode on');
        sketch.mode = 'molecule';
        expect(sketch.mode, 'molecule', 'Query mode off');
        expect(sketch.cxsmiles, 'CCO', 'the drawing as it was');
        expect(address(), addressOf('CCO'), 'the address as it was');
      } finally {
        close(app);
      }
    }, OPENS_ONE);

  // Main's F1 (run19): a query's link is its CXSMARTS where the SMARTS in it is Crux's own `smarts`, so what SMARTS
  // cannot hold travels with the query. An alias does not: Crux Sketch reads no CXSMARTS atom label, and RDKit's
  // molfile holds none (the lane's handoff, [OWNER]). Each query opened from its link and from the link it is
  // rewritten as.
  test('a query keeps through its link what SMARTS cannot hold, an AND group and an attachment point: the link its ' +
    'CXSMARTS, reopened as drawn, in query mode, the host saying nothing of it', async () => {
    const cases: [string, string, string][] = [
      [QUERY_AND, QUERY_AND_LINK, 'STERAC1'],
      [QUERY_ATTACHED, QUERY_ATTACHED, 'ATTCHPT=1'],
    ];
    const problems: string[] = [];
    for (const [given, link, field] of cases) {
      for (const opened of new Set([given, link])) {
        let app: SketcherApp | null = null;
        try {
          app = await openAt(`${ADDRESS}?${LINK_PARAMETER}=${encodeURIComponent(opened)}`);
          const sketch = app.sketch!;
          if (sketch.mode !== 'query')
            problems.push(`${opened}: opened in ${sketch.mode} mode`);
          if (!sketch.molV3000.includes(field))
            problems.push(`${opened}: drawn without ${field}`);
          const said = app.view.root.querySelector('.chem-invalid-molecule-warning')?.textContent ?? '';
          if (said !== '')
            problems.push(`${opened}: the host says "${said}"`);
          await addressBecomes(addressOf(link), opened).catch(() => problems.push(
            `${opened}: the link is ${new URLSearchParams(location.search).get(LINK_PARAMETER)}, not ${link}`));
        } finally {
          close(app);
        }
      }
    }
    expect(problems.join('\n'), '', 'what each query kept through its link');
  }, {timeout: 180000});

  test('every change of the drawing writes the address through the view\'s path, in place: a molecule as its ' +
    'CXSMILES, a query as its SMARTS, an empty drawing as an empty parameter', async () => {
    let app: SketcherApp | null = null;
    try {
      app = await openAt(`${ADDRESS}?${LINK_PARAMETER}=CCO`);
      await addressBecomes(addressOf('CCO'), 'opened');
      const entries = history.length;
      const steps: [string, string][] = [['CCN', 'CCN'], [ENHANCED, ENHANCED], [QUERY, QUERY_SMARTS], ['', '']];
      for (const [set, link] of steps) {
        app.host.setMolecule(set);
        if (set === '')
          app.sketch!.setValue('');
        await addressBecomes(addressOf(link), `the drawing set to ${set || 'nothing'}`);
        expect(app.view.path.endsWith(linkPath(link)), true, `the view's path holds the link: ${app.view.path}`);
      }
      expect(history.length, entries, 'the address changed in place (history.replaceState): Back leaves the app');
    } finally {
      close(app);
    }
  }, OPENS_ONE);

  test('a malformed link: the app opens, Crux empty, the host says "Malformed molecule", the address as given',
    async () => {
      let app: SketcherApp | null = null;
      const given = `${ADDRESS}?${LINK_PARAMETER}=${encodeURIComponent('C1CC(')}`;
      try {
        app = await openAt(given);
        expect(app.sketch!.isEmpty, true, 'nothing drawn');
        await awaitCheck(() => (app!.view.root.querySelector('.chem-invalid-molecule-warning')?.textContent ?? '')
          .includes('Malformed molecule'), 'the host shows the molecule as malformed', 10000);
        expect(address(), given, 'the address as given until something is drawn');
        app.host.setMolecule('CC');
        await addressBecomes(addressOf('CC'), 'drawn after');
      } finally {
        close(app);
      }
    }, OPENS_ONE);

  test('a drawing too large for the address: the address drops the molecule and says so once; back under, the ' +
    'link again', async () => {
    let app: SketcherApp | null = null;
    const warnings: string[] = [];
    const warning = grok.shell.warning;
    // 500 ammonium ions: '[NH4+].' is 13 characters URL-encoded
    const large = Array(500).fill('[NH4+]').join('.');
    try {
      app = await openAt(`${ADDRESS}?${LINK_PARAMETER}=CC`);
      grok.shell.warning = (s: string | HTMLElement) => {
        warnings.push(String(s));
      };
      app.host.setMolecule(large);
      const length = encodeURIComponent(app.sketch!.cxsmiles).length;
      expect(length > MAX_LINK_LENGTH, true, `the link would be ${length} characters`);
      await addressBecomes(ADDRESS, 'too large');
      expect(JSON.stringify(warnings), JSON.stringify([tooLongForLink(length)]), 'said once');
      app.host.setMolecule(large + '.C');
      await awaitCheck(() => app!.sketch!.cxsmiles.endsWith('.C') || app!.sketch!.cxsmiles.startsWith('C.'),
        'the larger drawing', 10000);
      expect(warnings.length, 1, 'not said again while it stays too large');
      app.host.setMolecule('CCC');
      await addressBecomes(addressOf('CCC'), 'back under');
    } finally {
      grok.shell.warning = warning;
      close(app);
    }
  }, OPENS_ONE);

  test('Export: every format in both styles, from the drawing (the example)', async () => {
    let app: SketcherApp | null = null;
    try {
      app = await openAt(ADDRESS);
      const sketch = app.sketch!;
      const exported = app.exporter()!;
      for (const style of [ACS_1996, AS_DRAWN]) {
        for (const format of FORMATS) {
          const r = await exported.result(format, style, false);
          if ('error' in r)
            throw new Error(`${format} (${style}): ${r.error}`);
          expect(r.file, FORMAT_FILES[format].file);
          const text = typeof r.content === 'string' ? r.content : '';
          const what = `${format} (${style})`;
          switch (format) {
          case SVG:
            expect(/^<\?xml/.test(text) && text.includes('<svg'), true, what);
            expect(new RegExp(`<svg[^>]*\\swidth='[\\d.]+${style === ACS_1996 ? 'pt' : 'px'}'`).test(text), true,
              `${what}: sized in ${style === ACS_1996 ? 'points' : 'CSS px'}`);
            break;
          case PNG: {
            const bytes = new Uint8Array(await (r.content as Blob).arrayBuffer());
            expect(Array.from(bytes.slice(1, 4)).map((b) => String.fromCharCode(b)).join(''), 'PNG', what);
            const at = findChunk(bytes, 'pHYs');
            const perMetre = at < 0 ? 0 : new DataView(bytes.buffer, at + 8, 4).getUint32(0);
            expect(Math.round(perMetre * 0.0254), PNG_DPI, `${what}: its pHYs says ${PNG_DPI} dpi`);
            break;
          }
          case CDXML:
            expect(text.includes('<CDXML'), true, what);
            break;
          case MOL_V2000:
            expect(text.includes('V2000') && text.includes('M  END'), true, what);
            break;
          case MOL_V3000:
            expect(text.includes('V3000') && text.includes('M  END'), true, what);
            break;
          case SDF:
            expect(text.includes('M  END') && text.endsWith('$$$$\n'), true, what);
            break;
          case SMILES:
            expect(text, sketch.plainSmiles, what);
            expect(canonical(text), canonical(EXAMPLE), what);
            break;
          case CXSMILES:
            expect(text, EXAMPLE, what);
            break;
          case SMARTS:
            expect(text, sketch.smarts, what);
            break;
          case INCHI: {
            const inchi = rdkit.get_mol(EXAMPLE);
            try {
              expect(text, inchi.get_inchi(), `${what}: as RDKit writes it`);
            } finally {
              inchi.delete();
            }
            break;
          }
          case INCHIKEY: {
            const mol = rdkit.get_mol(EXAMPLE);
            try {
              expect(text, rdkit.get_inchikey_for_inchi(mol.get_inchi()), `${what}: as RDKit writes it`);
            } finally {
              mol.delete();
            }
            expect(text, EXAMPLE_INCHIKEY, `${what}: atorvastatin's`);
            break;
          }
          }
        }
      }
    } finally {
      close(app);
    }
  }, EXPORTS);

  test('Export styles: ACS 1996 lays the molecule out afresh and draws it black at 14.4 pt bonds; As drawn is the ' +
    'canvas\'s drawing, its coordinates and colours', async () => {
    let app: SketcherApp | null = null;
    try {
      app = await openAt(ADDRESS);
      const sketch = app.sketch!;
      // drawn as no layout would draw it: 3-aminopropanol's five atoms on one vertical line
      app.host.setMolFile(STRAIGHT_UP);
      await awaitCheck(() => sketch.cxsmiles === 'NCCCO', 'the straight drawing', 10000);
      const exported = app.exporter()!;
      const acs = await exported.cdxml(ACS_1996);
      const drawn = await exported.cdxml(AS_DRAWN);
      expect(!!acs && !!drawn && acs !== drawn, true, 'the two CDXML pages differ: one laid out afresh');
      expect(drawn, sketch.cdxml, 'As drawn: the drawing\'s own CDXML');
      // ACS: the CDXML Crux writes for the molecule read from its text, laid out by Crux
      const copy = await exported.sketcherOf(ACS_1996);
      expect(copy !== null && copy !== sketch && canonical(copy!.cxsmiles) === canonical(sketch.cxsmiles), true,
        'the ACS copy is the same molecule');
      const acsSvg = (await exported.picture(ACS_1996, false))!.svg;
      const drawnSvg = (await exported.picture(AS_DRAWN, false))!.svg;
      const colours = (svg: string) => new Set(Array.from(svg.matchAll(/(?:fill|stroke)[:=]['"]?\s*(#[0-9A-Fa-f]{6})/g),
        (m) => m[1].toUpperCase()));
      const inked = (svg: string) => [...colours(svg)].filter((c) => c !== '#FFFFFF');
      expect(inked(acsSvg).every((c) => c === '#000000'), true, `ACS 1996 is black: ${inked(acsSvg)}`);
      expect(inked(drawnSvg).some((c) => c !== '#000000'), true, `As drawn has the atoms' colours: ${inked(drawnSvg)}`);
    } finally {
      close(app);
    }
  }, EXPORTS);

  test('Export of a query: SMARTS, the MOL files and the pictures; SMILES, CXSMILES and InChI say why they have none',
    async () => {
      let app: SketcherApp | null = null;
      try {
        app = await openAt(`${ADDRESS}?${LINK_PARAMETER}=${encodeURIComponent(QUERY)}`);
        const exported = app.exporter()!;
        for (const format of FORMATS) {
          const r = await exported.result(format, ACS_1996, false);
          if (format === SMILES || format === CXSMILES)
            expect('error' in r && r.error, QUERY_HAS_NO_SMILES, format);
          else if (format === INCHI || format === INCHIKEY)
            expect('error' in r && r.error.startsWith(`${format}: InChI cannot hold this drawing`), true, format);
          else
            expect('error' in r ? r.error : '', '', `${format} is written`);
        }
        const smarts = await exported.result(SMARTS, ACS_1996, false);
        expect('content' in smarts && smarts.content, QUERY_SMARTS);
      } finally {
        close(app);
      }
    }, EXPORTS);

  test('the ribbon\'s Export... opens the dialog: style and format chosen as a person does, Copy copies what it ' +
    'shows and leaves the dialog open, Save downloads and closes it', async () => {
    let app: SketcherApp | null = null;
    const download = DG.Utils.download;
    const writeText = navigator.clipboard.writeText;
    const downloads: [string, BlobPart, string | undefined][] = [];
    const copied: string[] = [];
    const dialog = () => DG.Dialog.getOpenDialogs().find((d) => d.title === EXPORT_TITLE);
    const open = async (sketcherApp: SketcherApp) => {
      const button = sketcherApp.view.getRibbonPanels().flat().find((b) => b.textContent === EXPORT_BUTTON)!;
      (button.tagName === 'BUTTON' ? button : button.querySelector('button')!).click();
      await awaitCheck(() => !!dialog(), 'the Export dialog opens', 5000);
      const d = dialog()!;
      const preview = d.root.querySelector(`.${PREVIEW_CLASS}`) as HTMLElement;
      const ready = () => preview.dataset.ready === 'true';
      await awaitCheck(() => ready() && !!preview.querySelector('img'), 'the SVG previewed', 10000);
      return {d, inputs: moleculeExportOf(d.root)!.inputs, preview, ready};
    };
    try {
      app = await openAt(ADDRESS);
      DG.Utils.download = (name: string, content: BlobPart, type?: string) => {
        downloads.push([name, content, type]);
      };
      Object.defineProperty(navigator.clipboard, 'writeText', {configurable: true, value: async (s: string) => {
        copied.push(s);
      }});
      const {d, inputs, preview, ready} = await open(app);
      expect(inputs.style.value, ACS_1996, 'ACS 1996 first, as in the reactions\' panel');
      expect(inputs.format.value, SVG);
      inputs.format.value = MOL_V3000;
      await awaitCheck(() => ready() && !!preview.querySelector(`.${TEXT_CLASS}`), 'the molfile shown as text', 10000);
      expect(inputs.style.enabled, false, 'the style applies to SVG, PNG and CDXML alone');
      const shown = (preview.querySelector(`.${TEXT_CLASS}`) as HTMLTextAreaElement).value;
      expect(shown, app.sketch!.molV3000, 'the text shown is the drawing\'s MOL V3000');
      d.getButton(COPY).click();
      await awaitCheck(() => copied.length === 1, 'Copy copies', 5000);
      expect(copied[0], shown, 'what it copies is what it shows');
      expect(dialog() === d, true, 'Copy leaves the dialog open');
      inputs.format.value = CDXML;
      inputs.style.value = AS_DRAWN;
      await awaitCheck(() => ready() && preview.dataset.style === AS_DRAWN && !!preview.querySelector('img'),
        'CDXML as drawn previewed', 10000);
      expect((d.root.querySelector(`.${NOTE_CLASS}`)?.textContent ?? '').includes('as drawn'), true, 'the note says');
      d.getButton(SAVE).click();
      await awaitCheck(() => downloads.length === 1, 'Save downloads the CDXML', 5000);
      expect(JSON.stringify([downloads[0][0], downloads[0][2]]),
        JSON.stringify([FORMAT_FILES[CDXML].file, FORMAT_FILES[CDXML].type]));
      expect(String(downloads[0][1]), app.sketch!.cdxml, 'the drawing\'s own CDXML');
      await awaitCheck(() => !dialog(), 'Save closes the dialog, as the reactions\' panel does', 5000);
      const again = await open(app);
      expect(again.inputs.style.value, ACS_1996, 'opened again: ACS 1996');
      expect(again.inputs.format.value, SVG, 'opened again: SVG');
      again.d.getButton(SAVE).click();
      await awaitCheck(() => downloads.length === 2, 'Save downloads the SVG', 5000);
      expect(downloads[1][0], FORMAT_FILES[SVG].file);
      expect(STYLED_FORMATS.length, 3);
    } finally {
      DG.Utils.download = download;
      Object.defineProperty(navigator.clipboard, 'writeText', {configurable: true, value: writeText});
      dialog()?.close();
      close(app);
    }
  }, EXPORTS);

  test('help: the view and the Export dialog point at Chem\'s Sketcher page, which has their sections', async () => {
    const response = await fetch(sketcherHelpUrl());
    expect(response.ok, true, `${sketcherHelpUrl()} is served from Chem's files`);
    const page = await response.text();
    expect(page.includes('# Sketcher'), true, 'the page');
    for (const heading of ['## The link', '## Export'])
      expect(page.includes(heading), true, `its section "${heading}"`);
    expect(sketcherHelpUrl('export').endsWith(`${SKETCHER_HELP_FILE}#export`), true, 'the dialog\'s section');
  });
});

/** The offset of the first `type` chunk of PNG `bytes` (its length field), or -1. */
function findChunk(bytes: Uint8Array, type: string): number {
  let at = 8;
  const view = new DataView(bytes.buffer);
  while (at + 8 <= bytes.length) {
    const length = view.getUint32(at);
    const name = String.fromCharCode(...Array.from(bytes.slice(at + 4, at + 8)));
    if (name === type)
      return at;
    at += 12 + length;
  }
  return -1;
}
