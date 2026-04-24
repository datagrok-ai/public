/* eslint-disable max-lines-per-function */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';

import {
  CHEM_ENUM_MAX_RESULTS,
  ChemEnumCore,
  ChemEnumMode,
  ChemEnumModes,
  ChemEnumRGroup,
  enumerateRaw,
  enumerateSampleRaw,
  extractRNumbers,
  makeCore,
  makeRGroup,
  validateParams,
} from './pt-chem-enum';
import {_package} from '../package';
import {defaultErrorHandler} from '../utils/err-info';

const DIALOG_TITLE = 'PolyTool Chem Enumeration';

const MODE_TOOLTIPS: Record<ChemEnumMode, string> = {
  [ChemEnumModes.Zip]: 'Zip: every R-group list must have the same length N. The i-th result uses the i-th entry from every list. Produces N results per core.',
  [ChemEnumModes.Cartesian]: 'Cartesian: every combination of entries across R-group lists. Produces ∏|Rᵢ| results per core.',
};

// Small thumbnail in cards; big render in hover tooltip.
const CARD_W = 110;
const CARD_H = 104;
const THUMB_W = 100;
const THUMB_H = 72;
const TOOLTIP_W = 320;
const TOOLTIP_H = 260;

// ─── Card primitives ────────────────────────────────────────────────────────

interface CardOpts {
  smiles: string;
  subtitle: string;
  error?: string;
  onEdit?: () => void;
  onRemove?: () => void;
}

/** Draws a molecule into a fixed-size host, constraining SVG dimensions. */
function drawMolInto(host: HTMLElement, smi: string, w: number, h: number): void {
  ui.empty(host);
  try {
    const el = grok.chem.drawMolecule(smi, w, h);
    el.style.width = `${w}px`;
    el.style.height = `${h}px`;
    el.style.maxWidth = `${w}px`;
    el.style.maxHeight = `${h}px`;
    el.style.display = 'block';
    host.appendChild(el);
  } catch {
    host.appendChild(ui.divText('—', {style: {color: 'var(--grey-4)'}}));
  }
}

function buildCard(opts: CardOpts): HTMLElement {
  const thumbHost = ui.div([], {style: {
    width: `${THUMB_W}px`, height: `${THUMB_H}px`,
    display: 'flex', alignItems: 'center', justifyContent: 'center',
    background: 'transparent', overflow: 'hidden', flex: '0 0 auto',
  }});
  if (opts.smiles && !opts.error) drawMolInto(thumbHost, opts.smiles, THUMB_W, THUMB_H);
  else thumbHost.appendChild(ui.divText('—', {style: {color: 'var(--grey-4)'}}));

  const subtitleEl = ui.divText(opts.subtitle, {style: {
    fontSize: '10px', color: opts.error ? 'var(--red-3)' : 'var(--grey-5)',
    textAlign: 'center', padding: '2px 4px', whiteSpace: 'nowrap',
    overflow: 'hidden', textOverflow: 'ellipsis',
  }});

  const card = ui.divV([thumbHost, subtitleEl], {style: {
    width: `${CARD_W}px`, minWidth: `${CARD_W}px`, height: `${CARD_H}px`,
    border: `2px solid ${opts.error ? 'var(--red-3)' : 'var(--grey-2)'}`,
    borderRadius: '4px', position: 'relative', background: 'var(--white)',
    boxSizing: 'border-box', margin: '0 4px 0 0', overflow: 'hidden', flex: '0 0 auto',
  }});

  // Large preview on hover — only when the card is valid. Wrap in a function for typing.
  if (opts.smiles && !opts.error) {
    let bigHost: HTMLElement | null = null;
    ui.tooltip.bind(card, () => {
      if (!bigHost) {
        bigHost = ui.div([], {style: {width: `${TOOLTIP_W}px`, height: `${TOOLTIP_H}px`}});
        drawMolInto(bigHost, opts.smiles, TOOLTIP_W, TOOLTIP_H);
      }
      return bigHost;
    });
  } else if (opts.error) {
    ui.tooltip.bind(card, opts.error);
  }

  // Hover-revealed action icons (edit + delete).
  const actions = ui.divH([], {style: {
    position: 'absolute', top: '2px', right: '2px', gap: '2px',
    background: 'rgba(255,255,255,0.85)', borderRadius: '4px',
    padding: '1px 2px', display: 'none',
  }});
  if (opts.onEdit) {
    const editBtn = ui.icons.edit((e: MouseEvent) => { e.stopPropagation(); opts.onEdit!(); }, 'Edit');
    actions.appendChild(editBtn);
  }
  if (opts.onRemove) {
    const delBtn = ui.icons.delete((e: MouseEvent) => { e.stopPropagation(); opts.onRemove!(); }, 'Remove');
    actions.appendChild(delBtn);
  }
  card.addEventListener('mouseenter', () => { actions.style.display = 'flex'; });
  card.addEventListener('mouseleave', () => { actions.style.display = 'none'; });
  if (opts.onEdit || opts.onRemove) card.appendChild(actions);

  return card;
}

// ─── Draw core / R-group (single dialog with sketcher + R# below) ───────────

async function smilesFromSketcher(sk: DG.chem.Sketcher, rdkit: RDModule): Promise<string | null> {
  const mol = sk.getMolFile();
  if (!mol || mol.trim() === '') return null;
  if (!DG.chem.isMolBlock(mol)) return mol.trim();
  try {
    const smi = await grok.functions.call('Chem:convertMolNotation', {
      molecule: mol,
      sourceNotation: DG.chem.Notation.MolBlock,
      targetNotation: DG.chem.Notation.Smiles,
    });
    return (smi ?? '').toString().trim() || null;
  } catch {
    return null;
  }
}

/**
 * Draw-a-core dialog. OK is enabled iff the sketched molecule has ≥1 R-label.
 * `initialMolfile` preloads the sketcher for editing an existing core.
 */
async function openCoreSketchDialog(rdkit: RDModule, initialMolfile?: string): Promise<ChemEnumCore | null> {
  return new Promise((resolve) => {
    const sk = new DG.chem.Sketcher(DG.chem.SKETCHER_MODE.INPLACE);
    sk.syncCurrentObject = false;
    if (initialMolfile) sk.setMolFile(initialMolfile);

    const status = ui.divText('', {style: {fontSize: '11px', padding: '4px 0', minHeight: '18px'}});
    let currentSmiles: string | null = null;
    let okBtn: HTMLButtonElement | null = null;

    const revalidate = async () => {
      currentSmiles = await smilesFromSketcher(sk, rdkit);
      const rs = currentSmiles ? extractRNumbers(currentSmiles) : [];
      if (!currentSmiles) {
        status.innerText = 'Draw a molecule with at least one R-group label.';
        status.style.color = 'var(--grey-4)';
      } else if (rs.length === 0) {
        status.innerText = 'No R-group detected — add e.g. [*:1] to mark an attachment point.';
        status.style.color = 'var(--red-3)';
      } else {
        status.innerText = `Detected R-groups: ${rs.map((n) => 'R' + n).join(', ')}.`;
        status.style.color = 'var(--green-2)';
      }
      if (okBtn) okBtn.disabled = !(currentSmiles && rs.length > 0);
    };

    sk.onChanged.subscribe(() => { revalidate(); });

    const body = ui.divV([sk.root, status]);
    const dialog = ui.dialog({title: initialMolfile ? 'Edit Core' : 'Draw Core'})
      .add(body)
      .onOK(() => resolve(currentSmiles ? makeCore(currentSmiles, '', rdkit) : null))
      .onCancel(() => resolve(null));
    dialog.show({resizable: true});
    okBtn = dialog.getButton('OK') as HTMLButtonElement;
    okBtn.disabled = true;
    revalidate();
  });
}

/**
 * Draw-an-R-group dialog. R# is auto-inferred from the sketched SMILES.
 * OK is enabled iff the molecule has exactly one R-label AND R# is set.
 */
async function openRGroupSketchDialog(rdkit: RDModule, initialMolfile?: string, initialRNum?: number): Promise<ChemEnumRGroup | null> {
  return new Promise((resolve) => {
    const sk = new DG.chem.Sketcher(DG.chem.SKETCHER_MODE.INPLACE);
    sk.syncCurrentObject = false;
    if (initialMolfile) sk.setMolFile(initialMolfile);

    const rNumberInput = ui.input.int('Target R#', {value: initialRNum, min: 1});
    const status = ui.divText('', {style: {fontSize: '11px', padding: '4px 0', minHeight: '18px'}});
    let currentSmiles: string | null = null;
    let detectedRs: number[] = [];
    let okBtn: HTMLButtonElement | null = null;
    let userTouchedR = false;
    rNumberInput.onChanged.subscribe(() => { userTouchedR = true; updateOk(); });

    const updateOk = () => {
      const n = rNumberInput.value;
      const ok = currentSmiles != null && n != null && n >= 1 && detectedRs.length === 1;
      if (okBtn) okBtn.disabled = !ok;
    };

    const revalidate = async () => {
      currentSmiles = await smilesFromSketcher(sk, rdkit);
      detectedRs = currentSmiles ? extractRNumbers(currentSmiles) : [];
      if (!currentSmiles) {
        status.innerText = 'Draw an R-group with exactly one attachment point.';
        status.style.color = 'var(--grey-4)';
      } else if (detectedRs.length === 0) {
        status.innerText = 'No R-group detected — add [*:N] to mark the attachment point.';
        status.style.color = 'var(--red-3)';
      } else if (detectedRs.length > 1) {
        status.innerText = `R-group must contain exactly one attachment point. Found ${detectedRs.length}: ${detectedRs.map((n) => 'R' + n).join(', ')}.`;
        status.style.color = 'var(--red-3)';
      } else {
        if (!userTouchedR || rNumberInput.value == null) rNumberInput.value = detectedRs[0];
        status.innerText = `Detected R${detectedRs[0]}. Target R# is used for joining; change it to remap.`;
        status.style.color = 'var(--green-2)';
      }
      updateOk();
    };

    sk.onChanged.subscribe(() => { revalidate(); });

    const body = ui.divV([sk.root, rNumberInput.root, status]);
    const dialog = ui.dialog({title: initialMolfile ? 'Edit R-Group' : 'Draw R-Group'})
      .add(body)
      .onOK(() => {
        const n = rNumberInput.value;
        if (currentSmiles && n != null) resolve(makeRGroup(currentSmiles, n, '', rdkit));
        else resolve(null);
      })
      .onCancel(() => resolve(null));
    dialog.show({resizable: true});
    okBtn = dialog.getButton('OK') as HTMLButtonElement;
    okBtn.disabled = true;
    revalidate();
  });
}

// ─── Import wizard ──────────────────────────────────────────────────────────

function inferRNumberFromColumnName(name: string): number | null {
  const m = name.match(/r[\s_\-:]*(\d+)/i);
  return m ? parseInt(m[1], 10) : null;
}

interface ImportResult {
  cores?: ChemEnumCore[];
  rGroups?: ChemEnumRGroup[];
}

async function openImportWizard(
  kind: 'cores' | 'rgroups', rdkit: RDModule, defaultDedup: boolean,
): Promise<ImportResult | null> {
  return new Promise((resolve) => {
    const isPickable = (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE;
    const tableInput = ui.input.table('Table', {value: grok.shell.t ?? undefined});
    const columnInput = ui.input.choice<string>('Column', {items: [] as string[], nullable: true});
    const rNumberInput = kind === 'rgroups' ?
      ui.input.int('Target R#', {value: 1, min: 1}) : null;
    const dedupInput = ui.input.bool('Remove duplicates', {value: defaultDedup});
    ui.tooltip.bind(dedupInput.input,
      kind === 'cores' ?
        'When checked, identical cores in the selected column are collapsed into one entry.' :
        'When checked, identical R-groups are collapsed. Leave off for Zip mode if your lists are intentionally aligned by position.');

    const currentColumn = (): DG.Column | null => {
      const t = tableInput.value;
      const n = columnInput.value;
      return t && n ? t.col(n) : null;
    };

    const repopulateColumns = (autoPick: boolean) => {
      const t = tableInput.value;
      if (!t) { columnInput.items = []; columnInput.value = null; return; }
      const cols = t.columns.toList().filter(isPickable);
      columnInput.items = cols.map((c) => c.name);
      if (autoPick) {
        const mol = cols.find((c) => c.semType === DG.SEMTYPE.MOLECULE);
        columnInput.value = (mol ?? cols[0])?.name ?? null;
      }
    };

    const previewHost = ui.div([], {style: {height: `${CARD_H + 26}px`, overflow: 'hidden', padding: '4px'}});
    const countText = ui.divText('', {style: {fontSize: '11px', color: 'var(--grey-5)', margin: '4px 0'}});
    const badge = makeErrorBadge();

    let built: ImportResult = {};
    let okButton: HTMLButtonElement | null = null;
    let vv: DG.VirtualView | null = null;

    const rebuildPreview = () => {
      badge.setErrors([]);
      built = {};
      const col = currentColumn();
      if (!col) {
        countText.innerText = '';
        if (okButton) okButton.disabled = true;
        if (vv) vv.setData(0, () => ui.div());
        return;
      }

      const values: string[] = [];
      for (let i = 0; i < col.length; i++) {
        if (col.isNone(i)) continue;
        const v = col.get(i);
        if (v == null || String(v).trim() === '') continue;
        values.push(String(v));
      }

      const dedup = dedupInput.value;
      const items: { smi: string, err?: string, subtitle: string }[] = [];
      const errs: string[] = [];
      if (kind === 'cores') {
        const parsed = values.map((v) => makeCore(v, '', rdkit));
        let cores: ChemEnumCore[] = parsed;
        let dupCount = 0;
        if (dedup) {
          const seen = new Set<string>();
          cores = [];
          for (const c of parsed) {
            const k = c.error ? `err:${c.originalSmiles}` : c.smiles;
            if (seen.has(k)) { dupCount++; continue; }
            seen.add(k);
            cores.push(c);
          }
        }
        built.cores = cores;
        cores.forEach((c, i) => {
          if (c.error) errs.push(`Core ${i + 1}: ${c.error}`);
          items.push({
            smi: c.error ? '' : c.smiles, err: c.error,
            subtitle: c.error ? 'invalid' : c.rNumbers.map((n) => 'R' + n).join(', '),
          });
        });
        countText.innerText = `${cores.length} core${cores.length === 1 ? '' : 's'}` +
          (dedup && dupCount ? ` (${dupCount} duplicate${dupCount === 1 ? '' : 's'} skipped)` : '') + '.';
      } else {
        const n = rNumberInput!.value ?? 1;
        const parsed = values.map((v) => makeRGroup(v, n, '', rdkit));
        let rGroups: ChemEnumRGroup[] = parsed;
        let dupCount = 0;
        if (dedup) {
          const seen = new Set<string>();
          rGroups = [];
          for (const rg of parsed) {
            const k = rg.error ? `err:${rg.originalSmiles}` : rg.smiles;
            if (seen.has(k)) { dupCount++; continue; }
            seen.add(k);
            rGroups.push(rg);
          }
        }
        built.rGroups = rGroups;
        rGroups.forEach((rg, i) => {
          if (rg.error) errs.push(`R-group ${i + 1}: ${rg.error}`);
          items.push({
            smi: rg.error ? '' : rg.smiles, err: rg.error,
            subtitle: rg.error ? 'invalid' : `R${rg.rNumber}${rg.sourceRNumber != null && rg.sourceRNumber !== rg.rNumber ? ` (from R${rg.sourceRNumber})` : ''}`,
          });
        });
        countText.innerText = `${rGroups.length} R-group${rGroups.length === 1 ? '' : 's'} for R${n}` +
          (dedup && dupCount ? ` (${dupCount} duplicate${dupCount === 1 ? '' : 's'} skipped)` : '') + '.';
      }

      badge.setErrors(errs);
      if (!vv) {
        vv = ui.virtualView(items.length, (i) => buildCard({smiles: items[i].smi, subtitle: items[i].subtitle, error: items[i].err}), false, 1);
        vv.root.style.height = `${CARD_H + 12}px`;
        vv.root.style.width = '100%';
        previewHost.appendChild(vv.root);
      } else {
        vv.setData(items.length, (i) => buildCard({smiles: items[i].smi, subtitle: items[i].subtitle, error: items[i].err}));
      }
      if (okButton) okButton.disabled = values.length === 0;
    };

    tableInput.onChanged.subscribe(async () => {
      const t = tableInput.value;
      if (t) { await t.meta.detectSemanticTypes(); await grok.data.detectSemanticTypes(t); }
      repopulateColumns(true);
      rebuildPreview();
    });
    columnInput.onChanged.subscribe(() => {
      if (kind === 'rgroups' && columnInput.value) {
        const inferred = inferRNumberFromColumnName(columnInput.value);
        if (inferred != null) rNumberInput!.value = inferred;
      }
      rebuildPreview();
    });
    if (rNumberInput) rNumberInput.onChanged.subscribe(() => rebuildPreview());
    dedupInput.onChanged.subscribe(() => rebuildPreview());

    const body = ui.divV([
      tableInput.root,
      columnInput.root,
      ...(rNumberInput ? [rNumberInput.root] : []),
      dedupInput.root,
      ui.divH([countText, badge.root], {style: {alignItems: 'center', gap: '8px'}}),
      previewHost,
    ]);

    const dialog = ui.dialog({title: kind === 'cores' ? 'Import Cores' : 'Import R-Groups'})
      .add(body)
      .onOK(() => resolve(built))
      .onCancel(() => resolve(null));
    dialog.show({resizable: true, width: 640});
    okButton = dialog.getButton('OK') as HTMLButtonElement;
    okButton.disabled = true;

    if (tableInput.value) {
      tableInput.value.meta.detectSemanticTypes().then(() => {
        repopulateColumns(true);
        rebuildPreview();
      });
    }
  });
}

// ─── Error badge ────────────────────────────────────────────────────────────

function makeErrorBadge(): { root: HTMLElement, setErrors: (errs: string[]) => void } {
  const root = ui.div([], {style: {
    display: 'none', alignItems: 'center', gap: '4px',
    padding: '2px 8px', borderRadius: '10px',
    background: 'var(--red-2)', color: 'var(--red-3)',
    fontSize: '11px', fontWeight: '600', cursor: 'help',
  }});
  return {
    root,
    setErrors(errs: string[]) {
      if (errs.length === 0) { root.style.display = 'none'; return; }
      root.style.display = 'inline-flex';
      root.innerText = `⚠ ${errs.length} issue${errs.length === 1 ? '' : 's'}`;
      ui.tooltip.bind(root, errs.join('\n'));
    },
  };
}

// ─── Main dialog ────────────────────────────────────────────────────────────

interface ChemEnumDialogState {
  cores: ChemEnumCore[];
  rGroupsByNum: Map<number, ChemEnumRGroup[]>;
  mode: ChemEnumMode;
  appendToTable: DG.DataFrame | null;
}

/** Extracts a molfile string from a cell whose column is of MOLECULE semtype. */
async function cellToMolfile(cell: DG.Cell): Promise<string | null> {
  try {
    if (!cell || cell.rowIndex < 0 || cell.column?.semType !== DG.SEMTYPE.MOLECULE) return null;
    const raw = cell.value;
    if (!raw) return null;
    if (DG.chem.isMolBlock(raw)) return raw;
    return (await grok.functions.call('Chem:convertMolNotation', {
      molecule: raw,
      sourceNotation: cell.column.getTag(DG.TAGS.UNITS) ?? DG.chem.Notation.Unknown,
      targetNotation: DG.chem.Notation.MolBlock,
    })) ?? null;
  } catch {
    return null;
  }
}

export async function polyToolEnumerateChemUI(cell?: DG.Cell): Promise<void> {
  await _package.initPromise;

  try {
    const rdkit = await getRdKitModule();

    // If launched from a molecule cell context menu, preload that molecule as the first core.
    let preloadCore: ChemEnumCore | null = null;
    if (cell) {
      const mol = await cellToMolfile(cell);
      if (mol) {
        const smi = (await grok.functions.call('Chem:convertMolNotation', {
          molecule: mol,
          sourceNotation: DG.chem.Notation.MolBlock,
          targetNotation: DG.chem.Notation.Smiles,
        })) as string;
        if (smi) {
          const c = makeCore(smi, '', rdkit);
          if (!c.error) preloadCore = c;
          else grok.shell.info('Selected molecule has no R-group — add at least one R-label to use it as a core.');
        }
      }
    }

    const panel = buildChemEnumPanel(rdkit, preloadCore);
    const dialog = ui.dialog({title: DIALOG_TITLE})
      .add(panel.root)
      .onOK(async () => { await panel.execute(); });
    panel.bindActionButton(dialog.getButton('OK') as HTMLButtonElement);
    dialog.show({resizable: true, width: 960});
  } catch (err: any) {
    defaultErrorHandler(err);
  }
}

/** Opens the chem enumeration UI as a top-level view (app entry point). */
export async function polyToolEnumerateChemApp(): Promise<DG.View | null> {
  await _package.initPromise;
  try {
    const rdkit = await getRdKitModule();
    const panel = buildChemEnumPanel(rdkit, null);
    const runBtn = ui.bigButton('Enumerate', async () => { await panel.execute(); });
    panel.bindActionButton(runBtn as HTMLButtonElement);

    const view = DG.View.create();
    view.name = DIALOG_TITLE;
    view.box = true;
    view.root.appendChild(ui.divV([panel.root, ui.div([runBtn], {style: {padding: '8px 4px 0'}})],
      {style: {height: '100%', width: '100%', padding: '8px'}}));
    return view;
  } catch (err: any) {
    defaultErrorHandler(err);
    return null;
  }
}

interface ChemEnumPanel {
  root: HTMLElement;
  state: ChemEnumDialogState;
  execute: () => Promise<void>;
  /** Bind the action (OK / Enumerate) button so it tracks validation state. */
  bindActionButton: (btn: HTMLButtonElement) => void;
}

function buildChemEnumPanel(rdkit: RDModule, preloadCore: ChemEnumCore | null): ChemEnumPanel {
  const state: ChemEnumDialogState = {
    cores: preloadCore ? [preloadCore] : [],
    rGroupsByNum: new Map(),
    mode: ChemEnumModes.Cartesian,
    appendToTable: null,
  };

  // ── Cores: single-row horizontal virtualView ───────────────────────────────
  const ROW_H = CARD_H + 16; // card height + horizontal scrollbar breathing room
  const coresEmpty = ui.divText('No cores — draw or import at least one.', {style: {color: 'var(--grey-4)', padding: '20px 12px', fontSize: '12px'}});
  const coresVvHost = ui.div([], {style: {width: '100%', height: `${ROW_H}px`, overflow: 'hidden'}});
  let coresVv: DG.VirtualView | null = null;

  const coresRenderer = (i: number): HTMLElement => {
    const c = state.cores[i];
    return buildCard({
      smiles: c.error ? '' : c.smiles,
      subtitle: c.error ? 'invalid' : `core ${i + 1} · ${c.rNumbers.map((n) => 'R' + n).join(', ')}`,
      error: c.error,
      onEdit: async () => {
        const edited = await openCoreSketchDialog(rdkit, c.smiles);
        if (edited) { state.cores[i] = edited; refresh(); }
      },
      onRemove: () => { state.cores.splice(i, 1); refresh(); },
    });
  };

  const redrawCores = () => {
    if (state.cores.length === 0) {
      if (coresVvHost.firstChild) ui.empty(coresVvHost);
      coresVv = null;
      if (coresEmpty.parentElement !== coresVvHost.parentElement)
        coresVvHost.parentElement?.insertBefore(coresEmpty, coresVvHost.nextSibling);

      coresVvHost.style.display = 'none';
      coresEmpty.style.display = 'block';
      return;
    }
    coresVvHost.style.display = 'block';
    coresEmpty.style.display = 'none';
    if (!coresVv) {
      coresVv = ui.virtualView(state.cores.length, coresRenderer, false, 1);
      applyHorizontalRowStyle(coresVv.root);
      coresVvHost.appendChild(coresVv.root);
    } else {
      coresVv.setData(state.cores.length, coresRenderer);
    }
  };

  /** VirtualView defaults to a vertical-scroll viewport; clamp it to horizontal-only for our rows. */
  function applyHorizontalRowStyle(vvRoot: HTMLElement) {
    vvRoot.style.width = '100%';
    vvRoot.style.height = `${ROW_H}px`;
    vvRoot.style.setProperty('overflow-y', 'hidden', 'important');
    vvRoot.style.setProperty('overflow-x', 'auto', 'important');
    const viewport = vvRoot.firstElementChild as HTMLElement | null;
    if (viewport) {
      viewport.style.setProperty('overflow-y', 'hidden', 'important');
      viewport.style.setProperty('overflow-x', 'auto', 'important');
      viewport.style.height = `${ROW_H}px`;
    }
  }

  // ── R-Groups: capped-height column, one horizontal virtualView per R# ─────
  const rGroupsHost = ui.div([], {style: {
    width: '100%', padding: '2px 0',
    flex: '1 1 auto', minHeight: '0',
    overflowY: 'auto', overflowX: 'hidden',
  }});
  const rGroupsEmpty = ui.divText('No R-groups — draw or import for each R number used by your cores.', {style: {color: 'var(--grey-4)', padding: '12px', fontSize: '12px'}});

  const rGroupsRenderers = new Map<number, {row: HTMLElement, vv: DG.VirtualView, header: HTMLElement}>();

  const redrawRGroups = () => {
    ui.empty(rGroupsHost);
    rGroupsRenderers.clear();

    if (state.rGroupsByNum.size === 0) { rGroupsHost.appendChild(rGroupsEmpty); return; }

    const sortedNums = [...state.rGroupsByNum.keys()].sort((a, b) => a - b);
    for (const n of sortedNums) {
      const list = state.rGroupsByNum.get(n)!;
      const label = ui.divText(`R${n} (${list.length})`, {style: {fontWeight: '600', fontSize: '12px', color: 'var(--grey-6)', alignSelf: 'center'}});
      const clearBtn = ui.button('Remove all', () => { state.rGroupsByNum.delete(n); refresh(); });
      clearBtn.style.marginLeft = '4px';
      ui.tooltip.bind(clearBtn, `Remove all R${n} groups`);
      const header = ui.divH([label, clearBtn], {style: {alignItems: 'center', justifyContent: 'flex-start', gap: '0', margin: '6px 0 2px'}});
      const renderer = (i: number): HTMLElement => {
        const rg = list[i];
        const remap = rg.sourceRNumber != null && rg.sourceRNumber !== rg.rNumber ? ` (from R${rg.sourceRNumber})` : '';
        return buildCard({
          smiles: rg.error ? '' : rg.smiles,
          subtitle: rg.error ? 'invalid' : `r group ${i + 1} · R${rg.rNumber}${remap}`,
          error: rg.error,
          onEdit: async () => {
            const edited = await openRGroupSketchDialog(rdkit, rg.smiles, rg.rNumber);
            if (!edited) return;
            // If the R# changed, move between lists.
            if (edited.rNumber !== rg.rNumber) {
              list.splice(i, 1);
              if (list.length === 0) state.rGroupsByNum.delete(rg.rNumber);
              const target = state.rGroupsByNum.get(edited.rNumber) ?? [];
              target.push(edited);
              state.rGroupsByNum.set(edited.rNumber, target);
            } else {
              list[i] = edited;
            }
            refresh();
          },
          onRemove: () => {
            list.splice(i, 1);
            if (list.length === 0) state.rGroupsByNum.delete(n);
            refresh();
          },
        });
      };
      const vv = ui.virtualView(list.length, renderer, false, 1);
      applyHorizontalRowStyle(vv.root);
      const row = ui.divV([header, vv.root]);
      rGroupsHost.appendChild(row);
      rGroupsRenderers.set(n, {row, vv, header});
    }
  };

  // ── Preview (up to 12 random samples, wrapped flex) ───────────────────────
  const PREVIEW_COUNT = 12;
  const previewHost = ui.div([], {style: {
    width: '100%', display: 'flex', flexWrap: 'wrap', gap: '6px',
    alignContent: 'flex-start', padding: '2px 0', maxHeight: '450px', overflow: 'scroll'
  }});

  const redrawPreview = () => {
    ui.empty(previewHost);
    // Uncanonicalized but parseable — no sync RDKit work during enumeration;
    // drawMolecule parses each visible SMILES lazily when the card renders.
    const samples = enumerateSampleRaw(
      {cores: state.cores, rGroups: state.rGroupsByNum, mode: state.mode, maxResults: Math.min(CHEM_ENUM_MAX_RESULTS, 1000)},
      PREVIEW_COUNT);
    if (samples.length === 0) {
      previewHost.appendChild(ui.divText('Preview will appear once cores and R-groups are valid.', {style: {color: 'var(--grey-4)', padding: '20px 12px', fontSize: '12px'}}));
      return;
    }
    samples.forEach((s: {smiles: string; rGroupSmilesByNum: Map<number, string>}, i: number) => {
      const rTag = [...s.rGroupSmilesByNum.keys()].sort((a, b) => a - b).map((n) => 'R' + n).join(',');
      previewHost.appendChild(buildCard({smiles: s.smiles, subtitle: `sample ${i + 1} · ${rTag}`}));
    });
  };

  // ── Status line: count text + error badge ──────────────────────────────────
  const countText = ui.divText('', {style: {fontSize: '12px', color: 'var(--grey-6)'}});
  const errorBadge = makeErrorBadge();

  const modeInput = ui.input.choice<ChemEnumMode>('Enumerator type', {
    value: state.mode,
    items: Object.values(ChemEnumModes),
    onValueChanged: (v) => { state.mode = v!; updateModeTooltip(); refresh(); },
  }) as DG.ChoiceInput<ChemEnumMode>;
  const updateModeTooltip = () => ui.tooltip.bind(modeInput.input, MODE_TOOLTIPS[state.mode] ?? '');
  updateModeTooltip();

  const appendToTableInput = ui.input.table('Append to table', {
    items: grok.shell.tables, nullable: true,
    onValueChanged: (v) => { state.appendToTable = v; },
  });

  let okButton: HTMLButtonElement | null = null;

  /** Collects positional error messages (no reference to internal ids). */
  const collectPositionalErrors = (): string[] => {
    const errs: string[] = [];
    state.cores.forEach((c, i) => { if (c.error) errs.push(`Core ${i + 1}: ${c.error}`); });
    const sortedNums = [...state.rGroupsByNum.keys()].sort((a, b) => a - b);
    for (const n of sortedNums) {
      const list = state.rGroupsByNum.get(n)!;
      list.forEach((rg, i) => { if (rg.error) errs.push(`R${n} group ${i + 1}: ${rg.error}`); });
    }
    const v = validateParams({cores: state.cores, rGroups: state.rGroupsByNum, mode: state.mode});
    // Strip duplicates — core/rg errors are already listed above; add only global ones.
    for (const m of v.errors) {
      if (/^Core "/.test(m)) continue; // id-based core messages — already represented positionally
      errs.push(m);
    }
    return errs;
  };

  // eslint-disable-next-line prefer-const
  let coresClearBtn: HTMLButtonElement | undefined;
  // eslint-disable-next-line prefer-const
  let rGroupsClearBtn: HTMLButtonElement | undefined;

  const refresh = () => {
    redrawCores();
    redrawRGroups();
    const v = validateParams({cores: state.cores, rGroups: state.rGroupsByNum, mode: state.mode});
    countText.innerText = v.predictedCount > 0 ?
      `${v.predictedCount.toLocaleString()} molecule${v.predictedCount === 1 ? '' : 's'} will be generated` : '';
    errorBadge.setErrors(collectPositionalErrors());
    if (okButton) okButton.disabled = !v.ok;
    if (coresClearBtn) coresClearBtn.disabled = state.cores.length === 0;
    if (rGroupsClearBtn) rGroupsClearBtn.disabled = state.rGroupsByNum.size === 0;
    redrawPreview();
  };

  // ── Add/import handlers (dedup now lives in the import wizard, gated by its checkbox) ──
  const addCoreFromSketcher = async () => {
    const c = await openCoreSketchDialog(rdkit);
    if (!c) return;
    state.cores.push(c);
    refresh();
  };
  const addCoresFromImport = async () => {
    const res = await openImportWizard('cores', rdkit, state.mode === ChemEnumModes.Cartesian);
    if (!res?.cores) return;
    state.cores.push(...res.cores);
    refresh();
  };
  const addRGroupFromSketcher = async () => {
    const rg = await openRGroupSketchDialog(rdkit);
    if (!rg) return;
    const list = state.rGroupsByNum.get(rg.rNumber) ?? [];
    list.push(rg);
    state.rGroupsByNum.set(rg.rNumber, list);
    refresh();
  };
  const addRGroupsFromImport = async () => {
    const res = await openImportWizard('rgroups', rdkit, state.mode === ChemEnumModes.Cartesian);
    if (!res?.rGroups) return;
    for (const rg of res.rGroups) {
      const list = state.rGroupsByNum.get(rg.rNumber) ?? [];
      list.push(rg);
      state.rGroupsByNum.set(rg.rNumber, list);
    }
    refresh();
  };

  const sectionHeader = (
    label: string,
    onDraw?: () => void, onImport?: () => void, onClear?: () => void,
  ): {root: HTMLElement, clearBtn?: HTMLButtonElement} => {
    const parts: HTMLElement[] = [
      ui.divText(label, {style: {fontWeight: '600', fontSize: '12px', color: 'var(--grey-6)', alignSelf: 'center', minWidth: '60px'}}),
    ];
    if (onDraw) {
      const drawBtn = ui.button('+ Draw', onDraw);
      drawBtn.style.marginLeft = '4px';
      ui.tooltip.bind(drawBtn, `${label}: open sketcher`);
      parts.push(drawBtn);
    }
    if (onImport) {
      const importBtn = ui.button('↓ Import…', onImport);
      importBtn.style.marginLeft = '4px';
      ui.tooltip.bind(importBtn, `${label}: pick a table + column`);
      parts.push(importBtn);
    }
    let clearBtn: HTMLButtonElement | undefined;
    if (onClear) {
      clearBtn = ui.button('Remove all', onClear) as HTMLButtonElement;
      clearBtn.style.marginLeft = '4px';
      ui.tooltip.bind(clearBtn, `Remove all ${label.toLowerCase()}`);
      parts.push(clearBtn);
    }
    const root = ui.divH(parts, {style: {
      alignItems: 'center', justifyContent: 'flex-start',
      gap: '0', margin: '0 0 2px', flex: '0 0 auto',
    }});
    return {root, clearBtn};
  };

  const sectionStyle = {
    display: 'flex', flexDirection: 'column',
    minHeight: '150px', padding: '4px 0',
  } as const;

  // ── Layout: horizontal split — cores + r-groups (left 60%) | preview (right 40%) ──
  const coresHeader = sectionHeader(
    'Cores', addCoreFromSketcher, addCoresFromImport,
    () => { state.cores.splice(0, state.cores.length); refresh(); },
  );
  coresClearBtn = coresHeader.clearBtn;
  const coresSection = ui.divV([
    coresHeader.root,
    ui.div([coresVvHost, coresEmpty], {style: {padding: '0'}}),
  ], {style: sectionStyle});

  const rGroupsHeader = sectionHeader(
    'R-Groups', addRGroupFromSketcher, addRGroupsFromImport,
    () => { state.rGroupsByNum.clear(); refresh(); },
  );
  rGroupsClearBtn = rGroupsHeader.clearBtn;
  const rGroupsSection = ui.divV([
    rGroupsHeader.root,
    rGroupsHost,
  ], {style: {
    display: 'flex', flexDirection: 'column',
    maxHeight: '300px', padding: '4px 0',
  }});

  const previewSection = ui.divV([
    sectionHeader(`Preview`).root,
    previewHost,
  ], {style: {
    ...sectionStyle,
    width: '100%', height: '100%',
    overflowY: 'auto', overflowX: 'hidden',
  }});

  const leftColumn = ui.divV([coresSection, rGroupsSection], {style: {
    flex: '0 0 60%', width: '60%',
    display: 'flex', flexDirection: 'column', gap: '4px',
    paddingRight: '8px', borderRight: '1px solid var(--grey-2)',
  }});

  const rightColumn = ui.divV([previewSection], {style: {
    flex: '0 0 40%', width: '40%',
    display: 'flex', flexDirection: 'column',
    paddingLeft: '8px',
  }});

  const split = ui.divH([leftColumn, rightColumn], {style: {
    width: '100%', alignItems: 'stretch', gap: '0',
  }});

  const statusLine = ui.divH([
    modeInput.root,
    countText,
    errorBadge.root,
  ], {style: {alignItems: 'center', gap: '16px', padding: '4px 0'}});

  const footer = ui.div([appendToTableInput.root], {style: {padding: '2px 0'}});

  const body = ui.divV([
    split, statusLine, footer,
  ], {style: {
    width: '100%', padding: '4px',
    display: 'flex', flexDirection: 'column', gap: '4px',
  }});

  refresh();
  return {
    root: body,
    state,
    execute: () => executeEnumeration(state, rdkit),
    bindActionButton: (btn) => { okButton = btn; refresh(); },
  };
}

// ─── Execution ──────────────────────────────────────────────────────────────

async function executeEnumeration(state: ChemEnumDialogState, _rdkit: RDModule): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('Enumerating...');
  try {
    // Stage 1 — string-level join, no RDKit calls. Produces parseable but uncanonical SMILES.
    pi.update(10, 'Building molecules...');
    const results = enumerateRaw({cores: state.cores, rGroups: state.rGroupsByNum, mode: state.mode});
    if (!results) { grok.shell.warning('Enumeration failed — check validation messages in the dialog.'); return; }
    if (results.length === 0) { grok.shell.warning('No molecules produced.'); return; }

    const rNumbersUsed = new Set<number>();
    for (const r of results) for (const n of r.rGroupSmilesByNum.keys()) rNumbersUsed.add(n);
    const sortedRs = [...rNumbersUsed].sort((a, b) => a - b);

    const smilesCol = DG.Column.fromStrings('Enumerated', results.map((r) => r.smiles));
    smilesCol.semType = DG.SEMTYPE.MOLECULE;
    const coreCol = DG.Column.fromStrings('Core', results.map((r) => r.coreSmiles));
    coreCol.semType = DG.SEMTYPE.MOLECULE;
    const rCols = sortedRs.map((n) =>
      DG.Column.fromStrings(`R${n}`, results.map((r) => r.rGroupSmilesByNum.get(n) ?? '')));
    for (const c of rCols) c.semType = DG.SEMTYPE.MOLECULE;

    const df = DG.DataFrame.fromColumns([smilesCol, coreCol, ...rCols]);
    df.name = 'Chem Enumeration';

    // Stage 2 — canonicalize the whole Enumerated column in parallel via Chem workers.
    pi.update(40, `Canonicalizing ${results.length.toLocaleString()} molecule(s)...`);
    try {
      await grok.functions.call('Chem:convertNotation', {
        data: df,
        molecules: smilesCol,
        targetNotation: DG.chem.Notation.Smiles,
        overwrite: true,
        join: false,
        kekulize: false,
      });
    } catch (err: any) {
      // Canonicalization is a nice-to-have; the uncanonical SMILES are still valid output.
      _package.logger.warning(`Canonicalization skipped: ${err?.message ?? err}`);
    }

    pi.update(90, 'Finalizing...');
    await grok.data.detectSemanticTypes(df);

    if (state.appendToTable) {
      state.appendToTable.append(df, true);
      await state.appendToTable.meta.detectSemanticTypes();
    } else {
      grok.shell.addTableView(df);
    }
  } catch (err: any) {
    defaultErrorHandler(err);
  } finally {
    pi.close();
  }
}
