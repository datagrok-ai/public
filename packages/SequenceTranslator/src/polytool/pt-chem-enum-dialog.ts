/* eslint-disable max-lines */
/* eslint-disable max-lines-per-function */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';

import {
  addRGroupsFromSmiles,
  buildExportColumns,
  CHEM_ENUM_MAX_RESULTS,
  ChemEnumCore,
  ChemEnumMode,
  ChemEnumModes,
  ChemEnumRGroup,
  copyRGroupList,
  enumerateRaw,
  enumerateSampleRaw,
  extractRNumbers,
  BUILTIN_R_GROUP_TEMPLATES,
  isValidTemplateSmiles,
  makeCore,
  makeRGroup,
  normalizeRLabels,
  parseRGroupTemplates,
  pickDefaultTargetR,
  RGroupTemplate,
  RGroupTemplateItem,
  rGroupTargetWarnings,
  uniqueKeepMask,
  validateParams,
} from './pt-chem-enum';
import {_package} from '../package';
import {defaultErrorHandler} from '../utils/err-info';

const DIALOG_TITLE = 'Markush Enumerator';

const DEFAULT_TABLE_NAME = 'Markush enumeration';
const DEFAULT_REMOVE_DUPLICATES = true;

/** Single source of truth for the result-table name: a trimmed user value, or the default when blank/missing. */
function resolveTableName(raw: string | undefined): string {
  return raw?.trim() || DEFAULT_TABLE_NAME;
}

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
  onDuplicate?: () => void;
  onRemove?: () => void;
}

// Popular multi symbol single atoms for quick lookup in card builder
const SINGLE_ATOM_SYMBOLS_LOOKUP = new Set([
  'Cl', 'Br', 'Al', 'Si', 'Li', 'Na', 'Mg', 'Ca', 'Ti', 'At', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Kr', 'Rb',
  'Au', 'Ag', 'Pt', 'Pb', 'Sn', 'Sb', 'Te', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho']);

/** Draws a molecule into a fixed-size host, constraining SVG dimensions. */
function drawMolInto(host: HTMLElement, smi: string, w: number, h: number): void {
  ui.empty(host);
  try {
    const correctedSmi = smi.length === 1 || SINGLE_ATOM_SYMBOLS_LOOKUP.has(smi) ? `[${smi}]` : smi;
    const el = grok.chem.drawMolecule(correctedSmi, w, h);
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
  if (opts.smiles && !opts.error)
    drawMolInto(thumbHost, opts.smiles, THUMB_W, THUMB_H);
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
    editBtn.style.color = 'var(--blue-3)';
    actions.appendChild(editBtn);
  }
  if (opts.onDuplicate) {
    const dupBtn = ui.icons.copy((e: MouseEvent) => { e.stopPropagation(); opts.onDuplicate!(); }, 'Duplicate');
    dupBtn.style.color = 'var(--blue-3)';
    actions.appendChild(dupBtn);
  }
  if (opts.onRemove) {
    const delBtn = ui.icons.delete((e: MouseEvent) => { e.stopPropagation(); opts.onRemove!(); }, 'Remove');
    delBtn.style.color = 'var(--red-3)';
    actions.appendChild(delBtn);
  }
  card.addEventListener('mouseenter', () => { actions.style.display = 'flex'; });
  card.addEventListener('mouseleave', () => { actions.style.display = 'none'; });
  if (opts.onEdit || opts.onDuplicate || opts.onRemove) card.appendChild(actions);

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
    let isSingleAtom = false;
    let okBtn: HTMLButtonElement | null = null;
    let userTouchedR = false;
    rNumberInput.onChanged.subscribe(() => { userTouchedR = true; updateOk(); });

    const updateOk = () => {
      const n = rNumberInput.value;
      const ok = currentSmiles != null && n != null && n >= 1 &&
        (detectedRs.length === 1 || isSingleAtom);
      if (okBtn) okBtn.disabled = !ok;
    };

    const revalidate = async () => {
      currentSmiles = await smilesFromSketcher(sk, rdkit);
      detectedRs = currentSmiles ? extractRNumbers(currentSmiles) : [];
      // Probe single-atom mode only when there's no R-label — keeps the
      // common labeled path free of an extra RDKit parse.
      isSingleAtom = false;
      if (currentSmiles && detectedRs.length === 0) {
        const probe = makeRGroup(currentSmiles, rNumberInput.value ?? 1, '', rdkit);
        isSingleAtom = !!probe.isSingleAtom;
      }
      if (!currentSmiles) {
        status.innerText = 'Draw an R-group with one attachment point, or a single atom (e.g. N, O, Cl).';
        status.style.color = 'var(--grey-4)';
      } else if (detectedRs.length === 0 && !isSingleAtom) {
        status.innerText = 'No R-group detected — add [*:N] to mark the attachment point, or draw a single atom (e.g. N, O, Cl).';
        status.style.color = 'var(--red-3)';
      } else if (detectedRs.length === 0 && isSingleAtom) {
        const targetN = rNumberInput.value;
        status.innerText = `Single atom — will be substituted into [*:${targetN ?? '?'}] in the core.`;
        status.style.color = 'var(--green-2)';
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
    if (rNumberInput) {
      ui.tooltip.bind(rNumberInput.input,
        'Target R-number for these R-groups. Single-atom rows (no [*:N] label, e.g. just N or O) get substituted into the core\'s R# slot directly.');
    }
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
          const subtitle = rg.error ? 'invalid' :
            rg.isSingleAtom ? `R${rg.rNumber} · atom` :
              `R${rg.rNumber}${rg.sourceRNumber != null && rg.sourceRNumber !== rg.rNumber ? ` (from R${rg.sourceRNumber})` : ''}`;
          items.push({smi: rg.error ? '' : rg.smiles, err: rg.error, subtitle});
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

// ─── Copy R-group list to another slot ──────────────────────────────────────

/** R#s referenced by at least one valid (non-errored) core — feeds the default target and warnings. */
function collectCoreRNumbers(state: ChemEnumDialogState): Set<number> {
  const nums = new Set<number>();
  for (const c of state.cores) {
    if (c.error) continue;
    for (const rn of c.rNumbers) nums.add(rn);
  }
  return nums;
}


/**
 * Opens a dialog to copy R`srcN`'s substituent list into another R-slot — core-aware default
 * target, Append/Replace, advisory warnings. The copy + re-labeling is done by `copyRGroupList`.
 */
function openCopyToRGroupDialog(
  srcN: number,
  state: ChemEnumDialogState,
  rdkit: RDModule,
  refresh: () => void,
): void {
  const srcList = state.rGroupsByNum.get(srcN);
  if (!srcList || srcList.length === 0) {
    grok.shell.info(`R${srcN} has no substituents to copy.`);
    return;
  }

  // Snapshot (the dialog is modal over a static panel).
  const coreRNumbers = collectCoreRNumbers(state);
  const invalidCount = srcList.filter((rg) => rg.error != null).length;

  const defaultTarget = pickDefaultTargetR(new Set(state.rGroupsByNum.keys()), coreRNumbers, srcN);

  const targetInput = ui.input.int('Target R#', {value: defaultTarget, min: 1});
  const mergeInput = ui.input.choice<string>('Merge policy', {items: ['Append', 'Replace'], value: 'Append'});
  ui.tooltip.bind(mergeInput.input,
    'Append: add the copied substituents after any existing entries in the target slot.\n' +
    'Replace: clear the target slot first, then populate it with the copies.');

  const status = ui.divText('', {style: {fontSize: '11px', padding: '4px 0', minHeight: '18px'}});
  let okBtn: HTMLButtonElement | null = null;

  const fail = (m: string) => {
    status.innerText = m;
    status.style.color = 'var(--red-3)';
    if (okBtn) okBtn.disabled = true;
  };

  const updateStatus = () => {
    const t = targetInput.value;
    if (t == null || !Number.isInteger(t) || t < 1) return fail('Enter a valid R-group number (≥ 1).');
    if (t === srcN) return fail(`Target R# must differ from the source (R${srcN}).`);
    const existing = state.rGroupsByNum.get(t)?.length ?? 0;
    const replace = mergeInput.value === 'Replace';
    const noun = srcList.length === 1 ? 'substituent' : 'substituents';
    const msg = existing > 0 ?
      `Copy ${srcList.length} ${noun} R${srcN} → R${t} (${replace ? 'replacing' : 'appending to'} ${existing} existing).` :
      `Copy ${srcList.length} ${noun} R${srcN} → R${t} (new slot).`;
    const warnings = rGroupTargetWarnings(t, invalidCount, coreRNumbers);
    if (warnings.length > 0) {
      status.innerText = `${msg} ⚠ ${warnings.join('; ')}.`;
      status.style.color = 'var(--orange-2)';
    } else {
      status.innerText = msg;
      status.style.color = 'var(--green-2)';
    }
    if (okBtn) okBtn.disabled = false;
  };

  targetInput.onChanged.subscribe(updateStatus);
  mergeInput.onChanged.subscribe(updateStatus);

  const dialog = ui.dialog({title: `Copy R${srcN} to…`})
    .add(ui.divV([targetInput.root, mergeInput.root, status]))
    .onOK(() => {
      const t = targetInput.value;
      if (t == null || !Number.isInteger(t) || t < 1 || t === srcN) return;
      copyRGroupList(state.rGroupsByNum, srcN, t, mergeInput.value === 'Replace' ? 'replace' : 'append', rdkit);
      refresh();
    });
  dialog.show({resizable: false, width: 360});
  okBtn = dialog.getButton('OK') as HTMLButtonElement;
  updateStatus();
}

/**
 * Opens a dialog to insert a ready-made R-group template (e.g. the alkyl series) into an R-slot,
 * with a target R# and Append/Replace. The substituents are built + re-labeled by `addRGroupsFromSmiles`.
 */
// Loads every *.json in files/enumeration/r-group-templates/ once (memoized) and merges them, so the
// shipped catalogue and any client-dropped template files all show up in the picker. Client files are
// listed first and our default (DEFAULT_TEMPLATES_FILE) last, so the built-in sets sit at the bottom.
// A malformed file is skipped, not fatal; if nothing loads, falls back to BUILTIN (never empty).
const TEMPLATES_DIR = 'enumeration/r-group-templates';
const DEFAULT_TEMPLATES_FILE = 'r-group-templates.json';
let templatesPromise: Promise<RGroupTemplate[]> | null = null;
function loadRGroupTemplates(rdkit: RDModule): Promise<RGroupTemplate[]> {
  return templatesPromise ??= (async () => {
    const all: RGroupTemplate[] = [];
    try {
      const files = (await _package.files.list(TEMPLATES_DIR)).filter((f) => f.extension.toLowerCase() === 'json');
      files.sort((a, b) =>
        (a.fileName === DEFAULT_TEMPLATES_FILE ? 1 : 0) - (b.fileName === DEFAULT_TEMPLATES_FILE ? 1 : 0) ||
        a.fileName.localeCompare(b.fileName));
      for (const f of files) {
        try {
          for (const t of parseRGroupTemplates(await _package.files.readAsText(f))) {
            const items = t.items.filter((it) => isValidTemplateSmiles(it.smiles, rdkit));
            if (items.length > 0) all.push({...t, items});
          }
        } catch (e) { _package.logger.warning(`Skipped R-group template file ${f.fileName}: ${e}`); }
      }
    } catch { /* fall back to the built-in set below */ }
    return all.length > 0 ? all : BUILTIN_R_GROUP_TEMPLATES;
  })();
}

function openAddTemplateDialog(
  templates: RGroupTemplate[], state: ChemEnumDialogState, rdkit: RDModule, refresh: () => void,
): void {
  const coreRNumbers = collectCoreRNumbers(state);

  const templateInput = ui.input.choice<string>('Template',
    {items: templates.map((t) => t.name), value: templates[0].name});
  const defaultTarget = pickDefaultTargetR(new Set(state.rGroupsByNum.keys()), coreRNumbers);
  const targetInput = ui.input.int('Target R#', {value: defaultTarget, min: 1});
  const mergeInput = ui.input.choice<string>('Merge policy', {items: ['Append', 'Replace'], value: 'Append'});
  ui.tooltip.bind(mergeInput.input,
    'Append: add the template after any existing entries.\nReplace: overwrite the target slot.');

  const status = ui.divText('', {style: {fontSize: '11px', padding: '4px 0', minHeight: '18px'}});
  let okBtn: HTMLButtonElement | null = null;
  const fail = (m: string) => { status.innerText = m; status.style.color = 'var(--red-3)'; if (okBtn) okBtn.disabled = true; };

  // `working` is the curated subset that will be inserted. It starts as the whole template; each
  // preview card has a hover-revealed ✕ that drops its substituent. Switching templates resets it.
  let working: RGroupTemplateItem[] = [];
  const resetWorking = () => {
    const tmpl = templates.find((x) => x.name === templateInput.value);
    working = tmpl ? [...tmpl.items] : [];
  };

  // Thumbnail grid of the working set; hover a card to enlarge or ✕ to remove it before inserting.
  const previewHost = ui.div([], {style: {display: 'flex', flexWrap: 'wrap', overflowY: 'auto',
    gap: '4px', padding: '4px 0', maxWidth: '100%', maxHeight: `${CARD_H * 2 + 16}px`, minHeight: `${CARD_H}px`}});
  const renderPreview = () => {
    ui.empty(previewHost);
    for (const item of working) {
      // Show a generic `*` attachment (not `*:1`) so the preview reads the same for any target R#.
      const display = item.smiles.replace(/\[\*:\d+\]/g, '[*]');
      // Remove only this card's node (template items are unique) — no redraw of the survivors.
      const card = buildCard({smiles: display, subtitle: item.label ?? display, onRemove: () => {
        const i = working.indexOf(item);
        if (i >= 0) working.splice(i, 1);
        card.remove();
        updateStatus();
      }});
      previewHost.appendChild(card);
    }
  };

  const updateStatus = () => {
    const t = targetInput.value;
    if (t == null || !Number.isInteger(t) || t < 1) return fail('Enter a valid R-group number (≥ 1).');
    if (working.length === 0) return fail('No substituents selected — keep at least one.');
    const existing = state.rGroupsByNum.get(t)?.length ?? 0;
    const replace = mergeInput.value === 'Replace';
    const msg = existing > 0 ?
      `Add ${working.length} substituents to R${t} (${replace ? 'replacing' : 'appending to'} ${existing} existing).` :
      `Add ${working.length} substituents to R${t} (new slot).`;
    const warnings = rGroupTargetWarnings(t, 0, coreRNumbers);
    status.innerText = warnings.length > 0 ? `${msg} ⚠ ${warnings.join('; ')}.` : msg;
    status.style.color = warnings.length > 0 ? 'var(--orange-2)' : 'var(--green-2)';
    if (okBtn) okBtn.disabled = false;
  };

  templateInput.onChanged.subscribe(() => { resetWorking(); renderPreview(); updateStatus(); });
  targetInput.onChanged.subscribe(updateStatus);
  mergeInput.onChanged.subscribe(updateStatus);

  const dialog = ui.dialog({title: 'Add R-group template'})
    .add(ui.divV([templateInput.root, previewHost, targetInput.root, mergeInput.root, status]))
    .onOK(() => {
      const t = targetInput.value;
      if (t == null || !Number.isInteger(t) || t < 1 || working.length === 0) return;
      const smiles = working.map((it) => it.smiles);
      addRGroupsFromSmiles(state.rGroupsByNum, smiles, t, mergeInput.value === 'Replace' ? 'replace' : 'append', rdkit);
      refresh();
    });
  dialog.show({resizable: true, width: 520});
  okBtn = dialog.getButton('OK') as HTMLButtonElement;
  resetWorking();
  renderPreview();
  updateStatus();
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
  removeDuplicates: boolean;
  tableName: string;
}

// ─── History (localStorage-backed, shared by dialog + app) ──────────────────

const HISTORY_KEY = 'd4-markush-enumeration-history';
const HISTORY_MAX = 10;

interface ChemEnumHistoryEntry {
  /** ISO timestamp of when the enumeration was recorded. */
  date: string;
  /** One-line summary: `"6 cores · R1:6 · R2:4"`. */
  summary: string;
  mode: ChemEnumMode;
  /** Original (unnormalized) SMILES for each core. */
  cores: string[];
  /** R-number → list of original R-group SMILES. JSON object so we can round-trip. */
  rGroups: { [rNumber: string]: string[] };
  /** Optional for backward-compatibility with entries stored before this option existed. */
  removeDuplicates?: boolean;
  /** Optional for backward-compatibility with entries stored before this option existed. */
  tableName?: string;
}

function readHistory(): ChemEnumHistoryEntry[] {
  try {
    const raw = localStorage.getItem(HISTORY_KEY);
    if (!raw) return [];
    const arr = JSON.parse(raw);
    return Array.isArray(arr) ? arr : [];
  } catch {
    return [];
  }
}

function writeHistory(entries: ChemEnumHistoryEntry[]): void {
  try {
    localStorage.setItem(HISTORY_KEY, JSON.stringify(entries));
  } catch {
    // Quota exceeded or storage unavailable — silently skip; history is best-effort.
  }
}

function summarizeState(state: ChemEnumDialogState): string {
  const partsR: string[] = [];
  const sortedNums = [...state.rGroupsByNum.keys()].sort((a, b) => a - b);
  for (const n of sortedNums) partsR.push(`R${n}:${state.rGroupsByNum.get(n)!.length}`);
  return `${state.cores.length} core${state.cores.length === 1 ? '' : 's'}` +
    (partsR.length ? ' · ' + partsR.join(' · ') : '');
}

function buildHistoryEntry(state: ChemEnumDialogState): ChemEnumHistoryEntry {
  const rGroups: { [n: string]: string[] } = {};
  for (const [n, list] of state.rGroupsByNum)
    rGroups[String(n)] = list.map((rg) => rg.originalSmiles);
  return {
    date: new Date().toISOString(),
    summary: summarizeState(state),
    mode: state.mode,
    cores: state.cores.map((c) => c.originalSmiles),
    rGroups,
    removeDuplicates: state.removeDuplicates,
    tableName: resolveTableName(state.tableName),
  };
}

/** Prepends `state` as a new entry, keeps at most {@link HISTORY_MAX}. */
function recordHistory(state: ChemEnumDialogState): void {
  const entry = buildHistoryEntry(state);
  const prev = readHistory();
  writeHistory([entry, ...prev].slice(0, HISTORY_MAX));
}

/** Reparses SMILES through `makeCore` / `makeRGroup` so validation runs with the current rdkit. */
function applyHistoryEntry(entry: ChemEnumHistoryEntry, state: ChemEnumDialogState, rdkit: RDModule): void {
  state.cores = entry.cores.map((smi) => makeCore(smi, '', rdkit));
  state.rGroupsByNum = new Map();
  for (const [nStr, list] of Object.entries(entry.rGroups ?? {})) {
    const num = parseInt(nStr, 10);
    if (!Number.isFinite(num)) continue;
    state.rGroupsByNum.set(num, list.map((smi) => makeRGroup(smi, num, '', rdkit)));
  }
  if (entry.mode === ChemEnumModes.Zip || entry.mode === ChemEnumModes.Cartesian)
    state.mode = entry.mode;
  state.removeDuplicates = entry.removeDuplicates ?? DEFAULT_REMOVE_DUPLICATES;
  state.tableName = resolveTableName(entry.tableName);
}

function formatHistoryDate(iso: string): string {
  try {
    const d = new Date(iso);
    return `${d.toLocaleDateString()} ${d.toLocaleTimeString([], {hour: '2-digit', minute: '2-digit'})}`;
  } catch {
    return iso;
  }
}

/** Pops a context menu listing history entries. `onPick` applies + refreshes. */
function showHistoryMenu(onPick: (entry: ChemEnumHistoryEntry) => void): void {
  const entries = readHistory();
  const menu = DG.Menu.popup();
  if (entries.length === 0) {
    menu.item('No entries', () => {});
  } else {
    for (const entry of entries) {
      const label = `${formatHistoryDate(entry.date)} — ${entry.summary}`;
      menu.item(label, () => onPick(entry));
    }
  }
  menu.show();
}

// ─── CSV preload (app mode only, when history is empty) ─────────────────────

/**
 * Seeds `state` from the shipped `files/enumeration/{cores,rgroups}.csv`. The cores CSV has a
 * single "Core" column; the R-groups CSV has one column per R-number (`R1`, `R2`, …).
 * Returns whether anything was loaded.
 */
async function preloadFromFiles(state: ChemEnumDialogState, rdkit: RDModule): Promise<boolean> {
  try {
    const coresText = await _package.files.readAsText('enumeration/chem_enum_cores.csv');
    const coresDf = DG.DataFrame.fromCsv(coresText);
    const coreCol = coresDf.col('Core') ?? coresDf.columns.byIndex(0);
    if (coreCol) {
      for (let i = 0; i < coreCol.length; i++) {
        if (coreCol.isNone(i)) continue;
        const v = String(coreCol.get(i) ?? '').trim();
        if (!v) continue;
        state.cores.push(makeCore(v, '', rdkit));
      }
    }

    const rgText = await _package.files.readAsText('enumeration/chem_enum_rgroups.csv');
    const rgDf = DG.DataFrame.fromCsv(rgText);
    for (const col of rgDf.columns.toList()) {
      const m = col.name.match(/^r\s*(\d+)$/i);
      if (!m) continue;
      const num = parseInt(m[1], 10);
      const list: ChemEnumRGroup[] = [];
      for (let i = 0; i < col.length; i++) {
        if (col.isNone(i)) continue;
        const v = String(col.get(i) ?? '').trim();
        if (!v) continue;
        list.push(makeRGroup(v, num, '', rdkit));
      }
      if (list.length > 0) state.rGroupsByNum.set(num, list);
    }
    return state.cores.length > 0 || state.rGroupsByNum.size > 0;
  } catch {
    return false;
  }
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
      .addButton('History', () => {
        showHistoryMenu((entry) => {
          applyHistoryEntry(entry, panel.state, rdkit);
          panel.refresh();
        });
      })
      .onOK(async () => { await panel.execute(); });
    panel.bindActionButton(dialog.getButton('OK') as HTMLButtonElement);
    dialog.show({resizable: true, width: 960});
    const historyButton = dialog.getButton('History') as HTMLButtonElement;
    if (historyButton) {
      historyButton.style.order = '1';
      historyButton.style.marginRight = 'auto';
      ui.tooltip.bind(historyButton, 'View and apply past enumerations');
      historyButton.innerHTML = ui.iconFA('history', () => {}).outerHTML;
    }
  } catch (err: any) {
    defaultErrorHandler(err);
  }
}

/** Opens the chem enumeration UI as a top-level view (app entry point). */
export async function polyToolEnumerateChemApp(): Promise<DG.View | null> {
  await _package.initPromise;
  try {
    const rdkit = await getRdKitModule();
    const panel = buildChemEnumPanel(rdkit, null, 'app');
    const runBtn = ui.button('Enumerate', async () => { await panel.execute(); });
    panel.bindActionButton(runBtn as HTMLButtonElement);
    panel.appActionHost?.appendChild(runBtn);

    // Seed the panel: most-recent history wins; if there is none, fall back to the shipped
    // demo CSVs so the app never opens to an empty screen. Both paths mutate panel.state and
    // panel.refresh() picks the changes up.
    const history = readHistory();
    if (history.length > 0) {
      applyHistoryEntry(history[0], panel.state, rdkit);
      panel.refresh();
    } else {
      preloadFromFiles(panel.state, rdkit).then((loaded) => {
        if (loaded) panel.refresh();
      });
    }

    const view = DG.View.create();
    view.name = DIALOG_TITLE;
    view.box = true;
    view.root.appendChild(ui.div([panel.root],
      {style: {height: '100%', width: '100%', padding: '8px', boxSizing: 'border-box'}}));
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
  /** Re-render after mutating `state` from outside (e.g. applying a history entry). */
  refresh: () => void;
  /**
   * Bottom-right "controls" cell host for the caller's Enumerate button.
   * Only populated in app layout; undefined in dialog layout.
   */
  appActionHost?: HTMLElement;
}

type ChemEnumLayout = 'dialog' | 'app';

function buildChemEnumPanel(
  rdkit: RDModule, preloadCore: ChemEnumCore | null, layout: ChemEnumLayout = 'dialog',
): ChemEnumPanel {
  const state: ChemEnumDialogState = {
    cores: preloadCore ? [preloadCore] : [],
    rGroupsByNum: new Map(),
    mode: ChemEnumModes.Cartesian,
    appendToTable: null,
    removeDuplicates: DEFAULT_REMOVE_DUPLICATES,
    tableName: DEFAULT_TABLE_NAME,
  };

  // ── Cores: single-row horizontal virtualView (dialog) or wrapped flow (app) ─
  const ROW_H = CARD_H + 16; // card height + horizontal scrollbar breathing room
  const coresEmpty = ui.divText('No cores — draw or import at least one.', {style: {color: 'var(--grey-4)', padding: '20px 12px', fontSize: '12px'}});
  // App layout puts cores into the top-left grid cell — fill height, scroll vertically as cards wrap.
  // Dialog layout keeps the original single-row horizontal scroller.
  const coresVvHost = layout === 'app' ?
    ui.div([], {style: {width: '100%', flex: '1 1 auto', minHeight: '0', overflowY: 'auto', overflowX: 'hidden'}}) :
    ui.div([], {style: {width: '100%', height: `${ROW_H}px`, overflow: 'hidden'}});
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
      onDuplicate: async () => {
        const dup = await openCoreSketchDialog(rdkit, c.smiles);
        if (dup) { state.cores.push(dup); refresh(); }
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
    if (layout === 'app') {
      // Cores are realistically dozens at most — render all directly in a wrapped flow so
      // the top-left cell uses its full area instead of constraining cards to one row.
      ui.empty(coresVvHost);
      coresVv = null;
      const wrap = ui.div([], {style: {
        display: 'flex', flexWrap: 'wrap', gap: '4px',
        width: '100%', alignContent: 'flex-start', padding: '2px 0',
      }});
      for (let i = 0; i < state.cores.length; i++)
        wrap.appendChild(coresRenderer(i));
      coresVvHost.appendChild(wrap);
      return;
    }
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
      const copyBtn = ui.button('Copy to…', () => openCopyToRGroupDialog(n, state, rdkit, refresh));
      copyBtn.style.marginLeft = '4px';
      ui.tooltip.bind(copyBtn, `Copy all R${n} substituents into another R-group slot`);
      const clearBtn = ui.button('Remove all', () => { state.rGroupsByNum.delete(n); refresh(); });
      clearBtn.style.marginLeft = '4px';
      ui.tooltip.bind(clearBtn, `Remove all R${n} groups`);
      const header = ui.divH([label, copyBtn, clearBtn], {style: {alignItems: 'center', justifyContent: 'flex-start', gap: '0', margin: '6px 0 2px'}});
      const renderer = (i: number): HTMLElement => {
        const rg = list[i];
        const subtitle = rg.error ? 'invalid' :
          rg.isSingleAtom ? `r group ${i + 1} · R${rg.rNumber} · atom` :
            `r group ${i + 1} · R${rg.rNumber}`;
        return buildCard({
          smiles: rg.error ? '' : rg.smiles,
          subtitle,
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
          onDuplicate: async () => {
            const dup = await openRGroupSketchDialog(rdkit, rg.smiles, rg.rNumber);
            if (!dup) return;
            const target = state.rGroupsByNum.get(dup.rNumber) ?? [];
            target.push(dup);
            state.rGroupsByNum.set(dup.rNumber, target);
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
  // Dialog caps preview height at 450px (fits next to the cores/r-groups column).
  // App layout drops the cap so preview fills its bottom-left cell.
  const previewHost = ui.div([], {style: layout === 'app' ? {
    width: '100%', height: '100%', display: 'flex', flexWrap: 'wrap', gap: '6px',
    alignContent: 'flex-start', padding: '2px 0', flex: '1 1 auto', minHeight: '0', overflow: 'auto',
  } : {
    width: '100%', display: 'flex', flexWrap: 'wrap', gap: '6px',
    alignContent: 'flex-start', padding: '2px 0', maxHeight: '450px', overflow: 'scroll',
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

  const removeDuplicatesInput = ui.input.bool('Remove duplicates', {
    value: state.removeDuplicates,
    onValueChanged: (v) => { state.removeDuplicates = !!v; },
  });
  ui.tooltip.bind(removeDuplicatesInput.input,
    'Removes duplicate molecules from the final enumerated output (compared by canonical SMILES). ' +
    'This is distinct from the import wizard\'s dedup, which collapses identical input cores/R-groups.');

  const tableNameInput = ui.input.string('Table name', {
    value: state.tableName,
    onValueChanged: (v) => { state.tableName = v ?? ''; },
  });
  ui.tooltip.bind(tableNameInput.input,
    'Name of the result table created by the enumeration (ignored when appending to an existing table).');

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
    // Re-sync the option inputs after external state mutation (e.g. applying a history entry).
    // Their onValueChanged only writes state (no refresh()), so this can't re-enter; modeInput is
    // left untouched on purpose — setting its value would re-enter refresh() and loop.
    removeDuplicatesInput.value = state.removeDuplicates;
    tableNameInput.value = state.tableName;
    redrawPreview();
  };

  // ── Add/import handlers. Input dedup (identical cores/R-groups) is gated by the import
  //    wizard's checkbox; output dedup (identical products) is the "Remove duplicates" option. ──
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
  const addRGroupsFromTemplate = () =>
    loadRGroupTemplates(rdkit).then((templates) => openAddTemplateDialog(templates, state, rdkit, refresh));

  // ── CSV export — one combined file. First column "Core", then one column per R# (R1, R2, …).
  // The columns have different lengths (cores vs each R# list), so shorter ones are padded with
  // empty cells. Not the tidiest table, but it carries the whole setup in a single CSV.
  const exportAllCsv = () => {
    const cols = buildExportColumns(state.cores, state.rGroupsByNum);
    if (cols.length === 0) { grok.shell.info('Nothing to export.'); return; }
    const df = DG.DataFrame.fromColumns(cols.map((c) => DG.Column.fromStrings(c.name, c.values)));
    DG.Utils.download('markush-enumeration.csv', df.toCsv());
  };

  const sectionHeader = (
    label: string,
    onDraw?: () => void, onImport?: () => void, onClear?: () => void, onTemplate?: () => void,
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
    if (onTemplate) {
      // Icon + label, matching "+ Draw" / "↓ Import…": leading swatchbook glyph at the same size as the +/↓.
      const tmplIcon = ui.iconFA('swatchbook');
      tmplIcon.style.fontSize = 'inherit';
      tmplIcon.style.marginRight = '4px';
      const templateBtn = ui.button([tmplIcon, 'Templates'], onTemplate, `${label}: insert a ready-made set`);
      templateBtn.style.marginLeft = '4px';
      parts.push(templateBtn);
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

  // ── Section headers (shared between layouts) ──────────────────────────────
  const coresHeader = sectionHeader(
    'Cores', addCoreFromSketcher, addCoresFromImport,
    () => { state.cores.splice(0, state.cores.length); refresh(); },
  );
  coresClearBtn = coresHeader.clearBtn;

  // Single combined-export icon — `arrow-to-bottom`, the same glyph the grid's download button uses.
  const exportBtn = ui.iconFA('arrow-to-bottom', exportAllCsv, 'Download cores + R-groups as one CSV');

  const rGroupsHeader = sectionHeader(
    'R-Groups', addRGroupFromSketcher, addRGroupsFromImport,
    () => { state.rGroupsByNum.clear(); refresh(); },
    addRGroupsFromTemplate,
  );
  rGroupsClearBtn = rGroupsHeader.clearBtn;

  let body: HTMLElement;
  let appActionHost: HTMLElement | undefined;

  if (layout === 'app') {
    // ── App: two flex columns. Left = cores + preview. Right = r-groups (grows) + controls (compact). ──
    const cellBaseStyle = {
      display: 'flex', flexDirection: 'column',
      minHeight: '0', minWidth: '0',
      padding: '8px', boxSizing: 'border-box',
      border: '1px solid var(--grey-2)', borderRadius: '4px',
      background: 'var(--white)', overflow: 'hidden',
    } as const;
    const growCellStyle = {...cellBaseStyle, flex: '1 1 0'} as const;

    countText.style.fontSize = '11px';
    countText.style.color = 'var(--grey-5)';

    const coresCell = ui.divV([
      coresHeader.root,
      coresVvHost,
      coresEmpty,
    ], {style: growCellStyle});

    const rGroupsCell = ui.divV([
      rGroupsHeader.root,
      rGroupsHost,
    ], {style: growCellStyle});

    const previewCell = ui.divV([
      sectionHeader('Preview').root,
      previewHost,
    ], {style: growCellStyle});

    // Run header carries the live status — count + issues — inline on the right.
    const runLabel = ui.divText('Run', {style: {
      fontWeight: '600', fontSize: '12px', color: 'var(--grey-6)', alignSelf: 'center', minWidth: '40px',
    }});
    const runHeaderStats = ui.divH([countText, errorBadge.root], {style: {
      alignItems: 'center', gap: '8px', marginLeft: 'auto',
    }});
    const runHeader = ui.divH([runLabel, runHeaderStats], {style: {
      alignItems: 'center', gap: '8px', margin: '0 0 6px',
      flex: '0 0 auto', width: '100%',
    }});

    // History icon — sits next to Enumerate. Reads localStorage on click and pops a menu.
    const historyBtn = ui.iconFA('history', () => {
      showHistoryMenu((entry) => {
        applyHistoryEntry(entry, state, rdkit);
        refresh();
      });
    }, 'Show recent enumerations');
    historyBtn.style.fontSize = '16px';
    historyBtn.style.padding = '4px 6px';
    historyBtn.style.cursor = 'pointer';
    historyBtn.style.color = 'var(--blue-3)';

    const clearAllBtn = ui.button('Clear all', () => {
      state.cores.splice(0, state.cores.length);
      state.rGroupsByNum.clear();
      refresh();
    });
    ui.tooltip.bind(clearAllBtn, 'Remove all cores and R-groups');

    appActionHost = ui.div([], {style: {flex: '0 0 auto'}});
    const actionRow = ui.divH([historyBtn, appActionHost], {style: {
      alignItems: 'center', gap: '8px', marginTop: '8px', flex: '0 0 auto',
    }});

    const controlsCell = ui.divV([
      runHeader,
      modeInput.root,
      appendToTableInput.root,
      removeDuplicatesInput.root,
      tableNameInput.root,
      ui.divH([clearAllBtn, exportBtn], {style: {gap: '8px', alignItems: 'center'}}),
      actionRow,
    ], {style: {
      ...cellBaseStyle, flex: '0 0 auto', overflow: 'visible',
    }});

    const leftColumn = ui.divV([coresCell, previewCell], {style: {
      flex: '1 1 50%', minWidth: '0',
      display: 'flex', flexDirection: 'column', gap: '8px',
    }});

    const rightColumn = ui.divV([rGroupsCell, controlsCell], {style: {
      flex: '1 1 50%', minWidth: '0',
      display: 'flex', flexDirection: 'column', gap: '8px',
    }});

    body = ui.divH([leftColumn, rightColumn], {style: {
      width: '100%', height: '100%',
      display: 'flex', flexDirection: 'row',
      gap: '8px',
      padding: '4px', boxSizing: 'border-box',
    }});
  } else {
    // ── Dialog: horizontal split — cores + r-groups (left 60%) | preview (right 40%) ──
    const coresSection = ui.divV([
      coresHeader.root,
      ui.div([coresVvHost, coresEmpty], {style: {padding: '0'}}),
    ], {style: sectionStyle});

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

    const footer = ui.divV([
      ui.divH([appendToTableInput.root, exportBtn],
        {style: {alignItems: 'center', gap: '8px'}}),
      removeDuplicatesInput.root,
      tableNameInput.root,
    ], {style: {padding: '2px 0'}});

    body = ui.divV([
      split, statusLine, footer,
    ], {style: {
      width: '100%', padding: '4px',
      display: 'flex', flexDirection: 'column', gap: '4px',
    }});
  }

  refresh();
  return {
    root: body,
    state,
    execute: () => executeEnumeration(state, rdkit),
    bindActionButton: (btn) => { okButton = btn; refresh(); },
    refresh,
    appActionHost,
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

    // Record this configuration in history before any post-processing — if the user closes
    // mid-canonicalization, they can still recover the inputs that produced something.
    recordHistory(state);

    const rNumbersUsed = new Set<number>();
    for (const r of results) for (const n of r.rGroupSmilesByNum.keys()) rNumbersUsed.add(n);
    const sortedRs = [...rNumbersUsed].sort((a, b) => a - b);

    const smilesCol = DG.Column.fromStrings('Enumerated', results.map((r) => r.smiles));
    smilesCol.semType = DG.SEMTYPE.MOLECULE;
    const coreCol = DG.Column.fromStrings('Core', results.map((r) => normalizeRLabels(r.coreSmiles ?? '')));
    coreCol.semType = DG.SEMTYPE.MOLECULE;
    const rCols = sortedRs.map((n) =>
      DG.Column.fromStrings(`R${n}`, results.map((r) => r.rGroupSmilesByNum.get(n) ?? '')));
    for (const c of rCols) c.semType = DG.SEMTYPE.MOLECULE;

    let df = DG.DataFrame.fromColumns([smilesCol, coreCol, ...rCols]);

    // Stage 2 — canonicalize the whole Enumerated column in parallel via Chem workers.
    pi.update(40, `Canonicalizing ${results.length.toLocaleString()} molecule(s)...`);
    let canonical: string[] | null = null;
    try {
      const res: DG.Column = await grok.functions.call('Chem:convertNotation', {
        data: df,
        molecules: smilesCol,
        targetNotation: DG.chem.Notation.Smiles,
        overwrite: false,
        join: false,
        kekulize: false,
      });
      // in older version of the chem, overwrite is super slow, it has been updated but we can do it like this here
      const resArr = res.toList();
      smilesCol.init((i) => resArr[i]);
      smilesCol.meta.units = DG.chem.Notation.Smiles;
      canonical = resArr;
    } catch (err: any) {
      // Canonicalization is a nice-to-have; the uncanonical SMILES are still valid output.
      _package.logger.warning(`Canonicalization skipped: ${err?.message ?? err}`);
    }

    // Dedup by canonical SMILES — reuse the array materialized during canonicalization
    // (falling back to the column only when canonicalization was skipped).
    if (state.removeDuplicates) {
      const mask = uniqueKeepMask(canonical ?? smilesCol.toList());
      const keep = DG.BitSet.create(df.rowCount, (i) => mask[i]);
      const removed = df.rowCount - keep.trueCount;
      if (removed > 0) {
        df = df.clone(keep);
        grok.shell.info(`Removed ${removed} duplicate molecule${removed === 1 ? '' : 's'}.`);
      }
    }

    pi.update(90, 'Finalizing...');
    await grok.data.detectSemanticTypes(df);

    if (state.appendToTable) {
      state.appendToTable.append(df, true);
      await state.appendToTable.meta.detectSemanticTypes();
    } else {
      // Set name on the final frame: `df` may have been reassigned by `clone()` above,
      // and `name` is a dedicated Dart property that `clone(saveTags)` does not carry.
      df.name = resolveTableName(state.tableName);
      grok.shell.addTableView(df);
    }
  } catch (err: any) {
    defaultErrorHandler(err);
  } finally {
    pi.close();
  }
}
