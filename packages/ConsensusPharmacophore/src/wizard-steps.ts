/*
 * Per-step content builders for the Consensus Pharmacophore wizard.
 *
 * Each `buildStepN(ctx)` returns an HTMLElement that the wizard shell mounts
 * into its left content panel. Step 1 holds the input form (lifted from the
 * old `buildInputPanel` in `orchestrator.ts`). Steps 2-5 show the active
 * stage's result data as compact HTML tables.
 *
 * Step builders are pure rendering functions — they do NOT run stages.
 * Stage execution is owned by `WizardShell.runCurrentStep()` which calls
 * the orchestrator's existing `previewFetch / previewAlign / ...` methods.
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DEFAULT_OPTIONS, PipelineOptions, PocketMethod, PocketRep} from './orchestrator-types';
import {FAMILY_CODES, FAMILY_MAP, resolveFamily} from './family-map';
import {loadDemoPdbIds, parsePdbIds, parsePdbIdsDetailed,
  SUMMARY_FAMILY_COL} from './orchestrator';
import {DroppedPdb, fetchPdbSummaries, PdbSummary, searchPdbsByLigand} from './rcsb-client';
import {fetchUniProtById, getPdbIdsForAccession, isUniProtAccession,
  searchUniProt, UniProtHit} from './uniprot-client';
import {isValidRcsbId, looksLikePdb, registerLocalPdb} from './local-pdb-store';

// ---------------------------------------------------------------------------
// Step context — passed from the wizard shell to each step builder
// ---------------------------------------------------------------------------

/**
 * Per-step context object. The shell constructs one of these each time
 * `goToStep()` runs and hands it to the step builder. Step builders read
 * the current state and register callbacks for input changes.
 */
export interface StepContext {
  // ---- read-only state ----
  /** Latest parsed PDB IDs from the textarea. */
  pdbIds: string[];
  /** Latest pipeline options (resolution, pocket, consensus knobs). */
  options: PipelineOptions;
  /** Whether the user came back here AFTER editing inputs further downstream. */
  isStale: boolean;
  /** Stage outputs (may be null when the corresponding stage has not run). */
  pdbQcData: DG.DataFrame | null;
  /** Per-PDB drop records from Stage 1 (resolution / method / ligand filters).
   *  Empty until Stage 1 has run. Used by Step 1 to render the "Dropped" section. */
  pdbQcDropped: DroppedPdb[];
  alignedData: DG.DataFrame | null;
  pocketAtoms: DG.DataFrame | null;
  summaryData: DG.DataFrame | null;
  consensusData: DG.DataFrame | null;
  /** Live count of interactions / PDBs the current consensus exclusions leave
   *  in (after both per-PDB and per-interaction filtering). */
  consensusUsed?: {interactions: number; pdbs: number};

  // ---- callbacks into the shell ----
  /** Fired when the user edits the PDB textarea or runs the demo loader. */
  onPdbInputChanged: (ids: string[], demoLoaded?: boolean) => void;
  /** Fired when any input in the Advanced Options accordion changes. */
  onOptionsChanged: (next: PipelineOptions) => void;
  /** Re-render the current step's content panel. Used by interactions that
   *  change the step's own state in a way that can't be reflected via direct
   *  DOM mutation — e.g. a Step 2 row click that needs to update both the
   *  dropdown selection AND the table row highlight AND the hint text.
   *  Must be lightweight: re-renders only the step content (left column),
   *  not the persistent Mol* viewer or detail panel. */
  refreshStep?: () => void;
  /** Fired when the user clicks a PDB row in the Step 2 RMSD table to
   *  isolate / release that PDB in the Mol* viewer. Currently unused — the
   *  row click was repurposed to "set as alignment reference" (see buildStep2),
   *  but the callback + orchestrator.isolatePdb() are kept for a possible
   *  future re-introduction (e.g. via an eye-icon column). Pass null to
   *  release the isolation. */
  onIsolatePdb?: (pdbId: string | null) => void;
  /** Currently-isolated PDB id, for visual highlighting of the row. Unused
   *  in the current row-click-sets-reference UX; see onIsolatePdb above. */
  isolatedPdb?: string | null;
  /** Step 4 "Use" checkbox handler — hide (`hidden = true`) or show this PDB's
   *  protein + interaction overlays in Mol*. Complements onOptionsChanged,
   *  which updates the consensus exclusion; this drives the 3D visibility. */
  onTogglePdbVisibility?: (pdbId: string, hidden: boolean) => void;

  // ---- one-time seed values ----
  /** Initial PDB IDs (for Step 1's textarea population). */
  seedPdbIds?: string[];
  /** Whether the demo label should be shown above the textarea on first paint. */
  demoLoaded?: boolean;
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/** Family color chip — a small colored dot inline before a family name. */
function familyChip(code: string): HTMLElement {
  const fam = resolveFamily(code);
  const dot = document.createElement('span');
  dot.className = 'cp-family-chip';
  dot.style.background = fam.hexColor;
  return dot;
}

function statItem(value: string | number, label: string): HTMLElement {
  return ui.divV([
    ui.divText(String(value), 'cp-wizard-stat-val'),
    ui.divText(label, 'cp-wizard-stat-label'),
  ], 'cp-wizard-stat-item');
}

function statStrip(items: Array<{value: string | number; label: string}>): HTMLElement {
  return ui.divH(items.map((it) => statItem(it.value, it.label)), 'cp-wizard-stat-strip');
}

function panelTitle(title: string, description?: string): HTMLElement {
  const els: HTMLElement[] = [ui.divText(title, 'cp-step-panel-title')];
  if (description) els.push(ui.divText(description, 'cp-step-panel-description'));
  return ui.divV(els);
}

function placeholderText(text: string): HTMLElement {
  const el = ui.divText(text);
  el.style.fontStyle = 'italic';
  el.style.color = '#888';
  el.style.padding = '20px 0';
  el.style.textAlign = 'center';
  el.style.fontSize = '12px';
  return el;
}

/** Build a compact `<table>` from a DataFrame given column names and row limit.
 *  The optional `format` callback receives the cell value AND the row index
 *  so format functions can synthesize whole-row content (e.g. the Step 4
 *  "Families" column renders a chip row built from per-family columns). */
function dataFrameTable(
  df: DG.DataFrame,
  columns: ReadonlyArray<{
    key: string;
    header: string;
    format?: (v: any, rowIndex: number, df: DG.DataFrame) => string | HTMLElement;
  }>,
  maxRows = 50,
): HTMLElement {
  const table = document.createElement('table');
  table.className = 'cp-step-table';
  const thead = document.createElement('thead');
  const tr = document.createElement('tr');
  for (const c of columns) {
    const th = document.createElement('th');
    th.textContent = c.header;
    tr.append(th);
  }
  thead.append(tr);
  table.append(thead);

  const tbody = document.createElement('tbody');
  const n = Math.min(df.rowCount, maxRows);
  for (let i = 0; i < n; i++) {
    const row = document.createElement('tr');
    for (const c of columns) {
      const td = document.createElement('td');
      const col = df.columns.byName(c.key);
      // null when column is missing — format() can then distinguish missing
      // from an empty/zero value and show "–" appropriately.
      const v = col ? col.get(i) : null;
      if (c.format) {
        const formatted = c.format(v, i, df);
        if (formatted instanceof HTMLElement) td.append(formatted);
        else td.textContent = formatted;
      } else {
        td.textContent = (v == null) ? '–' : String(v);
      }
      row.append(td);
    }
    tbody.append(row);
  }
  table.append(tbody);

  if (df.rowCount > maxRows) {
    const tfoot = document.createElement('tfoot');
    const tr2 = document.createElement('tr');
    const td = document.createElement('td');
    td.colSpan = columns.length;
    td.textContent = `... +${df.rowCount - maxRows} more rows`;
    td.style.textAlign = 'center';
    td.style.color = '#888';
    td.style.fontStyle = 'italic';
    td.style.padding = '4px';
    tr2.append(td);
    tfoot.append(tr2);
    table.append(tfoot);
  }

  return table;
}

// ---------------------------------------------------------------------------
// Per-step option sections
//
// Each option section builds its own input controls from `ctx.options`,
// wires the onChanged handlers to call `ctx.onOptionsChanged` with a
// merged PipelineOptions (its own fields + everything else from ctx.options).
// This keeps each step's options visible right where they apply (QC on
// Step 1, Pocket on Step 3, Consensus on Step 5) instead of bunching them
// all under Step 1's "Advanced options" accordion.
// ---------------------------------------------------------------------------

/** Stage 1 (Fetch) filters: resolution, experimental methods, ligand size. */
function buildQcOptionsSection(ctx: StepContext): HTMLElement {
  const o = ctx.options;
  const resolutionInput = ui.input.float('Max resolution (Å)', {value: o.maxResolution});
  resolutionInput.root.title = 'Drop PDB entries with resolution worse than this cutoff. ' +
    'Lower numbers = sharper structures. 2.5 Å is the typical cut for druggable kinase data.';

  // Experimental-method toggles. Three traditional experimental methods are
  // on by default; AI-predicted models are off (RCSB's main entry index
  // doesn't include them anyway).
  const allowXrayInput = ui.input.bool('X-ray', {value: o.allowXray});
  allowXrayInput.root.title = 'Include X-ray crystal structures (the historical gold ' +
    'standard for small-molecule ligand positions).';
  const allowNmrInput = ui.input.bool('NMR', {value: o.allowNmr});
  allowNmrInput.root.title = 'Include NMR ensembles (solution / solid-state). Ligand ' +
    'positions are averaged over the ensemble, so they can be softer than X-ray data.';
  const allowCryoEmInput = ui.input.bool('Cryo-EM', {value: o.allowCryoEm});
  allowCryoEmInput.root.title = 'Include cryo-EM single-particle reconstructions. Modern ' +
    'cryo-EM at 2–3 Å is comparable to X-ray for ligand positioning.';
  // No "AI predicted" toggle: UniProt PDB cross-refs never include AlphaFold
  // models, and predicted apo structures have no bound ligand to contribute to
  // a ligand-based consensus pharmacophore.

  const mwInput = ui.input.float('Min ligand MW (Da)', {value: o.minLigandMw});
  mwInput.root.title = 'Drop ligands lighter than this MW. Filters out water, ions, ' +
    'cryoprotectants (glycerol, PEG, etc.) so only drug-sized small molecules are kept. ' +
    'The cofactor deny-list runs separately.';

  const allInputs = [resolutionInput, allowXrayInput, allowNmrInput,
    allowCryoEmInput, mwInput];
  for (const inp of allInputs) {
    (inp as any).onChanged.subscribe(() => ctx.onOptionsChanged({
      ...ctx.options,
      maxResolution: resolutionInput.value ?? DEFAULT_OPTIONS.maxResolution,
      allowXray: allowXrayInput.value ?? DEFAULT_OPTIONS.allowXray,
      allowNmr: allowNmrInput.value ?? DEFAULT_OPTIONS.allowNmr,
      allowCryoEm: allowCryoEmInput.value ?? DEFAULT_OPTIONS.allowCryoEm,
      minLigandMw: mwInput.value ?? DEFAULT_OPTIONS.minLigandMw,
    }));
  }

  // Small label above the 4 method toggles so the user knows what they
  // collectively control.
  const methodsLabel = ui.divText('Source');
  methodsLabel.style.fontSize = '11px';
  methodsLabel.style.color = '#666';
  methodsLabel.style.marginTop = '8px';
  methodsLabel.style.marginBottom = '2px';

  // Lay the 4 method toggles out horizontally as an even 4-column grid so
  // each toggle takes exactly one quarter of the row's width regardless of
  // its label length (X-ray=5 chars vs AI predicted=12 chars).
  const methodsRow = ui.div([
    allowXrayInput.root,
    allowNmrInput.root,
    allowCryoEmInput.root,
  ]);
  methodsRow.style.display = 'grid';
  methodsRow.style.gridTemplateColumns = 'repeat(3, 1fr)';
  methodsRow.style.gap = '8px';
  methodsRow.style.width = '100%';
  methodsRow.style.alignItems = 'center';
  // Override the input wrappers so the checkbox + label sit inline within
  // each grid cell, with no margin pushing them off-cell.
  for (const inp of [allowXrayInput, allowNmrInput, allowCryoEmInput]) {
    (inp.root as HTMLElement).style.display = 'inline-flex';
    (inp.root as HTMLElement).style.alignItems = 'center';
    (inp.root as HTMLElement).style.margin = '0';
    (inp.root as HTMLElement).style.minWidth = '0';
  }

  // Resolution + Min ligand MW share one row (both are short numeric inputs).
  const numericRow = ui.divH([resolutionInput.root, mwInput.root], 'cp-filter-numeric-row');
  numericRow.style.gap = '12px';
  resolutionInput.root.style.flex = '1';
  mwInput.root.style.flex = '1';

  return ui.divV([
    numericRow,
    methodsLabel,
    methodsRow,
  ]);
}

// ---------------------------------------------------------------------------
// Step 1 — Results summary (after Stage 1 runs)
// ---------------------------------------------------------------------------

/** Classify the RCSB `experimental_method` string into a short label
 *  ("X-ray" / "NMR" / "Cryo-EM" / "AI" / "Other"). Mirrors the substring
 *  rules in rcsb-client.ts so users see a consistent label across the UI. */
function shortMethod(method: string): string {
  const m = (method ?? '').toLowerCase();
  if (m.includes('x-ray')) return 'X-ray';
  if (m.includes('nmr')) return 'NMR';
  if (m.includes('electron microscopy') || m.includes('cryo-em')) return 'Cryo-EM';
  if (m.includes('alphafold') || m.includes('predicted') ||
      m.includes('computational') || m.includes('computed')) return 'AI';
  return method || '–';
}

/** Build a compact "what got accepted by Stage 1" table for the Step 1
 *  panel. One row per accepted PDB with resolution, method, and the
 *  ligands that passed the cofactor + MW filters. If `dropped` is
 *  non-empty, also renders a "Dropped" section underneath with the
 *  per-input reason (resolution > cutoff, method filter, etc.). */
function buildPdbQcSummary(
  df: DG.DataFrame, requestedIds: string[], dropped: DroppedPdb[],
  ctx: StepContext,
): HTMLElement {
  // Group pdb_qc rows by pdb_id. Each row in pdb_qc is one (PDB, ligand)
  // pair, so per-PDB info repeats across rows.
  type AccLigand = {comp_id: string; mw: number | null};
  type AccPdb = {pdb_id: string; resolution: number | null;
    method: string; ligands: AccLigand[]};
  const acc = new Map<string, AccPdb>();
  const pdbCol = df.col('pdb_id');
  const resCol = df.col('resolution');
  const methodCol = df.col('experimental_method');
  const ligCol = df.col('ligand_comp_id');
  const mwCol = df.col('ligand_formula_weight');
  if (pdbCol) {
    for (let i = 0; i < df.rowCount; i++) {
      const id = String(pdbCol.get(i) ?? '');
      if (!id) continue;
      let row = acc.get(id);
      if (!row) {
        const r = resCol ? resCol.get(i) : null;
        row = {
          pdb_id: id,
          resolution: r == null || r === '' ? null : Number(r),
          method: methodCol ? String(methodCol.get(i) ?? '') : '',
          ligands: [],
        };
        acc.set(id, row);
      }
      const lig = ligCol ? String(ligCol.get(i) ?? '').trim() : '';
      if (lig && !row.ligands.find((l) => l.comp_id === lig)) {
        const mwRaw = mwCol ? mwCol.get(i) : null;
        const mw = mwRaw == null || mwRaw === '' || !Number.isFinite(Number(mwRaw)) ?
          null : Number(mwRaw);
        row.ligands.push({comp_id: lig, mw});
      }
    }
  }

  const accepted = Array.from(acc.values()).sort((a, b) =>
    a.pdb_id < b.pdb_id ? -1 : (a.pdb_id > b.pdb_id ? 1 : 0));
  const acceptedIds = new Set(accepted.map((r) => r.pdb_id.toUpperCase()));
  const droppedIds = requestedIds
    .map((id) => id.toUpperCase())
    .filter((id) => !acceptedIds.has(id));

  // Build the summary table. New leading "Use" column: un-tick to drop a
  // PDB from the pipeline (Stage 2 onward) AND hide it in Mol*. Row click
  // (not on the checkbox) isolates that PDB in the viewer — mirrors Step 4.
  const excludedSet = new Set(
    (ctx.options.excludedInputPdbs ?? []).map((s) => s.toUpperCase()));
  const isolatedUpper = (ctx.isolatedPdb ?? '').toUpperCase();
  const table = document.createElement('table');
  table.className = 'cp-step-table';
  const thead = document.createElement('thead');
  const headerRow = document.createElement('tr');
  for (const h of ['Use', 'PDB', 'Resolution', 'Method', 'Ligand(s)']) {
    const th = document.createElement('th');
    th.textContent = h;
    headerRow.append(th);
  }
  thead.append(headerRow);
  table.append(thead);
  const tbody = document.createElement('tbody');
  for (const row of accepted) {
    const idUp = row.pdb_id.toUpperCase();
    const tr = document.createElement('tr');
    // Use the same hover/isolated highlight classes as Step 2 / Step 4 rows.
    tr.classList.add('cp-rmsd-row');
    if (isolatedUpper && idUp === isolatedUpper)
      tr.classList.add('cp-rmsd-row-isolated');
    tr.style.cursor = 'pointer';
    tr.title = isolatedUpper === idUp
      ? `Click again to show all PDBs.`
      : `Click to isolate ${row.pdb_id} (hides other proteins).`;

    // Use checkbox — drives excludedInputPdbs + Mol* visibility. stopPropagation
    // so ticking the box doesn't also trigger row-isolate.
    const tdUse = document.createElement('td');
    const cb = document.createElement('input');
    cb.type = 'checkbox';
    cb.checked = !excludedSet.has(idUp);
    cb.title = 'Include this PDB in alignment + downstream pipeline (and show it in 3D).';
    cb.style.cursor = 'pointer';
    cb.addEventListener('click', (e) => e.stopPropagation());
    cb.addEventListener('change', () => {
      const next = new Set(
        (ctx.options.excludedInputPdbs ?? []).map((s) => s.toUpperCase()));
      if (cb.checked) next.delete(idUp); else next.add(idUp);
      ctx.onOptionsChanged({...ctx.options, excludedInputPdbs: [...next]});
      // Hide / show this PDB in Mol* immediately (the persistent viewer carries
      // the aligned structures from the prior step).
      ctx.onTogglePdbVisibility?.(idUp, !cb.checked);
      ctx.refreshStep?.();
    });
    tdUse.append(cb);

    const tdId = document.createElement('td');
    tdId.textContent = row.pdb_id;
    tdId.style.fontFamily = 'monospace';
    tdId.style.fontWeight = '600';
    const tdRes = document.createElement('td');
    tdRes.textContent = row.resolution == null ? '–' :
      `${row.resolution.toFixed(2)} Å`;
    const tdMethod = document.createElement('td');
    tdMethod.textContent = shortMethod(row.method);
    tdMethod.title = row.method || '';
    const tdLig = document.createElement('td');
    if (row.ligands.length === 0) {
      tdLig.textContent = '–';
    } else {
      // Render as "FMM 550 Da, W2R 412 Da" — comp_id monospace, MW small grey.
      row.ligands.forEach((lig, i) => {
        if (i > 0) {
          const sep = document.createElement('span');
          sep.textContent = ', ';
          tdLig.append(sep);
        }
        const compEl = document.createElement('span');
        compEl.textContent = lig.comp_id;
        compEl.style.fontFamily = 'monospace';
        compEl.style.fontSize = '10px';
        tdLig.append(compEl);
        if (lig.mw != null) {
          const mwEl = document.createElement('span');
          mwEl.textContent = ` ${lig.mw.toFixed(0)} Da`;
          mwEl.style.color = '#888';
          mwEl.style.fontSize = '10px';
          tdLig.append(mwEl);
        }
      });
    }
    tr.append(tdUse, tdId, tdRes, tdMethod, tdLig);
    // Row click → isolate / release (Mol* show-only-it). Re-clicking an
    // already-isolated row releases isolation. stopPropagation on the checkbox
    // above keeps the two interactions independent.
    tr.addEventListener('click', () => {
      if (!ctx.onIsolatePdb) return;
      ctx.onIsolatePdb(isolatedUpper === idUp ? null : idUp);
    });
    tbody.append(tr);
  }
  table.append(tbody);

  // Wrapper with a header line + the outcome counts + (optional) dropped list.
  const els: HTMLElement[] = [];
  const heading = ui.h3('Accepted PDBs');
  els.push(heading);

  const outcome = ui.divText(
    `${accepted.length} of ${requestedIds.length} input(s) passed the filters.`);
  outcome.style.fontSize = '11px';
  outcome.style.color = accepted.length === requestedIds.length ? '#388E3C' : '#888';
  outcome.style.marginBottom = '6px';
  els.push(outcome);

  els.push(table);

  // Dropped section — only render if Stage 1 actually produced drop records.
  // (`droppedIds` from set-diff would catch RCSB returning fewer rows than
  // input, but `dropped[]` carries the human-readable per-PDB reason.)
  if (dropped.length > 0) {
    els.push(buildDroppedPdbSection(dropped));
  } else if (droppedIds.length > 0) {
    // Fallback when Stage 1 wasn't the source of truth (e.g. legacy cache).
    const drop = ui.divText(
      `Dropped: ${droppedIds.join(', ')} (didn't pass filters).`);
    drop.style.fontSize = '11px';
    drop.style.color = '#B26A00';
    drop.style.marginTop = '4px';
    els.push(drop);
  }

  return ui.divV(els);
}

/** Build the per-PDB "Dropped" table for Step 1. Each row shows the PDB
 *  ID, the reason it was filtered out, and (when known) its resolution
 *  and experimental method as reported by RCSB. */
function buildDroppedPdbSection(dropped: DroppedPdb[]): HTMLElement {
  const heading = ui.h3('Dropped PDBs');
  const blurb = ui.divText(
    `${dropped.length} input(s) didn't pass the filters. Adjust filter ` +
    'settings above and re-run to include them.');
  blurb.style.fontSize = '11px';
  blurb.style.color = '#B26A00';
  blurb.style.marginBottom = '6px';

  const table = document.createElement('table');
  table.className = 'cp-step-table';
  const thead = document.createElement('thead');
  const headerRow = document.createElement('tr');
  for (const h of ['PDB', 'Resolution', 'Method', 'Reason']) {
    const th = document.createElement('th');
    th.textContent = h;
    headerRow.append(th);
  }
  thead.append(headerRow);
  table.append(thead);
  const tbody = document.createElement('tbody');
  for (const d of dropped) {
    const tr = document.createElement('tr');
    const tdId = document.createElement('td');
    tdId.textContent = d.pdb_id;
    tdId.style.fontFamily = 'monospace';
    tdId.style.fontWeight = '600';
    tdId.style.color = '#B26A00';
    const tdRes = document.createElement('td');
    tdRes.textContent = d.resolution == null ? '–' : `${d.resolution.toFixed(2)} Å`;
    const tdMethod = document.createElement('td');
    tdMethod.textContent = d.method ? shortMethod(d.method) : '–';
    tdMethod.title = d.method || '';
    const tdReason = document.createElement('td');
    tdReason.textContent = d.reason;
    tdReason.style.fontSize = '10px';
    tdReason.style.color = '#666';
    tr.append(tdId, tdRes, tdMethod, tdReason);
    tbody.append(tr);
  }
  table.append(tbody);

  return ui.divV([heading, blurb, table]);
}

/** Stage 2a (Align) options. Currently just the reference-PDB override —
 *  if the user doesn't pick one, the Python script auto-selects the
 *  highest-resolution structure. Same field (`refPdbId`) is also used by
 *  Stage 5b for the final consensus rendering frame. */
function buildAlignOptionsSection(ctx: StepContext): HTMLElement {
  // Build the choice list from the ACCEPTED PDBs (Stage 1 QC output), NOT the
  // raw textarea input — so structures dropped by the filters (e.g. 1M14 over
  // the resolution cap) don't appear as pickable references. Each option shows
  // the PDB's resolution so the user can choose a reference informedly (Auto
  // picks the highest-resolution one). Falls back to the raw input list (no
  // resolution) if QC data isn't available yet.
  const AUTO_LABEL = 'Auto (highest resolution)';
  const acceptedList: Array<{id: string; res: number | null}> = [];
  if (ctx.pdbQcData) {
    const qPid = ctx.pdbQcData.col('pdb_id');
    const qRes = ctx.pdbQcData.col('resolution');
    if (qPid) {
      for (let i = 0; i < ctx.pdbQcData.rowCount; i++) {
        const id = String(qPid.get(i) ?? '').trim();
        if (!id) continue;
        const r = qRes ? Number(qRes.get(i)) : NaN;
        acceptedList.push({id, res: Number.isFinite(r) ? r : null});
      }
    }
  }
  const refList = acceptedList.length > 0
    ? acceptedList : ctx.pdbIds.map((id) => ({id, res: null as number | null}));
  // Each dropdown option's DISPLAY label encodes the resolution, but the VALUE
  // we act on is the bare PDB id — keep maps in both directions.
  const labelFor = (e: {id: string; res: number | null}): string =>
    e.res == null ? e.id : `${e.id} · ${e.res.toFixed(2)} Å`;
  const labelToId = new Map<string, string>();
  const idToLabel = new Map<string, string>();
  for (const e of refList) {
    const lbl = labelFor(e);
    labelToId.set(lbl, e.id);
    idToLabel.set(e.id.toUpperCase(), lbl);
  }
  const items = [AUTO_LABEL, ...refList.map(labelFor)];
  const current = (ctx.options.refPdbId ?? '').trim();
  const initial = (current && idToLabel.has(current.toUpperCase()))
    ? idToLabel.get(current.toUpperCase())! : AUTO_LABEL;

  // Resolved reference of the displayed alignment (★ marker). When the user
  // picks this PDB explicitly, the Stage 2a output is identical to "Auto"
  // — so we normalize to undefined to keep the fingerprint stable and
  // avoid spurious stale state.
  const resolvedRefUpper = (() => {
    const af = ctx.alignedData;
    const col = af?.col('ref_pdb_id');
    if (!af || !col || af.rowCount === 0) return '';
    return String(col.get(0) ?? '').toUpperCase();
  })();

  const applyRefChoice = (pickedLabel: string | null | undefined): void => {
    const lbl = String(pickedLabel ?? AUTO_LABEL);
    // Map the display label ("1XKK · 2.40 Å") back to the bare PDB id. Fallback
    // to the label itself for the raw-input path where label === id.
    const pickedId = lbl === AUTO_LABEL ? undefined : (labelToId.get(lbl) ?? lbl);
    let refPdbId: string | undefined = pickedId;
    if (refPdbId && resolvedRefUpper && refPdbId.toUpperCase() === resolvedRefUpper) {
      // Picked PDB equals the resolved reference → semantic Auto.
      refPdbId = undefined;
    }
    ctx.onOptionsChanged({...ctx.options, refPdbId});
    // Also isolate the picked PDB in Mol* so the user can see the single
    // protein they're about to align to. Auto / re-click → null → show all.
    const isolateTarget = (lbl === AUTO_LABEL) ? null : (pickedId ?? null);
    ctx.onIsolatePdb?.(isolateTarget);
    ctx.refreshStep?.();
  };

  const refInput = ui.input.choice('Reference', {
    value: initial,
    items,
    // `onValueChanged` is the documented Datagrok API and fires reliably for
    // choice inputs. The legacy `.onChanged.subscribe()` path was unreliable
    // here — empirically it only fired when both 'input' and 'change' events
    // were dispatched on the underlying <select>, but native user clicks on a
    // <select> only fire 'change'. The user's dropdown pick was therefore
    // silently swallowed, the Stage 2a cache wasn't invalidated, and the
    // Run-step click re-rendered the previous alignment.
    onValueChanged: (v: string) => applyRefChoice(v),
  });
  refInput.root.title = 'PDB used as the alignment template. Every other ' +
    'structure is rotated + translated onto this one, so its RMSD is 0 by ' +
    'construction. Auto = use the highest-resolution entry.';

  // Belt-and-braces: also listen to the raw <select>'s 'change' event so
  // that even if Datagrok's InputBase wrapper drifts in a future version we
  // still capture the user's pick. Idempotent — applyRefChoice is safe to
  // call multiple times with the same value (options are compared by the
  // fingerprint, not identity).
  const selectEl = refInput.root.querySelector('select') as HTMLSelectElement | null;
  if (selectEl) selectEl.addEventListener('change', () => applyRefChoice(selectEl.value));

  return ui.divV([refInput.root]);
}

/** Stage 3 (Pocket) options: method, radius, representation. */
function buildPocketOptionsSection(ctx: StepContext): HTMLElement {
  const o = ctx.options;

  // Builds the next options snapshot from the live input values — used by
  // both the InputBase onValueChanged callbacks and the backup raw-select
  // change listeners (see buildAlignOptionsSection for the rationale).
  const emit = (): void => ctx.onOptionsChanged({
    ...ctx.options,
    pocketMethod: (pocketMethodInput.value as PocketMethod) ?? DEFAULT_OPTIONS.pocketMethod,
    pocketRadius: pocketRadiusInput.value ?? DEFAULT_OPTIONS.pocketRadius,
    pocketRep: (pocketRepInput.value as PocketRep) ?? DEFAULT_OPTIONS.pocketRep,
  });

  const pocketMethodInput = ui.input.choice<PocketMethod>('Pocket method', {
    value: o.pocketMethod,
    items: ['cutoff', 'dbscan'] as PocketMethod[],
    onValueChanged: () => emit(),
  });
  pocketMethodInput.root.title = '"cutoff" keeps protein atoms within pocket-radius of any drug ligand. ' +
    '"dbscan" clusters the union point cloud and keeps the largest non-noise cluster.';
  const pocketRadiusInput = ui.input.float('Pocket radius (Å)', {
    value: o.pocketRadius,
    onValueChanged: () => emit(),
  });
  pocketRadiusInput.root.title = 'Distance cutoff in Å. Cutoff path only.';
  const pocketRepInput = ui.input.choice<PocketRep>('Pocket rep', {
    value: o.pocketRep,
    items: ['spacefill', 'gaussian-surface'] as PocketRep[],
    onValueChanged: () => emit(),
  });
  pocketRepInput.root.title = 'Mol* binding-site rendering style.';

  // Backup listeners on the raw <select>s in case the InputBase wrapper's
  // onValueChanged drifts in a future version (see buildAlignOptionsSection).
  for (const inp of [pocketMethodInput, pocketRepInput]) {
    const sel = inp.root.querySelector('select') as HTMLSelectElement | null;
    if (sel) sel.addEventListener('change', () => emit());
  }

  return ui.divV([pocketMethodInput.root, pocketRadiusInput.root, pocketRepInput.root]);
}

/** Stage 5a (Consensus) options + output reference PDB. */
function buildConsensusOptionsSection(ctx: StepContext): HTMLElement {
  const o = ctx.options;
  const kqInput = ui.input.int('Avg features per cluster (kq)',
    {value: o.kq, min: 1, max: 30});
  kqInput.root.title = 'k-means n_clusters = ceil(n_features / kq). Lower = more, smaller clusters; ' +
    'higher = fewer, fatter clusters. Default 7 matches TeachOpenCADD T009.';
  const minFracInput = ui.input.float('Min ligand fraction',
    {value: o.minClusterSizeFraction, min: 0.0, max: 1.0});
  minFracInput.root.title = 'A cluster is kept only if at least this fraction of distinct ligands ' +
    'contributed a feature to it. Drop to 0.5 to surface less-conserved hotspots.';
  const topClusterInput = ui.input.int('Max clusters per family',
    {value: o.topClusterNumber, min: 1, max: 10});
  topClusterInput.root.title = 'After filtering by min ligand fraction, keep at most N clusters per ' +
    'family (sorted by number of contributing ligands, then by cluster size).';
  // Reference PDB now lives on Step 2 (Align). It's used both there (as
  // alignment template) AND by Stage 5b (rendering frame), so we don't
  // duplicate the input on Step 5 — just keep ctx.options.refPdbId as the
  // single source of truth.
  for (const inp of [kqInput, minFracInput, topClusterInput]) {
    (inp as any).onChanged.subscribe(() => ctx.onOptionsChanged({
      ...ctx.options,
      kq: kqInput.value ?? DEFAULT_OPTIONS.kq,
      minClusterSizeFraction: minFracInput.value ?? DEFAULT_OPTIONS.minClusterSizeFraction,
      topClusterNumber: topClusterInput.value ?? DEFAULT_OPTIONS.topClusterNumber,
    }));
  }
  return ui.divV([kqInput.root, minFracInput.root, topClusterInput.root]);
}

// ---------------------------------------------------------------------------
// Step 1 — Fetch PDBs (input form)
// ---------------------------------------------------------------------------

/** Center a checkbox in its table cell. */
function wrapCheckbox(cb: HTMLInputElement): HTMLElement {
  const d = ui.div([cb]);
  d.style.textAlign = 'center';
  return d;
}

/** A text cell with an optional hover tooltip (for long titles / ligand names). */
function cellWithTitle(text: string, title: string): HTMLElement {
  const el = ui.divText(text);
  if (title) el.title = title;
  return el;
}

/** Drop stereo markers (@, @@, cis/trans /\) so the depiction renders flat —
 *  no wedges, R/S, "abs", "this enantiomer" or "unknown chirality" clutter. */
function stripStereo(smiles: string): string {
  return smiles.replace(/@+/g, '').replace(/[/\\]/g, '');
}

/** Render the bound ligand(s) as small flat 2D structures with comp-id
 *  captions, so the user can eyeball them against the query molecule. Falls
 *  back to the comp-id text when a SMILES is missing or won't render. Capped
 *  at 3. */
function ligandCell(ligands: {compId: string; name: string; smiles: string}[]): HTMLElement {
  if (!ligands.length) return ui.divText('—');
  const items = ligands.slice(0, 3).map((l) => {
    let pic: HTMLElement;
    if (l.smiles) {
      try {
        pic = DG.chem.svgMol(stripStereo(l.smiles), 96, 70, {
          suppressChiralText: true,
          suppressCIPParity: true,
          suppressESR: true,
        });
      } catch { pic = ui.divText(l.compId); }
    } else {
      pic = ui.divText(l.compId);
    }
    const cap = ui.divText(l.compId, 'cp-ligand-cap');
    cap.title = l.name || l.compId;
    return ui.divV([pic, cap], 'cp-ligand-struct');
  });
  const row = ui.divH(items, 'cp-ligand-struct-row');
  if (ligands.length > 3)
    row.appendChild(ui.divText(`+${ligands.length - 3}`, 'cp-ligand-cap'));
  return row;
}

/** A small RCSB structure-snapshot thumbnail for a PDB id (the same kind of
 *  preview image Datagrok's PDB cell renderer shows). Hidden if the image
 *  fails to load (e.g. an entry with no precomputed assembly image). */
function pdbSnapshot(pdbId: string): HTMLElement {
  const img = document.createElement('img');
  img.className = 'cp-pdb-snapshot';
  img.loading = 'lazy';
  img.alt = pdbId;
  img.title = pdbId;
  img.src = `https://cdn.rcsb.org/images/structures/${pdbId.toLowerCase()}_assembly-1.jpeg`;
  img.addEventListener('error', () => { img.style.visibility = 'hidden'; });
  return img;
}

/** Whether a candidate passes the current resolution / experimental-method
 *  filters (the same gates Stage 1 QC applies). Min-ligand-MW is ligand-
 *  specific and stays a QC-time filter. */
function passesFilters(s: PdbSummary, o: PipelineOptions): boolean {
  // No resolution reported (e.g. NMR) — don't pre-select it.
  if (s.resolution == null) return false;
  if (s.resolution > o.maxResolution) return false;
  const m = (s.method || '').toLowerCase();
  if (m.includes('x-ray')) return o.allowXray;
  if (m.includes('nmr')) return o.allowNmr;
  if (m.includes('electron microscopy') || m.includes('cryo-em')) return o.allowCryoEm;
  if (m.includes('alphafold') || m.includes('predicted') ||
      m.includes('computational') || m.includes('computed')) return o.allowAlphaFold;
  return true; // unknown / missing method — keep it
}

/** A "Filters" block (reuses the Stage-1 QC option inputs) for embedding in the
 *  search dialogs. */
function buildFiltersBlock(ctx: StepContext): HTMLElement {
  const h = ui.divText('Filters', 'cp-search-filters-title');
  h.style.fontWeight = '600';
  h.style.marginTop = '4px';
  return ui.divV([h, buildQcOptionsSection(ctx)], 'cp-search-filters');
}

/** Build a results table (snapshot · PDB · title · res · method · ligand ·
 *  similarity, with a checkbox) for `summaries`. Registers each checkbox into
 *  the shared `checks` map (does NOT clear it) and sets the initial tick state
 *  to `allChecked`. Returns the table container; its own "Select all" toggles
 *  only this table's rows. */
function buildResultsTable(
  summaries: PdbSummary[],
  scoreById: Map<string, number>,
  checks: Map<string, HTMLInputElement>,
  onSelectionChange: (() => void) | undefined,
  allChecked: boolean,
): HTMLElement {
  const showSim = scoreById.size > 0;
  const local: HTMLInputElement[] = [];

  const selectAll = document.createElement('input');
  selectAll.type = 'checkbox';
  selectAll.checked = allChecked;
  selectAll.addEventListener('change', () => {
    for (const c of local) c.checked = selectAll.checked;
    onSelectionChange?.();
  });
  const selectAllRow = ui.divH([wrapCheckbox(selectAll), ui.divText('Select all')]);
  selectAllRow.style.alignItems = 'center';
  selectAllRow.style.gap = '6px';
  selectAllRow.style.margin = '4px 0';

  const table = document.createElement('table');
  table.className = 'cp-step-table cp-ligand-picker-table';
  const thead = document.createElement('thead');
  const htr = document.createElement('tr');
  const headers = ['', '', 'PDB', 'Title', 'Res (Å)', 'Method', 'Ligand(s)'];
  if (showSim) headers.push('Similarity');
  for (const h of headers) {
    const th = document.createElement('th');
    th.textContent = h;
    htr.appendChild(th);
  }
  thead.appendChild(htr);
  table.appendChild(thead);

  const tbody = document.createElement('tbody');
  for (const s of summaries) {
    const tr = document.createElement('tr');
    const cb = document.createElement('input');
    cb.type = 'checkbox';
    cb.checked = allChecked;
    cb.addEventListener('change', () => onSelectionChange?.());
    checks.set(s.id, cb);
    local.push(cb);

    const score = scoreById.get(s.id);
    const cells: HTMLElement[] = [
      wrapCheckbox(cb),
      pdbSnapshot(s.id),
      ui.divText(s.id),
      cellWithTitle(s.title || '—', s.title),
      ui.divText(s.resolution == null ? '—' : s.resolution.toFixed(2)),
      ui.divText(s.method || '—'),
      ligandCell(s.ligands),
    ];
    if (showSim) cells.push(ui.divText(score == null ? '—' : score.toFixed(2)));
    for (const c of cells) {
      const td = document.createElement('td');
      td.appendChild(c);
      tr.appendChild(td);
    }
    // Clicking the row toggles its checkbox (except the box and the snapshot).
    tr.addEventListener('click', (e) => {
      const t = e.target as HTMLElement;
      if (t !== cb && t.tagName !== 'IMG') {
        cb.checked = !cb.checked;
        onSelectionChange?.();
      }
    });
    tbody.appendChild(tr);
  }
  table.appendChild(tbody);

  return ui.divV([selectAllRow, ui.div([table], 'cp-ligand-picker-scroll')]);
}

/** Render results into `host`: structures that pass the filters and those that
 *  don't are split into "Passed" / "Excluded" tabs (passed pre-ticked, excluded
 *  un-ticked) so it's easy to see and adjust. Repopulates the shared `checks`. */
function renderLigandResults(
  host: HTMLElement,
  summaries: PdbSummary[],
  scoreById: Map<string, number>,
  checks: Map<string, HTMLInputElement>,
  onSelectionChange?: () => void,
  passes?: (s: PdbSummary) => boolean,
): void {
  checks.clear();
  ui.empty(host);
  if (!summaries.length) {
    host.appendChild(ui.divText('No matching PDB structures.', 'cp-ligand-empty'));
    onSelectionChange?.();
    return;
  }
  const pass = passes ?? ((): boolean => true);
  const passed = summaries.filter(pass);
  const failed = summaries.filter((s) => !pass(s));

  if (passed.length && failed.length) {
    // Lightweight tab toggle (full-width) instead of ui.tabControl, which
    // constrains the pane width and clipped the table.
    const passedEl = buildResultsTable(passed, scoreById, checks, onSelectionChange, true);
    const failedEl = buildResultsTable(failed, scoreById, checks, onSelectionChange, false);
    failedEl.style.display = 'none';
    const passedTab = ui.divText(`Passed (${passed.length})`,
      'cp-result-tab cp-result-tab-active');
    const excludedTab = ui.divText(`Excluded (${failed.length})`, 'cp-result-tab');
    const activate = (showPassed: boolean): void => {
      passedEl.style.display = showPassed ? '' : 'none';
      failedEl.style.display = showPassed ? 'none' : '';
      passedTab.classList.toggle('cp-result-tab-active', showPassed);
      excludedTab.classList.toggle('cp-result-tab-active', !showPassed);
    };
    passedTab.onclick = (): void => activate(true);
    excludedTab.onclick = (): void => activate(false);
    host.appendChild(ui.divV([
      ui.divH([passedTab, excludedTab], 'cp-result-tabbar'),
      passedEl,
      failedEl,
    ]));
  } else {
    host.appendChild(buildResultsTable(
      passed.length ? passed : failed, scoreById, checks, onSelectionChange,
      passed.length > 0));
  }
  onSelectionChange?.();
}

/**
 * Ligand-search dialog: a Chem sketcher to draw/paste the query molecule, a
 * Search button that queries RCSB for PDBs with a similar bound ligand, and a
 * results table (structure snapshot, metadata, similarity, checkbox). On OK,
 * `onConfirm` receives the ticked PDB ids.
 */
/** Organism-chooser label: always shows the PDB count in brackets (incl. 0). */
function uniProtLabel(h: UniProtHit): string {
  return `${h.organism || '?'} · ${h.geneName || '(no gene)'} · ` +
    `${h.accession} (${h.pdbCount} PDB)`;
}

/**
 * Ligand-search dialog: a Chem sketcher to draw/paste the query molecule, a
 * Search button that queries RCSB for PDBs with a similar bound ligand, and a
 * results table (snapshot / metadata / similarity / checkbox). On OK,
 * `onConfirm(ids)` runs with the ticked PDB ids.
 */
function openLigandSearchDialog(
  onConfirm: (ids: string[]) => void, ctx: StepContext,
): void {
  const molInput = ui.input.string('Molecule', {value: '', nullable: true});
  const molInputEl = molInput.root.querySelector('input') as HTMLInputElement | null;
  molInputEl?.setAttribute('placeholder', 'SMILES — or click Sketch to draw');
  if (molInputEl) molInputEl.style.minWidth = '280px';
  molInput.root.style.flex = '1';

  const status = ui.divText('Paste a SMILES or sketch a molecule, then Search.',
    'cp-ligand-search-status');
  status.style.minHeight = '16px';
  status.style.color = '#777';
  status.style.fontSize = '12px';

  const resultsHost = ui.div([], 'cp-ligand-results-host');
  const checks = new Map<string, HTMLInputElement>();
  let summaries: PdbSummary[] = [];
  const someChecked = (): boolean => [...checks.values()].some((c) => c.checked);
  const selInfo = ui.divText('', 'cp-sel-count');
  let okBtn: HTMLButtonElement | null = null;
  const refreshOk = (): void => {
    const n = [...checks.values()].filter((c) => c.checked).length;
    if (okBtn) okBtn.disabled = n === 0;
    selInfo.textContent = n ? `${n} structure(s) selected` : '';
  };

  // Sketching opens the Chem sketcher in its own (correctly-sized) dialog and
  // writes the result back into the SMILES field — embedding the sketcher
  // inline fought the layout, so we keep the search dialog clean.
  const sketchBtn = ui.button('Sketch', () => {
    const sk = new DG.chem.Sketcher();
    const cur = (molInput.value ?? '').trim();
    if (cur) sk.setSmiles(cur);
    ui.dialog('Sketch molecule')
      .add(sk.root)
      .onOK(() => { const s = sk.getSmiles(); if (s) molInput.value = s; })
      .show();
  }, 'Draw the query molecule with the Chem sketcher.');

  const lSearch = ui.button('Search', async () => {
    const smiles = (molInput.value ?? '').trim();
    if (!smiles) { grok.shell.warning('Enter a SMILES or sketch a molecule first.'); return; }
    lSearch.disabled = true;
    status.textContent = 'Searching RCSB for similar bound ligands...';
    ui.empty(resultsHost);
    resultsHost.appendChild(ui.loader());
    try {
      const hits = await searchPdbsByLigand(smiles);
      if (!hits.length) {
        summaries = [];
        renderLigandResults(resultsHost, [], new Map(), checks, refreshOk);
        status.textContent = 'No PDBs with a similar bound ligand.';
        return;
      }
      status.textContent = `Found ${hits.length} — loading details...`;
      summaries = await fetchPdbSummaries(hits.map((h) => h.id), smiles);
      const scoreById = new Map(hits.map((h) => [h.id, h.score]));
      const rank = new Map(hits.map((h, i) => [h.id, i]));
      summaries.sort((a, b) => (rank.get(a.id) ?? 0) - (rank.get(b.id) ?? 0));
      const pass = (s: PdbSummary): boolean => passesFilters(s, ctx.options);
      renderLigandResults(resultsHost, summaries, scoreById, checks, refreshOk, pass);
      const nPass = summaries.filter(pass).length;
      status.textContent =
        `${summaries.length} candidate(s); ${nPass} pass the filters (pre-ticked). ` +
        'Adjust ticks, then OK.';
    } catch (e: any) {
      const msg = e?.message ?? String(e);
      status.textContent = `Search failed: ${msg.slice(0, 120)}`;
      grok.shell.error(`Ligand search failed: ${msg}`);
      ui.empty(resultsHost);
    } finally {
      lSearch.disabled = false;
    }
  }, 'Search RCSB for PDB structures whose bound ligand is similar to this molecule.');
  lSearch.classList.add('cp-demo-hero-btn');

  molInputEl?.addEventListener('keydown', (e) => {
    if (e.key === 'Enter') {
      e.preventDefault();
      e.stopPropagation();
      if (!lSearch.disabled) lSearch.click();
    }
  });

  const inputRow = ui.divH([molInput.root, sketchBtn, lSearch], 'cp-target-lookup-row');
  inputRow.style.alignItems = 'flex-end';
  inputRow.style.gap = '6px';

  const body = ui.divV([inputRow, buildFiltersBlock(ctx), status, selInfo, resultsHost],
    'cp-search-dialog-body');
  body.style.width = '880px';

  const dlg = ui.dialog('Search by ligand');
  dlg.add(body);
  dlg.onOK(() => {
    if (!someChecked()) return;
    onConfirm(summaries.map((s) => s.id).filter((id) => checks.get(id)?.checked));
  });
  dlg.show();
  okBtn = dlg.getButton('OK');
  refreshOk(); // disabled until a search produces (checked) results
}

/**
 * Target-search dialog: type a gene name, protein name, or UniProt accession,
 * search, pick the species inline (when a gene matches several orthologs), and
 * the cross-referenced PDBs. On OK, `onConfirm(ids, accession)` runs.
 */
function openTargetSearchDialog(
  onConfirm: (ids: string[], accession: string) => void, ctx: StepContext,
): void {
  const input = ui.input.string('Target', {value: '', nullable: true});
  const inputEl = input.root.querySelector('input') as HTMLInputElement | null;
  inputEl?.setAttribute('placeholder', 'UniProt ID or name, e.g. P00533, EGFR');
  if (inputEl) inputEl.style.minWidth = '280px';
  input.root.style.flex = '1';

  const status = ui.divText('Enter a UniProt accession or a gene/protein name, then Search.',
    'cp-ligand-search-status');
  status.style.minHeight = '16px';
  status.style.color = '#777';
  status.style.fontSize = '12px';

  const organismHost = ui.div([], 'cp-target-organism-host');
  organismHost.style.display = 'none';

  const resultsHost = ui.div([], 'cp-ligand-results-host');
  const checks = new Map<string, HTMLInputElement>();
  let summaries: PdbSummary[] = [];
  let accession = '';
  const someChecked = (): boolean => [...checks.values()].some((c) => c.checked);
  const selInfo = ui.divText('', 'cp-sel-count');
  let okBtn: HTMLButtonElement | null = null;
  const refreshOk = (): void => {
    const n = [...checks.values()].filter((c) => c.checked).length;
    if (okBtn) okBtn.disabled = n === 0;
    selInfo.textContent = n ? `${n} structure(s) selected` : '';
  };

  const loadAndRender = async (hit: UniProtHit): Promise<void> => {
    accession = hit.accession;
    status.textContent = `Fetching PDB cross-references for ${hit.accession}...`;
    ui.empty(resultsHost);
    resultsHost.appendChild(ui.loader());
    try {
      const pdbIds = await getPdbIdsForAccession(hit.accession);
      if (!pdbIds.length) {
        summaries = [];
        renderLigandResults(resultsHost, [], new Map(), checks, refreshOk);
        status.textContent = `${hit.accession} (${hit.organism}) has no PDB cross-references.`;
        return;
      }
      const MAX = 50;
      const capped = pdbIds.slice(0, MAX);
      status.textContent = `${hit.accession} · ${pdbIds.length} PDBs — loading details...`;
      summaries = await fetchPdbSummaries(capped);
      const pass = (s: PdbSummary): boolean => passesFilters(s, ctx.options);
      renderLigandResults(resultsHost, summaries, new Map(), checks, refreshOk, pass);
      const note = pdbIds.length > MAX ? ` (first ${MAX} of ${pdbIds.length})` : '';
      const nPass = summaries.filter(pass).length;
      status.textContent =
        `${hit.accession} (${hit.organism}) — ${summaries.length} structure(s)${note}; ` +
        `${nPass} pass the filters (pre-ticked). Adjust ticks, then OK.`;
    } catch (e: any) {
      const msg = e?.message ?? String(e);
      status.textContent = `Lookup failed: ${msg.slice(0, 120)}`;
      grok.shell.error(`Target lookup failed: ${msg}`);
      ui.empty(resultsHost);
    }
  };

  const showOrganismChooser = (hits: UniProtHit[]): void => {
    ui.empty(organismHost);
    const labels = hits.map(uniProtLabel);
    const byLabel = new Map<string, UniProtHit>(hits.map((h, i) => [labels[i], h]));
    const choice = ui.input.choice<string>('Organism', {value: labels[0], items: labels});
    const tip = ui.divText('', 'cp-target-picker-tooltip');
    tip.style.fontSize = '11px';
    tip.style.color = '#555';
    tip.style.marginTop = '4px';
    const update = (): void => {
      tip.textContent = byLabel.get((choice.value ?? '') as string)?.proteinName ?? '';
    };
    update();
    choice.onChanged.subscribe(() => {
      update();
      const h = byLabel.get((choice.value ?? '') as string);
      if (h) void loadAndRender(h);
    });
    organismHost.appendChild(ui.divV([
      ui.divText(`Found ${hits.length} SwissProt entries — species matters ` +
        '(PDB structures differ across orthologs).', 'cp-target-picker-hint'),
      choice.root,
      tip,
    ]));
    organismHost.style.display = '';
  };

  const tSearch = ui.button('Search', async () => {
    const raw = (input.value ?? '').trim();
    if (!raw) { grok.shell.warning('Enter a UniProt accession or a gene/protein name.'); return; }
    // Join a stray space between a gene name and its trailing number so
    // "TLR 8" / "CDK 2" resolve like "TLR8" / "CDK2".
    const query = raw.replace(/([A-Za-z])\s+(\d)/g, '$1$2');
    tSearch.disabled = true;
    ui.empty(organismHost);
    organismHost.style.display = 'none';
    ui.empty(resultsHost);
    status.textContent = 'Searching UniProt...';
    try {
      if (isUniProtAccession(query)) {
        const hit = await fetchUniProtById(query);
        if (!hit) { status.textContent = `UniProt accession ${query} not found.`; return; }
        await loadAndRender(hit);
        return;
      }
      const hits = await searchUniProt(query);
      if (!hits.length) {
        status.textContent = `No SwissProt entries matched "${query}". Try a UniProt accession.`;
        return;
      }
      if (hits.length === 1) { await loadAndRender(hits[0]); return; }
      showOrganismChooser(hits);
      await loadAndRender(hits[0]);
    } catch (e: any) {
      const msg = e?.message ?? String(e);
      status.textContent = `Lookup failed: ${msg.slice(0, 120)}`;
      grok.shell.error(`Target lookup failed: ${msg}`);
    } finally {
      tSearch.disabled = false;
    }
  }, 'Resolve a UniProt accession or gene/protein name to PDB structures.');
  tSearch.classList.add('cp-demo-hero-btn');
  // Enter searches — and must NOT bubble to the dialog (which would treat it as
  // OK and close the dialog).
  inputEl?.addEventListener('keydown', (e) => {
    if (e.key === 'Enter') {
      e.preventDefault();
      e.stopPropagation();
      if (!tSearch.disabled) tSearch.click();
    }
  });

  const inputRow = ui.divH([input.root, tSearch], 'cp-target-lookup-row');
  inputRow.style.alignItems = 'flex-end';
  inputRow.style.gap = '6px';

  const body = ui.divV([inputRow, organismHost, buildFiltersBlock(ctx), status, selInfo,
    resultsHost], 'cp-search-dialog-body');
  body.style.width = '880px';

  const dlg = ui.dialog('Search by target');
  dlg.add(body);
  dlg.onOK(() => {
    if (!someChecked()) return;
    onConfirm(summaries.map((s) => s.id).filter((id) => checks.get(id)?.checked), accession);
  });
  dlg.show();
  okBtn = dlg.getButton('OK');
  refreshOk(); // disabled until a search produces (checked) results
}

export function buildStep1(ctx: StepContext): HTMLElement {
  const initial = (ctx.seedPdbIds ?? ctx.pdbIds ?? []).join('\n');

  // Textarea (empty label — the section header serves as the label).
  const pdbInput = ui.input.textArea('', {
    value: initial,
    size: {width: 380, height: 120},
  });
  (pdbInput.root.querySelector('textarea') as HTMLTextAreaElement | null)?.setAttribute(
    'placeholder', 'Type or paste PDB IDs, one per line (e.g. 3W2S, 4WKQ, 5HG8)');

  // Demo label (visible only if seed PDBs were loaded from the demo).
  const demoLabel = ui.label('Demo: 5 active-state EGFR kinase structures',
    {classes: 'cp-demo-label'});
  demoLabel.style.display = ctx.demoLoaded ? 'block' : 'none';

  // Status line (accepted/rejected counts).
  const pdbStatus = ui.divText('', 'cp-pdb-status');
  pdbStatus.style.fontSize = '11px';
  pdbStatus.style.color = '#777';
  pdbStatus.style.marginTop = '2px';
  pdbStatus.style.minHeight = '14px';

  const refreshPdbStatus = (): void => {
    const {accepted, rejected} = parsePdbIdsDetailed(pdbInput.value ?? '');
    if (accepted.length === 0 && rejected.length === 0) {
      pdbStatus.textContent = '';
      pdbStatus.style.color = '#777';
      return;
    }
    let msg = `${accepted.length} PDB ID${accepted.length === 1 ? '' : 's'} accepted`;
    if (rejected.length > 0) {
      const shown = rejected.slice(0, 3).join(', ');
      const more = rejected.length > 3 ? `, +${rejected.length - 3} more` : '';
      msg += ` · ${rejected.length} ignored (${shown}${more})`;
      pdbStatus.style.color = '#B26A00';
    } else {
      pdbStatus.style.color = '#4CAF50';
    }
    pdbStatus.textContent = msg;
  };

  pdbInput.onChanged.subscribe(() => {
    demoLabel.style.display = 'none';
    refreshPdbStatus();
    ctx.onPdbInputChanged(parsePdbIds(pdbInput.value ?? ''), false);
  });

  // --- Add local .pdb files via the OS file picker --------------------------
  // RCSB-coded files (e.g. 1XKK.pdb) are added as plain IDs; custom files are
  // registered in the local store under a sanitized name and routed around
  // RCSB (see local-pdb-store.ts).
  //
  // A file picker is used rather than OS drag-and-drop: Datagrok installs a
  // window-level file-drop handler (the full-screen "Drop your CSV files to
  // open" overlay) that intercepts dragged files before they reach this
  // textarea, so a package-local drop target can't win that race.
  const appendIds = (ids: string[]): void => {
    if (!ids.length) return;
    const merged = parsePdbIds(pdbInput.value ?? '');
    for (const id of ids) if (!merged.includes(id)) merged.push(id);
    pdbInput.value = merged.join('\n');
    demoLabel.style.display = 'none';
    refreshPdbStatus();
    ctx.onPdbInputChanged(parsePdbIds(pdbInput.value ?? ''), false);
  };

  const handleFiles = async (files: File[]): Promise<void> => {
    const added: string[] = [];
    const skipped: string[] = [];
    for (const file of files) {
      try {
        const text = await file.text();
        if (!looksLikePdb(text)) { skipped.push(file.name); continue; }
        const stem = file.name.replace(/\.[^.]*$/, '');
        added.push(isValidRcsbId(stem) ? stem.toUpperCase() : registerLocalPdb(file.name, text).id);
      } catch {
        skipped.push(file.name);
      }
    }
    appendIds(added);
    if (added.length)
      grok.shell.info(`Added ${added.length} structure(s): ${added.join(', ')}.`);
    if (skipped.length)
      grok.shell.warning(`Ignored ${skipped.length} file(s) with no PDB ` +
        `coordinates: ${skipped.join(', ')}.`);
  };

  const fileInput = document.createElement('input');
  fileInput.type = 'file';
  fileInput.accept = '.pdb,.ent,.pdb1';
  fileInput.multiple = true;
  fileInput.style.display = 'none';
  fileInput.addEventListener('change', () => {
    const files = fileInput.files ? Array.from(fileInput.files) : [];
    if (files.length) void handleFiles(files);
    fileInput.value = ''; // allow re-picking the same file
  });

  const addFileBtn = ui.button('Add .pdb file(s)', () => fileInput.click(),
    'Add one or more .pdb files from your computer. An RCSB-coded file (e.g. ' +
    '1XKK.pdb) adds that ID; a custom structure is registered locally and ' +
    'named after the file.');
  addFileBtn.classList.add('cp-add-pdb-btn');
  const addFileRow = ui.divH([addFileBtn, fileInput], 'cp-add-pdb-row');

  // Demo hero button. The demo set is 6 PDB IDs: five active-state EGFR
  // kinase structures that pass QC + 1M14 (2.6 Å), included deliberately to
  // demonstrate the "Dropped" section (it fails the default 2.5 Å resolution
  // filter). Hence "6 PDBs" — the older "5 PDBs" label undercounted.
  const demoHero = ui.bigButton('Load EGFR demo (6 PDBs)', async () => {
    const ids = await loadDemoPdbIds();
    pdbInput.value = ids.join('\n');
    demoLabel.style.display = 'block';
    refreshPdbStatus();
    ctx.onPdbInputChanged(ids, true);
  });
  demoHero.classList.add('cp-demo-hero-btn');
  demoHero.title = '1XKK / 3W2S / 4WKQ / 5HG5 / 5HG8 — five active-state EGFR ' +
    'kinase structures bound to known inhibitors, plus 1M14 (2.6 Å) which is ' +
    'dropped by the default resolution filter to demo the "Dropped" section.';

  // Search — two buttons (target + ligand) under one "Search structures"
  // heading; each opens its own dialog.
  const searchStatus = ui.divText('', 'cp-target-lookup-status');
  searchStatus.style.fontSize = '11px';
  searchStatus.style.color = '#777';
  searchStatus.style.marginTop = '4px';
  searchStatus.style.minHeight = '14px';
  const note = (added: number, suffix: string): void => {
    searchStatus.textContent = added ?
      `Added ${added} structure(s)${suffix}.` : 'No structures selected.';
    if (added)
      grok.shell.info(`Search: added ${added} PDB ID${added === 1 ? '' : 's'}.`);
  };
  const targetBtn = ui.button('Search by target', () => {
    openTargetSearchDialog((selected, accession) => {
      if (selected.length) appendIds(selected);
      note(selected.length, accession ? ` for ${accession}` : '');
    }, ctx);
  }, 'Find PDB structures by gene name, protein name, or UniProt accession.');
  targetBtn.classList.add('cp-demo-hero-btn');
  const ligandBtn = ui.button('Search by ligand', () => {
    openLigandSearchDialog((selected) => {
      if (selected.length) appendIds(selected);
      note(selected.length, ' from ligand search');
    }, ctx);
  }, 'Draw or paste a ligand and find PDB structures with a similar bound ligand.');
  ligandBtn.classList.add('cp-demo-hero-btn');
  const searchRow = ui.divV([targetBtn, ligandBtn], 'cp-search-btn-row');

  // Stage-1 filters now live inside the search dialogs (target / ligand),
  // where the candidate list is filtered live and the PDBs are picked.

  refreshPdbStatus();

  // After Stage 1 runs, show the per-PDB results summary so the user can
  // see which inputs passed the filters and at what resolution / method.
  // The Dropped section underneath surfaces the ones that didn't make it,
  // each with a per-PDB reason from Stage 1.
  const resultsSection = ctx.pdbQcData ?
    buildPdbQcSummary(ctx.pdbQcData, ctx.pdbIds, ctx.pdbQcDropped, ctx) :
    null;

  return ui.divV([
    panelTitle('1. Fetch PDBs',
      'Search by target or ligand, paste PDB IDs, or add local .pdb files — ' +
      'then click "Run step" in the footer to fetch + QC. New here? Load the ' +
      'EGFR demo at the bottom.'),
    ui.h3('Search structures'),
    searchRow,
    searchStatus,
    ui.h3('PDB IDs'),
    pdbInput.root,
    addFileRow,
    pdbStatus,
    ...(resultsSection ? [resultsSection] : []),
    ui.h3('Demo'),
    demoLabel,
    demoHero,
  ], 'cp-step-panel');
}

// ---------------------------------------------------------------------------
// Step 2 — Align (RMSD table)
// ---------------------------------------------------------------------------

export function buildStep2(ctx: StepContext): HTMLElement {
  const els: HTMLElement[] = [
    panelTitle('2. Align',
      'Stage 2a Kabsch superimposition. Each PDB is rotated + translated into a shared reference frame.'),
    // Reference-PDB selector lives here so users can decide which structure
    // is the template (RMSD = 0 row) before clicking Run step. Empty list
    // (no PDBs yet) → just show Auto.
    ui.h3('Alignment options'),
    buildAlignOptionsSection(ctx),
  ];

  if (!ctx.alignedData) {
    els.push(placeholderText('Click "Run step" in the footer to compute alignment.'));
    return ui.divV(els, 'cp-step-panel');
  }

  const df = ctx.alignedData;
  // Stage 2a output uses `rmsd_to_ref`; mode is "Cα" for 2a (global) vs
  // "Pocket Cα" if we ended up with the 2b refinement.
  const rmsdCol = df.columns.byName('rmsd_to_ref');
  const refPdbCol = df.columns.byName('ref_pdb_id');
  const pdbIdCol = df.columns.byName('pdb_id');

  // The reference PDB is the one named in `ref_pdb_id` (same value on every
  // row). It's also the row with RMSD = 0 by construction.
  let referencePdb = '–';
  if (refPdbCol && df.rowCount > 0)
    referencePdb = String(refPdbCol.get(0) ?? '–');

  // Skip the reference row (rmsd = 0) and NaN values when computing avg.
  let sum = 0;
  let n = 0;
  if (rmsdCol) {
    for (let i = 0; i < df.rowCount; i++) {
      const v = Number(rmsdCol.get(i));
      if (Number.isFinite(v) && v > 0) { sum += v; n += 1; }
    }
  }
  const avgRmsd = n > 0 ? (sum / n).toFixed(2) : '–';
  const mode = df.name === 'aligned_structures (pass 2)' ? 'Pocket Cα' : 'Cα';

  els.push(statStrip([
    {value: df.rowCount, label: 'Structures'},
    {value: `${avgRmsd} Å`, label: 'Avg RMSD'},
    {value: referencePdb, label: 'Reference'},
    {value: mode, label: 'Mode'},
  ]));

  // Compact RMSD table. The reference row gets a ★ marker so the user can
  // tell why its RMSD is 0.
  //
  // Row-click UX: clicking a row sets that PDB as the alignment reference
  // for the next run (same as picking it from the dropdown above). Re-clicking
  // the currently-selected reference reverts to Auto. The picked row gets a
  // highlight so the user sees their selection before clicking Run step.
  // The ★ marker still shows the reference of the LAST completed run — when
  // it differs from the user's new pick, the stale banner explains why the
  // displayed RMSDs no longer match.
  const pendingRef = (ctx.options.refPdbId ?? '').toUpperCase();
  // Per-PDB resolution (Å) from the Stage 1 QC data. Auto-reference picks the
  // highest-resolution entry, so surfacing resolution here lets the user make
  // an informed manual reference choice. '–' when QC data is unavailable.
  const resByPdb = new Map<string, number>();
  if (ctx.pdbQcData) {
    const qPid = ctx.pdbQcData.col('pdb_id');
    const qRes = ctx.pdbQcData.col('resolution');
    if (qPid && qRes) {
      for (let i = 0; i < ctx.pdbQcData.rowCount; i++) {
        const id = String(qPid.get(i) ?? '').toUpperCase();
        const r = Number(qRes.get(i));
        if (id && Number.isFinite(r)) resByPdb.set(id, r);
      }
    }
  }
  const cols: Array<{key: string; header: string;
    format?: (v: any, rowIndex: number, df: DG.DataFrame) => string | HTMLElement}> = [
    {key: 'pdb_id', header: 'PDB',
      format: (v, rowIndex, theDf) => {
        const id = String(v ?? '');
        const isRef = pdbIdCol && refPdbCol &&
          String(pdbIdCol.get(rowIndex)) === String(refPdbCol.get(rowIndex));
        if (!isRef) return id;
        const wrap = document.createElement('span');
        wrap.textContent = id;
        const star = document.createElement('span');
        star.textContent = ' ★';
        star.title = 'Reference structure of the last completed alignment — ' +
          'all RMSDs are computed against this PDB.';
        star.style.color = '#FFB300';
        star.style.fontWeight = '700';
        wrap.append(star);
        return wrap;
      }},
    {key: 'pdb_id', header: 'Resolution (Å)',
      format: (_v, rowIndex) => {
        const id = pdbIdCol ? String(pdbIdCol.get(rowIndex) ?? '').toUpperCase() : '';
        const r = resByPdb.get(id);
        return r == null ? '–' : r.toFixed(2);
      }},
    {key: 'rmsd_to_ref', header: 'RMSD (Å)',
      format: (v) => (v == null || v === '' || !Number.isFinite(Number(v)) ?
        '–' : Number(v).toFixed(2))},
  ];
  const tableEl = dataFrameTable(df, cols, 20);
  const isolatedUpper = (ctx.isolatedPdb ?? '').toUpperCase();
  const refUpper = referencePdb.toUpperCase();
  const rows = tableEl.querySelectorAll<HTMLTableRowElement>('tbody > tr');
  rows.forEach((tr, i) => {
    const id = pdbIdCol ? String(pdbIdCol.get(i) ?? '') : '';
    if (!id) return;
    const idUpper = id.toUpperCase();
    tr.classList.add('cp-rmsd-row');
    // Highlight the row that is currently ISOLATED in the viewer (one PDB
    // shown, the rest hidden). We deliberately do NOT highlight on the pending
    // reference choice: a *committed* reference is shown by the ★ marker, and a
    // *pending* (not-yet-run) reference is always also isolated — both the
    // row-click and the dropdown set isolation alongside the ref — so
    // `isolatedUpper` already covers the pre-run selection. Keying the
    // highlight on isolation alone means a completed "Run step" (which clears
    // isolation) visually DESELECTS the row, leaving just the ★ on the new
    // reference instead of a stuck blue highlight.
    const selected = isolatedUpper && idUpper === isolatedUpper;
    if (selected) tr.classList.add('cp-rmsd-row-isolated');
    tr.title = `Click to align to ${id} (sets the reference for the next run).`;
    tr.addEventListener('click', () => {
      const idUpper = id.toUpperCase();
      const resolvedRefUpper = referencePdb.toUpperCase();
      const alreadyPending = pendingRef && idUpper === pendingRef;

      // The reference choice:
      //  - If we're toggling off the currently-pending ref → Auto.
      //  - If the clicked PDB IS the resolved reference of the displayed
      //    data → Auto. Setting refPdbId='5HG8' explicitly when 5HG8 is
      //    already the auto-picked ref produces an identical Stage 2a
      //    output, so we treat it as a no-op for the alignment and avoid
      //    spurious "stale" state.
      //  - Otherwise → set as new explicit reference.
      let nextRefPdbId: string | undefined;
      if (alreadyPending) nextRefPdbId = undefined;
      else if (idUpper === resolvedRefUpper) nextRefPdbId = undefined;
      else nextRefPdbId = idUpper;

      // The isolation choice is independent of the ref choice: toggle off
      // if the user clicks the row that's already isolated, else isolate
      // to the clicked PDB. This decouples the two interactions so the
      // user can still visually focus on the resolved-ref PDB.
      const currentIsolated = (ctx.isolatedPdb ?? '').toUpperCase();
      const nextIsolatedPdb = currentIsolated === idUpper ? null : idUpper;

      ctx.onOptionsChanged({...ctx.options, refPdbId: nextRefPdbId});
      ctx.onIsolatePdb?.(nextIsolatedPdb);
      // Re-render the step so the dropdown selection, row highlight, and
      // hint text all reflect the new pending reference. The Mol* viewer +
      // detail panel are NOT re-rendered (different DOM regions).
      ctx.refreshStep?.();
    });
  });

  // Hint above the table. Four states the user can be in:
  //   1. No selection (no pending ref, nothing isolated): default invitation.
  //   2. Pending ref differs from the resolved ref ('★') → STALE (hot color),
  //      explains the next "Run step" will re-align.
  //   3. Pending ref equals resolved ref → explicit lock; click again for Auto.
  //   4. Isolation only (pending ref is Auto, viewer shows just one PDB) →
  //      explains how to clear, and that the ref is unchanged.
  let hintText: string;
  let hintHot = false;
  if (pendingRef && pendingRef !== refUpper) {
    hintText = `Pending reference: ${pendingRef} — click "Run step" in the footer to re-align. The ★ still marks the reference of the last completed run (${referencePdb}).`;
    hintHot = true;
  } else if (pendingRef) {
    hintText = `Reference: ${pendingRef}. Click any other row to switch, or click ${pendingRef} again to revert to Auto.`;
  } else if (isolatedUpper && isolatedUpper === refUpper) {
    hintText = `Showing only ${isolatedUpper} (the current alignment reference). Click ${isolatedUpper} again to show all, or click another row to switch the reference.`;
  } else if (isolatedUpper) {
    hintText = `Showing only ${isolatedUpper}. Click ${isolatedUpper} again to show all 5, or click another row to switch.`;
  } else {
    hintText = 'Click any row to set that PDB as the alignment reference (or use the dropdown above). Auto picks the highest-resolution entry.';
  }
  const hint = ui.divText(hintText, 'cp-rmsd-hint');
  hint.style.fontSize = '11px';
  hint.style.color = hintHot ? '#B26A00' : '#888';
  hint.style.fontStyle = 'italic';
  hint.style.marginBottom = '4px';
  els.push(hint);
  els.push(tableEl);

  return ui.divV(els, 'cp-step-panel');
}

// ---------------------------------------------------------------------------
// Step 3 — Pocket
// ---------------------------------------------------------------------------

export function buildStep3(ctx: StepContext): HTMLElement {
  const els: HTMLElement[] = [
    panelTitle('3. Pocket',
      `Stage 3 isolates binding-site residues within ${ctx.options.pocketRadius.toFixed(1)} Å of any drug ligand. ` +
      'A mid-pipeline picker will appear if the structures split into multiple conformational clusters.'),
    // Pocket options live on this step so users can tweak the cutoff /
    // method / representation right where they take effect — kept always
    // visible so the user can change a knob and re-run without hunting
    // for a collapsed accordion.
    ui.h3('Pocket options'),
    buildPocketOptionsSection(ctx),
  ];

  if (!ctx.pocketAtoms) {
    els.push(placeholderText('Click "Run step" in the footer to extract the pocket.'));
    return ui.divV(els, 'cp-step-panel');
  }

  const df = ctx.pocketAtoms;
  // Pocket atom counts per PDB
  const pdbCol = df.columns.byName('pdb_id');
  const counts = new Map<string, number>();
  if (pdbCol) {
    for (let i = 0; i < df.rowCount; i++) {
      const id = String(pdbCol.get(i));
      counts.set(id, (counts.get(id) ?? 0) + 1);
    }
  }

  els.push(statStrip([
    {value: counts.size, label: 'PDBs'},
    {value: df.rowCount, label: 'Pocket atoms'},
    {value: `${ctx.options.pocketRadius.toFixed(1)} Å`, label: 'Cutoff'},
  ]));

  // Build summary table: PDB | atom count
  const summaryRows = Array.from(counts.entries());
  const table = document.createElement('table');
  table.className = 'cp-step-table';
  const thead = document.createElement('thead');
  const headerRow = document.createElement('tr');
  for (const h of ['PDB', 'Pocket atoms']) {
    const th = document.createElement('th');
    th.textContent = h;
    headerRow.append(th);
  }
  thead.append(headerRow);
  table.append(thead);
  const tbody = document.createElement('tbody');
  // Row-click isolation (same UX as Step 2 / Step 4): clicking a PDB row hides
  // every other PDB's protein AND pocket overlay; clicking the isolated row
  // again shows all. Pure visualization toggle — no Python re-run.
  const isolatedUpper = (ctx.isolatedPdb ?? '').toUpperCase();
  for (const [id, n] of summaryRows) {
    const idUpper = id.toUpperCase();
    const tr = document.createElement('tr');
    tr.classList.add('cp-rmsd-row');
    if (isolatedUpper && idUpper === isolatedUpper)
      tr.classList.add('cp-rmsd-row-isolated');
    tr.title = isolatedUpper === idUpper
      ? 'Click again to show all PDBs.'
      : `Click to isolate ${id} (hides other proteins + their pocket atoms).`;
    const tdPdb = document.createElement('td');
    tdPdb.textContent = id;
    const tdN = document.createElement('td');
    tdN.textContent = String(n);
    tr.append(tdPdb, tdN);
    tr.addEventListener('click', () => {
      if (!ctx.onIsolatePdb) return;
      ctx.onIsolatePdb(isolatedUpper === idUpper ? null : idUpper);
    });
    tbody.append(tr);
  }
  table.append(tbody);

  const hint = ui.divText(
    isolatedUpper
      ? `Isolated: ${isolatedUpper} — click ${isolatedUpper} again to show all.`
      : 'Click any row to isolate that PDB in the viewer (hides the others + their pocket atoms).',
    'cp-rmsd-hint');
  hint.style.fontSize = '11px';
  hint.style.color = isolatedUpper ? '#1976D2' : '#888';
  hint.style.fontStyle = 'italic';
  hint.style.marginBottom = '4px';

  els.push(hint, table);

  return ui.divV(els, 'cp-step-panel');
}

// ---------------------------------------------------------------------------
// Step 4 — Features (ProLIF interaction summary)
// ---------------------------------------------------------------------------

export function buildStep4(ctx: StepContext): HTMLElement {
  const els: HTMLElement[] = [
    panelTitle('4. Features',
      'Stage 4 (ProLIF) extracts protein-ligand interactions per PDB. ' +
      'Each row in the right panel is one detected interaction. Click a PDB · ligand row for details.'),
  ];

  if (!ctx.summaryData) {
    els.push(placeholderText('Click "Run step" in the footer to extract features.\n(~15–30 s per PDB)'));
    return ui.divV(els, 'cp-step-panel');
  }

  const df = ctx.summaryData;
  // Build chip totals per row (D×n A×n ...).
  els.push(statStrip([
    {value: df.rowCount, label: 'PDB · Ligand'},
    {value: getTotalInteractions(df), label: 'Interactions'},
    {value: ctx.pdbIds.length, label: 'PDBs'},
  ]));

  // Precompute which family columns exist on this summary df so we don't
  // repeat the lookup once per row.
  const familyCols = FAMILY_CODES
    .map((code) => ({code, col: df.columns.byName(SUMMARY_FAMILY_COL[code])}))
    .filter((p) => p.col != null);

  const cols: Array<{key: string; header: string;
    format?: (v: any, rowIndex: number, df: DG.DataFrame) => string | HTMLElement}> = [
    {
      // "Use" checkbox — include/exclude this whole PDB from the Step 5
      // consensus. Reads/writes options.consensusExcludedPdbs. stopPropagation
      // so ticking the box doesn't also trigger the row's Mol* isolation.
      key: 'pdb_id', header: 'Use',
      format: (_v, rowIndex, theDf) => {
        const id = String(theDf.col('pdb_id')?.get(rowIndex) ?? '').toUpperCase();
        const excluded = new Set((ctx.options.consensusExcludedPdbs ?? []).map((s) => s.toUpperCase()));
        const cb = document.createElement('input');
        cb.type = 'checkbox';
        cb.checked = !excluded.has(id);
        cb.title = 'Include this PDB in the Step 5 consensus (and show it in the 3D viewer).';
        cb.style.cursor = 'pointer';
        cb.addEventListener('click', (e) => e.stopPropagation());
        cb.addEventListener('change', () => {
          const next = new Set((ctx.options.consensusExcludedPdbs ?? []).map((s) => s.toUpperCase()));
          if (cb.checked) next.delete(id); else next.add(id);
          ctx.onOptionsChanged({...ctx.options, consensusExcludedPdbs: [...next]});
          // Un-ticking a PDB also hides it in Mol* (and re-ticking shows it).
          ctx.onTogglePdbVisibility?.(id, !cb.checked);
          ctx.refreshStep?.();
        });
        return cb;
      }},
    {key: 'pdb_id', header: 'PDB'},
    {key: 'ligand_comp_id', header: 'Ligand'},
    {key: 'n_total', header: 'Total'},
    {
      // Synthetic "Families" column — the cell renders the family chips
      // (D×n A×n …) for THIS row. We use n_total as the key so the cell
      // gets *some* value passed in, but `format` ignores it and reads
      // family-specific columns from the row directly.
      key: 'n_total',
      header: 'Families',
      format: (_v, rowIndex) => buildFamilyChipsCell(familyCols, rowIndex),
    },
  ];

  // Row-click isolation: clicking a PDB row hides the other PDBs (proteins
  // AND their chain-F interaction overlays) in Mol*. Re-clicking the
  // already-isolated row releases isolation (back to all PDBs). Wiring is
  // analogous to Step 2's row-click-sets-reference, but here it's a pure
  // visualization toggle — no Python re-run, no stale state.
  const isolatedUpper = (ctx.isolatedPdb ?? '').toUpperCase();
  const pdbIdCol = df.columns.byName('pdb_id');
  const tableEl = dataFrameTable(df, cols, 20);
  const rows = tableEl.querySelectorAll<HTMLTableRowElement>('tbody > tr');
  rows.forEach((tr, i) => {
    const id = pdbIdCol ? String(pdbIdCol.get(i) ?? '').toUpperCase() : '';
    if (!id) return;
    tr.classList.add('cp-rmsd-row');  // reuses Step 2's hover style
    if (isolatedUpper && id === isolatedUpper)
      tr.classList.add('cp-rmsd-row-isolated');
    tr.title = isolatedUpper === id
      ? `Click again to show all PDBs.`
      : `Click to isolate ${id} (hides other proteins + their interactions).`;
    tr.addEventListener('click', () => {
      if (!ctx.onIsolatePdb) return;
      // Drive the right-hand per-interaction detail panel to THIS PDB. The
      // panel subscribes to the summary DataFrame's onCurrentRowChanged (see
      // orchestrator.setupPdbDetailPanel), so setting currentRowIdx here is
      // what makes "click a PDB row" update the interaction list — otherwise
      // the panel stays stuck on row 0 and only the Mol* isolation changes.
      df.currentRowIdx = i;
      const next = isolatedUpper === id ? null : id;
      ctx.onIsolatePdb(next);
    });
  });

  const hint = ui.divText(
    isolatedUpper
      ? `Isolated: ${isolatedUpper} — click ${isolatedUpper} again to show all.`
      : 'Click any row to isolate that PDB in the viewer (hides the others + their interactions).',
    'cp-rmsd-hint');
  hint.style.fontSize = '11px';
  hint.style.color = isolatedUpper ? '#1976D2' : '#888';
  hint.style.fontStyle = 'italic';
  hint.style.marginBottom = '4px';

  // Consensus-selection summary: how many interactions / PDBs are currently
  // ticked to feed Step 5 (after BOTH the per-PDB "Use" boxes and the
  // per-interaction boxes in the right detail panel). Updates live on toggle.
  const used = ctx.consensusUsed ?? {interactions: 0, pdbs: 0};
  const consensusNote = ui.divText(
    `Step 5 consensus will use ${used.interactions} interaction(s) across ` +
    `${used.pdbs} of ${df.rowCount} PDB(s) — untick "Use" here, or an ` +
    `interaction in the right panel, to drop it.`,
    'cp-rmsd-hint');
  consensusNote.style.fontSize = '11px';
  consensusNote.style.color = '#555';
  consensusNote.style.marginBottom = '4px';

  els.push(hint, consensusNote, tableEl);

  return ui.divV(els, 'cp-step-panel');
}

/** Build a single-cell HTMLElement showing the family chip row
 *  (D×1 A×2 …) for the given row of a pdb_interaction_summary df. */
function buildFamilyChipsCell(
  familyCols: Array<{code: string; col: DG.Column | null}>,
  rowIndex: number,
): HTMLElement {
  const wrap = document.createElement('span');
  wrap.style.display = 'inline-flex';
  wrap.style.flexWrap = 'wrap';
  wrap.style.gap = '6px';
  wrap.style.alignItems = 'center';
  let anyShown = false;
  for (const {code, col} of familyCols) {
    const n = Number(col!.get(rowIndex) ?? 0);
    if (n <= 0) continue;
    anyShown = true;
    const chip = document.createElement('span');
    chip.style.display = 'inline-flex';
    chip.style.alignItems = 'center';
    chip.style.gap = '2px';
    chip.append(familyChip(code));
    const txt = document.createElement('span');
    txt.textContent = `${code}×${n}`;
    txt.style.color = '#555';
    chip.append(txt);
    wrap.append(chip);
  }
  if (!anyShown) {
    const dash = document.createElement('span');
    dash.textContent = '–';
    dash.style.color = '#aaa';
    wrap.append(dash);
  }
  return wrap;
}

function getTotalInteractions(summary: DG.DataFrame): number {
  const col = summary.columns.byName('n_total');
  if (!col) return 0;
  let total = 0;
  for (let i = 0; i < summary.rowCount; i++) total += Number(col.get(i) ?? 0);
  return total;
}

// ---------------------------------------------------------------------------
// Step 5 — Consensus
// ---------------------------------------------------------------------------

export function buildStep5(ctx: StepContext): HTMLElement {
  const els: HTMLElement[] = [
    panelTitle('5. Consensus',
      'Stage 5a k-means clustering per family. Each row is a consensus pharmacophore point ' +
      '(chain P in Mol*, sphere radius = max spread, B-factor = ligand frequency).'),
    // Consensus knobs + reference-PDB override — always visible so the user
    // can adjust kq / min-fraction and re-run without re-expanding.
    ui.h3('Consensus options'),
    buildConsensusOptionsSection(ctx),
  ];

  if (!ctx.consensusData) {
    els.push(placeholderText('Click "Run step" in the footer to compute the consensus pharmacophore.'));
    return ui.divV(els, 'cp-step-panel');
  }

  const df = ctx.consensusData;

  els.push(buildCompletionBanner(df.rowCount));

  els.push(statStrip([
    {value: df.rowCount, label: 'Consensus points'},
    {value: countFamilies(df), label: 'Families'},
    {value: ctx.pdbIds.length, label: 'PDBs'},
  ]));

  // Provenance note: if the user excluded PDBs in Step 4, show what actually
  // fed this consensus so the result is traceable. (Silent when nothing's
  // excluded — the stat strip already covers the all-in case.)
  const excludedPdbs = new Set((ctx.options.consensusExcludedPdbs ?? []).map((s) => s.toUpperCase()));
  const nExclInts = (ctx.options.consensusExcludedInteractions ?? []).length;
  if ((excludedPdbs.size > 0 || nExclInts > 0) && ctx.consensusUsed) {
    const totPdbs = ctx.summaryData?.rowCount ?? ctx.consensusUsed.pdbs;
    const what: string[] = [];
    if (excludedPdbs.size > 0) what.push(`PDB(s): ${[...excludedPdbs].join(', ')}`);
    if (nExclInts > 0) what.push(`${nExclInts} interaction(s)`);
    const note = ui.divText(
      `Built from ${ctx.consensusUsed.interactions} interaction(s) across ` +
      `${ctx.consensusUsed.pdbs} of ${totPdbs} PDB(s) — excluded in Step 4: ${what.join('; ')}.`,
      'cp-rmsd-hint');
    note.style.fontSize = '11px';
    note.style.color = '#B26A00';
    note.style.fontStyle = 'italic';
    note.style.marginBottom = '6px';
    els.push(note);
  }

  const cols: Array<{key: string; header: string;
    format?: (v: any, rowIndex: number, df: DG.DataFrame) => string | HTMLElement}> = [
    {
      key: 'family',
      header: 'Family',
      format: (v) => {
        const wrap = document.createElement('span');
        wrap.append(familyChip(String(v)));
        const txt = document.createElement('span');
        const fam = resolveFamily(String(v));
        txt.textContent = ` ${fam.name}`;
        wrap.append(txt);
        return wrap;
      },
    },
    {key: 'x', header: 'x', format: (v) => Number(v).toFixed(2)},
    {key: 'y', header: 'y', format: (v) => Number(v).toFixed(2)},
    {key: 'z', header: 'z', format: (v) => Number(v).toFixed(2)},
    {
      key: 'frequency',
      header: 'Freq',
      format: (v) => (v == null ? '–' : (Number(v) * 100).toFixed(0) + '%'),
    },
  ];
  if (df.columns.byName('cluster_radius_a')) {
    cols.push({
      key: 'cluster_radius_a',
      header: 'r (Å)',
      format: (v) => (v == null ? '–' : Number(v).toFixed(2)),
    });
  }
  els.push(dataFrameTable(df, cols, 50));

  return ui.divV(els, 'cp-step-panel');
}

function countFamilies(consensus: DG.DataFrame): number {
  const col = consensus.columns.byName('family');
  if (!col) return 0;
  const s = new Set<string>();
  for (let i = 0; i < consensus.rowCount; i++) s.add(String(col.get(i)));
  return s.size;
}

function buildCompletionBanner(n: number): HTMLElement {
  const icon = ui.divText('✓', 'cp-wizard-completion-banner-icon');
  const titleEl = ui.divText('Consensus pharmacophore ready',
    'cp-wizard-completion-banner-title');
  const body = ui.divText(`${n} consensus point${n === 1 ? '' : 's'} identified. ` +
    'Inspect spheres in the Mol* viewer (chain P) or export below.',
    'cp-wizard-completion-banner-body');
  const text = ui.divV([titleEl, body]);
  return ui.divH([icon, text], 'cp-wizard-completion-banner');
}

