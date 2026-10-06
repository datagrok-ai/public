/* eslint-disable max-len */
/**
 * AI functions for the PolyTool enumerators (HELM Enumerator dialog, Markush
 * Enumerator dialog and app view).
 *
 * The vocabulary is registered ONCE per session through `grok.functions.register`
 * (real platform Funcs, no package.ts entries). Every function takes the generic
 * `widget` argument the AI assistant injects (the targeted open dialog, or the
 * current view when no widget ref is passed) and resolves the live enumerator
 * behind it at call time. The overrides are attached to the wrapper the shell
 * actually reports (`DG.toJs(dart)`) — a dialog/view built with `new Dialog(...)`
 * / `DG.View.create()` is a different JS object than the one
 * `Dialog.getOpenDialogs()` / `grok.shell.v` later return.
 */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

import {PolyToolEnumeratorType, PolyToolEnumeratorTypes} from './types';
import {ChemEnumModes, ChemEnumMode, makeCore, addRGroupsFromSmiles, validateParams} from './pt-chem-enum';
import type {ChemEnumDialogState} from './pt-chem-enum-dialog';

const AI_NAMESPACE = 'SequenceTranslatorAI';

const _funcs = new Map<string, DG.Func>();

/** Registers one vocabulary Func on first use; later calls return the cached one. */
function aiFunc(name: string, signature: string, description: string, run: (...args: any[]) => any): DG.Func {
  let f = _funcs.get(name);
  if (!f) {
    f = grok.functions.register({signature, run, isAsync: true, namespace: AI_NAMESPACE,
      tags: 'st-ai-function', options: {description}});
    _funcs.set(name, f);
  }
  return f;
}

/** The wrapper instance the shell reports for a dialog or view — overrides live there. */
function shellWrapperOf(widget: {dart?: any}): any {
  return (widget?.dart ? DG.toJs(widget.dart) : null) ?? widget;
}

/** Standard functions of a Dart-backed dialog/widget (getDialogInfo, setInput, clickButton, ...). */
function dartFunctionsOf(w: any): DG.Func[] {
  try {
    return DG.DartWidget.prototype.getFunctions.call(w) ?? [];
  } catch (_) {
    return [];
  }
}

// ─────────────────────── HELM Enumerator (dialog) ───────────────────────

/** One breadth placeholder row as the AI passes it. */
export type BreadthPlaceholderInput = {start: number, end: number, monomers: string[] | string};

/** What the dialog closure hands over so the shared vocabulary can drive it. */
export interface HelmEnumeratorAiContext {
  getState(): {[key: string]: unknown};
  /** null when there is no validation problem. */
  validationError(): string | null;
  setMacromolecule(helm: string): Promise<void>;
  setPlaceholders(placeholders: {[position: string]: string[] | string}): void;
  setBreadthPlaceholders(placeholders: BreadthPlaceholderInput[]): void;
  setOption(name: string, value: unknown): boolean;
  execute(): Promise<DG.DataFrame | null>;
}

const helmCtxByWrapper = new WeakMap<object, HelmEnumeratorAiContext>();

function helmCtxOf(widget: any): HelmEnumeratorAiContext {
  const ctx = (widget != null ? helmCtxByWrapper.get(widget) : null) ??
    (widget?.jsView != null ? helmCtxByWrapper.get(widget.jsView) : null);
  if (!ctx)
    throw new Error('Target the open HELM Enumerator dialog (pass its widget ref from list_view_widgets)');
  return ctx;
}

function helmEnumeratorFunctions(): DG.Func[] {
  return [
    aiFunc('getHelmEnumeratorState', 'dynamic getHelmEnumeratorState(widget widget)',
      'The HELM Enumerator state - the source macromolecule (HELM), placeholders, enumerator type, options and the current validation problem. Call this first',
      async (widget: any) => helmCtxOf(widget).getState()),

    aiFunc('setHelmEnumeratorMacromolecule', 'dynamic setHelmEnumeratorMacromolecule(widget widget, string helm)',
      'Set the source macromolecule of the HELM Enumerator as a HELM string, e.g. PEPTIDE1{A.C.D}$$$$',
      async (widget: any, helm: string) => {
        const ctx = helmCtxOf(widget);
        if (!helm?.trim())
          return {success: false, error: 'helm is required'};
        await ctx.setMacromolecule(helm.trim());
        return {success: true};
      }),

    aiFunc('setHelmEnumeratorPlaceholders', 'dynamic setHelmEnumeratorPlaceholders(widget widget, map placeholders, list breadthPlaceholders)',
      'Replace the enumeration placeholders. placeholders maps a 1-based monomer position to a list of substitute monomer symbols, e.g. position 3 to A, C, D. breadthPlaceholders is an optional list of objects with start, end and monomers - each substitutes every position of the range. Pass at least one of the two',
      async (widget: any, placeholders?: {[position: string]: string[] | string}, breadthPlaceholders?: BreadthPlaceholderInput[]) => {
        const ctx = helmCtxOf(widget);
        const hasPoint = placeholders && Object.keys(placeholders).length > 0;
        const hasBreadth = Array.isArray(breadthPlaceholders) && breadthPlaceholders.length > 0;
        if (!hasPoint && !hasBreadth)
          return {success: false, error: 'Pass placeholders (position to monomers map) and/or breadthPlaceholders'};
        if (hasPoint) {
          const badPos = Object.keys(placeholders!).find((p) => !(parseInt(p) >= 1));
          if (badPos != null)
            return {success: false, error: `Invalid position '${badPos}' — positions are 1-based integers`};
          ctx.setPlaceholders(placeholders!);
        }
        if (hasBreadth) {
          const bad = breadthPlaceholders!.find((b) => !(b?.start >= 1) || !(b?.end >= b?.start));
          if (bad)
            return {success: false, error: 'Each breadth placeholder needs start >= 1 and end >= start'};
          ctx.setBreadthPlaceholders(breadthPlaceholders!);
        }
        return {success: true, ...(ctx.validationError() ? {validationError: ctx.validationError()} : {})};
      }),

    aiFunc('setHelmEnumeratorOptions', 'dynamic setHelmEnumeratorOptions(widget widget, map options)',
      'Set HELM Enumerator options. Keys: enumeratorType (single, parallel or matrix), keepOriginal, toAtomicLevel, generateHelm, chiralityEngine, highlightMonomers (booleans), trivialName (string)',
      async (widget: any, options: {[key: string]: unknown}) => {
        const ctx = helmCtxOf(widget);
        if (!options || !Object.keys(options).length)
          return {success: false, error: 'Pass options as a map, e.g. enumeratorType, toAtomicLevel'};
        if (options.enumeratorType != null &&
          !Object.values(PolyToolEnumeratorTypes).includes(options.enumeratorType as PolyToolEnumeratorType))
          return {success: false, error: `enumeratorType must be one of: ${Object.values(PolyToolEnumeratorTypes).join(', ')}`};
        const unknown = Object.keys(options).filter((k) => !ctx.setOption(k, options[k]));
        return {success: true, ...(unknown.length ? {unknownOptions: unknown} : {})};
      }),

    aiFunc('runHelmEnumerator', 'dynamic runHelmEnumerator(widget widget)',
      'Run the HELM enumeration - builds the table of enumerated sequences (a new table view, or appended to the chosen table). The dialog stays open',
      async (widget: any) => {
        const ctx = helmCtxOf(widget);
        const problem = ctx.validationError();
        if (problem)
          return {success: false, error: problem};
        const df = await ctx.execute();
        if (!df)
          return {success: false, error: 'Enumeration produced nothing — check the macromolecule and placeholders'};
        return {success: true, rows: df.rowCount, table: df.name || 'the result table',
          note: 'The enumerated sequences are in the result table'};
      }),
  ];
}

/** Called by the dialog builder: wires the shared vocabulary to this dialog instance. */
export function attachHelmEnumeratorAi(dialog: DG.Dialog, ctx: HelmEnumeratorAiContext): void {
  const w = shellWrapperOf(dialog);
  helmCtxByWrapper.set(w, ctx);
  if (w !== dialog)
    helmCtxByWrapper.set(dialog, ctx);
  w.aiDescription = 'HELM Enumerator — enumerates a macromolecule (HELM) by substituting monomers at chosen ' +
    'positions, producing a table of variants. Functions: getHelmEnumeratorState (call first), ' +
    'setHelmEnumeratorMacromolecule, setHelmEnumeratorPlaceholders, setHelmEnumeratorOptions, runHelmEnumerator. ' +
    'The generic dialog functions (setInput, clickButton) also work, but the enumerator functions are preferred.';
  w.getFunctions = () => [...helmEnumeratorFunctions(), ...dartFunctionsOf(w)];
}

// ─────────────────────── Markush Enumerator (dialog + app view) ───────────────────────

/** The slice of ChemEnumPanel the vocabulary needs (the panel interface itself is module-private). */
export interface MarkushAiContext {
  state: ChemEnumDialogState;
  rdkit: RDModule;
  refresh: () => void;
  execute: () => Promise<void>;
}

const markushCtxByWrapper = new WeakMap<object, MarkushAiContext>();

function markushCtxOf(widget: any): MarkushAiContext {
  const ctx = (widget != null ? markushCtxByWrapper.get(widget) : null) ??
    (widget?.jsView != null ? markushCtxByWrapper.get(widget.jsView) : null);
  if (!ctx)
    throw new Error('Target the open Markush Enumerator (the dialog widget ref, or its app view)');
  return ctx;
}

function markushValidation(ctx: MarkushAiContext) {
  return validateParams({cores: ctx.state.cores, rGroups: ctx.state.rGroupsByNum, mode: ctx.state.mode});
}

function markushEnumeratorFunctions(): DG.Func[] {
  return [
    aiFunc('getMarkushEnumeratorState', 'dynamic getMarkushEnumeratorState(widget widget)',
      'The Markush Enumerator state - the cores (SMILES with R-labels), the R-group substituents per R number, the mode and output options, plus validation errors and the predicted result count. Call this first',
      async (widget: any) => {
        const ctx = markushCtxOf(widget);
        const validation = markushValidation(ctx);
        return {
          cores: ctx.state.cores.map((c) => ({smiles: c.smiles, rNumbers: c.rNumbers, ...(c.error ? {error: c.error} : {})})),
          rGroups: Object.fromEntries([...ctx.state.rGroupsByNum.entries()].map(([n, list]) =>
            [`R${n}`, list.map((g) => ({smiles: g.smiles, ...(g.error ? {error: g.error} : {})}))])),
          mode: ctx.state.mode,
          removeDuplicates: ctx.state.removeDuplicates,
          tableName: ctx.state.tableName,
          appendToTable: ctx.state.appendToTable?.name ?? null,
          valid: validation.ok,
          ...(validation.errors.length ? {validationErrors: validation.errors} : {}),
          predictedCount: validation.predictedCount,
        };
      }),

    aiFunc('addMarkushCores', 'dynamic addMarkushCores(widget widget, list smiles)',
      'Add core structures to the Markush Enumerator - each a SMILES with at least one R-label, e.g. [*:1]c1ccccc1',
      async (widget: any, smiles: string[]) => {
        const ctx = markushCtxOf(widget);
        if (!Array.isArray(smiles) || !smiles.length)
          return {success: false, error: 'Pass smiles as a non-empty list'};
        const invalid: {smiles: string, error: string}[] = [];
        let added = 0;
        for (const smi of smiles) {
          const core = makeCore(String(smi ?? ''), '', ctx.rdkit);
          if (core.error)
            invalid.push({smiles: String(smi ?? ''), error: core.error});
          else {
            ctx.state.cores.push(core);
            added++;
          }
        }
        ctx.refresh();
        return {success: added > 0, added, ...(invalid.length ? {invalid} : {}),
          ...(added === 0 ? {error: 'No valid cores — each needs at least one R-label like [*:1]'} : {})};
      }),

    aiFunc('addMarkushRGroups', 'dynamic addMarkushRGroups(widget widget, int rNumber, list smiles, string mode)',
      'Add R-group substituents for one R number - each a SMILES with exactly one R-label, or a single atom like N, O, Cl. mode is optional: append (default) or replace',
      async (widget: any, rNumber: number, smiles: string[], mode?: string) => {
        const ctx = markushCtxOf(widget);
        if (!(rNumber >= 1))
          return {success: false, error: 'rNumber must be a positive integer'};
        if (!Array.isArray(smiles) || !smiles.length)
          return {success: false, error: 'Pass smiles as a non-empty list'};
        const m = (mode ?? 'append').toLowerCase();
        if (m !== 'append' && m !== 'replace')
          return {success: false, error: "mode must be 'append' or 'replace'"};
        const added = addRGroupsFromSmiles(ctx.state.rGroupsByNum, smiles.map((s) => String(s ?? '')),
          Math.round(rNumber), m as 'append' | 'replace', ctx.rdkit);
        ctx.refresh();
        const invalid = (ctx.state.rGroupsByNum.get(Math.round(rNumber)) ?? []).filter((g) => g.error)
          .map((g) => ({smiles: g.originalSmiles, error: g.error!}));
        return {success: true, rNumber: Math.round(rNumber), added, ...(invalid.length ? {invalid} : {})};
      }),

    aiFunc('clearMarkushItems', 'dynamic clearMarkushItems(widget widget, string kind, int rNumber)',
      'Remove Markush Enumerator inputs. kind is cores, rgroups or all. rNumber is optional - clears only that R-group slot',
      async (widget: any, kind: string, rNumber?: number) => {
        const ctx = markushCtxOf(widget);
        const k = (kind ?? '').toLowerCase();
        if (!['cores', 'rgroups', 'all'].includes(k))
          return {success: false, error: "kind must be 'cores', 'rgroups' or 'all'"};
        if (k === 'cores' || k === 'all')
          ctx.state.cores.length = 0;
        if (k === 'rgroups' || k === 'all') {
          if (rNumber != null && k === 'rgroups') {
            if (!ctx.state.rGroupsByNum.delete(Math.round(rNumber)))
              return {success: false, error: `No R-group slot R${rNumber}`};
          } else
            ctx.state.rGroupsByNum.clear();
        }
        ctx.refresh();
        return {success: true};
      }),

    aiFunc('setMarkushOptions', 'dynamic setMarkushOptions(widget widget, map options)',
      'Set Markush Enumerator options. Keys: mode (Zip or Cartesian), removeDuplicates (boolean), tableName (string)',
      async (widget: any, options: {[key: string]: unknown}) => {
        const ctx = markushCtxOf(widget);
        if (!options || !Object.keys(options).length)
          return {success: false, error: 'Pass options as a map, e.g. mode, removeDuplicates, tableName'};
        const unknown: string[] = [];
        for (const [key, value] of Object.entries(options)) {
          if (key === 'mode') {
            if (!Object.values(ChemEnumModes).includes(value as ChemEnumMode))
              return {success: false, error: `mode must be one of: ${Object.values(ChemEnumModes).join(', ')}`};
            ctx.state.mode = value as ChemEnumMode;
          } else if (key === 'removeDuplicates')
            ctx.state.removeDuplicates = !!value;
          else if (key === 'tableName')
            ctx.state.tableName = String(value ?? '');
          else
            unknown.push(key);
        }
        ctx.refresh();
        return {success: true, ...(unknown.length ? {unknownOptions: unknown} : {})};
      }),

    aiFunc('runMarkushEnumeration', 'dynamic runMarkushEnumeration(widget widget)',
      'Run the Markush enumeration - assembles every core and R-group combination into a molecule table (a new table view, or appended to the chosen table)',
      async (widget: any) => {
        const ctx = markushCtxOf(widget);
        const validation = markushValidation(ctx);
        if (!validation.ok)
          return {success: false, error: validation.errors.join(' ')};
        await ctx.execute();
        return {success: true, predictedCount: validation.predictedCount,
          note: `The enumerated molecules are in the '${ctx.state.tableName}' table`};
      }),
  ];
}

const MARKUSH_AI_DESCRIPTION = 'Markush Enumerator — enumerates a Markush structure: cores (SMILES with ' +
  'R-labels) combined with per-R-number substituent lists into a molecule table (Zip pairs the lists, ' +
  'Cartesian builds all combinations). Functions: getMarkushEnumeratorState (call first), addMarkushCores, ' +
  'addMarkushRGroups, clearMarkushItems, setMarkushOptions, runMarkushEnumeration.';

/** Called by the Markush dialog builder. */
export function attachMarkushDialogAi(dialog: DG.Dialog, ctx: MarkushAiContext): void {
  const w = shellWrapperOf(dialog);
  markushCtxByWrapper.set(w, ctx);
  if (w !== dialog)
    markushCtxByWrapper.set(dialog, ctx);
  w.aiDescription = MARKUSH_AI_DESCRIPTION;
  w.getFunctions = () => [...markushEnumeratorFunctions(), ...dartFunctionsOf(w)];
}

/** Called by the Markush app entry point. */
export function attachMarkushViewAi(view: DG.ViewBase, ctx: MarkushAiContext): void {
  const w = shellWrapperOf(view);
  markushCtxByWrapper.set(w, ctx);
  if (w !== view)
    markushCtxByWrapper.set(view, ctx);
  w.aiDescription = MARKUSH_AI_DESCRIPTION +
    ' Act through the view functions (search list_view_functions with "markush").';
  w.getFunctions = () => markushEnumeratorFunctions();
}
