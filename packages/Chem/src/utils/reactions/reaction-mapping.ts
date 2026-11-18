/**
 * Atom mapping for the Reaction panel and the sketchers' Auto-map, never run over a column and never written into a
 * cell unasked: RXNMapper in Chem's Crux file (`initRxnMapper`, `chemMapReaction`, since crux-core f290a6d33), run on a
 * Crux worker so that the UI thread never waits for it, its weights (upstream's `pytorch_model.bin`, 3.2 MB, MIT:
 * files/rxnmapper; its notice NOTICE.rxnmapper beside this file) a file of this package read on the first mapping. The
 * panel offers mapping when Chem's Crux file exports the mapper.
 */
import * as grok from 'datagrok-api/grok';
import {cruxModuleSync} from '../../engine/crux-module';
import {getCruxService} from '../../crux/crux-searches';

/** The weights, a file of this package (NOTICE.rxnmapper beside this file: where it comes from, its licence), by its
 * full path: package.ts's `_package` is the platform's in the package's bundle only, not in its test bundle. */
const WEIGHTS_FILE = 'System:AppData/Chem/rxnmapper/pytorch_model.bin';
/** How long one mapping may take on a worker (13.5 ms a reaction median, 127 ms the most in crux-core's set; the first
 * one loads the weights). */
const MAP_TIMEOUT_MS = 30000;

/** The mapper's calls in Chem's Crux file: the weights handed over once (again: a no-op), then one reaction mapped,
 * as JSON `{mapped_rxn, confidence}` or `{error}` (`canonicalize`: the reaction's SMILES canonical first, the default,
 * as upstream's). */
export interface RxnMapper {
  initRxnMapper(weights: Uint8Array): void;
  chemMapReaction(reaction: string, canonicalize?: boolean): string;
}

/** A mapped reaction: its reaction SMILES with atom maps, and the mapper's confidence (0..1). */
export interface MappedReaction {
  mapped: string;
  confidence: number;
}

let weights: Promise<Uint8Array> | null = null;

/** Whether Chem's Crux file maps reactions. */
export function mapperExported(): boolean {
  const crux = cruxModuleSync();
  return !!crux && 'initRxnMapper' in crux && 'chemMapReaction' in crux;
}

/** The mapper ready to map: its weights read (once a page) and Chem's Crux workers started (once a page), what the
 * first mapping of a page waits for before its reaction is mapped; the Reaction panel's 5 s limit counts from after it
 * (a cold page's start is no mapping's time). */
export async function mapperReady(): Promise<void> {
  weights ??= grok.dapi.files.readAsBytes(WEIGHTS_FILE).catch((e) => {
    weights = null;
    throw e;
  });
  await weights;
  await getCruxService();
}

/** `reaction` with its atoms mapped, on a Crux worker (crux/jobs/reactions.ts); the mapper's error thrown. The agents
 * stay agents ({@link agentsBackOverTheArrow}). `timeoutMs` is the job's own limit: a job still running past it (and
 * CruxService's margin) has its worker replaced, which is how a mapping given up is cancelled (the Reaction panel's
 * 5 s). */
export async function mapReaction(reaction: string, timeoutMs: number = MAP_TIMEOUT_MS): Promise<MappedReaction> {
  weights ??= grok.dapi.files.readAsBytes(WEIGHTS_FILE).catch((e) => {
    weights = null;
    throw e;
  });
  const bytes = await weights;
  const result = await (await getCruxService())
    .idleWorkerJob<MappedReaction>('map reaction', [reaction, bytes], timeoutMs);
  return {...result, mapped: agentsBackOverTheArrow(result.mapped)};
}

/** Whether a SMILES component has a mapped atom (`[...:n]`, n > 0). */
const hasMappedAtom = (component: string): boolean => /:[1-9]\d*\]/.test(component);

/**
 * Whether reaction value `text` carries an atom mapping of its own: any mapped atom (map number > 0), in a top-level
 * bracket atom of its SMILES or SMARTS (its CX extension aside; a recursive SMARTS in a bracket atom kept whole), or
 * in an atom line of its RXN block (V2000: columns 61-63; V3000: the field after the coordinates). The Reaction panel
 * draws such a value's own mapping and never maps it (automatic mapping must not replace a mapping the value already
 * has).
 */
export function mappedInTheValue(text: string): boolean {
  if (text.trimStart().startsWith('$RXN')) {
    let v3000Atoms = false;
    return text.split('\n').some((line) => {
      if (/^M {2}V30 BEGIN ATOM/.test(line))
        v3000Atoms = true;
      else if (/^M {2}V30 END ATOM/.test(line))
        v3000Atoms = false;
      else if (v3000Atoms) { // M  V30 index type x y z aamap ...
        const fields = line.trim().split(/\s+/);
        return fields.length > 7 && +fields[7] > 0;
      } else // a V2000 atom line: x, y and z in 10 columns each, a space, the symbol
        return /^[ \d.+-]{30} [A-Za-z*]/.test(line) && +(line.substring(60, 63).trim() || '0') > 0;
      return false;
    });
  }
  const body = text.replace(/\s+\|[^|]*\|\s*$/, '');
  let depth = 0;
  let start = -1;
  for (let i = 0; i < body.length; i++) {
    if (body[i] === '[') {
      if (depth++ === 0)
        start = i + 1;
    } else if (body[i] === ']' && depth > 0 && --depth === 0) {
      const map = /:(\d+)$/.exec(body.slice(start, i));
      if (map && +map[1] > 0)
        return true;
    }
  }
  return false;
}

/** A CX section that holds fragment groups only (`f:0.1.2,4.5`): components grouped as one, by their numbers. */
const FRAGMENT_GROUPS = /^f:\d+(?:\.\d+)*(?:,\d+(?:\.\d+)*)*$/;

/**
 * RXNMapper writes every component left of the arrow as a reactant, the agents among them (upstream's convention: a
 * left-hand component none of whose atoms is mapped is a reagent). Those components go back over the arrow, after any
 * agents the mapped text already has, so that a mapped step draws its agents where the unmapped one drew them (else
 * H2O, H2, Pd and EtOH are drawn as reactants with plus signs once mapped). A reaction whose
 * left side has no mapped atom at all is returned as it is. The CX section's fragment groups (`|f:0.1.2|`: a salt as
 * one agent), which number the components in order, are numbered anew with them (the Isotopes example's K2CO3 once
 * kept `f:0.1.2` after the move, a group over the phenol and the iodomethane, which Crux refused, so the app's pane
 * wrote the mapped reaction as text); a CX section with anything else in it, its atoms numbered, is left as it is, the
 * agents with it.
 */
export function agentsBackOverTheArrow(mapped: string): string {
  const cx = /^(.*?)\s+\|([^|]*)\|\s*$/s.exec(mapped);
  const [body, groups] = cx ? [cx[1], cx[2]] : [mapped, null];
  if (groups !== null && !FRAGMENT_GROUPS.test(groups))
    return mapped;
  const parts = body.split('>');
  if (parts.length !== 3)
    return mapped;
  const [left, middle, right] = parts.map((p) => p.split('.').filter((c) => c));
  const kept = [...left.keys()].filter((i) => hasMappedAtom(left[i]));
  if (!kept.length || kept.length === left.length)
    return mapped;
  const moved = [...left.keys()].filter((i) => !hasMappedAtom(left[i]));
  // the components' new order: the mapped reactants, the agents (the text's, then the moved ones), the products (their
  // numbers as they were: as many components come before them)
  const order = [...kept, ...middle.map((_, i) => left.length + i), ...moved];
  const renumber = new Map(order.map((old, i) => [old, i]));
  const text = `${kept.map((i) => left[i]).join('.')}>${[...middle, ...moved.map((i) => left[i])].join('.')}>` +
    right.join('.');
  if (groups === null)
    return text;
  const renumbered = groups.slice(2).split(',').map((g) => g.split('.').map((n) => renumber.get(+n) ?? +n)
    .sort((a, b) => a - b).join('.')).join(',');
  return `${text} |f:${renumbered}|`;
}
