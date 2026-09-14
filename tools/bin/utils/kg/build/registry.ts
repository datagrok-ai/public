/// The ordered extractors `grok kg build` runs (build-plan.md WO-1) and the context they share.
import {TypeSystem} from '../types';
import {HomeSet} from '../homes';
import {Emitter} from './emitter';
import {homesExtractor} from './extract/homes';
import {packagesExtractor} from './extract/ts/packages';
import {functionsExtractor} from './extract/ts/functions';

export type Mode = 'full' | 'public';

export interface BuildContext {
  system: TypeSystem;
  /** Absolute paths. */
  kgRoot: string;
  repoRoot: string;
  mode: Mode;
  backlogDir?: string;
  /** The home documents, when the caller already loaded them for check. */
  homes?: HomeSet;
}

export interface Extractor {
  name: string;
  layer: string;
  modes: Mode[];
  run(ctx: BuildContext, emitter: Emitter): void | Promise<void>;
}

export const EXTRACTORS: Extractor[] = [homesExtractor, packagesExtractor, functionsExtractor];

/** The extractors for [mode], narrowed by `--only`; names that match nothing are returned for the caller to refuse. */
export function selectExtractors(mode: Mode, only?: string[]): {selected: Extractor[], unknown: string[]} {
  const known = EXTRACTORS.filter((e) => e.modes.includes(mode));
  if (!only) return {selected: known, unknown: []};
  return {selected: known.filter((e) => only.includes(e.name)), unknown: only.filter((n) => !EXTRACTORS.some((e) => e.name === n))};
}

/** Runs the extractors in order; an extractor that reports no status is `ok`, one that throws is `missing`. */
export async function runExtractors(extractors: Extractor[], ctx: BuildContext, emitter: Emitter): Promise<void> {
  for (const extractor of extractors) {
    try {
      await extractor.run(ctx, emitter);
      emitter.sources[extractor.name] ??= 'ok';
    } catch (e: any) {
      emitter.source(extractor.name, 'missing');
      emitter.problem('extractor_errors', `${extractor.name}: ${e.stack ?? e.message}`);
      console.error(`${extractor.name}: ${e.message}`);
    }
  }
}
