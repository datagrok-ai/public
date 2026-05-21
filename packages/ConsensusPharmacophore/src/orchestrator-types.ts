/*
 * Pipeline-wide option types. Lives in its own module so rcsb-client.ts,
 * future Python-wrapper helpers, and tests can import without pulling in the
 * orchestrator's runtime code.
 *
 * Blueprint reference - section 0 (PipelineOptions), Appendix B.
 */

export type PocketMethod = 'cutoff' | 'dbscan';
export type PocketRep = 'spacefill' | 'gaussian-surface';

export interface PipelineOptions {
  maxResolution: number;       // A; default 2.5
  requireXray: boolean;        // default true
  minLigandMw: number;         // Da; default 100
  pocketMethod: PocketMethod;  // default 'cutoff'
  pocketRadius: number;        // A; default 5.0
  pocketRep: PocketRep;        // default 'spacefill' — viz only
  refPdbId?: string;           // optional reference overlay
}

export const DEFAULT_OPTIONS: PipelineOptions = {
  maxResolution: 2.5,
  requireXray: true,
  minLigandMw: 100,
  pocketMethod: 'cutoff',
  pocketRadius: 5.0,
  pocketRep: 'spacefill',
};
