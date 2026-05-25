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
  // Stage 5a consensus k-means knobs. Pinned defaults match the values that
  // were hard-coded in the Python script before they were exposed to the UI.
  kq: number;                  // Avg features per cluster (n_clusters = ceil(n_features/kq)); default 7
  minClusterSizeFraction: number; // 0.0-1.0; cluster kept only if >= this fraction of ligands contribute; default 0.75
  topClusterNumber: number;    // Max clusters retained per family; default 4
}

export const DEFAULT_OPTIONS: PipelineOptions = {
  maxResolution: 2.5,
  requireXray: true,
  minLigandMw: 100,
  pocketMethod: 'cutoff',
  pocketRadius: 5.0,
  pocketRep: 'spacefill',
  kq: 7,
  minClusterSizeFraction: 0.75,
  topClusterNumber: 4,
};
