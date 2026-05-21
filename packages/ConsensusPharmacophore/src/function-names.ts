/*
 * String constants for grok.functions.call(...) targets. Keeps the orchestrator
 * from silently failing when a Python script's name directive drifts from the
 * TS-side call site. Update both files together.
 *
 * Blueprint reference - section 0 (orchestration shape), Phase 6.
 */

export const FN = {
  STAGE2A: 'ConsensusPharmacophore:PharmacophoreStage2AlignPockets',
  STAGE2B: 'ConsensusPharmacophore:PharmacophoreStage2bAlignPocketsPocket',
  STAGE3:  'ConsensusPharmacophore:PharmacophoreStage3IsolatePocket',
  STAGE4:  'ConsensusPharmacophore:PharmacophoreStage4ExtractFeatures',
  STAGE5A: 'ConsensusPharmacophore:PharmacophoreStage5aConsensusKmeans',
} as const;
