import * as DG from 'datagrok-api/dg';

export type MultipleSequenceAlignmentUIOptions = {
  col?: DG.Column<string> | null;
  clustersCol?: DG.Column | null;
  kalign?: { gapOpen?: number; gapExtend?: number; terminalGap?: number };
  /** Name of the MSA engine function to use for non-canonical sequences */
  engine?: string;
  /** Engine-specific parameters (e.g. method, gapOpen, gapExtend) */
  engineParams?: Record<string, any>;
};
