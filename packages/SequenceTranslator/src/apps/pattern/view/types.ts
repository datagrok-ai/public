import * as DG from 'datagrok-api/dg';

import {STRANDS, STRAND_ENDS, TERMINI} from '../model/const';

export type BooleanInput = DG.InputBase<boolean | null>;
export type StringInput = DG.InputBase<string | null>;
export type NumberInput = DG.InputBase<number | null>;

export type Position = { x: number, y: number };

export type StrandEndToNumberMap = Record<typeof STRAND_ENDS[number], number>;
export type StrandToNumberMap = Record<typeof STRANDS[number], number>;
export type StrandEndToSVGElementsMap = Record<typeof STRAND_ENDS[number], SVGElement | null>;
export type TerminusToSVGElementMap = Record<typeof TERMINI[number], SVGElement | null>
