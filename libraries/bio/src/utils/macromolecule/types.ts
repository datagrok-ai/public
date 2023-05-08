import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ALPHABET} from './consts';
import {UnitsHandler} from '../units-handler';

export type SeqColStats = { freq: MonomerFreqs, sameLength: boolean }
export type SplitterFunc = (seq: string) => string[];
export type MonomerFreqs = { [m: string]: number };

/** Alphabet candidate type */
export type CandidateType = [string, Set<string>, number];
/** Alphabet candidate similarity type */
export type CandidateSimType = [string, Set<string>, number, MonomerFreqs, number];
