import {ALIGNMENT, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

export type DetectorTestData = { [testName: string]: { csv: string, neg?: string[], pos?: { [colName: string]: PosCol } } };

export type DfReaderFunc = () => Promise<DG.DataFrame>;

export class PosCol {
  constructor(
    public readonly units: NOTATION,
    public readonly aligned: ALIGNMENT | null,
    public readonly alphabet: string | null,
    public readonly alphabetSize: number,
    public readonly alphabetIsMultichar?: boolean,
    public readonly separator?: string,
  ) { };
}
