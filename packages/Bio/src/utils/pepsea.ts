/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {NOTATION, TAGS as bioTAGS, ALIGNMENT, ALPHABET} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';

import {checkForSingleSeqClusters} from './multiple-sequence-alignment';
import * as C from './constants';
import {_package} from '../package';

export const pepseaMethods = ['mafft --auto', 'mafft', 'linsi', 'ginsi', 'einsi', 'fftns', 'fftnsi', 'nwns', 'nwnsi'];

const ALIGNMENT_META_KEYS = ['AlignedSeq', 'AlignedSubpeptide', 'HELM', 'ID', 'PolymerID'];

type PepseaResponse = {
  Alignment: {
    PolymerID: string; AlignedSubpeptide: string; HELM: string;
    ID: string; AlignedSeq: string; [key: string]: string;
  }[];
  AlignmentScore: {[key: string]: number | null};
};

type PepseaBodyUnit = {ID: string; HELM: string};

export const Pepsea = new class {
  public readonly dcName: string = 'bio';

  public async getDockerContainer(): Promise<DG.DockerContainer> {
    return await grok.dapi.docker.dockerContainers.filter(this.dcName).first();
  }
}();


/** Aligns all sequences in the column using PepSeA Docker container.
 * Does not handle clustering - aligns all rows as a single group.
 * Used by the registered sequenceMSA function. */
export async function alignWithPepsea(
  srcCol: DG.Column<string>,
  method: string = 'mafft --auto',
  gapOpen: number = 1.53,
  gapExtend: number = 0,
): Promise<DG.Column<string>> {
  const container = await Pepsea.getDockerContainer();
  const rowCount = srcCol.length;

  const body: PepseaBodyUnit[] = [];
  for (let i = 0; i < rowCount; i++) {
    const seq = srcCol.get(i);
    if (seq)
      body.push({ID: i.toString(), HELM: seq});
  }

  const response = await requestAlignedObjects(container.id, body, method, gapOpen, gapExtend);
  const aligned = parseAlignmentResponse(response, rowCount);

  const colName = srcCol.dataFrame?.columns?.getUnusedName(`msa(${srcCol.name})`) ?? `msa(${srcCol.name})`;
  return createPepseaResultColumn(colName, aligned);
}


/** Aligns sequences with PepSeA, supporting per-cluster alignment.
 * Used by tests and legacy code paths. */
export async function runPepsea(
  table: DG.DataFrame, srcCol: DG.Column<string>, unUsedName: string,
  method: typeof pepseaMethods[number] = 'ginsi', gapOpen: number = 1.53, gapExtend: number = 0.0,
  clustersCol: DG.Column<string | number> | null = null, logger?: ILogger, onlySelected: boolean = false,
): Promise<DG.Column<string>> {
  const container = await Pepsea.getDockerContainer();
  const rowCount = srcCol.length;

  clustersCol ??= DG.Column.int('Clusters', rowCount).init(0);
  if (clustersCol.type !== DG.COLUMN_TYPE.STRING)
    clustersCol = clustersCol.convertTo(DG.TYPE.STRING);

  const categories = clustersCol.categories;
  const data = clustersCol.getRawData();
  const bodies: PepseaBodyUnit[][] = new Array(categories.length);
  const clusterIndexes: number[][] = new Array(categories.length);

  const rows = onlySelected ? selectedRows(table.selection) : allRows(rowCount);
  for (const rowIndex of rows) {
    const catIdx = data[rowIndex];
    if (!categories[catIdx]) continue;
    const helmSeq = srcCol.get(rowIndex);
    if (helmSeq) {
      (bodies[catIdx] ??= []).push({ID: rowIndex.toString(), HELM: helmSeq});
      (clusterIndexes[catIdx] ??= []).push(rowIndex);
    }
  }
  checkForSingleSeqClusters(clusterIndexes, categories);

  const alignedSequences: string[] = new Array(rowCount).fill(null);
  for (const body of bodies) {
    if (!body || body.length === 0) continue;
    const response = await requestAlignedObjects(container.id, body, method, gapOpen, gapExtend, logger);
    for (const alignment of response.Alignment)
      alignedSequences[parseInt(alignment.ID)] = extractAlignedSequence(alignment);
  }

  return createPepseaResultColumn(unUsedName, alignedSequences);
}


// --- Helpers ---

function extractAlignedSequence(alignment: PepseaResponse['Alignment'][0]): string {
  return Object.entries(alignment)
    .filter(([key]) => !ALIGNMENT_META_KEYS.includes(key))
    .map(([, val]) => val !== '-' ? val : '')
    .join(C.PEPSEA.SEPARATOR);
}

function parseAlignmentResponse(response: PepseaResponse, rowCount: number): string[] {
  const aligned: string[] = new Array(rowCount).fill(null);
  for (const alignment of response.Alignment)
    aligned[parseInt(alignment.ID)] = extractAlignedSequence(alignment);
  return aligned;
}

function createPepseaResultColumn(name: string, sequences: string[]): DG.Column<string> {
  const col = DG.Column.fromStrings(name, sequences);
  col.meta.units = NOTATION.SEPARATOR;
  col.setTag(bioTAGS.separator, C.PEPSEA.SEPARATOR);
  col.setTag(bioTAGS.aligned, ALIGNMENT.SEQ_MSA);
  col.setTag(bioTAGS.alphabet, ALPHABET.UN);
  col.setTag(bioTAGS.alphabetIsMultichar, 'true');
  col.semType = DG.SEMTYPE.MACROMOLECULE;
  return col;
}

function* selectedRows(selection: DG.BitSet): Generator<number> {
  for (let i = -1; (i = selection.findNext(i, true)) !== -1;)
    yield i;
}

function* allRows(count: number): Generator<number> {
  for (let i = 0; i < count; i++)
    yield i;
}

async function requestAlignedObjects(
  dockerfileId: string, body: PepseaBodyUnit[], method: string,
  gapOpen: number, gapExtend: number, logger?: ILogger,
): Promise<PepseaResponse> {
  const params = {
    method: 'POST',
    headers: {'Accept': 'application/json', 'Content-Type': 'application/json'},
    body: JSON.stringify(body),
  };
  const path = `/align?method=${method}&gap_open=${gapOpen}&gap_extend=${gapExtend}`;

  const t1 = window.performance.now();
  // @ts-ignore
  const response: Response = await grok.dapi.docker.dockerContainers.fetchProxy(dockerfileId, path, params);
  const t2 = window.performance.now();
  _package.logger.debug(`Bio: requestAlignedObjects() ET: ${(t2 - t1)} ms`);

  const contentType = response.headers.get('content-type');
  const isJson = contentType === 'application/json';

  if (!response.ok) {
    if (isJson) {
      const json = await response.json();
      if (json['pepsea-error']) throw new Error(`PepSeA error: ${json['pepsea-error']}`);
      if (json['datagrok-error']) throw new Error(`Datagrok error: ${json['datagrok-error']}`);
      throw new Error(response.statusText);
    }
    const text = await response.text();
    throw new Error(`Error: ${text}`);
  }

  if (!isJson) {
    const text = await response.text();
    throw new Error(`Error: PepSeA expected JSON response, got '${text}'.`);
  }

  const responseObj = await response.json();
  if ('pepsea-stderr' in responseObj)
    logger?.warning(responseObj['pepsea-stderr'] as string);
  return responseObj as PepseaResponse;
}
