/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Subject} from 'rxjs';

import {testEvent} from '@datagrok-libraries/utils/src/test';
import {NOTATION, TAGS as bioTAGS, ALIGNMENT, ALPHABET} from '@datagrok-libraries/bio/src/utils/macromolecule';
import { fetchWrapper } from '@datagrok-libraries/utils/src/fetch-utils';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';

import {checkForSingleSeqClusters} from './multiple-sequence-alignment';
import * as C from './constants';

import {_package} from '../package';

export const Pepsea = new class {
  public readonly dcName: string = 'bio';

  public async getDockerContainer(): Promise<DG.DockerContainer> {
    return await grok.dapi.docker.dockerContainers.filter(this.dcName).first();
  }
}();

export const pepseaMethods = ['mafft --auto', 'mafft', 'linsi', 'ginsi', 'einsi', 'fftns', 'fftnsi', 'nwns', 'nwnsi'];
const alignmentObjectMetaKeys = ['AlignedSeq', 'AlignedSubpeptide', 'HELM', 'ID', 'PolymerID'];
type PepseaResponse = {
  Alignment: {
    PolymerID: string, AlignedSubpeptide: string, HELM: string, ID: string, AlignedSeq: string, [key: string]: string,
  }[],
  AlignmentScore: { [key: string]: number | null },
};
type PepseaBodyUnit = { ID: string, HELM: string };

/** Gets the column containing MSA sequences produced by the 'PepSeA' tool from the {@link srcCol} column.
 * Does not add the result column to the dataframe of {@link srcCol}.
 * @async
 * @param {DG.Column} srcCol - The column containing the sequences to be aligned.
 * @param {string} unUsedName - The name of the result column.
 * @param {string} method - The method used for alignment.
 * @param {number} gapOpen - The gap open penalty.
 * @param {number} gapExtend - The gap extension penalty.
 * @param {DG.Column} clustersCol - The column containing the clusters of the sequences.
 * @param logger {ILogger} Logger
 */
export async function runPepsea(srcCol: DG.Column<string>, unUsedName: string,
  method: typeof pepseaMethods[number] = 'ginsi', gapOpen: number = 1.53, gapExtend: number = 0.0,
  clustersCol: DG.Column<string | number> | null = null, logger?: ILogger
): Promise<DG.Column<string>> {
  const pepseaContainer = await Pepsea.getDockerContainer();
  const peptideCount = srcCol.length;
  clustersCol ??= DG.Column.int('Clusters', peptideCount).init(0);
  if (clustersCol.type != DG.COLUMN_TYPE.STRING)
    clustersCol = clustersCol.convertTo(DG.TYPE.STRING);

  const clustersColCategories = clustersCol.categories;
  const clustersColData = clustersCol.getRawData();
  const bodies: PepseaBodyUnit[][] = new Array(clustersColCategories.length);
  const clusterIndexes: number[][] = new Array(clustersColCategories.length);

  // Grouping data by clusters
  for (let rowIndex = 0; rowIndex < peptideCount; ++rowIndex) {
    const clusterCategoryIdx = clustersColData[rowIndex];
    const cluster = clustersColCategories[clusterCategoryIdx];
    if (cluster === '')
      continue;

    const clusterId = clustersColCategories.indexOf(cluster);
    const helmSeq = srcCol.get(rowIndex);
    if (helmSeq) {
      (bodies[clusterId] ??= []).push({ID: rowIndex.toString(), HELM: helmSeq});
      (clusterIndexes[clusterCategoryIdx] ??= []).push(rowIndex);
    }
  }
  checkForSingleSeqClusters(clusterIndexes, clustersColCategories);

  const alignedSequences: string[] = new Array(peptideCount);
  for (const body of bodies) { // getting aligned sequences for each cluster
    const alignedObject = await requestAlignedObjects(pepseaContainer.id, body, method, gapOpen, gapExtend, logger);
    const alignments = alignedObject.Alignment;

    for (const alignment of alignments) { // filling alignedSequencesCol
      alignedSequences[parseInt(alignment.ID)] = Object.entries(alignment)
        .filter((v) => !alignmentObjectMetaKeys.includes(v[0]))
        .map((v) => v[1] !== '-' ? v[1] : '')
        .join(C.PEPSEA.SEPARATOR);
    }
  }

  const alignedSequencesCol: DG.Column<string> = DG.Column.fromStrings(unUsedName, alignedSequences);
  alignedSequencesCol.meta.units = NOTATION.SEPARATOR;
  alignedSequencesCol.setTag(bioTAGS.separator, C.PEPSEA.SEPARATOR);
  alignedSequencesCol.setTag(bioTAGS.aligned, ALIGNMENT.SEQ_MSA);
  alignedSequencesCol.setTag(bioTAGS.alphabet, ALPHABET.UN);
  alignedSequencesCol.setTag(bioTAGS.alphabetIsMultichar, 'true');
  alignedSequencesCol.semType = DG.SEMTYPE.MACROMOLECULE;

  return alignedSequencesCol;
}

async function requestAlignedObjects(
  dockerfileId: string, body: PepseaBodyUnit[], method: string, gapOpen: number, gapExtend: number, logger?: ILogger
): Promise<PepseaResponse> {
  const params = {
    method: 'POST',
    headers: {'Accept': 'application/json', 'Content-Type': 'application/json'},
    body: JSON.stringify(body),
  };
  const path = `/align?method=${method}&gap_open=${gapOpen}&gap_extend=${gapExtend}`;
  let responseObj: any;
  if ('fetchProxy' in grok.dapi.docker.dockerContainers) {
    // new dockerContainers API
    const t1: number = window.performance.now();
    // @ts-ignore
    const response: Response = await grok.dapi.docker.dockerContainers.fetchProxy(dockerfileId, path, params);
    const t2: number = window.performance.now();
    _package.logger.debug(`Bio: requestAlignedObjects() dockerContainers.fetchProxy(), ET: ${(t2 - t1)} ms`);
    const responseContentType = response.headers.get('content-type');
    const isJson: boolean = responseContentType === 'application/json';
    if (!response.ok && isJson) {
      const responseJson = await response.json();
      const pepseaErrorMsg = responseJson['pepsea-error'];
      if (!!pepseaErrorMsg)
        throw new Error(`PepSeA error: ${pepseaErrorMsg}`);

      const datagrokErrorMsg = responseJson['datagrok-error'];
      if (!!datagrokErrorMsg)
        throw new Error(`Datagrok error: ${datagrokErrorMsg}`);

      throw new Error(response.statusText);
    } else if (!response.ok && !isJson) {
      const responseStr = await response.text();
      throw new Error(`Error: ${responseStr}`);
    } else if (!isJson) {
      const responseStr = await response.text();
      throw new Error(`Error: PepSeA expected JSON response, got '${responseStr}'.`);
    }
    responseObj = await response.json();
  } else {
    // @ts-ignore
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(dockerfileId, path, params)!;
    const responseStr = await response.text();
    if (!responseStr)
      throw new Error('Empty response');
    responseObj = JSON.parse(responseStr);

    const pepseaErrorMsg = responseObj['pepsea-error'];
    if (!!pepseaErrorMsg)
      throw new Error(`PepSeA error: ${pepseaErrorMsg}`);

    const datagrokErrorMsg = responseObj['datagrok-error'];
    if (!!datagrokErrorMsg)
      throw new Error(`Datagrok error: ${datagrokErrorMsg}`);
  }
  // Check for pepsea stderr output
  if ('pepsea-stderr' in responseObj) {
    const pepseaStdErr: string = responseObj['pepsea-stderr'] as string;
    logger?.warning(pepseaStdErr);
  }
  return responseObj as PepseaResponse;
}
