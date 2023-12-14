/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Subject} from 'rxjs';

import {testEvent} from '@datagrok-libraries/utils/src/test';
import {NOTATION, TAGS as bioTAGS, ALIGNMENT, ALPHABET} from '@datagrok-libraries/bio/src/utils/macromolecule';

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
 */
export async function runPepsea(srcCol: DG.Column<string>, unUsedName: string,
  method: typeof pepseaMethods[number] = 'ginsi', gapOpen: number = 1.53, gapExtend: number = 0.0,
  clustersCol: DG.Column<string | number> | null = null,
): Promise<DG.Column<string> | null> {
  const pepseaContainer = await Pepsea.getDockerContainer();
  if (pepseaContainer.status !== 'started' && pepseaContainer.status !== 'checking') {
    grok.log.warning('PepSeA container has not started yet');
    return null;
  }


  const peptideCount = srcCol.length;
  clustersCol ??= DG.Column.int('Clusters', peptideCount).init(0);
  if (clustersCol.type != DG.COLUMN_TYPE.STRING)
    clustersCol = clustersCol.convertTo(DG.TYPE.STRING);

  const clusters = clustersCol.categories;
  const bodies: PepseaBodyUnit[][] = new Array(clusters.length);

  // Grouping data by clusters
  for (let rowIndex = 0; rowIndex < peptideCount; ++rowIndex) {
    const cluster = clustersCol.get(rowIndex) as string;
    if (cluster === '')
      continue;

    const clusterId = clusters.indexOf(cluster);
    const helmSeq = srcCol.get(rowIndex);
    if (helmSeq)
      (bodies[clusterId] ??= []).push({ID: rowIndex.toString(), HELM: helmSeq});
  }

  const alignedSequences: string[] = new Array(peptideCount);
  for (const body of bodies) { // getting aligned sequences for each cluster
    const alignedObject = await requestAlignedObjects(pepseaContainer.id, body, method, gapOpen, gapExtend);
    const alignments = alignedObject.Alignment;

    for (const alignment of alignments) { // filling alignedSequencesCol
      alignedSequences[parseInt(alignment.ID)] = Object.entries(alignment)
        .filter((v) => !alignmentObjectMetaKeys.includes(v[0]))
        .map((v) => v[1] !== '-' ? v[1] : '')
        .join(C.PEPSEA.SEPARATOR);
    }
  }

  const alignedSequencesCol: DG.Column<string> = DG.Column.fromStrings(unUsedName, alignedSequences);
  alignedSequencesCol.setTag(DG.TAGS.UNITS, NOTATION.SEPARATOR);
  alignedSequencesCol.setTag(bioTAGS.separator, C.PEPSEA.SEPARATOR);
  alignedSequencesCol.setTag(bioTAGS.aligned, ALIGNMENT.SEQ_MSA);
  alignedSequencesCol.setTag(bioTAGS.alphabet, ALPHABET.UN);
  alignedSequencesCol.setTag(bioTAGS.alphabetIsMultichar, 'true');
  alignedSequencesCol.semType = DG.SEMTYPE.MACROMOLECULE;

  return alignedSequencesCol;
}

async function requestAlignedObjects(dockerfileId: string, body: PepseaBodyUnit[], method: string, gapOpen: number,
  gapExtend: number): Promise<PepseaResponse> {
  const params = {
    method: 'POST',
    headers: {'Accept': 'application/json', 'Content-Type': 'application/json'},
    body: JSON.stringify(body),
  };
  const path = `/align?method=${method}&gap_open=${gapOpen}&gap_extend=${gapExtend}`;
  const response = await grok.dapi.docker.dockerContainers.request(dockerfileId, path, params);
  return JSON.parse(response ?? '{}');
}
