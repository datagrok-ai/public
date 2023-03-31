/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {NOTATION, TAGS as bioTAGS, ALIGNMENT, ALPHABET} from '@datagrok-libraries/bio/src/utils/macromolecule';
import * as C from './constants';

export const pepseaMethods = ['mafft --auto', 'mafft', 'linsi', 'ginsi', 'einsi', 'fftns', 'fftnsi', 'nwns', 'nwnsi'];
const alignmentObjectMetaKeys = ['AlignedSeq', 'AlignedSubpeptide', 'HELM', 'ID', 'PolymerID'];
type PepseaRepsonse = {
  Alignment: {
    PolymerID: string, AlignedSubpeptide: string, HELM: string, ID: string, AlignedSeq: string, [key: string]: string,
  }[],
  AlignmentScore: {[key: string]: number | null},
};
type PepseaBodyUnit = {ID: string, HELM: string};

export async function runPepsea(srcCol: DG.Column<string>, unUsedName: string,
  method: typeof pepseaMethods[number] = 'ginsi', gapOpen: number = 1.53, gapExtend: number = 0.0,
  clustersCol: DG.Column<string | number> | null = null,
  ): Promise<DG.Column<string>> {
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

  //@ts-ignore: this is a temporary workaround for the issue with docker containers. This will be fixed in 1.14.0
  const pepseaContainer = await (grok.dapi.docker !== undefined ? grok.dapi.docker.dockerContainers : grok.dapi.dockerfiles).filter('bio').first();
  const alignedSequences: string[] = new Array(peptideCount);
  for (const body of bodies) { // getting aligned sequences for each cluster
    const alignedObject = await requestAlignedObjects(pepseaContainer.id, body, method, gapOpen, gapExtend);
    const alignments = alignedObject.Alignment;

    for (const alignment of alignments) {  // filling alignedSequencesCol
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
  alignedSequencesCol.semType = DG.SEMTYPE.MACROMOLECULE;

  return alignedSequencesCol;
}

async function requestAlignedObjects(dockerfileId: string, body: PepseaBodyUnit[], method: string, gapOpen: number,
  gapExtend: number): Promise<PepseaRepsonse> {
  const params = {
    method: 'POST',
    headers: {'Accept': 'application/json', 'Content-Type': 'application/json'},
    body: JSON.stringify(body),
  };
  const path = `/align?method=${method}&gap_open=${gapOpen}&gap_extend=${gapExtend}`;
  //@ts-ignore: this is a temporary workaround for the issue with docker containers
  const response = await (grok.dapi.docker !== undefined ? grok.dapi.docker.dockerContainers : grok.dapi.dockerfiles)
    .request(dockerfileId, path, params);
  return JSON.parse(response ?? '{}');
}
