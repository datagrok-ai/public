/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {NOTATION, TAGS as bioTAGS, ALIGNMENT, ALPHABET} from '@datagrok-libraries/bio/src/utils/macromolecule';

export const pepseaMethods = ['mafft --auto', 'mafft', 'linsi', 'ginsi', 'einsi', 'fftns', 'fftnsi', 'nwns', 'nwnsi'];
const alignmentObjectMetaKeys = ['AlignedSeq', 'AlignedSubpeptide', 'HELM', 'ID', 'PolymerID'];
type PepseaRepsonse = {
  Alignment: {
    PolymerID: string, AlignedSubpeptide: string, HELM: string, ID: string, AlignedSeq: string, [key: string]: string,
  }[],
  AlignmentScore: {[key: string]: number | null},
};
type PepseaBodyUnit = {ID: string, HELM: string};

export function pepseaDialog(): void {
  const table = grok.shell.t;
  const colInput = ui.columnInput('Sequences', table, table.columns.bySemType('Macromolecule'),
    (col: DG.Column<string>) => {
      if (col.getTag(DG.TAGS.UNITS) != 'helm' || col.semType != DG.SEMTYPE.MACROMOLECULE)
        grok.shell.info('Sequence column must contain sequences in HELM notation!');
    });
  colInput.setTooltip('Sequences column to use for alignment');
  const methodInput = ui.choiceInput('Method', 'ginsi', pepseaMethods);
  methodInput.setTooltip('Alignment method');
  const gapOpenInput = ui.floatInput('Gap open', 1.53);
  gapOpenInput.setTooltip('Gap opening penalty at group-to-group alignment');
  const gapExtendInput = ui.floatInput('Gap extend', 0.0);
  gapExtendInput.setTooltip('Gap extension penalty to skip the alignment');
  const clusterColInput = ui.columnInput('Clusters', table, null);
  clusterColInput.setTooltip('Clusters column to perform in-cluster alignment');

  ui.dialog('PepSeA Multiple Sequence Alignment')
    .add(colInput)
    .add(methodInput)
    .add(gapOpenInput)
    .add(gapExtendInput)
    .add(clusterColInput)
    .onOK(async () => {
      const progress = DG.TaskBarProgressIndicator.create('Performing MSA...');
      try {
        const msaCol = await runPepsea(colInput.value!, methodInput.stringValue, gapOpenInput.value,
          gapExtendInput.value, clusterColInput.value);
        table.columns.add(msaCol);

        // This call is required to enable cell renderer activation
        await grok.data.detectSemanticTypes(table);
      } catch (e) {
        grok.shell.error('PepseaMsaError: Could not perform alignment. See console for details.');
        console.error(e);
      }
      progress.close();
    })
    .show();
}

export async function runPepsea(col: DG.Column<string>, method: typeof pepseaMethods[number] = 'ginsi',
  gapOpen: number | null = 1.53, gapExtend: number | null = 0.0, clustersCol: DG.Column<string | number> | null = null,
  ): Promise<DG.Column<string>> {
  const peptideCount = col.length;
  gapOpen ??= 1.53;
  gapExtend ??= 0.0;
  clustersCol ??= DG.Column.int('Clusters', peptideCount).init(0);
  clustersCol = (clustersCol.type !== DG.TYPE.STRING ? clustersCol.convertTo(DG.TYPE.STRING) :
    clustersCol) as DG.Column<string>;

  const clusters = clustersCol.categories;
  const bodies: PepseaBodyUnit[][] = new Array(clusters.length);

  // Grouping data by clusters
  for (let rowIndex = 0; rowIndex < peptideCount; ++rowIndex) {
    const cluster = clustersCol.get(rowIndex) as string;
    if (cluster == '')
      continue;

    const clusterId = clusters.indexOf(cluster);
    const helmSeq = col.get(rowIndex);
    if (helmSeq)
      (bodies[clusterId] ??= []).push({ID: rowIndex.toString(), HELM: helmSeq});
  }

  const dockerfileId = (await grok.dapi.dockerfiles.filter('bio').first()).id;

  const alignedSequences: string[] = new Array(peptideCount);
  for (const body of bodies) { // getting aligned sequences for each cluster
    const alignedObject = await requestAlignedObjects(dockerfileId, body, method, gapOpen, gapExtend);
    const alignments = alignedObject.Alignment;

    for (const alignment of alignments) {  // filling alignedSequencesCol
      alignedSequences[parseInt(alignment.ID)] = Object.entries(alignment)
        .filter((v) => !alignmentObjectMetaKeys.includes(v[0]))
        .map((v) => v[1])
        .join('.');
    }
  }

  const newColName = col.dataFrame.columns.getUnusedName(`msa(${col.name})`);
  const alignedSequencesCol: DG.Column<string> = DG.Column.fromStrings(newColName, alignedSequences);
  alignedSequencesCol.setTag(DG.TAGS.UNITS, NOTATION.SEPARATOR);
  alignedSequencesCol.setTag(bioTAGS.aligned, ALIGNMENT.SEQ_MSA);
  alignedSequencesCol.setTag(bioTAGS.alphabet, ALPHABET.PT);
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
  const response = await grok.dapi.dockerfiles.request(dockerfileId, path, params);
  return JSON.parse(response ?? '{}');
}
