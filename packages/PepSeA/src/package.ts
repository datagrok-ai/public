/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

const methods = ['mafft --auto', 'mafft', 'linsi', 'ginsi', 'einsi', 'fftns', 'fftnsi', 'nwns', 'nwnsi'];
type PepseaRepsonse = {
  Alignment: {
    PolymerID: string, AlignedSubpeptide: string, HELM: string, ID: string, AlignedSeq: string, [key: string]: string,
  }[],
  AlignmentScore: {[key: string]: number | null},
};
type PepseaBodyUnit = {ID: string, HELM: string};

//top-menu: Bio | PepSeA MSA...
export function pepseaMSA(): void {
  const table = grok.shell.t;
  const colInput = ui.columnInput('Sequences', table, table.columns.bySemType('Macromolecule'),
    (col: DG.Column<string>) => {
      if (col.getTag(DG.TAGS.UNITS) != 'helm' || col.semType != DG.SEMTYPE.MACROMOLECULE)
        grok.shell.info('Sequence column must contain sequences in HELM notation!');
    });
  const methodInput = ui.choiceInput('Method', 'ginsi', methods);
  const gapOpenInput = ui.floatInput('Gap open', 1.53);
  const gapExtendInput = ui.floatInput('Gap extend', 0.0);
  const clusterColInput = ui.columnInput('Clusters', table, table.columns.byIndex(1));
  const alignByClusterInput = ui.boolInput('Use clusters?', false,
    () => clusterColInput.root.hidden = !alignByClusterInput.value!);
  alignByClusterInput.fireChanged();

  ui.dialog('PepSeA Multiple Sequence Alignment')
    .add(colInput)
    .add(methodInput)
    .add(gapOpenInput)
    .add(gapExtendInput)
    .add(alignByClusterInput)
    .add(clusterColInput)
    .onOK(async () => {
      const progress = DG.TaskBarProgressIndicator.create('Performing MSA...');
      try {
        await perfromPepseaMSA(colInput.value!, methodInput.stringValue, gapOpenInput.value, gapExtendInput.value,
          alignByClusterInput.value ? clusterColInput.value : null);
      } catch (e) {
        grok.shell.error('PepseaMsaError: Could not perform alignment. See console for details.');
        console.error(e);
      }
      progress.close();
    })
    .show();
}

async function perfromPepseaMSA(col: DG.Column<string>, method: string, gapOpen: number | null,
  gapExtend: number | null, clustersCol: DG.Column<string | number> | null): Promise<void> {
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
    if (cluster === '')
      continue;

    const clusterId = clusters.indexOf(cluster);
    const helmSeq = col.get(rowIndex);
    if (helmSeq)
      (bodies[clusterId] ??= []).push({ID: rowIndex.toString(), HELM: helmSeq});
  }

  const dockerfileId = (await grok.dapi.dockerfiles.filter('pepsea').first()).id;
  const alignedSequencesCol = DG.Column.string('Aligned', peptideCount);

  for (const body of bodies) { // getting aligned sequences for each cluster
    const alignedObject = await requestAlignedObjects(dockerfileId, body, method, gapOpen, gapExtend);
    const alignments = alignedObject.Alignment;

    for (const alignment of alignments) // filling alignedSequencesCol
      alignedSequencesCol.set(parseInt(alignment.ID), alignment.AlignedSubpeptide);
  }

  const semType = await grok.functions.call('Bio:detectMacromolecule', {col: alignedSequencesCol}) as string;
  if (semType)
    alignedSequencesCol.semType = semType;

  col.dataFrame.columns.add(alignedSequencesCol);
}

async function requestAlignedObjects(dockerfileId: string, body: PepseaBodyUnit[], method: string,
  gapOpen: number | null, gapExtend: number | null): Promise<PepseaRepsonse> {
  const params = {
    method: 'POST',
    headers: {'Accept': 'application/json', 'Content-Type': 'application/json'},
    body: JSON.stringify(body),
  };
  const path = `/align?method=${method}&gap_open=${gapOpen}&gap_extend=${gapExtend}`;
  const response = await grok.dapi.dockerfiles.request(dockerfileId, path, params);
  return JSON.parse(response ?? '{}');
}
