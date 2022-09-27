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
type pepseaResonse = {
  Alignment: {
    PolymerID: string, AlignedSubpeptide: string, HELM: string, ID: string, AlignedSeq: string, [key: string]: string,
  }[],
  AlignmentScore: {[key: string]: number | null},
};

//top-menu: Bio | PepSeA MSA...
export function pepseaMSA(): void {
  const table = grok.shell.t;
  const colInput = ui.columnInput('Sequences', table, null, (col: DG.Column<string>) => {
    if (col.getTag(DG.TAGS.UNITS) != 'helm' || col.semType != DG.SEMTYPE.MACROMOLECULE)
      grok.shell.info('Sequence column must have be of HELM notation!')
  });
  const methodInput = ui.choiceInput('Method', 'ginsi', methods);
  const gapOpenInput = ui.floatInput('Gap open', 1.53);
  const gapExtendInput = ui.floatInput('Gap extend', 0.0);

  ui.dialog('PepSeA Multiple Sequence Alignment')
    .add(colInput)
    .add(methodInput)
    .add(gapOpenInput)
    .add(gapExtendInput)
    .onOK(async () => {
      const progress = DG.TaskBarProgressIndicator.create('Performing MSA...');
      await perfromPepseaMSA(colInput.value!, methodInput.stringValue, gapOpenInput.value, gapExtendInput.value);
      progress.close();
    })
    .show();
}

async function perfromPepseaMSA(col: DG.Column<string>, method: string, gapOpen: number | null,
  gapExtend: number | null): Promise<void> {
  gapOpen ??= 1.53;
  gapExtend ??= 0.0;

  const body: {ID: string, HELM: string}[] = [];
  for (let colIndex = 0; colIndex < col.length; ++colIndex) {
    const helmSeq = col.get(colIndex);
    if (helmSeq)
      body.push({ID: colIndex.toString(), HELM: helmSeq});
  }

  const dockerfileId = (await grok.dapi.dockerfiles.filter('pepsea').first()).id;

  const params = {
    method: 'POST',
    headers: {'Accept': 'application/json', 'Content-Type': 'application/json'},
    body: JSON.stringify(body),
  };
  const response = await grok.dapi.dockerfiles.request(dockerfileId, `/align?method=${method}&gap_open=${gapOpen}&gap_extend=${gapExtend}`, params);
  const alignedObjects: pepseaResonse = JSON.parse(response ?? '{}');
  const alignments = alignedObjects.Alignment;

  const alignedSeqCol = col.dataFrame.columns.addNewString('Aligned');
  for (const alignment of alignments)
    alignedSeqCol.set(parseInt(alignment.ID), alignment.AlignedSeq);

  await grok.data.detectSemanticTypes(col.dataFrame);

  grok.shell.info(`Alignment score: ${Object.values(alignedObjects.AlignmentScore)[0]}`);
}
