import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as NGL from 'NGL';

export async function viewNglUI(fileContent: string): Promise<void> {
  const view = grok.shell.newView('NGL');
  const host = ui.div([], 'd4-ngl-viewer');
  const stage = new NGL.Stage(host);
  const blob = new Blob([fileContent], {type: 'application/octet-binary'});
  await stage.loadFile(blob, {defaultRepresentation: true});
  handleResize(host, stage);
}

/**
 * @param {any} file
 * @return {DG.View}
 */
export function previewNglUI(file: DG.FileInfo): { view: DG.View, loadingPromise: Promise<void> } {
  const view = DG.View.create();
  const host = ui.div([], 'd4-ngl-viewer');
  const stage = new NGL.Stage(host);

  const loadingPromise = new Promise<void>(async (resolve, reject) => {
    try {
      const data = await file.readAsBytes();
      const blob = new Blob([data], {type: 'application/octet-binary'});
      await stage.loadFile(blob, {defaultRepresentation: true, ext: file.extension});
    } catch (err: any) {
      reject(err);
    }
  });

  handleResize(host, stage);
  view.append(host);
  return {view, loadingPromise};
}

export function nglWidgetUI(pdbId: string): DG.Widget {
  const host = ui.div([], {classes: 'd4-ngl-viewer', style: {width: '100%'}});
  const stage = new NGL.Stage(host);

  //const pdbIdPath: strint = `rcsb://${pdbId}`;
  // Link `rcsb://${pdbId}` causes CORS error
  const pdbIdPath: string = `https://files.rcsb.org/download/${pdbId}.cif`;

  stage.loadFile(pdbIdPath, {defaultRepresentation: true})
    .then(() => {
      stage.compList[0].autoView();
      handleResize(host, stage);
    })
    .catch((reason) => {
      grok.shell.warning(`Can\'t load structure for PDB ID '${pdbId}' error: ${reason}`);
    });
  return DG.Widget.fromRoot(host);
}

function handleResize(host: HTMLDivElement, stage: NGL.Stage) {
  const canvas = stage.viewer.renderer.domElement;

  function resize() {
    canvas!.width = Math.floor(canvas!.clientWidth * window.devicePixelRatio);
    canvas!.height = Math.floor(canvas!.clientHeight * window.devicePixelRatio);
    stage.handleResize();
  }

  ui.onSizeChanged(host).subscribe((_) => resize());
  resize();
}
