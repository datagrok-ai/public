import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as NGL from 'NGL';

/** @returns {DG.View} */
export function nglViewUI(file: any): DG.View {
  const view = DG.View.create();
  const host = ui.div([], 'd4-ngl-viewer');
  const stage = new NGL.Stage(host);
  handleResize(host, stage);

  function loadBytes(bytes: any) {
    const blob = new Blob([bytes], {type: 'application/octet-binary'});
    stage.loadFile(blob, {defaultRepresentation: true, ext: file.extension});
  }

  file
    .readAsBytes()
    .then(loadBytes);

  view.append(host);
  return view;
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

function handleResize(host: HTMLDivElement, stage: any) {
  const canvas = host.querySelector('canvas');

  function resize() {
    canvas!.width = Math.floor(canvas!.clientWidth * window.devicePixelRatio);
    canvas!.height = Math.floor(canvas!.clientHeight * window.devicePixelRatio);
    stage.handleResize();
  }

  ui.onSizeChanged(host).subscribe((_) => resize());
  resize();
}
