import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';

let sketcherId = 0;

export class OpenChemLibSketcher extends DG.chem.SketcherBase {
  _sketcher: any;

  constructor() {
    super();

    let id = `ocl-sketcher-${sketcherId++}`;
    let host = ui.div([], {
      id: id,
      style: {width: '600px', height: '600px'}
    });
    this.root.appendChild(host);

    let that = this;
    setTimeout(function() {
      that._sketcher = OCL.StructureEditor.createSVGEditor(id, 1);
      that._sketcher.setChangeListenerCallback(function() {
        that.onChanged.next(null);
      });
    }, 100);
  }
}
