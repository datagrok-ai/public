import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as chem from 'datagrok-api/src/chem';

let sketcherId = 0;

export class OpenChemLibSketcher extends chem.SketcherBase {
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
      // @ts-ignore
      that._sketcher = OCL.StructureEditor.createSVGEditor(id, true, 1);
      that._sketcher.setChangeListenerCallback(function() {
        that.onChanged.next(null);
      });
    }, 100);
  }
}