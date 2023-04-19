import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {RcsbFvStructureInterface} from '../../libs/rcsb-saguaro-3d/RcsbFvStructure/RcsbFvStructure';
import {CustomViewInterface} from '../../libs/rcsb-saguaro-3d/RcsbFvSequence/SequenceViews/CustomView';
import {RcsbFv3DAbstractInterface, RcsbFv3dBase} from './molstar-twin-viewer-base';

export interface RcsbFv3DCustomInterface extends RcsbFv3DAbstractInterface {
  structurePanelConfig: RcsbFvStructureInterface;
  sequencePanelConfig: {
    config: CustomViewInterface;
    title?: string;
    subtitle?: string;
  };
}


export class RcsbFv3DCustom extends RcsbFv3dBase {

  constructor(config?: RcsbFv3DCustomInterface) {
    super(config);
  }

  init(config: RcsbFv3DCustomInterface) {
    this.elementId = config.elementId ?? 'RcsbFv3D_mainDiv_' + Math.random().toString(36).substr(2);
    this.structureConfig = config.structurePanelConfig;
    this.sequenceConfig = {
      ...config.sequencePanelConfig,
      type: 'custom'
    };
    this.cssConfig = config.cssConfig;
  }
}

export class BiostructureAndTrackViewer extends DG.JsViewer {

}
