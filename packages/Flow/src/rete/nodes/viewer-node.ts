/** Viewer nodes — build a Datagrok viewer from a wired table via `table.plot.fromType(type, {})`
 *  + `setOptions(look)`. NB `fromType` is async and passing a look straight to it throws. */

import {ClassicPreset} from 'rete';
import {FlowNode} from '../scheme';
import {getSocket} from '../sockets';
import {categoricalColor, CAT} from '../../types/type-map';

const COLOR_VIEWER = categoricalColor(CAT.cyan);

/** One option surfaced in the property panel; `kind` drives the editor. */
export interface ViewerOption {
  /** The look key written into `viewerLook` (e.g. `xColumnName`, `title`). */
  key: string;
  label: string;
  kind: 'column' | 'string';
}

export interface ViewerSpec {
  /** The DG viewer type passed to `plot.fromType` (a `DG.VIEWER` value). */
  type: string;
  /** Node title / browser label. */
  label: string;
  options: ViewerOption[];
}

// "Chart title", not "Title" — the panel header already has the node's Title field.
const TITLE: ViewerOption = {key: 'title', label: 'Chart title', kind: 'string'};

/** Curated core viewers with a few verified options; the full look is reached via "Edit settings". */
export const CORE_VIEWER_SPECS: ViewerSpec[] = [
  {type: 'Scatter plot', label: 'Scatter Plot', options: [
    {key: 'xColumnName', label: 'X', kind: 'column'},
    {key: 'yColumnName', label: 'Y', kind: 'column'},
    {key: 'colorColumnName', label: 'Color', kind: 'column'},
    {key: 'sizeColumnName', label: 'Size', kind: 'column'},
    TITLE,
  ]},
  {type: 'Histogram', label: 'Histogram', options: [
    {key: 'valueColumnName', label: 'Value', kind: 'column'}, TITLE,
  ]},
  {type: 'Line chart', label: 'Line Chart', options: [
    {key: 'xColumnName', label: 'X', kind: 'column'}, TITLE,
  ]},
  {type: 'Bar chart', label: 'Bar Chart', options: [
    {key: 'splitColumnName', label: 'Split by', kind: 'column'}, TITLE,
  ]},
  {type: 'Pie chart', label: 'Pie Chart', options: [TITLE]},
  {type: 'Box plot', label: 'Box Plot', options: [
    {key: 'valueColumnName', label: 'Value', kind: 'column'}, TITLE,
  ]},
  {type: 'Heat map', label: 'Heat Map', options: [TITLE]},
  {type: 'Grid', label: 'Grid', options: [TITLE]},
  {type: 'Trellis plot', label: 'Trellis Plot', options: [TITLE]},
  {type: 'Network diagram', label: 'Network Diagram', options: [TITLE]},
];

export const VIEWER_TYPE_PREFIX = 'Viewers/';

export function genericViewerSpec(type: string, label: string): ViewerSpec {
  return {type, label, options: [TITLE]};
}

export class ViewerNode extends FlowNode {
  constructor(spec: ViewerSpec) {
    super(spec.label);
    this.dgNodeType = 'utility';
    this.dgOutputType = 'viewer';
    // `viewerOptionSpecs` is serialized so the panel can render the fields after a reload.
    this.properties = {
      viewerType: spec.type,
      viewerLook: {} as Record<string, unknown>,
      viewerOptionSpecs: spec.options,
    };
    (this as unknown as {color: string}).color = COLOR_VIEWER;
    this.addInput('table', new ClassicPreset.Input(getSocket('dataframe'), 'table'));
    // Each column option is also a connectable input socket; a connected column wins at compile time.
    for (const o of spec.options) {
      if (o.kind === 'column')
        this.addInput(o.key, new ClassicPreset.Input(getSocket('column'), o.label));
    }
    this.addOutput('viewer', new ClassicPreset.Output(getSocket('viewer'), 'viewer'));
    this.requiredInputs = ['table'];
  }
}
