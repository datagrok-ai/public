/** Viewer nodes — manually-built nodes that create a Datagrok viewer from a
 *  table, bypassing the default viewer *functions* (which need a TableView
 *  lifecycle). At run time a viewer node emits
 *    `let v = await table.plot.fromType('<Type>', {}); v.setOptions(<look>);`
 *  producing a live `DG.Viewer` whose `.root` renders in the bottom preview.
 *
 *  A node exposes a curated handful of the most-needed options (column choices,
 *  title) in the property panel; the **full** look is editable via the preview's
 *  "Edit settings" button (which does `grok.shell.o = viewer` and writes the
 *  changed `getOptions().look` — minus the `#type` tag — back into the node).
 *
 *  Empirically grounded (see the retired viewer-diag): `df.plot.fromType` is
 *  async and returns the viewer; `setOptions(look)` applies a look keyed by
 *  `xColumnName`/`title`/… ; passing a look straight to `fromType` throws. */

import {ClassicPreset} from 'rete';
import {FlowNode} from '../scheme';
import {getSocket} from '../sockets';
import {categoricalColor, CAT} from '../../types/type-map';

const COLOR_VIEWER = categoricalColor(CAT.cyan);

/** One option surfaced as an editable field in the property panel. `kind`
 *  drives the editor (a `column` field is a column-name text box). */
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

const TITLE: ViewerOption = {key: 'title', label: 'Title', kind: 'string'};

/** Curated core viewers (from the `DataFrame.plot` namespace) with a few
 *  verified, high-value options exposed. Everything else is reached via the
 *  live "Edit settings" editor, so a sparse list here is fine. */
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

/** The registered-type-name prefix for viewer nodes. */
export const VIEWER_TYPE_PREFIX = 'Viewers/';

/** A generic spec for a (non-core) discovered viewer — title only; full
 *  settings via the live editor. */
export function genericViewerSpec(type: string, label: string): ViewerSpec {
  return {type, label, options: [TITLE]};
}

export class ViewerNode extends FlowNode {
  constructor(spec: ViewerSpec) {
    super(spec.label);
    this.dgNodeType = 'utility';
    this.dgOutputType = 'viewer';
    // `viewerType` marks this node for the viewer emit path; `viewerLook` holds
    // the accumulated options (look keys, minus `#type`); `viewerOptionSpecs`
    // lets the panel render the exposed fields after a reload (serialized).
    this.properties = {
      viewerType: spec.type,
      viewerLook: {} as Record<string, unknown>,
      viewerOptionSpecs: spec.options,
    };
    (this as unknown as {color: string}).color = COLOR_VIEWER;
    this.addInput('table', new ClassicPreset.Input(getSocket('dataframe'), 'table'));
    // Each `column` option is also a connectable input socket (keyed by its look
    // key, e.g. `xColumnName`): wire a column in, or type a name in the panel.
    // A connected column wins at compile time (its `.name` feeds `setOptions`).
    for (const o of spec.options) {
      if (o.kind === 'column')
        this.addInput(o.key, new ClassicPreset.Input(getSocket('column'), o.label));
    }
    this.addOutput('viewer', new ClassicPreset.Output(getSocket('viewer'), 'viewer'));
    this.requiredInputs = ['table'];
  }
}
