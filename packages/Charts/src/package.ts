import {_MultiPlotViewer} from './package.g';
import {_WordCloudViewer} from './package.g';
import {_TreeViewer} from './package.g';
import {_SurfacePlot} from './package.g';
import {_SunburstViewer} from './package.g';
import {_SankeyViewer} from './package.g';
import {_RadarViewer} from './package.g';
import {_GroupAnalysisViewer} from './package.g';
import {_GlobeViewer} from './package.g';
import {_ChordViewer} from './package.g';
import * as DG from 'datagrok-api/dg';
import {FlagCellRenderer} from './renderers/flag-cell-renderer';
import {viewerDemo} from './demos/demo';
import * as grok from 'datagrok-api/grok';
import {TimelinesViewer} from './viewers/timelines/timelines-viewer';

export * from './package.g';
export const _package = new DG.Package();

export class PackageFunctions {
  @grok.decorators.func({
    name: 'Timelines',
    description: 'Creates a timelines viewer',
    meta: {
      showInGallery: 'false',
      icon: 'icons/timelines-viewer.svg',
      role: 'viewer',
    },
    outputs: [{name: 'result', type: 'viewer'}],
  })
  static timelinesViewer(): TimelinesViewer {
    return new TimelinesViewer();
  }

  @grok.decorators.func({
    outputs: [{name: 'result', type: 'grid_cell_renderer'}],
    meta: {role: 'cellRenderer'},
  })
  static flagCellRenderer() {
    return new FlagCellRenderer();
  }

  @grok.decorators.demo({
    description: 'Chord viewer visualizes weighted relationships between several entities',
    meta: { demoPath: 'Visualization | General | Chord', demoWait: '4000' },
  })
  static async chordViewerDemo() {
    await viewerDemo('Chord');
  }

  // Globe demo is commented out in original, so skip

  @grok.decorators.demo({
    description: 'Radar viewer is used on multivariate data to plot groups of values over several common variables',
    meta: { demoPath: 'Visualization | General | Radar', demoWait: '4000' },
  })
  static async radarViewerDemo() {
    await viewerDemo('Radar');
  }

  @grok.decorators.demo({
    description: 'Sankey viewer depicts a flow from one set of values to another',
    meta: { demoPath: 'Visualization | General | Sankey', demoWait: '4000' },
  })
  static async sankeyViewerDemo() {
    await viewerDemo('Sankey');
  }

  @grok.decorators.demo({
    description: 'Sunburst viewer displays hierarchical data',
    meta: { demoPath: 'Visualization | General | Sunburst', demoWait: '4000' },
  })
  static async sunburstViewerDemo() {
    await viewerDemo('Sunburst');
  }

  @grok.decorators.demo({
    description: 'Surface plot viewer displays a set of three-dimensional data as a mesh surface',
    meta: { demoPath: 'Visualization | General | Surface Plot', demoWait: '4000' },
  })
  static async surfacePlotDemo() {
    await viewerDemo('SurfacePlot');
  }

  @grok.decorators.demo({
    description: 'Timelines viewer displays the flow of events over time',
    meta: { demoPath: 'Visualization | General | Timelines', demoWait: '4000' },
  })
  static async timelinesViewerDemo() {
    await viewerDemo('Timelines', {lineWidth: 4, markerPosition: 'above main line'});
  }

  @grok.decorators.demo({
    description: 'Tree viewer visualizes hierarchical data by categories',
    meta: { demoPath: 'Visualization | Data Flow and Hierarchy | Tree', demoWait: '4000' },
  })
  static async treeViewerDemo() {
    await viewerDemo('Tree', {left: '40px', right: '75px'});
  }

  @grok.decorators.demo({
    description: 'Word Cloud viewer visualizes unstructured text data',
    meta: { demoPath: 'Visualization | General | Word Cloud', demoWait: '4000' },
  })
  static async wordCloudViewerDemo() {
    await viewerDemo(DG.VIEWER.WORD_CLOUD, {drawOutOfBound: false});
  }
}

export {_ChordViewer};
export {_GlobeViewer};
export {_GroupAnalysisViewer};
export {_RadarViewer};
export {_SankeyViewer};
export {_SunburstViewer};
export {_SurfacePlot};
export {_TreeViewer};
export {_WordCloudViewer};
export {_MultiPlotViewer};
