import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function flagCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/charts:FlagCellRenderer', {});
  }

  //Chord viewer visualizes weighted relationships between several entities
  export async function chordViewerDemo(): Promise<any> {
    return await grok.functions.call('@datagrok/charts:ChordViewerDemo', {});
  }

  //Radar viewer is used on multivariate data to plot groups of values over several common variables
  export async function radarViewerDemo(): Promise<any> {
    return await grok.functions.call('@datagrok/charts:RadarViewerDemo', {});
  }

  //Sankey viewer depicts a flow from one set of values to another
  export async function sankeyViewerDemo(): Promise<any> {
    return await grok.functions.call('@datagrok/charts:SankeyViewerDemo', {});
  }

  //Sunburst viewer displays hierarchical data
  export async function sunburstViewerDemo(): Promise<any> {
    return await grok.functions.call('@datagrok/charts:SunburstViewerDemo', {});
  }

  //Surface plot viewer displays a set of three-dimensional data as a mesh surface
  export async function surfacePlotDemo(): Promise<any> {
    return await grok.functions.call('@datagrok/charts:SurfacePlotDemo', {});
  }

  //Timelines viewer displays the flow of events over time
  export async function timelinesViewerDemo(): Promise<any> {
    return await grok.functions.call('@datagrok/charts:TimelinesViewerDemo', {});
  }

  //Tree viewer visualizes hierarchical data by categories
  export async function treeViewerDemo(): Promise<any> {
    return await grok.functions.call('@datagrok/charts:TreeViewerDemo', {});
  }

  //Word Cloud viewer visualizes unstructured text data
  export async function wordCloudViewerDemo(): Promise<any> {
    return await grok.functions.call('@datagrok/charts:WordCloudViewerDemo', {});
  }

  //Creates a chord viewer
  export async function chordViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/charts:ChordViewer', {});
  }

  //Creates a globe viewer
  export async function globeViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/charts:GlobeViewer', {});
  }

  //Creates a group analysis viewer
  export async function groupAnalysisViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/charts:GroupAnalysisViewer', {});
  }

  //Creates a multiplot viewer
  export async function multiPlotViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/charts:MultiPlotViewer', {});
  }

  //Creates a radar viewer
  export async function radarViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/charts:RadarViewer', {});
  }

  //Creates a sankey viewer
  export async function sankeyViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/charts:SankeyViewer', {});
  }

  //Creates a sunburst viewer
  export async function sunburstViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/charts:SunburstViewer', {});
  }

  //Creates a surface plot viewer
  export async function surfacePlot(): Promise<any> {
    return await grok.functions.call('@datagrok/charts:SurfacePlot', {});
  }

  //Creates a timelines viewer
  export async function timelinesViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/charts:TimelinesViewer', {});
  }

  //Creates a tree viewer
  export async function treeViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/charts:TreeViewer', {});
  }

  //Creates a word cloud viewer
  export async function wordCloudViewer(): Promise<any> {
    return await grok.functions.call('@datagrok/charts:WordCloudViewer', {});
  }
}
