import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';


export class Demo {
  static async chordViewerDemo() {
    const df = await grok.data.getDemoTable('energy_uk.csv');
    const tableView = grok.shell.addTableView(df);
    tableView.addViewer('ChordViewer');
  }

  static async globeViewerDemo() {
    const df = await grok.data.getDemoTable('geo/earthquakes.csv');
    const tableView = grok.shell.addTableView(df);
    // timeout is for the table view to fully load, activate detectors and then start the globe viewer
    setTimeout(() => {
      tableView.addViewer('GlobeViewer');
    }, 300);
  }

  static async groupAnalysisViewerDemo() {
    const df = await grok.data.getDemoTable('demog.csv');
    const tableView = grok.shell.addTableView(df);
    tableView.addViewer('GroupAnalysisViewer');
  }

  static async radarViewerDemo() {
    const df = await grok.data.getDemoTable('demog.csv');
    const tableView = grok.shell.addTableView(df);
    tableView.addViewer('RadarViewer');
  }

  static async sankeyViewerDemo() {
    const df = await grok.data.getDemoTable('energy_uk.csv');
    const tableView = grok.shell.addTableView(df);
    tableView.addViewer('SankeyViewer');
  }

  static async sunburstViewerDemo() {
    const df = await grok.data.getDemoTable('demog.csv');
    const tableView = grok.shell.addTableView(df);
    tableView.addViewer('SunburstViewer');
  }

  static async surfacePlotDemo() {
    const df = await grok.data.getDemoTable('demog.csv');
    const tableView = grok.shell.addTableView(df);
    tableView.addViewer('SurfacePlot');
  }

  static timelinesViewerDemo() {
    const df = DG.DataFrame.fromCsv(
      `USUBJID, AESTDY, AEENDY, SEX, AGE
      s1, 10/02/2018, 10/09/2018, F, 48
      s2, 10/04/2018, 10/07/2018, M, 51
      s3, 10/02/2018, 10/05/2018, F, 39
      s4, 10/07/2018, 10/08/2018, M, 43`);
    const tableView = grok.shell.addTableView(df);
    tableView.addViewer('TimelinesViewer');
  }

  static async treeViewerDemo() {
    const df = await grok.data.getDemoTable('demog.csv');
    const tableView = grok.shell.addTableView(df);
    tableView.addViewer('TreeViewer');
  }

  static async wordCloudViewerDemo() {
    const df = await grok.data.getDemoTable('word_cloud.csv');
    const tableView = grok.shell.addTableView(df);
    tableView.addViewer('WordCloudViewer');
  }
}
