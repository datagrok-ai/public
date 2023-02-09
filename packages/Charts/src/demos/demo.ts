import * as grok from 'datagrok-api/grok';
import {energyUK, demog, adverseEvents, earthquakes} from './demo-data';

export class Demo {
  static chordViewerDemo() {
    const tableView = grok.shell.addTableView(energyUK);
    tableView.addViewer('ChordViewer');
  }

  static globeViewerDemo() {
    const tableView = grok.shell.addTableView(earthquakes);
    // timeout is for the table view to fully load and then start the globe viewer
    setTimeout(() => {
      tableView.addViewer('GlobeViewer');
    }, 300);
  }

  static groupAnalysisViewerDemo() {
    const tableView = grok.shell.addTableView(demog);
    tableView.addViewer('GroupAnalysisViewer');
  }

  static radarViewerDemo() {
    const tableView = grok.shell.addTableView(demog);
    tableView.addViewer('RadarViewer');
  }

  static sankeyViewerDemo() {
    const tableView = grok.shell.addTableView(energyUK);
    tableView.addViewer('SankeyViewer');
  }

  static sunburstViewerDemo() {
    const tableView = grok.shell.addTableView(demog);
    tableView.addViewer('SunburstViewer');
  }

  static surfacePlotDemo() {
    const tableView = grok.shell.addTableView(demog);
    tableView.addViewer('SurfacePlot');
  }

  static timelinesViewerDemo() {
    const tableView = grok.shell.addTableView(adverseEvents);
    tableView.addViewer('TimelinesViewer');
  }

  static treeViewerDemo() {
    const tableView = grok.shell.addTableView(demog);
    tableView.addViewer('TreeViewer');
  }

  static wordCloudViewerDemo() {
    const tableView = grok.shell.addTableView(energyUK);
    tableView.addViewer('WordCloudViewer');
  }
}
