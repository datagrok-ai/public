import * as grok from 'datagrok-api/grok';

import {energyUK, demog, adverseEvents, earthquakes} from './demo-data';

export class Demo {
  static chordViewerDemo() {
    const tableView = grok.shell.addTableView(energyUK);
    tableView.addViewer('ChordViewer');
  }

  static globeViewerDemo() {
    const tableView = grok.shell.addTableView(earthquakes);
    tableView.addViewer('GlobeViewer');
  }

  static groupAnalysisViewerDemo() {

  }

  static radarViewerDemo() {

  }

  static sankeyViewerDemo() {

  }

  static sunburstViewerDemo() {

  }

  static surfacePlotViewerDemo() {

  }

  static timelinesViewerDemo() {

  }

  static treeViewerDemo() {

  }

  static wordCloudViewerDemo() {

  }
}
