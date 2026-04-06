import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PrismFile, PrismSheet, PrismAnalysis} from './prism-types';
import {prismSheetToDataFrame, prismAnalysisToDataFrame} from './prism-to-dataframe';
import {createGraphViewer} from './prism-graph-viewer';


/** Builds a tree-based preview for a parsed .prism file. */
export function buildPrismView(prismFile: PrismFile): HTMLElement {
  const tree = DG.TreeViewGroup.tree();
  const contentPanel = ui.div([], {style: {flex: '1', overflow: 'auto', padding: '8px'}});

  function showContent(element: HTMLElement): void {
    contentPanel.innerHTML = '';
    contentPanel.appendChild(element);
  }

  // Data sheets group
  if (prismFile.sheets.length > 0) {
    const dataGroup = tree.group('Data Sheets', null, true);
    for (const sheet of prismFile.sheets) {
      const label = sheet.title || sheet.id;
      const formatTag = sheet.dataFormat ? ` (${sheet.dataFormat})` : '';
      const node = dataGroup.item(label + formatTag, sheet);
      node.onSelected.subscribe(() => showSheetContent(sheet, contentPanel));
    }
  }

  // Analyses group
  if (prismFile.analyses.length > 0) {
    const analysesGroup = tree.group('Analyses', null, true);
    for (const analysis of prismFile.analyses) {
      const node = analysesGroup.item(analysis.title || analysis.id, analysis);
      node.onSelected.subscribe(() => showAnalysisContent(analysis, contentPanel));
    }
  }

  // Show first sheet by default
  if (prismFile.sheets.length > 0)
    showSheetContent(prismFile.sheets[0], contentPanel);

  tree.root.style.minWidth = '200px';
  tree.root.style.maxWidth = '300px';
  tree.root.style.overflow = 'auto';
  tree.root.style.borderRight = '1px solid var(--grey-2)';

  return ui.splitH([tree.root, contentPanel], null, true);
}


/** Displays a data sheet in the content panel. */
function showSheetContent(sheet: PrismSheet, panel: HTMLElement): void {
  panel.innerHTML = '';

  const df = prismSheetToDataFrame(sheet);
  const grid = DG.Viewer.grid(df);
  grid.root.style.width = '100%';
  grid.root.style.height = '50%';
  panel.appendChild(grid.root);

  // Add a chart viewer for visual representation
  const viewer = createGraphViewer(df, sheet);
  if (viewer) {
    viewer.root.style.width = '100%';
    viewer.root.style.height = '50%';
    panel.appendChild(viewer.root);
  }
}


/** Displays analysis results in the content panel. */
function showAnalysisContent(analysis: PrismAnalysis, panel: HTMLElement): void {
  panel.innerHTML = '';

  const df = prismAnalysisToDataFrame(analysis);
  const grid = DG.Viewer.grid(df);
  grid.root.style.width = '100%';
  grid.root.style.height = '100%';
  panel.appendChild(grid.root);
}
