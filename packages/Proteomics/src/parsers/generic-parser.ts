import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {SEMTYPE} from '../utils/proteomics-types';
import {focusProtein} from '../panels/protein-focus';
import {
  log2TransformColumns,
  copyAsLog2Columns,
  addPrimaryColumnIfNeeded,
  detectLog2Status,
  detectDelimiter,
  autoSuggestProteinIdColumn,
  autoSuggestIntensityColumns,
  autoSuggestGeneNameColumn,
} from './shared-utils';
import {resolveGeneLabels} from '../utils/gene-label-resolver';

/** Opens a file picker and shows the column mapping dialog for generic matrix import. */
export function showGenericImportDialog(): void {
  DG.Utils.openFile({
    accept: '.csv,.tsv,.txt',
    open: (file: File) => {
      const reader = new FileReader();
      reader.onload = () => {
        try {
          const text = reader.result as string;
          const delimiter = detectDelimiter(text);
          const df = DG.DataFrame.fromCsv(text, {delimiter});
          showMappingDialog(df, file.name);
        } catch (e: any) {
          grok.shell.error('Failed to parse file: ' + e.message);
        }
      };
      reader.readAsText(file);
    },
  });
}

/** Shows the column mapping dialog for a parsed DataFrame.
 * @param df - Parsed DataFrame from the imported file
 * @param fileName - Original file name for naming the output DataFrame */
function showMappingDialog(df: DG.DataFrame, fileName: string): void {
  // Auto-suggest columns
  const suggestedProtein = autoSuggestProteinIdColumn(df);
  const suggestedGene = autoSuggestGeneNameColumn(df);
  const suggestedIntensityNames = autoSuggestIntensityColumns(df);

  // Protein ID input -- string columns only
  const proteinIdInput = ui.input.column('Protein ID', {
    table: df,
    filter: (c: DG.Column) => c.type === DG.COLUMN_TYPE.STRING,
    value: suggestedProtein ?? undefined,
  });

  // Gene Name input -- optional, string columns only
  const geneNameInput = ui.input.column('Gene Name (optional)', {
    table: df,
    filter: (c: DG.Column) => c.type === DG.COLUMN_TYPE.STRING,
    value: suggestedGene ?? undefined,
    nullable: true,
  });

  // Intensity columns multi-select — show all numeric columns, pre-select keyword matches
  const numericColNames = df.columns.toList()
    .filter((c) => c.type === DG.COLUMN_TYPE.FLOAT || c.type === DG.COLUMN_TYPE.INT || c.type === DG.COLUMN_TYPE.BIG_INT)
    .map((c) => c.name);
  const intensityColsInput = ui.input.columns('Intensity Columns', {
    table: df,
    available: numericColNames,
    value: suggestedIntensityNames.map((n) => df.col(n)!).filter((c) => c != null),
  });

  // Log2 detection on initial intensity columns
  const initialDetection = detectLog2Status(df, suggestedIntensityNames);

  // Log2 toggle -- ON means "needs transform" (raw intensities)
  const log2Toggle = ui.input.bool('Log2 Transform', {value: !initialDetection.isLog2});

  // Hint label for log2 detection message
  const hintDiv = ui.divText(initialDetection.message);
  hintDiv.style.cssText = 'font-style:italic;color:#888;margin-bottom:8px;font-size:12px;';

  // Preview container
  const previewContainer = document.createElement('div');
  previewContainer.style.cssText = 'border:1px solid #ddd;margin-top:8px;width:100%;';

  /** Updates the preview grid with currently selected columns. */
  function updatePreview(): void {
    const selectedNames: string[] = [];
    const proteinCol = proteinIdInput.value;
    if (proteinCol) selectedNames.push(proteinCol.name);
    const geneCol = geneNameInput.value;
    if (geneCol) selectedNames.push(geneCol.name);
    const intensityCols: DG.Column[] = intensityColsInput.value ?? [];
    const intensityNames = intensityCols.map((c) => c.name);
    selectedNames.push(...intensityNames);

    if (selectedNames.length === 0) {
      previewContainer.innerHTML = '<div style="padding:8px;color:#888">Select columns to preview</div>';
      return;
    }

    try {
      const mask = DG.BitSet.create(df.rowCount, (i) => i < 5);
      const previewDf = df.clone(mask, selectedNames);
      const grid = DG.Viewer.grid(previewDf);
      previewContainer.innerHTML = '';
      const countLabel = ui.divText(`Preview: ${selectedNames.length} columns, 5 rows`);
      countLabel.style.cssText = 'font-size:11px;color:#888;padding:2px 4px;';
      grid.root.style.cssText = 'width:100%;height:180px;';
      previewContainer.appendChild(countLabel);
      previewContainer.appendChild(grid.root);
    } catch {
      previewContainer.innerHTML = '<div style="padding:8px;color:#888">Unable to generate preview</div>';
    }
  }

  /** Updates log2 detection hint based on current intensity columns. */
  function updateLog2Detection(): void {
    const intensityCols: DG.Column[] = intensityColsInput.value ?? [];
    const names = intensityCols.map((c) => c.name);
    if (names.length === 0) {
      hintDiv.textContent = 'Select intensity columns';
      return;
    }
    const detection = detectLog2Status(df, names);
    hintDiv.textContent = detection.message;
    log2Toggle.value = !detection.isLog2;
  }

  // Wire up reactive updates
  proteinIdInput.onChanged.subscribe(() => updatePreview());
  geneNameInput.onChanged.subscribe(() => updatePreview());
  intensityColsInput.onChanged.subscribe(() => {
    updateLog2Detection();
    updatePreview();
  });

  // Initial preview
  updatePreview();

  // Build dialog
  ui.dialog('Import Generic Matrix')
    .add(proteinIdInput)
    .add(geneNameInput)
    .add(intensityColsInput)
    .add(log2Toggle)
    .add(ui.div([hintDiv]))
    .add(previewContainer)
    .onOK(async () => {
      const proteinCol = proteinIdInput.value;
      if (!proteinCol) {
        grok.shell.warning('Please select a Protein ID column');
        return;
      }

      const intensityCols: DG.Column[] = intensityColsInput.value ?? [];
      const intensityNames = intensityCols.map((c) => c.name);
      if (intensityNames.length === 0) {
        grok.shell.warning('Please select at least one intensity column');
        return;
      }

      // Assign semantic types
      proteinCol.semType = SEMTYPE.PROTEIN_ID;

      const geneCol = geneNameInput.value;
      if (geneCol)
        geneCol.semType = SEMTYPE.GENE_SYMBOL;

      // Log2 transform or copy
      if (log2Toggle.value)
        log2TransformColumns(df, intensityNames);
      else
        copyAsLog2Columns(df, intensityNames);

      // Primary columns (splits semicolon-delimited values)
      addPrimaryColumnIfNeeded(df, proteinCol.name, 'Primary Protein ID', SEMTYPE.PROTEIN_ID);
      if (geneCol)
        addPrimaryColumnIfNeeded(df, geneCol.name, 'Primary Gene Name', SEMTYPE.GENE_SYMBOL);

      // Source tag
      df.setTag('proteomics.source', 'generic');

      await resolveGeneLabels(df);

      // Name from file
      df.name = fileName.replace(/\.[^.]+$/, '');

      // Open in table view
      grok.shell.addTableView(df);
      grok.shell.info(`Imported ${df.rowCount} proteins from ${fileName}`);
      focusProtein(df);
    })
    .show();
}
