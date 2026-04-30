/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {HitDesignApp} from '../hit-design-app';
import {HitAppBase} from '../hit-app-base';
import {HitDesignMergeClashStrategy, HitDesignMergeConfig, HitDesignMergeMode} from '../types';
import {guessMoleculeColumn, guessVidColumn, isLocalUploadFileInfo, mergeIntoCampaign} from '../utils/merge-table';
import {_package} from '../../package';

// Mirrors the drag-n-drop path in HitDesignApp.handleUploadingDataframe: csv goes through
// DG.DataFrame.fromCsv, anything else is delegated to a registered file-handler function.
// HitAppBase.molFileExtReaders is statically filtered to only those extensions for which
// the platform actually exposes a reader, so this list reflects what we can really parse.
const ALLOWED_EXTENSIONS: string[] = ['csv', ...HitAppBase.molFileExtReaders.map((r) => r.ext)];

function getExtension(fi: DG.FileInfo): string {
  const ext = fi.extension?.toLowerCase();
  if (ext) return ext;
  const name = fi.name?.toLowerCase() ?? '';
  const dot = name.lastIndexOf('.');
  return dot >= 0 ? name.slice(dot + 1) : '';
}

async function readFileInfoToDataFrame(fi: DG.FileInfo): Promise<DG.DataFrame | null> {
  const ext = getExtension(fi);
  if (ext === 'csv') {
    const text = await fi.readAsString();
    return DG.DataFrame.fromCsv(text);
  }
  const reader = HitAppBase.molFileExtReaders.find((r) => r.ext === ext);
  if (!reader) return null;
  const bytes = await fi.readAsBytes();
  if (!bytes) return null;
  const bytesArg = reader.handlerFunc.inputs[0].name;
  const result: any = await reader.handlerFunc.apply({[bytesArg]: bytes});
  return Array.isArray(result) ? result[0] : result;
}

const MODE_LABELS: Record<HitDesignMergeMode, string> = {
  smiles: 'By canonical SMILES',
  vid: 'By V-iD',
};
const CLASH_LABELS: Record<HitDesignMergeClashStrategy, string> = {
  skip: 'Skip (keep campaign value)',
  overwrite: 'Overwrite (use incoming value)',
};

const NONE_OPTION = '(none)';

function setChoiceItems(input: DG.InputBase<string | null>, items: string[], desired?: string | null): void {
  (input as unknown as DG.ChoiceInput<string>).items = items;
  if (desired != null && items.includes(desired))
    input.value = desired;
  else if (items.length > 0)
    input.value = items[0];
  else
    input.value = null;
}

export async function openMergeTableDialog<T extends HitDesignApp>(app: T): Promise<void> {
  if (!app.dataFrame) {
    grok.shell.error('Cannot merge: campaign has no data.');
    return;
  }
  const saved: HitDesignMergeConfig | undefined = app.campaign?.mergeConfig;

  // Mutable dialog state
  let incomingDf: DG.DataFrame | null = null;
  let isLoadingFile = false;
  // Re-entrancy guard — onValueChanged on the file input fires for both user actions and
  // programmatic resets. We only want the validate-then-load behavior on user actions.
  let suppressFileChange = false;

  // ---- Inputs ----
  const fileInput = ui.input.file('Source file', {
    nullable: true,
    tooltipText: `Pick a file from your computer or from a Datagrok file share.\n` +
      `Allowed extensions: ${ALLOWED_EXTENSIONS.join(', ')}.\n` +
      `Files picked from a share are remembered with the campaign and reloaded automatically next time.`,
    onValueChanged: (v) => onFileChanged(v),
  });
  fileInput.input.setAttribute('placeholder', `Drop file or pick from share — ${ALLOWED_EXTENSIONS.join(', ')}`);

  const modeInput = ui.input.choice<string>('Merge by', {
    items: Object.values(MODE_LABELS),
    value: MODE_LABELS[saved?.mode ?? 'vid'],
    nullable: false,
    tooltipText: 'Join key. "Canonical SMILES" matches by canonicalized molecule structure; ' +
      '"V-iD" matches by the campaign\'s V-iD identifier.',
    onValueChanged: () => runValidation(),
  });

  const molColInput = ui.input.choice<string>('Molecule column (incoming)', {
    items: [],
    value: null as unknown as string,
    nullable: true,
    tooltipText: 'Column in the incoming table holding the molecules. ' +
      'Required for SMILES merge. In V-iD merge it is optional — if set, ' +
      'unmatched rows can be added to the campaign and registered with new V-iDs.',
    onValueChanged: () => runValidation(),
  });

  const vidColInput = ui.input.choice<string>('V-iD column (incoming)', {
    items: [],
    value: null as unknown as string,
    nullable: true,
    tooltipText: 'Column in the incoming table holding V-iD strings. Used as the join key when merging by V-iD.',
    onValueChanged: () => runValidation(),
  });

  const addRowsInput = ui.input.bool('Add unmatched rows', {
    value: saved?.addNewRows ?? true,
    tooltipText: 'When on, rows from the incoming table whose key is not present in the campaign are added as new rows. ' +
      'When off, only existing rows are updated. ' +
      'Disabled when no incoming molecule column is selected (nothing to register the new row with).',
    onValueChanged: () => runValidation(),
  });

  const runComputeInput = ui.input.bool('Run compute on new rows', {
    value: saved?.runComputeOnNewRows ?? true,
    tooltipText: 'Run the campaign\'s compute pipeline (descriptors, functions, scripts, queries) on the incoming molecules ' +
      'so new rows arrive with all calculated columns filled in. Disabled when no rows will be added.',
    onValueChanged: () => runValidation(),
  });

  const clashInput = ui.input.choice<string>('On clashing data', {
    items: Object.values(CLASH_LABELS),
    value: CLASH_LABELS[saved?.clashStrategy ?? 'overwrite'],
    nullable: false,
    tooltipText: 'How to resolve cells that have a value on both sides for a matched row. ' +
      '"Skip" keeps the existing campaign value; "Overwrite" replaces it with the incoming value. ' +
      'Empty campaign cells are always filled from the incoming table regardless of this setting.',
    onValueChanged: () => runValidation(),
  });

  const errorDiv = ui.divText('');
  errorDiv.style.color = 'var(--failure)';
  errorDiv.style.fontSize = '12px';
  errorDiv.style.marginTop = '8px';
  errorDiv.style.minHeight = '16px';

  const inputsRoot = ui.divV([
    fileInput.root,
    modeInput.root,
    molColInput.root,
    vidColInput.root,
    addRowsInput.root,
    runComputeInput.root,
    clashInput.root,
    errorDiv,
  ]);
  inputsRoot.style.minWidth = '420px';

  // ---- Helpers reading current input state ----
  const labelToMode = (label: string | null): HitDesignMergeMode =>
    label === MODE_LABELS.vid ? 'vid' : 'smiles';
  const labelToClash = (label: string | null): HitDesignMergeClashStrategy =>
    label === CLASH_LABELS.overwrite ? 'overwrite' : 'skip';
  const colChoiceValue = (input: DG.InputBase<string | null>): string | undefined => {
    const v = input.value;
    return v && v !== NONE_OPTION ? v : undefined;
  };

  // ---- Core validation / state machine ----
  function runValidation(): void {
    const mode = labelToMode(modeInput.value);
    const fileLoaded = !!incomingDf && !isLoadingFile;

    // Block / unblock inputs based on availability of incoming data + selected mode.
    modeInput.enabled = fileLoaded;
    molColInput.enabled = fileLoaded;
    vidColInput.enabled = fileLoaded && mode === 'vid';
    clashInput.enabled = fileLoaded;

    const molCol = colChoiceValue(molColInput);
    const vidCol = colChoiceValue(vidColInput);

    // Disable downstream toggles when their precondition isn't met, but DON'T clobber the
    // user's chosen value — the OK handler only honors them when they're effective, so the
    // value sticks across input changes (e.g. the user briefly unsets the molecule column).
    const canChooseAddRows = fileLoaded && !!molCol;
    addRowsInput.enabled = canChooseAddRows;

    const canRunCompute = canChooseAddRows && addRowsInput.value === true;
    runComputeInput.enabled = canRunCompute;

    let errorMessage = '';
    if (isLoadingFile)
      errorMessage = 'Reading file…';
    else if (!fileLoaded)
      errorMessage = 'Pick a source file to merge from.';
    else if (mode === 'smiles' && !molCol)
      errorMessage = 'Pick the molecule column from the incoming table.';
    else if (mode === 'vid' && !vidCol)
      errorMessage = 'Pick the V-iD column from the incoming table.';

    errorDiv.textContent = errorMessage;
    setOkEnabled(errorMessage === '');
  }

  // ---- File loading ----
  async function loadAndApplyFile(fi: DG.FileInfo, savedToApply: HitDesignMergeConfig | null): Promise<void> {
    isLoadingFile = true;
    incomingDf = null;
    runValidation();
    ui.setUpdateIndicator(inputsRoot, true, 'Reading file...');
    try {
      const df = await readFileInfoToDataFrame(fi);
      if (!df) throw new Error('Empty result');
      df.name = (fi.name ?? '').replace(/\.[^.]+$/, '');
      await df.meta.detectSemanticTypes();
      await grok.data.detectSemanticTypes(df);
      incomingDf = df;
    } catch (e) {
      _package.logger.error(e);
      grok.shell.error(`Failed to read "${fi.name}". Check that the file is accessible and well-formed.`);
      incomingDf = null;
    } finally {
      ui.setUpdateIndicator(inputsRoot, false);
      isLoadingFile = false;
    }
    populateColumnChoices(savedToApply);
    if (savedToApply)
      applySavedScalarChoices(savedToApply);
    runValidation();
  }

  function populateColumnChoices(savedToApply: HitDesignMergeConfig | null): void {
    if (!incomingDf) {
      setChoiceItems(molColInput, []);
      setChoiceItems(vidColInput, []);
      return;
    }
    const colNames = incomingDf.columns.toList().map((c) => c.name);

    const molDesired = savedToApply?.molColName ?? guessMoleculeColumn(incomingDf)?.name ?? null;
    const molItems = [NONE_OPTION, ...colNames];
    setChoiceItems(molColInput, molItems, molDesired ?? NONE_OPTION);

    const vidDesired = savedToApply?.vidColName ?? guessVidColumn(incomingDf)?.name ?? null;
    const vidItems = [NONE_OPTION, ...colNames];
    setChoiceItems(vidColInput, vidItems, vidDesired ?? NONE_OPTION);
  }

  function applySavedScalarChoices(savedToApply: HitDesignMergeConfig): void {
    modeInput.value = MODE_LABELS[savedToApply.mode];
    addRowsInput.value = savedToApply.addNewRows;
    runComputeInput.value = savedToApply.runComputeOnNewRows;
    clashInput.value = CLASH_LABELS[savedToApply.clashStrategy];
  }

  async function onFileChanged(fi: DG.FileInfo | null): Promise<void> {
    if (suppressFileChange) return;

    if (!fi) {
      incomingDf = null;
      populateColumnChoices(null);
      runValidation();
      return;
    }
    const ext = getExtension(fi);
    if (!ALLOWED_EXTENSIONS.includes(ext as typeof ALLOWED_EXTENSIONS[number])) {
      // Defer the clear to break the change → handler → change loop, and to make sure
      // we see the user's actual selection before reverting.
      setTimeout(() => {
        suppressFileChange = true;
        try {
          fileInput.value = null;
        } finally {
          suppressFileChange = false;
        }
        runValidation();
      }, 0);
      grok.shell.warning(`Unsupported file type ".${ext || '?'}". Allowed: ${ALLOWED_EXTENSIONS.join(', ')}.`);
      return;
    }
    await loadAndApplyFile(fi, null);
  }

  // ---- Dialog ----
  const dialog = ui.dialog('Merge with table');
  dialog.add(inputsRoot);
  dialog.onOK(async () => {
    await onOk();
  });

  function setOkEnabled(enabled: boolean): void {
    const okBtn = dialog.getButton('OK');
    if (!okBtn) return;
    if (enabled) {
      okBtn.removeAttribute('disabled');
      okBtn.classList.remove('disabled');
    } else {
      okBtn.setAttribute('disabled', 'true');
      okBtn.classList.add('disabled');
    }
  }

  async function onOk(): Promise<void> {
    if (!incomingDf) return;
    const mode = labelToMode(modeInput.value);
    const molCol = colChoiceValue(molColInput);
    const vidCol = colChoiceValue(vidColInput);
    const fi = fileInput.value;
    const filePath = fi && !isLocalUploadFileInfo(fi) ? fi.fullPath : undefined;

    const cfg: HitDesignMergeConfig = {
      mode,
      filePath,
      molColName: molCol,
      vidColName: mode === 'vid' ? vidCol : undefined,
      addNewRows: !!addRowsInput.value && !!molCol,
      runComputeOnNewRows: !!runComputeInput.value && !!addRowsInput.value && !!molCol,
      clashStrategy: labelToClash(clashInput.value),
    };

    const pg = DG.TaskBarProgressIndicator.create('Merging table into campaign...');
    try {
      await mergeIntoCampaign(app, incomingDf, cfg);
      if (app.campaign) {
        app.campaign.mergeConfig = cfg;
        await app.saveCampaign(false);
      }
      grok.shell.info('Table merged into campaign.');
    } catch (e) {
      _package.logger.error(e);
      grok.shell.error('Merge failed. See console for details.');
    } finally {
      pg.close();
    }
  }

  dialog.show();
  // Initial state: nothing loaded yet, all dependent inputs disabled.
  runValidation();

  // If the campaign already has a remembered file path, try to recreate the FileInfo and
  // walk through the same load → restore flow as a fresh user pick would.
  if (saved?.filePath) {
    try {
      const exists = await grok.dapi.files.exists(saved.filePath);
      if (!exists)
        grok.shell.warning(`Saved merge file not found or not accessible: ${saved.filePath}`);
      else {
        const restoredFi = DG.FileInfo.fromString(saved.filePath, '');
        suppressFileChange = true;
        try {
          fileInput.value = restoredFi;
        } finally {
          suppressFileChange = false;
        }
        await loadAndApplyFile(restoredFi, saved);
      }
    } catch (e) {
      _package.logger.error(e);
      grok.shell.warning('Could not restore previously saved merge file.');
    }
  }
}
