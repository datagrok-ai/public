/* Merge with table: a file from the user's computer or a share, merged into the campaign by canonical SMILES or by
   V-iD. OK stays off, saying why in the footer, until the dialog has what the merge needs; the dependent inputs wait
   for a file read. A file picked from a share is remembered with the campaign and read again when the dialog opens. */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {BoolInput, ChoiceInput, Dialog, Form, computed, signal} from '@datagrok-libraries/u2';
import {Platform, ViewHost, fileInput} from '@datagrok-libraries/u2/dg';
import type {HitDesignApp} from '../hit-design-app';
import {HitDesignMergeClashStrategy, HitDesignMergeConfig, HitDesignMergeMode} from '../types';
import {MERGE_ALLOWED_EXTENSIONS as ALLOWED_EXTENSIONS, getFileInfoExtension as getExtension, guessMoleculeColumn,
  guessVidColumn, isLocalUploadFileInfo, mergeIntoCampaign, readFileInfoToDataFrame} from '../utils/merge-table';
import {_package} from '../../package';

const MODE_LABELS: Record<HitDesignMergeMode, string> = {
  smiles: 'By canonical SMILES',
  vid: 'By V-iD',
};
const CLASH_LABELS: Record<HitDesignMergeClashStrategy, string> = {
  skip: 'Skip (keep campaign value)',
  overwrite: 'Overwrite (use incoming value)',
};

const NONE_OPTION = '(none)';

export async function openMergeTableDialog(app: HitDesignApp, host: ViewHost): Promise<void> {
  if (!app.dataFrame) {
    grok.shell.error('Cannot merge: campaign has no data.');
    return;
  }
  const saved: HitDesignMergeConfig | undefined = app.campaign?.mergeConfig;
  const dialog = new Dialog('Merge with table', {name: 'merge-table'});
  const load = dialog.build(() => {
    const incoming = signal<DG.DataFrame | null>(null);
    const reading = signal(false);
    const loaded = computed(() => incoming.value !== null && !reading.value);
    const listed = computed(() => incoming.value ? Platform.listedColumns(incoming.value) : []);
    const columns = computed(() => [NONE_OPTION, ...listed.value.map((c) => c.name)]);
    const picked = (input: ChoiceInput) => {
      const v = input.value;
      return v && v !== NONE_OPTION ? v : undefined;
    };

    const file = fileInput('Source file', {
      name: 'source-file', nullable: true, fill: true,
      placeholder: `Drop file or pick from share — ${ALLOWED_EXTENSIONS.join(', ')}`,
      tooltipText: `Pick a file from your computer or from a Datagrok file share.\n` +
        `Allowed extensions: ${ALLOWED_EXTENSIONS.join(', ')}.\n` +
        `Files picked from a share are remembered with the campaign and reloaded automatically next time.`,
      // Enter takes the typed path; it does not press OK
      onEnter: () => {},
      onChanged: (fi) => void read(fi, null),
    });
    file.addValidator((fi) => fi && !ALLOWED_EXTENSIONS.includes(getExtension(fi)) ?
      `Unsupported file type ".${getExtension(fi) || '?'}". Allowed: ${ALLOWED_EXTENSIONS.join(', ')}.` : null);
    const fileFirst = computed(() => reading.value ? 'Reading the source file…' :
      incoming.value === null ? 'Pick a source file first' : null);
    const mode = new ChoiceInput({label: 'Merge by', name: 'merge-by', items: Object.values(MODE_LABELS),
      value: MODE_LABELS[saved?.mode ?? 'vid'], nullable: false, disabled: fileFirst,
      tooltipText: 'Join key. "Canonical SMILES" matches by canonicalized molecule structure; ' +
        '"V-iD" matches by the campaign\'s V-iD identifier.'});
    const byVid = computed(() => mode.value === MODE_LABELS.vid);
    const molCol = new ChoiceInput({label: 'Molecule column (incoming)', name: 'molecule-column', items: columns,
      value: NONE_OPTION, nullable: true, disabled: fileFirst,
      tooltipText: 'Column in the incoming table holding the molecules. ' +
        'Required for SMILES merge. In V-iD merge it is optional — if set, ' +
        'unmatched rows can be added to the campaign and registered with new V-iDs.'});
    const vidCol = new ChoiceInput({label: 'V-iD column (incoming)', name: 'vid-column', items: columns,
      value: NONE_OPTION, nullable: true,
      disabled: computed(() => fileFirst.value ?? (byVid.value ? null : 'Used only when merging by V-iD')),
      tooltipText: 'Column in the incoming table holding V-iD strings. Used as the join key when merging by V-iD.'});
    const addOff = computed(() => fileFirst.value ??
      (picked(molCol) === undefined ? 'Pick the molecule column (incoming) first' : null));
    const addRows = new BoolInput({label: 'Add unmatched rows', name: 'add-rows', value: saved?.addNewRows ?? true,
      disabled: addOff,
      tooltipText: 'When on, rows from the incoming table whose key is not present in the campaign ' +
        'are added as new rows. ' +
        'When off, only existing rows are updated. ' +
        'Disabled when no incoming molecule column is selected (nothing to register the new row with).'});
    const runCompute = new BoolInput({label: 'Run compute on new rows', name: 'run-compute',
      value: saved?.runComputeOnNewRows ?? true,
      disabled: computed(() => addOff.value ?? (addRows.value ? null : 'Turn on Add unmatched rows first')),
      tooltipText: 'Run the campaign\'s compute pipeline (descriptors, functions, scripts, queries) ' +
        'on the incoming molecules ' +
        'so new rows arrive with all calculated columns filled in. Disabled when no rows will be added.'});
    const clash = new ChoiceInput({label: 'On clashing data', name: 'on-clash', items: Object.values(CLASH_LABELS),
      value: CLASH_LABELS[saved?.clashStrategy ?? 'overwrite'], nullable: false, disabled: fileFirst,
      tooltipText: 'How to resolve cells that have a value on both sides for a matched row. ' +
        '"Skip" keeps the existing campaign value; "Overwrite" replaces it with the incoming value. ' +
        'Empty campaign cells are always filled from the incoming table regardless of this setting.'});
    const autoMerge = new BoolInput({label: 'Auto-merge on campaign open', name: 'auto-merge',
      value: saved?.autoMergeOnOpen ?? false,
      disabled: computed(() => {
        const fi = file.value;
        return !fi ? 'Pick a source file first' : isLocalUploadFileInfo(fi) ?
          'Only a file picked from a Datagrok file share is merged again on open' : null;
      }),
      tooltipText: 'When on, this file is re-read and merged into the campaign every time the ' +
        'campaign is opened, using the configuration below. Compute is never re-run for this ' +
        'implicit pass and the campaign is not auto-saved. ' +
        'Only available for files picked from a Datagrok file share — local uploads cannot be ' +
        're-read on open.'});
    const form = new Form();
    for (const input of [file, mode, molCol, vidCol, addRows, runCompute, clash, autoMerge])
      form.add(input);
    dialog.add(form);

    dialog.okDisabled(computed(() => reading.value ? 'Reading file…' :
      !loaded.value ? 'Pick a source file to merge from.' :
        !byVid.value && picked(molCol) === undefined ? 'Pick the molecule column from the incoming table.' :
          byVid.value && picked(vidCol) === undefined ? 'Pick the V-iD column from the incoming table.' : null),
    {footer: true});

    // a file read after another was picked is dropped: the box holds the later one
    let turn = 0;
    async function read(fi: DG.FileInfo | null, preset: HitDesignMergeConfig | null): Promise<void> {
      const mine = ++turn;
      incoming.value = null;
      reading.value = false;
      if (!fi || !ALLOWED_EXTENSIONS.includes(getExtension(fi)))
        return;
      reading.value = true;
      let df: DG.DataFrame | null = null;
      try {
        df = await readFileInfoToDataFrame(fi);
        if (!df)
          throw new Error('Empty result');
        df.name = (fi.name ?? '').replace(/\.[^.]+$/, '');
        await df.meta.detectSemanticTypes();
        await grok.data.detectSemanticTypes(df);
      }
      catch (e) {
        _package.logger.error(e);
        grok.shell.error(`Failed to read "${fi.name}". Check that the file is accessible and well-formed.`);
        df = null;
      }
      if (mine !== turn)
        return;
      reading.value = false;
      incoming.value = df;
      if (!df)
        return;
      molCol.value = preset?.molColName ?? guessMoleculeColumn(listed.peek())?.name ?? NONE_OPTION;
      vidCol.value = preset?.vidColName ?? guessVidColumn(listed.peek())?.name ?? NONE_OPTION;
      if (preset) {
        mode.value = MODE_LABELS[preset.mode];
        addRows.value = preset.addNewRows;
        runCompute.value = preset.runComputeOnNewRows;
        clash.value = CLASH_LABELS[preset.clashStrategy];
        autoMerge.value = preset.autoMergeOnOpen ?? false;
      }
    }

    dialog.onOK(async () => {
      const df = incoming.peek();
      if (!df)
        return false;
      const molColName = picked(molCol);
      const fi = file.peek();
      const filePath = fi && !isLocalUploadFileInfo(fi) ? fi.fullPath : undefined;
      const add = addRows.peek() && !!molColName;
      const cfg: HitDesignMergeConfig = {
        mode: byVid.peek() ? 'vid' : 'smiles',
        filePath,
        molColName,
        vidColName: byVid.peek() ? picked(vidCol) : undefined,
        addNewRows: add,
        runComputeOnNewRows: add && runCompute.peek(),
        clashStrategy: clash.peek() === CLASH_LABELS.overwrite ? 'overwrite' : 'skip',
        // auto-merge needs a path to read again on open
        autoMergeOnOpen: autoMerge.peek() && !!filePath,
      };
      const progress = DG.TaskBarProgressIndicator.create('Merging table into campaign...');
      try {
        await mergeIntoCampaign(app, df, cfg);
        if (app.campaign) {
          app.campaign.mergeConfig = cfg;
          await app.saveCampaign(false);
        }
        grok.shell.info('Table merged into campaign.');
      }
      catch (e) {
        _package.logger.error(e);
        grok.shell.error('Merge failed. See console for details.');
        return false;
      } finally {
        progress.close();
      }
    });
    return {file, read};
  });
  host.show(dialog, {width: 560});

  if (!saved?.filePath)
    return;
  try {
    if (!await grok.dapi.files.exists(saved.filePath)) {
      grok.shell.warning(`Saved merge file not found or not accessible: ${saved.filePath}`);
      return;
    }
    // the user may have closed the dialog or picked a file meanwhile
    if (!dialog.isOpen.peek() || load.file.peek() !== null)
      return;
    const restored = DG.FileInfo.fromString(saved.filePath, '');
    // the dialog's own onChanged would read it without the remembered choices
    load.file.reset(restored);
    await load.read(restored, saved);
  }
  catch (e) {
    _package.logger.error(e);
    grok.shell.warning('Could not restore previously saved merge file.');
  }
}
