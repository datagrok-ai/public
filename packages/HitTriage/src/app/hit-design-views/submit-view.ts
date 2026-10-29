/* The Submit dialog of a Hit Design / PeptiHit campaign: what the campaign is (read-only rows, its name and folder
   editable in place), its status with the statuses other campaigns have, and the function it is submitted to. SAVE
   stores the status; SUBMIT runs the function on the campaign's rows, then saves. A campaign the user may only view
   offers neither, saying why. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ChoiceInput, Dialog, DynamicInput, Form, Input, LabelColumn, SuggestInput, TextInput, computed, div,
  iconButton, signal} from '@datagrok-libraries/u2';
import {ViewHost} from '@datagrok-libraries/u2/dg';
import type {HitDesignApp} from '../hit-design-app';
import {CampaignTableName, HitTriageSubmitTag, i18n} from '../consts';
import {checkFileExists, getFuncPackageNameSafe} from '../utils';
import {fieldText} from '../shared/campaign-catalog';

const VIEW_ONLY = 'You can view this campaign but not change it';
const NOT_SAVED = 'The campaign was not saved';

export class HitDesignSubmitView {
  constructor(private readonly app: HitDesignApp) {}

  /** The functions a campaign can be submitted to, by friendly name. */
  private static submitFunctions(): Map<string, DG.Func> {
    const funcs = new Set([...DG.Func.find({meta: {role: HitTriageSubmitTag}}),
      ...DG.Func.find({tags: [HitTriageSubmitTag]})]);
    return new Map([...funcs].map((f) => [f.friendlyName ?? f.name, f]));
  }

  render(): void {
    const app = this.app;
    const campaign = app.campaign;
    const canEdit = app.canEdit;
    const functions = HitDesignSubmitView.submitFunctions();
    const preset = app.submitParams;
    const presetName = [...functions].find(([, f]) => f.name === preset?.fName &&
      getFuncPackageNameSafe(f) === preset?.package)?.[0] ?? null;
    const dialog = new Dialog('Submit', {name: 'submit'});
    dialog.build(() => {
      const column = new LabelColumn();
      const summary = new Form({quiet: true, labelColumn: column});
      for (const input of this.summaryRows())
        summary.add(input);
      const inputs = new Form({labelColumn: column});
      const status = new SuggestInput({label: 'Status', name: 'status', value: campaign?.status ?? '',
        readOnly: !canEdit, minChars: 1, source: async (query) => {
          const q = query.toLowerCase();
          return app.existingStatuses.filter((s) => s?.toLowerCase().includes(q) && s.toLowerCase() !== q).slice(0, 5);
        }});
      const func = new ChoiceInput({label: 'Submit function', name: 'submit-function', items: [...functions.keys()],
        value: presetName, nullable: true, readOnly: !canEdit});
      inputs.add(status);
      inputs.add(func);
      dialog.add(div([summary.root, inputs.root]));

      const statusToSave = () => status.peek()?.trim() || campaign?.status || 'No Status';
      /** The chosen submit function into the campaign's template; answers the undo. */
      const keepSubmitFunction = () => {
        const template = app.campaign?.template;
        const before = template?.submit;
        const f = functions.get(func.peek() ?? '');
        if (template)
          template.submit = f ? {fName: f.name, package: getFuncPackageNameSafe(f)} : undefined;
        return () => template && (template.submit = before);
      };
      const chosen = computed(() => functions.get(func.value ?? ''));
      dialog.addButton('Save', async () => {
        const saving = app.campaign;
        if (!saving)
          return;
        const undo = keepSubmitFunction();
        const before = saving.status;
        saving.status = statusToSave();
        if (await app.saveCampaign() !== undefined)
          dialog.close();
        else {
          saving.status = before;
          undo();
        }
      }, {disabled: canEdit ? null : VIEW_ONLY});
      dialog.onOK(async () => {
        const undo = keepSubmitFunction();
        let submitted = false;
        try {
          submitted = await this.submit(statusToSave());
        } finally {
          // a submit function that threw leaves the template as it was; the dialog shows the error
          if (!submitted)
            undo();
        }
        return submitted;
      }, {text: 'Submit', primary: true});
      dialog.okDisabled(computed(() => !canEdit ? VIEW_ONLY :
        chosen.value === undefined ? 'Choose a submit function' : null));
      dialog.status(computed(() => chosen.value === undefined ? null :
        `Submits ${app.dataFrame?.filter.trueCount ?? 0} molecules to ${chosen.value.friendlyName}`));
    });
    ViewHost.of(app.designView!).show(dialog, {width: 640});
  }

  /** The summary rows: what the campaign is, its fields, its size and compute, its name and folder editable. */
  private summaryRows(): Input<any>[] {
    const app = this.app;
    const campaign = app.campaign;
    const template = app.template;
    const code = app.campaignId ?? campaign?.name ?? '';
    const text = (label: string, value: unknown) => new TextInput({label, value: String(value ?? ''), readOnly: true});
    const name = this.editable(text('Name', campaign?.friendlyName ?? campaign?.name ?? code), 'Edit Campaign Name',
      async (value) => {
        if (!value.trim())
          return 'Campaign name can not be empty';
        const saving = app.campaign!;
        const before = saving.friendlyName;
        saving.friendlyName = value;
        if (await app.saveCampaign(true) !== undefined)
          return null;
        saving.friendlyName = before;
        return NOT_SAVED;
      });
    const path = this.editable(text('File Path', app.campaignFolder), 'Edit file path', async (value) => {
      if (!await checkFileExists(value))
        return 'Given folder does not exist';
      const before = app.filePath;
      app.filePath = `${value.replace(/\/$/, '')}/${code}/${CampaignTableName}`;
      if (app.campaign)
        app.campaign.savePath = app.filePath;
      if (await app.saveCampaign(true) !== undefined)
        return null;
      app.filePath = before;
      if (app.campaign)
        app.campaign.savePath = before;
      return NOT_SAVED;
    });
    path.addOptions(iconButton('arrow-to-bottom', () => app.download(app.dataFrame!, code || 'Molecules'),
      {tooltip: i18n.download}));
    const rows: Input<any>[] = [name, text('Code', code), text('Template', template?.name ?? 'Molecules'), path];
    if (campaign?.authorUserFriendlyName)
      rows.push(text('Author', campaign.authorUserFriendlyName));
    if (campaign?.lastModifiedUserName)
      rows.push(text('Last Modified By', campaign.lastModifiedUserName));
    if (campaign?.createDate)
      rows.push(text('Create Date', campaign.createDate));
    for (const [key, value] of Object.entries(campaign?.campaignFields ?? app.campaignProps)) {
      const field = template?.campaignFields?.find((f) => f.name === key);
      if (field?.type === 'Molecule' && value) {
        rows.push(new DynamicInput<string>({label: key, value,
          render: (smiles) => smiles ? grok.chem.drawMolecule(smiles, 160, 100) : null}));
      }
      else
        rows.push(text(key, field?.type === 'Date' ? fieldText(value) : value));
    }
    rows.push(text('Number of Molecules', app.dataFrame?.rowCount ?? 0));
    const compute = template?.compute;
    const friendly = (f: {package: string, name: string}) =>
      DG.Func.find({package: f.package, name: f.name})[0]?.friendlyName ?? f.name;
    const methods = [compute?.descriptors.enabled ? 'Descriptors' : '', ...(compute?.functions ?? []).map(friendly)];
    rows.push(text('Enrichment Methods', methods.filter((m) => m.trim() !== '').join(', ')));
    return rows;
  }

  /** A read-only row with a pencil: it turns into an input with save and cancel (Enter and Escape alike); a save
   * that answers a problem keeps the input open, the problem under it. */
  private editable(input: TextInput, tooltip: string, save: (value: string) => Promise<string | null>): TextInput {
    const saved = signal(input.peek());
    const editing = signal(false);
    const problem = signal<{value: string, message: string} | null>(null);
    input.addValidator((value) => problem.value?.value === value ? problem.value.message : null);
    const done = (value: string) => {
      saved.value = value;
      problem.value = null;
      input.reset(value);
      editing.value = false;
    };
    const accept = async () => {
      const value = input.peek() ?? '';
      const error = await save(value);
      if (error === null)
        done(value);
      else {
        problem.value = {value, message: error};
        input.touch();
      }
    };
    const cancel = () => done(saved.peek() ?? '');
    const pencil = iconButton('pencil', () => {
      editing.value = true;
      input.focus();
    }, {tooltip});
    const acceptButton = iconButton('check', accept, {tooltip: 'Save'});
    const cancelButton = iconButton('times', cancel, {tooltip: 'Cancel'});
    for (const option of [pencil, acceptButton, cancelButton])
      input.addOptions(option);
    // while editing, Enter and Escape are the row's: the dialog would submit or close
    input.scope.listen(input.root, 'keydown', (e) => {
      const key = (e as KeyboardEvent).key;
      if (!editing.peek() || key !== 'Enter' && key !== 'Escape')
        return;
      e.preventDefault();
      e.stopPropagation();
      void (key === 'Enter' ? accept() : cancel());
    });
    input.effect(() => {
      const on = editing.value;
      input.readOnly = !on;
      pencil.hidden = on;
      acceptButton.hidden = !on;
      cancelButton.hidden = !on;
    });
    return input;
  }

  /** The AI's submit: the status given, else the campaign's, else "Submitted". */
  async submitWithStatus(status?: string): Promise<boolean> {
    return await this.submit(status?.trim() || this.app.campaign?.status || 'Submitted');
  }

  /** Runs the campaign's submit function on its rows (the filtered ones), then saves the campaign with `status`.
   * Answers whether both happened. */
  async submit(status: string = this.app.campaign?.status ?? 'No Status'): Promise<boolean> {
    const app = this.app;
    const params = app.submitParams;
    if (!params) {
      grok.shell.error('No submit function selected. Please select a function to submit the campaign.');
      return false;
    }
    const fn = DG.Func.find({name: params.fName, package: params.package})[0];
    if (!fn) {
      grok.shell.error(`Function ${params.fName} not found.`);
      return false;
    }
    const filtered = DG.DataFrame.fromCsv(app.dataFrame!.toCsv({filteredRowsOnly: true}));
    await fn.apply({df: filtered, molecules: app.molColName});
    const campaign = app.campaign!;
    const before = campaign.status;
    campaign.status = status;
    if (await app.saveCampaign() === undefined) {
      campaign.status = before;
      return false;
    }
    grok.shell.info('Submitted successfully.');
    return true;
  }
}
