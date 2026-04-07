/* eslint-disable max-len */
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {HitDesignApp} from '../hit-design-app';
import {_package} from '../../package';
import {HitDesignTemplate} from '../types';
import {HitBaseView} from '../base-view';
import {HitTriageSubmitTag} from '../consts';
import {getFuncPackageNameSafe} from '../utils';

export class HitDesignSubmitView extends HitBaseView<HitDesignTemplate, HitDesignApp> {
  private statusInput: DG.InputBase<string | undefined>;
  private statusSuggestionsMenu: DG.Menu;
  content: HTMLDivElement;
  constructor(app: HitDesignApp) {
    super(app);
    this.name = 'Submit';
    this.statusInput = ui.input.string('Status', {value: this.app.campaign?.status, nullable: false});
    this.statusSuggestionsMenu = DG.Menu.popup();
    this.statusInput.root.style.marginLeft = '12px';
    this.content = ui.div();

    this.statusInput.onChanged.subscribe(() => {
      this.statusSuggestionsMenu.clear();
      const status = (this.statusInput.value ?? '').toLowerCase();
      // eslint-disable-next-line max-len
      const similarStatuses = this.app.existingStatuses.filter((s) => s?.toLowerCase()?.includes(status) && s?.toLowerCase() !== status).filter((_, i) => i < 5);
      if (similarStatuses.length) {
        similarStatuses.forEach((s) => {
          this.statusSuggestionsMenu.item(s, () => {
            this.statusInput.value = s;
            this.statusSuggestionsMenu.root.remove();
            this.statusInput.root.focus();
          });
        });
        if (this.content.parentElement?.parentElement) {
          const xOffset = this.statusInput.root.offsetLeft + (this.statusInput.input?.offsetLeft ?? 0);
          const yOffset = this.statusInput.root.offsetTop + this.statusInput.root.offsetHeight + this.content.parentElement.offsetTop;
          this.statusSuggestionsMenu.show({element: this.content.parentElement.parentElement!, x: xOffset, y: yOffset});
        }
      } else
        this.statusSuggestionsMenu.root.remove();
    });
  }

  public getStatus() {
    return this.statusInput.value ?? this.app.campaign?.status ?? 'No Status';
  }

  render() {
    this.statusInput.value = this.app.campaign?.status ?? '';
    ui.empty(this.content);

    this.content = ui.divV([
      ui.h1('Summary'),
      ui.div([ui.tableFromMap(this.app.getSummary())]),
      this.statusInput.root,
    ]);


    const dlg = ui.dialog('Submit');
    dlg.add(this.content);
    dlg.addButton('Save', () => {
      this.app.campaign!.status = this.getStatus();
      this.app.saveCampaign();
      dlg.close();
    });
    //dlg.addButton('Submit', ()=>{this.submit(); dlg.close();}, undefined, 'Submit the campaign file to the specified function');
    dlg.show();

    const submitFunctions = Array.from(new Set(DG.Func.find({meta: {role: HitTriageSubmitTag}}).concat(DG.Func.find({tags: [HitTriageSubmitTag]}))));
    const submitFunctionMap = submitFunctions.reduce((acc, fn) => {
      acc[fn.friendlyName ?? fn.name] = fn;
      return acc;
    }, {} as Record<string, DG.Func>);
    if (submitFunctions.length > 0) {
      const chosenFunctionParams = this.app.submitParams;
      const foundFunction = chosenFunctionParams ? submitFunctions.find((fn) => fn.name === chosenFunctionParams.fName && getFuncPackageNameSafe(fn) == chosenFunctionParams.package) : undefined;
      const foundFunctionKey = foundFunction?.friendlyName ?? foundFunction?.name;
      const submitFunctionInput = ui.input.choice('Submit function', {
        value: foundFunctionKey,
        items: [...Object.keys(submitFunctionMap)],
        nullable: true,
      });

      submitFunctionInput.onChanged.subscribe(() => {
        const selectedFunction = submitFunctionMap[submitFunctionInput.value ?? ''];
        this.app.campaign!.template!.submit = selectedFunction ? {
          fName: selectedFunction.name,
          package: getFuncPackageNameSafe(selectedFunction),
        } : undefined;
        dlg.getButton('Submit')?.remove();
        if (selectedFunction) {
          dlg.addButton('Submit', () => {
            this.submit();
            dlg.close();
          }, undefined, 'Submit the campaign file to the specified function');
        }
      });
      submitFunctionInput.fireChanged();
      this.content.appendChild(submitFunctionInput.root);
    }
  }

  onActivated(): void {
  }

  async submit(): Promise<any> {
    const submitParams= this.app.submitParams;
    if (!submitParams) {
      grok.shell.error('No submit function selected. Please select a function to submit the campaign.');
      return;
    }

    const submitFn = DG.Func.find({name: submitParams.fName, package: submitParams.package})[0];
    if (!submitFn) {
      grok.shell.error(`Function ${submitParams.fName} not found.`);
      return;
    }
    const filteredDf = DG.DataFrame.fromCsv(this.app.dataFrame!.toCsv({filteredRowsOnly: true}));
    await submitFn.apply({df: filteredDf, molecules: this.app.molColName});
    this.app.campaign!.status = this.getStatus();
    this.app.saveCampaign();
    grok.shell.info('Submitted successfully.');
  }
}
