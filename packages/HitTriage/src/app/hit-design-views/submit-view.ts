/* eslint-disable max-len */
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {HitDesignApp} from '../hit-design-app';
import {_package} from '../../package';
import {HitDesignTemplate} from '../types';
import {HitBaseView} from '../base-view';

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

  render(): HTMLDivElement {
    this.statusInput.value = this.app.campaign?.status ?? '';
    ui.empty(this.content);

    this.content = ui.divV([
      ui.h1('Summary'),
      ui.div([ui.tableFromMap(this.app.getSummary())]),
      this.statusInput.root,
    ]);
    return this.content;
  }

  onActivated(): void {
    this.render();
  }

  async submit(): Promise<any> {
    const submitParams= this.app.template?.submit;
    if (!submitParams)
      return;
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
