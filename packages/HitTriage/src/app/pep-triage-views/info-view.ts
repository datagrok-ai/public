/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from '@datagrok-libraries/utils/src/u2';
import {_package} from '../../package';
import {PepTriageTemplate, HitTriageTemplate} from '../types';
import {PepTriageDataSourceTag, i18n} from '../consts';
import {InfoView} from '../hit-triage-views/info-view';
import {PepTriageApp} from '../pep-triage-app';
import {addBreadCrumbsToRibbons, popRibbonPannels} from '../utils';
import {newCampaignAccordeon} from '../accordeons/new-campaign-accordeon';
import {createPepTriageTemplateAccordeon} from '../accordeons/new-pep-triage-template-accordeon';

export class PepTriageInfoView extends InfoView {
  private pepTriageDataSourceFunctionsMap: {[key: string]: DG.Func | DG.DataQuery} = {};

  constructor(app: PepTriageApp) {
    super(app);
    this.name = 'PepTriage';
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp(_package.webRoot + 'README_HT.md');
  }

  override onActivated(): void {
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp(_package.webRoot + 'README_HT.md');
  }

  protected override getAppHeader() {
    return u2.appHeader({
      iconPath: _package.webRoot + '/images/icons/hit-triage-icon.png',
      learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/HitTriage/README_HT.md',
      description:
      '-  Configure your own workflow using the template editor\n' +
      '-  Load sequence datasets and auto-convert to molecules\n' +
      '-  Calculate different molecular properties\n' +
      '-  Filter sequences using different criteria\n' +
      '-  Save campaigns and continue any time from where you left off\n' +
      '-  Submit final selection to the function of your choice',
    });
  }

  protected override getTemplatesFolder(): string {
    return `${this.app.appName}/templates`;
  }

  protected override getDataSourceTag(): string {
    return PepTriageDataSourceTag;
  }

  protected override getDataSourceFunctionsMap(): {[key: string]: DG.Func | DG.DataQuery} {
    return this.pepTriageDataSourceFunctionsMap;
  }

  protected override async getNewCampaignAccordeon(template: HitTriageTemplate, campaignDetailsDiv: HTMLElement) {
    const ptTemplate = template as PepTriageTemplate;
    const seqColInput = ui.input.string('Sequence Column', {value: ptTemplate.sequenceColumnName ?? '', nullable: false});
    seqColInput.setTooltip('Name of the sequence column in the dataset');
    const molColInput = ui.input.string('Molecule Column', {value: ptTemplate.moleculeColumnName ?? '', nullable: true});
    molColInput.setTooltip('Name of the molecule column (leave empty to auto-convert from sequences)');

    ui.setUpdateIndicator(campaignDetailsDiv, true);
    const {root, promise, cancelPromise} = await newCampaignAccordeon(
      template, this.getDataSourceFunctionsMap(), this.getDataSourceTag(), {
        extraFormElements: [seqColInput.root, molColInput.root],
        getExtraData: () => ({
          sequenceColumnName: seqColInput.value?.trim() || ptTemplate.sequenceColumnName,
          moleculeColumnName: molColInput.value?.trim() || undefined,
        }),
      });
    ui.setUpdateIndicator(campaignDetailsDiv, false);
    promise.then(async (camp) => {
      // Override template column names with user-specified values
      const modifiedTemplate = {...ptTemplate,
        sequenceColumnName: camp.extraData?.sequenceColumnName ?? ptTemplate.sequenceColumnName,
        moleculeColumnName: camp.extraData?.moleculeColumnName ?? ptTemplate.moleculeColumnName,
      };
      this.app.dataFrame = camp.df;
      this.app._fileInputType = camp.type;
      await this.app.setTemplate(modifiedTemplate, undefined, undefined, undefined, camp.friendlyName);
      this.app.campaignProps = camp.campaignProps;
      await this.app.saveCampaign(undefined, false);
      if (modifiedTemplate.layoutViewState && this.app.campaign)
        this.app.campaign.layout = modifiedTemplate.layoutViewState;
    });

    cancelPromise.then(() => {
      this.init();
    });
    return root;
  }

  protected override async createNewTemplate(preset?: HitTriageTemplate): Promise<void> {
    const newTemplateAccordeon = await createPepTriageTemplateAccordeon(this.app, this.getDataSourceFunctionsMap(), preset as PepTriageTemplate);

    const newView = new DG.ViewBase();
    const curView = grok.shell.v;
    newView.name = 'New Template';
    newView.root.appendChild(newTemplateAccordeon.root);
    newView.parentCall = this.app.parentCall;
    grok.shell.addView(newView);
    newView.path = new URL(this.app.baseUrl).pathname + '/new-template';
    newView.parentView = curView;
    const {sub} = addBreadCrumbsToRibbons(grok.shell.v, this.app.appName, i18n.createNewTemplate, () => {
      grok.shell.v = curView;
      newView.close();
    });

    newTemplateAccordeon.template.then(async (t) => {
      await this.init(t as any);
      newView.close();
      popRibbonPannels(newView);
      grok.shell.v = curView;
      sub.unsubscribe();
    });
    newTemplateAccordeon.cancelPromise.then(() => {
      sub.unsubscribe();
      newView.close();
      grok.shell.v = curView;
    });
  }
}
