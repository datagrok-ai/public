/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {AppHeaderOptions} from '@datagrok-libraries/u2';
import type {TemplateFormContext} from '../shared/template-picker';
import {_package} from '../../package';
import {PepTriageTemplate, HitTriageTemplate} from '../types';
import {PepTriageDataSourceTag} from '../consts';
import {InfoView} from '../hit-triage-views/info-view';
import {PepTriageApp} from '../pep-triage-app';
import {INewCampaignResult, newCampaignAccordeon} from '../accordeons/new-campaign-accordeon';
import {createPepTriageTemplateAccordeon} from '../accordeons/new-pep-triage-template-accordeon';

export class PepTriageInfoView extends InfoView {
  private pepTriageDataSourceFunctionsMap: {[key: string]: DG.Func | DG.DataQuery} = {};

  constructor(app: PepTriageApp) {
    super(app);
    this.name = 'PepTriage';
  }

  protected override header(): AppHeaderOptions {
    return {
      icon: _package.webRoot + '/images/icons/hit-triage-icon.png',
      learnMore: 'https://github.com/datagrok-ai/public/blob/master/packages/HitTriage/README_HT.md',
      items: [
        'Configure your own workflow using the template editor',
        'Load sequence datasets and auto-convert to molecules',
        'Calculate different molecular properties',
        'Filter sequences using different criteria',
        'Save campaigns and continue any time from where you left off',
        'Submit final selection to the function of your choice',
      ],
    };
  }

  protected override getDataSourceTag(): string {
    return PepTriageDataSourceTag;
  }

  protected override getDataSourceFunctionsMap(): {[key: string]: DG.Func | DG.DataQuery} {
    return this.pepTriageDataSourceFunctionsMap;
  }

  protected override getNewCampaignForm(template: HitTriageTemplate, {owner, column}: TemplateFormContext): Promise<HTMLElement> {
    const ptTemplate = template as PepTriageTemplate;
    const seqColInput = ui.input.string('Sequence Column', {value: ptTemplate.sequenceColumnName ?? '', nullable: false});
    seqColInput.setTooltip('Name of the sequence column in the dataset');
    const molColInput = ui.input.string('Molecule Column', {value: ptTemplate.moleculeColumnName ?? '', nullable: true});
    molColInput.setTooltip('Name of the molecule column (leave empty to auto-convert from sequences)');
    return newCampaignAccordeon(template, this.getDataSourceFunctionsMap(), this.getDataSourceTag(), {
      owner,
      column,
      start: (camp) => this.startCampaign(template, camp),
      extraFormElements: [seqColInput.root, molColInput.root],
      getExtraData: () => ({
        sequenceColumnName: seqColInput.value?.trim() || ptTemplate.sequenceColumnName,
        moleculeColumnName: molColInput.value?.trim() || undefined,
      }),
    });
  }

  /** The template, with the sequence and molecule columns the user named. */
  protected override startCampaign(template: HitTriageTemplate, camp: INewCampaignResult) {
    const ptTemplate = template as PepTriageTemplate;
    return super.startCampaign({...ptTemplate,
      sequenceColumnName: camp.extraData?.sequenceColumnName ?? ptTemplate.sequenceColumnName,
      moleculeColumnName: camp.extraData?.moleculeColumnName ?? ptTemplate.moleculeColumnName,
    } as PepTriageTemplate, camp);
  }

  protected override async createNewTemplate(preset?: HitTriageTemplate): Promise<void> {
    ui.setUpdateIndicator(this.root, true);
    try {
      this.showTemplateEditor(await createPepTriageTemplateAccordeon(this.app, this.getDataSourceFunctionsMap(),
        preset as PepTriageTemplate));
    } finally {
      ui.setUpdateIndicator(this.root, false);
    }
  }
}
