import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from '@datagrok-libraries/utils/src/u2';
import {_package} from '../../package';
import {PeptiHitTemplate} from '../types';
import {newHitDesignCampaignAccordeon} from '../accordeons/new-hit-design-campaign-accordeon';
import {HitDesignInfoView} from '../hit-design-views/info-view';
import {PeptiHitApp} from '../pepti-hit-app';

export class PeptiHitInfoView extends HitDesignInfoView<PeptiHitTemplate, PeptiHitApp> {
  constructor(app: PeptiHitApp) {
    super(app);
    this.name = 'PeptiHit';
    grok.shell.windows.showHelp = true;
    grok.shell.windows.help.showHelp(_package.webRoot + 'README_HD.md');
  }

  override async getNewCampaignAccordeon(template: PeptiHitTemplate) {
    const {root, promise, cancelPromise} = newHitDesignCampaignAccordeon(template, true);
    promise.then(async (camp) => {
      this.app.dataFrame = camp.df;
      await this.app.setTemplate(template);
      this.app.campaignProps = camp.campaignProps;
      await this.app.saveCampaign(false);
      if (template.layoutViewState && this.app.campaign)
        this.app.campaign.layout = template.layoutViewState;
    });

    cancelPromise.then(() => {
      this.init();
    });
    return root;
  }

  override getAppHeader() {
    return u2.appHeader({
      iconPath: _package.webRoot + '/images/icons/pepti-hit-icon.png',
      learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/HitTriage/README_HD.md',
      description:
      '-  Configure your own workflow using the template editor\n' +
      '-  Sketch Helm sequences in the spreadsheet\n' +
      '-  Annotate and share ideas with the team\n' +
      '-  Convert to molecular form and calculate different properties\n' +
      '-  Save campaigns and continue from where you left off\n' +
      '-  Submit final selection to the function of your choice',
    });
  }
}
