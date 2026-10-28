import type {AppHeaderOptions} from '@datagrok-libraries/u2';
import {_package} from '../../package';
import {PeptiHitTemplate} from '../types';
import {HitDesignInfoView} from '../hit-design-views/info-view';
import {PeptiHitApp} from '../pepti-hit-app';

export class PeptiHitInfoView extends HitDesignInfoView<PeptiHitTemplate, PeptiHitApp> {
  constructor(app: PeptiHitApp) {
    super(app);
    this.name = 'PeptiHit';
  }

  protected override header(): AppHeaderOptions {
    return {
      icon: _package.webRoot + '/images/icons/pepti-hit-icon.png',
      learnMore: 'https://github.com/datagrok-ai/public/blob/master/packages/HitTriage/README_HD.md',
      items: [
        'Configure your own workflow using the template editor',
        'Sketch Helm sequences in the spreadsheet',
        'Annotate and share ideas with the team',
        'Convert to molecular form and calculate different properties',
        'Save campaigns and continue from where you left off',
        'Submit final selection to the function of your choice',
      ],
    };
  }
}
