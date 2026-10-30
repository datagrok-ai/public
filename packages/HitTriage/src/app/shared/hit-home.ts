/* The landing of the HitTriage apps, top to bottom in one scrolling column: the app header, "Continue Campaign" with
   the campaigns table and its icons, "New Campaign" with the template picker and its form. */
import {AppHeader, Control, Page, component, div, signal} from '@datagrok-libraries/u2';
import type {AppHeaderOptions} from '@datagrok-libraries/u2';
import {i18n} from '../consts';
import {CampaignCatalog, CampaignCatalogOptions, campaignCatalog} from './campaign-catalog';
import {CampaignTemplate, PICKER_PX, TemplatePicker, TemplatePickerOptions, templatePicker} from './template-picker';

export interface HitHomeOptions<T extends CampaignTemplate> {
  header: AppHeaderOptions;
  catalog: CampaignCatalogOptions;
  picker: TemplatePickerOptions<T>;
}

/** Under this width of the campaigns table, campaign names are cut shorter. */
const NARROW_PX = 820;

export interface HitHome {
  readonly catalog: CampaignCatalog;
  readonly picker: Control & TemplatePicker;
}

export function hitHome<T extends CampaignTemplate>(options: HitHomeOptions<T>): Control & HitHome {
  return component((self) => {
    const narrow = signal(false);
    const catalog = campaignCatalog({...options.catalog, narrow});
    const table = catalog.table.root;
    // a hidden landing (another view current) measures 0: it keeps the width it had
    const observer = new ResizeObserver(() => {
      if (table.clientWidth > 0)
        narrow.value = table.clientWidth < NARROW_PX;
    });
    observer.observe(table);
    self.own(() => observer.disconnect());
    const picker = templatePicker<T>(options.picker);
    const page = new Page({labelColumn: false, sections: [
      {title: i18n.continueCampaigns, actions: catalog.chrome(), content: catalog.table},
      {title: i18n.createNewCampaignHeader, content: picker, maxWidth: PICKER_PX},
    ]});
    // in a plain block the page is as tall as its content: the view scrolls the header away with the rest, as before
    const root = div([div([new AppHeader({name: 'app-header', ...options.header})], 'hit-triage-home-header'),
      div([page])], 'hit-triage-home');
    return {root, catalog, picker};
  }, {name: 'hit-home'});
}
