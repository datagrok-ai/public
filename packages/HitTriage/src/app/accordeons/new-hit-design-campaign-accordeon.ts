/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {CampaignFieldTypes, HitDesignTemplate} from '../types';
import {HitDesignMolColName, PeptiHitHelmColName, TileCategoriesColName, ViDColName, i18n} from '../consts';
import '../../../css/hit-triage.css';
import {getNewVid} from '../utils/calculate-single-cell';

type NewHitDesignCampaignRes = {
    df: DG.DataFrame,
    campaignProps: {[key: string]: any}
}

type HitDesignCampaignAccordeon = {
    promise: Promise<NewHitDesignCampaignRes>,
    root: HTMLElement,
    cancelPromise: Promise<void>
}

export function newHitDesignCampaignAccordeon(template: HitDesignTemplate, peptiHit = false): HitDesignCampaignAccordeon {
  const df = DG.DataFrame.create(1);
  if (peptiHit) {
    const pepCol = df.columns.addNew(PeptiHitHelmColName, DG.TYPE.STRING);
    pepCol.semType = DG.SEMTYPE.MACROMOLECULE;
    pepCol.meta.units = 'helm';
    pepCol.setTag('units', 'helm');
    pepCol.setTag('cell.renderer', 'helm');
    pepCol.setTag('.alphabetIsMultichar', 'true');
  }

  const molCol = df.columns.addNew(HitDesignMolColName, DG.TYPE.STRING);
  molCol.semType = DG.SEMTYPE.MOLECULE;
  molCol.meta.units = DG.UNITS.Molecule.SMILES;
  molCol.set(0, null);
  if (template?.stages?.length > 0) {
    const tileCategoryCol = df.columns.addNew(TileCategoriesColName, DG.TYPE.STRING);
    tileCategoryCol.set(0, template.stages[0]);
  }
  const vIdCol = df.columns.addNew(ViDColName, DG.TYPE.STRING);
  vIdCol.set(0, getNewVid(vIdCol));

  // campaign properties. each template might have number of additional fields that should
  // be filled by user for the campaign. they are cast into DG.Property objects and displayed as a form
  const campaignProps = template.campaignFields
    .map((field) => field.type === DG.SEMTYPE.MOLECULE ?
      ({...field, type: 'String', semtype: DG.SEMTYPE.MOLECULE}) : field)
    .map((field) =>
      DG.Property.fromOptions(
        {name: field.name, type: CampaignFieldTypes[field.type as keyof typeof CampaignFieldTypes],
          nullable: !field.required, ...(field.semtype ? {semType: field.semtype} : {})}));
  const campaignPropsObject: {[key: string]: any} = {};
  const campaignPropsForm = ui.input.form(campaignPropsObject, campaignProps);
  campaignPropsForm.classList.remove('ui-form');
  const form = ui.div([
    ...(template.campaignFields?.length > 0 ? []: []),
    campaignPropsForm,
  ]);
  const buttonsDiv = ui.buttonsInput([]); // div for create and cancel buttons
  form.appendChild(buttonsDiv);
  const okPromise = new Promise<NewHitDesignCampaignRes>((resolve) => {
    const startCampaignButton = ui.bigButton(i18n.StartCampaign,
      () => {
        for (const field of template.campaignFields) {
          if (field.required) {
            if (!campaignPropsObject[field.name] || campaignPropsObject[field.name] === '') {
              grok.shell.error(`Field '${field.name}' is required`);
              return;
            }
          }
        }
        resolve({df, campaignProps: campaignPropsObject});
      });
    buttonsDiv.appendChild(startCampaignButton);
  });
  const cancelPromise = new Promise<void>((resolve) => {
    const _cancelButton = ui.button(i18n.cancel, () => resolve());
    //buttonsDiv.appendChild(cancelButton);
  });

  return {promise: okPromise, root: form, cancelPromise};
}
