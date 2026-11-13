import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ObjectForm, propertyForm} from '@datagrok-libraries/u2/dg';

const SCORING = [
  DG.Property.fromOptions({name: 'method', caption: 'Method', type: DG.TYPE.STRING, choices: ['Lipinski', 'QED', 'MPO'],
    nullable: false}),
  DG.Property.fromOptions({name: 'cutoff', caption: 'Cutoff', type: DG.TYPE.FLOAT, min: 0, max: 1,
    description: 'Molecules scoring below it fail'}),
  DG.Property.fromOptions({name: 'normalize', caption: 'Normalize', type: DG.TYPE.BOOL}),
  DG.Property.fromOptions({name: 'column', caption: 'Result column', type: DG.TYPE.STRING, nullable: false,
    category: 'Output'}),
];

export function scoringSettings(settings: Record<string, unknown>): ObjectForm {
  return propertyForm(SCORING, settings,
    {onChanged: () => grok.userSettings.add('my-app', 'scoring', JSON.stringify(settings))});
}
