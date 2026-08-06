// Negative tsc fixture: every wrong-typed usage below must FAIL compilation.
import {sampleUi, sampleEventUi} from './src/generated/db-ui';

export async function badUi(): Promise<void> {
  await sampleUi.grid({query: {columns: ['nope']}});
  await sampleUi.grid({defaults: {nope: 1}});
  await sampleUi.grid({defaults: {count: 'not-a-number'}});
  await sampleUi.list({query: {expand: ['details:nope']}});
  await sampleUi.app({query: {filter: {property: 'nope', operator: '=', value: 1}}});
  // a column of ANOTHER table of the same schema is still a wrong column here
  await sampleUi.grid({defaults: {kind: 'created'}});
  await sampleEventUi.grid({defaults: {sample_id: 42}});
  sampleUi.row({nope: 1});
  sampleUi.query({columns: ['nope']});
  // choices are the generated literal union
  await sampleUi.grid({defaults: {status: 'bogus'}});
  // §1.7: a misspelled column in the form's seed values fails compile
  await sampleUi.form({values: {nmae: 'a'}});
  await sampleUi.formDialog({values: {count: 'not-a-number'}});
}
