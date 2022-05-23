import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { PersonView } from './views/person-view';
import { MeasurementView } from './views/measurement-view';
import { cohorts } from './cohorts';
import { TimelinesView } from './views/timelines';
import { convertColToString, joinCohorts } from './preprocessing/data-preparation';
import { PERSON_ID } from './constants';
import { PersonAchillesView } from './views/person-achilles';
import { ConditionOccurrenceAchillesView } from './views/cond-occurrence-achilles';
import { AchillesResultsView } from './views/achilles-results';

export let _package = new DG.Package();

//name: CDM
//tags: app
export async function CdmApp(): Promise<any> {

  const conn = grok.dapi.connections.list().then((res) => {
    const test = res;
  });

  let c: DG.FuncCall = grok.functions.getCurrentCall();
  function addView(view: DG.ViewBase): DG.ViewBase {
    view.box = true;
    view.parentCall = c;
    view.path = '/' + view.name;
    grok.shell.addView(view);
    return view;
  }

  await cohorts.initFromDatabase();

  const views = [];
 // views.push(<MeasurementView>addView(new MeasurementView('Measurement')));
 // views.push(<TimelinesView>addView(new TimelinesView('Timelines')));
 // views.push(<PersonAchillesView>addView(new PersonAchillesView('Person Achilles')));
 // views.push(<ConditionOccurrenceAchillesView>addView(new ConditionOccurrenceAchillesView('Condition Occurrence Achilles')));
  /* grok.data.query(`CDM:conditionOccurrence`).then(df => {
    convertColToString(df, PERSON_ID);
    joinCohorts(df);
    let tableView = DG.TableView.create(df, false);
    views.push(addView(tableView));
  });
  grok.data.query(`CDM:person`).then(df => {
    convertColToString(df, PERSON_ID);
    joinCohorts(df);
    let tableView = DG.TableView.create(df, false);
    views.push(addView(tableView));
  }); */
  // views.push(<PersonView>addView(new PersonView('Person')));
  views.push(<AchillesResultsView>addView(new AchillesResultsView('AchillesResults')));
}
