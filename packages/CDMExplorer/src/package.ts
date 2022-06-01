import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from "datagrok-api/ui";
import { DemographicsView } from './views/person-view';
import { MeasurementView } from './views/measurement-view';
import { cohorts } from './cohorts';
import { TimelinesView } from './views/timelines';
import { convertColToInt, convertColToString, dynamicComparedToBaseline, joinCohorts } from './preprocessing/data-preparation';
import { PERSON_ID, VIEWS } from './constants';
import { PersonAchillesView } from './views/person-achilles';
import { ConditionOccurrenceAchillesView } from './views/cond-occurrence-achilles';
import { AchillesResultsView } from './views/achilles-results';
import { addView, sortObject, updateDivInnerHTML } from './utils';
import { measurementTableView } from './views/measurements-table-view';
import { achillesResultsTableView } from './views/achilles-results-table-view';

export let _package = new DG.Package();
export let c: DG.FuncCall;

//name: CDM
//tags: app
export async function CdmApp(): Promise<any> {

  c = grok.functions.getCurrentCall();

  await cohorts.initFromDatabase();

 // views.push(<MeasurementView>addView(new MeasurementView('Measurement')));
 // views.push(<PersonAchillesView>addView(new PersonAchillesView('Person Achilles')));
/*    grok.data.query(`CDM:conditionOccurrence`).then(df => {
    convertColToString(df, PERSON_ID);
    joinCohorts(df);
    let tableView = DG.TableView.create(df, false);
    views.push(addView(tableView));
  }); */

 
   measurementTableView();
   achillesResultsTableView();
   VIEWS.push(<DemographicsView>addView(new DemographicsView('Demographics')));
   VIEWS.push(<TimelinesView>addView(new TimelinesView('Timelines')));
   //VIEWS.push(<ConditionOccurrenceAchillesView>addView(new ConditionOccurrenceAchillesView('Condition Occurrence Achilles')));

  let setObj = async (obj) => {
    grok.shell.o = await obj.propertyPanel();
  }

  grok.events.onCurrentViewChanged.subscribe((v) => {
    setTimeout(() => {
      let obj = VIEWS.find(it => it.name === grok.shell.v.name);
      if (obj) {
        if (obj.hasOwnProperty('loaded') && !obj.loaded) {
          obj.load();
        }
        if (obj.loaded) {
          setObj(obj);
        }
      }
    }, 100)
  });
}
