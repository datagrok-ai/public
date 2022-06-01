/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from "datagrok-api/ui";
import { PERSON_ID } from './constants';
import { convertColToString } from './preprocessing/data-preparation';


export class Cohorts {
    cohorts: DG.DataFrame = null;
    cohortsPivoted: DG.DataFrame;
    cohortFilters: DG.Viewer;

    async initFromDatabase(): Promise<void> {
        await grok.data.query(`CDM:cohorts`)
            .then(cohortsDf => {
                this.cohorts = cohortsDf;
                this.cohortsPivoted = this.cohorts
                    .groupBy([PERSON_ID])
                    .pivot('cohort_definition_id')
                    .first(PERSON_ID)
                    .aggregate();
                this.cohortsPivoted.columns.names().filter(it => it !== PERSON_ID).forEach(name => {
                    const cohortId = name.replace(` first(${PERSON_ID})`, '');
                    const cohortName = cohortsDf
                        .groupBy(['cohort_definition_id', 'cohort_definition_name'])
                        .where({ cohort_definition_id: `${cohortId}` })
                        .aggregate()
                        .get('cohort_definition_name', 0);
                    this.cohortsPivoted.col(name).name = cohortName;
                });
                this.cohortsPivoted.columns.remove('');
                this.convertCohortColsToBoolean();
                convertColToString(this.cohortsPivoted, PERSON_ID);
                this.createCohortsFilter();      
            });
    }

    private convertCohortColsToBoolean() {
        let cohortNames = this.cohortsPivoted.columns.names().filter(it => it !== PERSON_ID);
        cohortNames.forEach(cohortName => {
            let cohortCol = this.cohortsPivoted.col(cohortName);
            this.cohortsPivoted.columns.addNewBool(`${cohortName}_1`).init(i => !cohortCol.isNone(i));
            this.cohortsPivoted.columns.remove(cohortName);
            this.cohortsPivoted.col(`${cohortName}_1`).name = cohortName;
        })
    }

    private createCohortsFilter() {
        this.cohortFilters = DG.Viewer.fromType('Filters', this.cohortsPivoted, {
            'showContextMenu': false,
          });
            grok.shell.topMenu
              .group('Cohort')
              .item('Filter', () => {
                const dialog = ui.dialog({ title: '' })
                .add(ui.div(this.cohortFilters.root))
                .onOK(() => {})
                .show();
                dialog.root.addEventListener("mouseenter", (event) => {
                  this.cohortFilters.root.removeAttribute('data-widget');
                });
              });
    }

}


export let cohorts: Cohorts = new Cohorts();