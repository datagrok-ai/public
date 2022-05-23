/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { PERSON_ID } from './constants';
import { convertColToString } from './preprocessing/data-preparation';


export class Cohorts {
    cohorts: DG.DataFrame = null;
    cohortsPivoted: DG.DataFrame;

    async initFromDatabase(): Promise<void> {
        await grok.data.query(`CDM:cohorts`)
            .then(cohortsDf => {
                this.cohorts = cohortsDf;
                this.cohortsPivoted = this.cohorts
                    .groupBy([PERSON_ID])
                    .pivot('cohort_id')
                    .first(PERSON_ID)
                    .aggregate();
                this.cohortsPivoted.columns.names().filter(it => it !== PERSON_ID).forEach(name => {
                    this.cohortsPivoted.col(name).name = name.replace(` first(${PERSON_ID})`, '');
                })
                this.convertCohortColsToBoolean();
                convertColToString(this.cohortsPivoted, PERSON_ID);     
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

}


export let cohorts: Cohorts = new Cohorts();