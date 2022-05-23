import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { sortObject, updateDivInnerHTML } from '../utils';
import { convertColToInt, convertColToString, dynamicComparedToBaseline, joinCohorts } from '../preprocessing/data-preparation';
import { PERSON_ID } from '../constants';

export class MeasurementView extends DG.ViewBase {

    measurementTypes = {};

    constructor(name) {
        super({});
        this.name = name;
        this.createView();
    }

    async createView(): Promise<void> {
        grok.data.query(`CDM:measurementTypes`)
            .then(measurementTypesDf => {
               for (let i = 0; i < measurementTypesDf.rowCount; i++){
                this.measurementTypes[measurementTypesDf.get('measurement', i)] = measurementTypesDf.get('concept_id', i);
               }
               this.measurementTypes = sortObject(this.measurementTypes);
               let measurementChoices = ui.choiceInput('', Object.keys(this.measurementTypes)[0], Object.keys(this.measurementTypes));
               measurementChoices.onChanged((v) => {
                   grok.data.query(`CDM:measurementByConceptId`, {measurement_conc_id: this.measurementTypes[measurementChoices.value]})
                   .then(measurementsDf => {
                    measurementsDf.name = measurementChoices.value;
                    convertColToInt(measurementsDf, 'counter');
                    convertColToString(measurementsDf, PERSON_ID);
                    dynamicComparedToBaseline(measurementsDf, 'measurement_value', '1', 'counter', 'changes_from_bl');
                    joinCohorts(measurementsDf);
                    grok.shell.addTableView(measurementsDf);
                   })
               })
               measurementChoices.input.style.width = '150px';
               this.setRibbonPanels([
                [
                    measurementChoices.root
                ],
              ]); 

            });
    }

}