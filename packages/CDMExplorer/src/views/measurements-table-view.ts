import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from "datagrok-api/ui";
import { addView, sortObject, updateDivInnerHTML } from '../utils';
import { convertColToInt, convertColToString, dynamicComparedToBaseline, joinCohorts } from '../preprocessing/data-preparation';
import { PERSON_ID, VIEWS } from '../constants';
import { cohorts } from '../cohorts';

export function measurementTableView(): void {
    const measurementView = DG.TableView.create(DG.DataFrame.create());
    measurementView.name = 'Measurements';
    VIEWS.push(addView(measurementView));
    ui.setUpdateIndicator(measurementView.root, true);
    grok.data.query(`CDM:measurementTypes`)
      .then(measurementTypesDf => {
        let measurementTypes = {};
        for (let i = 0; i < measurementTypesDf.rowCount; i++) {
          measurementTypes[measurementTypesDf.get('measurement', i)] = measurementTypesDf.get('concept_id', i);
        }
        measurementTypes = sortObject(measurementTypes);
        let measurementsDf: DG.DataFrame = null;
        let selectedValues = [];
  
        let findCorrDiv = ui.div();
        let selectedCorrAggr = ['m_year', 'm_week'];
        const aggrCorrChoice = ui.choiceInput('', 'week', ['year', 'month', 'week']);
        aggrCorrChoice.onChanged(() => {
          selectedCorrAggr = aggrCorrChoice.value === 'year' ? ['m_year'] : 
          aggrCorrChoice.value === 'month' ? ['m_year', 'm_month'] :
          ['m_year', 'm_week'];
        });
  
        let findCorrButton = ui.button('Correlations', () => {
          ui.dialog({ title: 'Aggregation type' })
          .add(ui.div(aggrCorrChoice.root))
          .onOK(() => {
            const pivotedMeasurement = measurementsDf
            .groupBy([PERSON_ID].concat(selectedCorrAggr))
            .pivot('measurement')
            .med('measurement_value')
            .aggregate();
            selectedCorrAggr.forEach(col => pivotedMeasurement.columns.remove(col));
            grok.data.linkTables(cohorts.cohortsPivoted, pivotedMeasurement,
              [ PERSON_ID ], [ PERSON_ID ],
              [ DG.SYNC_TYPE.FILTER_TO_FILTER ]);
            const pivotedView = grok.shell.addTableView(pivotedMeasurement);
            pivotedView.name = `Correlations/aggr type ${aggrCorrChoice.value}/`;
            pivotedMeasurement.plot.fromType(DG.VIEWER.CORR_PLOT).then((v: any) => {
              pivotedView.addViewer(v);
            });
          })
          .show();
        });
  
        let selectMeasurementsButton = ui.button('Select measurements', () => {
          let measurementChoices = ui.multiChoiceInput('', selectedValues, Object.keys(measurementTypes));
          measurementChoices.onChanged((v) => {
            selectedValues = measurementChoices.value;
          });
          measurementChoices.input.style.maxWidth = '100%';
          measurementChoices.input.style.maxHeight = '100%';
          ui.dialog({ title: 'Select values' })
            .add(ui.div(measurementChoices.root, { style: { width: '400px', height: '300px' } }))
            .onOK(() => {
              ui.setUpdateIndicator(measurementView.root, true);
              Promise.all(selectedValues.map(it =>
                grok.data.query(`CDM:measurementByConceptId`,
                  { measurement_conc_id: measurementTypes[it] })))
                  .then(dfs => {
                    dfs.forEach(df => {
                      const concept_name = 
                        Object.keys(measurementTypes).find(key => measurementTypes[key] === df.get('measurement_concept_id', 0));
                      df.columns.addNewString('measurement').init((i)=> concept_name);
                      })
                    measurementsDf = dfs.reduce((previous, current) => previous.append(current));
                    measurementsDf.name = measurementChoices.value;
                    convertColToInt(measurementsDf, 'counter');
                    convertColToString(measurementsDf, PERSON_ID);
                    dynamicComparedToBaseline(measurementsDf, 'measurement_value', '1', 'counter', 'changes_from_bl');
                    joinCohorts(measurementsDf);
                    ui.setUpdateIndicator(measurementView.root, false);
                    measurementView.dataFrame = measurementsDf;
                    grok.data.linkTables(cohorts.cohortsPivoted, measurementView.dataFrame,
                      [ PERSON_ID ], [ PERSON_ID ],
                      [ DG.SYNC_TYPE.FILTER_TO_FILTER ]);
                    if(selectedValues.length > 1) {
                      updateDivInnerHTML(findCorrDiv, findCorrButton);
                    }
                  })
            })
            .show();
        });
  
        ui.setUpdateIndicator(measurementView.root, false);
        measurementView.setRibbonPanels([
          [
            selectMeasurementsButton
          ],
          [
            findCorrDiv
          ]
        ]);
      }); 
}