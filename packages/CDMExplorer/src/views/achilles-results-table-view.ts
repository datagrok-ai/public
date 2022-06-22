import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from "datagrok-api/ui";
import { addView, sortObject, updateDivInnerHTML } from '../utils';
import { convertColToInt, convertColToString, dynamicComparedToBaseline, joinCohorts } from '../preprocessing/data-preparation';
import { PERSON_ID, VIEWS } from '../constants';
import { cohorts } from '../cohorts';

export function achillesResultsTableView(): void {
    const stratumCount = 5;
    const achillesResultsView = DG.TableView.create(DG.DataFrame.create());
    achillesResultsView.name = 'Achilles Results';
    VIEWS.push(addView(achillesResultsView));
    ui.setUpdateIndicator(achillesResultsView.root, true);
    grok.data.query(`CDM:achillesAnalysis`)
            .then((analysisDf: DG.DataFrame) => {
               const analysesChoicesDiv = ui.div();
               const categoriesNames = analysisDf.col('category').categories;
               let categoriesChoices = ui.choiceInput('Category', '', categoriesNames);
               categoriesChoices.input.style.width = '150px';

               categoriesChoices.onChanged(() => {
               const analysesNames = analysisDf
                  .groupBy(['category', 'analysis_name'])
                  .where({category: `${categoriesChoices.value}`})
                  .aggregate()
                  .col('analysis_name').categories;
                let analysesChoices = ui.choiceInput('Analysis', '', analysesNames);
                analysesChoices.input.style.width = '150px';

                analysesChoices.onChanged(async (v) => {
                    ui.setUpdateIndicator(achillesResultsView.root, true);
                    const selectedAnalysisDf = analysisDf
                    .groupBy(analysisDf.columns.names())
                    .where({category: categoriesChoices.value, analysis_name: analysesChoices.value})
                    .aggregate();
                    const analysisId = selectedAnalysisDf.get('analysis_id', 0);
                    const resultsDf: DG.DataFrame = await grok.data.query(`CDM:achillesAnalysisByIds`, {analysis_id: analysisId });
                    let conceptIds = '';
                    let exisitngStratums = 0;
                    for (let i = 1; i <= stratumCount; i++){
                        if (selectedAnalysisDf.col(`stratum_${i}_name`).isNone(0)) {
                            resultsDf.columns.remove(`stratum_${i}`);
                        } else {
                            exisitngStratums += 1;
                            const stratumName = selectedAnalysisDf.get(`stratum_${i}_name`, 0);
                            if (stratumName.endsWith('_concept_id')) {
                                const stratumConcepts = resultsDf.col(`stratum_${i}`).categories;
                                conceptIds += stratumConcepts.join(', ');
                                resultsDf.columns.addNewInt(`stratum_${i}_num`).init((j) => parseInt(resultsDf.get(`stratum_${i}`, j)));
                            }
                        }
                    }
                    let conceptNamesDf: DG.DataFrame = null;
                    if (conceptIds)
                      conceptNamesDf = await grok.data.query(`CDM:conceptsByIds`, {ids: `in(${conceptIds})`});
                    for (let i = 1; i <= exisitngStratums; i++){
                        if (resultsDf.col(`stratum_${i}_num`)) {
                            resultsDf.join(conceptNamesDf,
                                [`stratum_${i}_num`], ['concept_id'],
                                resultsDf.columns.names(), ['concept_name'], DG.JOIN_TYPE.LEFT, true);
                            resultsDf.col('concept_name').name = selectedAnalysisDf.get(`stratum_${i}_name`, 0).replace('_concept_id', '');
                            resultsDf.columns.remove(`stratum_${i}_num`);
                            resultsDf.col(`stratum_${i}`).name = selectedAnalysisDf.get(`stratum_${i}_name`, 0);
                        } else {
                            resultsDf.col(`stratum_${i}`).name = selectedAnalysisDf.get(`stratum_${i}_name`, 0);
                        }
                    }
                    ui.setUpdateIndicator(achillesResultsView.root, false);
                    achillesResultsView.dataFrame = resultsDf;
                });
                updateDivInnerHTML(analysesChoicesDiv, analysesChoices.root);
               })
               ui.setUpdateIndicator(achillesResultsView.root, false);
               achillesResultsView.setRibbonPanels([
                [
                    categoriesChoices.root,
                    analysesChoicesDiv
                ],
              ]);
            });
}