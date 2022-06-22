import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { updateDivInnerHTML } from '../utils';

export class   AchillesResultsView extends DG.ViewBase {

    analysesNames: string [] = [];
    categoriesNames: string [] = [];
    resultsDiv = ui.box();
    analysesChoicesDiv = ui.div();
    stratumCount = 5;

    constructor(name) {
        super({});
        this.name = name;
        this.createView();
    }

    async createView(): Promise<void> {

        grok.data.query(`CDM:achillesAnalysis`)
            .then((analysisDf: DG.DataFrame) => {
               this.categoriesNames = analysisDf.col('category').categories;
               let categoriesChoices = ui.choiceInput('Category', '', this.categoriesNames);
               categoriesChoices.input.style.width = '150px';
               categoriesChoices.onChanged(() => {
                this.analysesNames = analysisDf
                  .groupBy(['category', 'analysis_name'])
                  .where({category: `${categoriesChoices.value}`})
                  .aggregate()
                  .col('analysis_name').categories;
                let analysesChoices = ui.choiceInput('Analysis', '', this.analysesNames);
                analysesChoices.input.style.width = '150px';
                analysesChoices.onChanged(async (v) => {

                    const selectedAnalysisDf = analysisDf
                    .groupBy(analysisDf.columns.names())
                    .where({category: categoriesChoices.value, analysis_name: analysesChoices.value})
                    .aggregate();

                    const analysisId = selectedAnalysisDf.get('analysis_id', 0);

                    const resultsDf: DG.DataFrame = await grok.data.query(`CDM:achillesAnalysisByIds`, {analysis_id: analysisId });

                    let conceptIds = '';
                    let exisitngStratums = 0;
                    for (let i = 1; i <= this.stratumCount; i++){
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
                      conceptNamesDf = await grok.data.query(`CDM:conceptsByIds`, {ids: conceptIds});
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
                    const view = grok.shell.addTableView(resultsDf);
                    view.name = `${analysesChoices.value}`;
                });
                updateDivInnerHTML(this.analysesChoicesDiv, analysesChoices.root);
               })

               this.setRibbonPanels([
                [
                    categoriesChoices.root,
                    this.analysesChoicesDiv
                ],
              ]); 
              this.root.append(this.resultsDiv);
            });
    }

}