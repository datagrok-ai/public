import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { updateDivInnerHTML } from '../utils';

export class AchillesResultsView extends DG.ViewBase {

    analysesNames: string [] = [];
    categoriesNames: string [] = [];
    resultsDiv = ui.box();
    analysesChoicesDiv = ui.div();

    constructor(name) {
        super({});
        this.name = name;
        this.createView();
    }

    async createView(): Promise<void> {

        grok.data.query(`CDM:achillesResults`)
            .then((resultsDf: DG.DataFrame) => {
               this.categoriesNames = resultsDf.col('category').categories;
               let categoriesChoices = ui.choiceInput('Category', '', this.categoriesNames);
               categoriesChoices.input.style.width = '150px';
               categoriesChoices.onChanged(() => {
                this.analysesNames = resultsDf.groupBy(['category', 'analysis_name']).aggregate().col('analysis_name').categories;
                let analysesChoices = ui.choiceInput('Analysis', '', this.analysesNames);
                analysesChoices.input.style.width = '150px';
                analysesChoices.onChanged((v) => {
                   let analysisDf = resultsDf.groupBy(resultsDf.columns.names())
                   .where({ analysis_name: analysesChoices.value })
                   .aggregate();
                   let formattedResultsDf = DG.DataFrame.create(analysisDf.rowCount);
                   const stratumName = this.getConcatenatedStratums(analysisDf, 0, true);
                   formattedResultsDf.columns.addNewString(stratumName).init((i) => this.getConcatenatedStratums(analysisDf, i));
                   formattedResultsDf.columns.addNewInt('count').init((i) => analysisDf.get('count_value', i));
                   /* if (stratumName.endsWith('_concept_id')) {
                    const concept_ids = formattedResultsDf.col(stratumName).categories.join(',');
                    grok.data.query(`CDM:conceptsByNumbers`, { ids: concept_ids }).then((conceptNamesDf) => {
                        formattedResultsDf.join(conceptNamesDf, 
                            [stratumName], ['concept_id'], formattedResultsDf.columns.names(), ['concept_name'], DG.JOIN_TYPE.LEFT, true);
                    })
                   } */
                   updateDivInnerHTML(this.resultsDiv, formattedResultsDf.plot.grid().root);
                   grok.shell.addTableView(formattedResultsDf);
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

    private getConcatenatedStratums(df: DG.DataFrame, idx: number, stratumName?: boolean) {
        let concatenated = '';
        for (let i = 1; i < 6; i++) {
            let colname = stratumName ? `stratum_${i}_name` : `stratum_${i}`; 
            let stratum = df.get(colname, idx);
            if (stratum) {
                concatenated = concatenated === '' ? stratum : concatenated += `|${stratum}`;
            } else
                break;
          }
          return concatenated;
    }

}