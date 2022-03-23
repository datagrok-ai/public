import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { addDataFromDmDomain } from '../data-preparation/utils';
import { ETHNIC, QS_CATEGORY, QS_RES, QS_RES_N, QS_SUB_CATEGORY, QS_TEST, RACE, SEX, VISIT_NUM } from '../constants/columns-constants';
import { updateDivInnerHTML } from '../utils/utils';
import { _package } from '../package';
import { ClinicalCaseViewBase } from '../model/ClinicalCaseViewBase';
import { TRT_ARM_FIELD, VIEWS_CONFIG } from '../views-config';


export class QuestionnaiesView extends ClinicalCaseViewBase {

  qsCategories = [];
  qsSubCategories = [];
  qsWithDm: DG.DataFrame;
  subCategoriesDiv = ui.div();
  categoriesChoice: DG.InputBase;
  subCategoriesChoice: DG.InputBase;
  selectedCategory = '';
  selectedSubCategory = '';
  questionsDiv = ui.box();
  questionsAcc: DG.Accordion;
  questions = [];
  questionsDf: DG.DataFrame;
  openedQuestionsDfs = {};
  switchButton: any;
  graphTypes = [ 'histogram', 'linechart' ];
  graphTypesDiv = ui.div();
  graphTypesChoice: DG.InputBase;
  selectedGraphType: string;
  splitByChoice: DG.InputBase;
  splitBy: string[];
  selectedSplitBy = '';

  constructor(name) {
    super({});
    this.name = name;
  }

  createView(): void {
    this.splitBy = [SEX, RACE, ETHNIC, [VIEWS_CONFIG[this.name][TRT_ARM_FIELD]]].filter(it => study.domains.dm.columns.names().includes(it)) as string[];
    this.selectedSplitBy = this.splitBy.length ? this.splitBy[0] : '';
    this.qsWithDm = addDataFromDmDomain(study.domains.qs, study.domains.dm, study.domains.qs.columns.names(), this.splitBy);
    this.qsCategories = this.qsWithDm.col(QS_CATEGORY).categories;
    this.selectedCategory = this.qsCategories[0];
    this.selectedGraphType = this.graphTypes[0];
    this.updateSubCategoriesChoice();
    this.updateQuestionsAcc();

    this.categoriesChoice = ui.choiceInput('Category', this.selectedCategory, this.qsCategories);
    this.categoriesChoice.onChanged((v) => {
      this.selectedCategory = this.categoriesChoice.value;
      this.updateSubCategoriesChoice();
      this.switchButton.value ? this.updateQuestionsAcc() : this.updateQuestionsPanel();
    });
    this.categoriesChoice.input.style.width = '100px';

    this.graphTypesChoice = ui.choiceInput('Graph type', this.selectedGraphType, this.graphTypes);
    this.graphTypesChoice.onChanged((v) => {
      this.selectedGraphType = this.graphTypesChoice.value;
      this.updateQuestionsPanel();
    });
    this.categoriesChoice.input.style.width = '100px';

    this.splitByChoice = ui.choiceInput('Split By', this.selectedSplitBy, this.splitBy);
    this.splitByChoice.onChanged((v) => {
      this.selectedSplitBy = this.splitByChoice.value;
      this.switchButton.value ? this.updateQuestionsAcc() : this.updateQuestionsPanel();
    });
    this.splitByChoice.input.style.width = '100px';

    this.switchButton = ui.switchInput('Accordion', true);
    this.switchButton.onChanged((v) => {
      if (this.switchButton.value) {
        this.updateQuestionsAcc();
        updateDivInnerHTML(this.graphTypesDiv, '');
      } else {
        this.updateQuestionsPanel();
        updateDivInnerHTML(this.graphTypesDiv, this.graphTypesChoice.root);

      }
    });

    this.setRibbonPanels(
      [[this.categoriesChoice.root], [this.subCategoriesDiv], [this.switchButton.root], [this.graphTypesDiv], [this.splitByChoice.root]]
    );

    this.root.className = 'grok-view ui-box';
    this.root.append(this.questionsDiv);

  }

  private updateSubCategories() {
    this.qsSubCategories = this.qsWithDm
      .groupBy([QS_CATEGORY, QS_SUB_CATEGORY])
      .where({ [QS_CATEGORY]: `${this.selectedCategory}` })
      .aggregate()
      .col(QS_SUB_CATEGORY)
      .categories;
  }

  private updateSubCategoriesChoice() {
    this.updateSubCategories();
    this.selectedSubCategory = this.qsSubCategories[0];
    this.subCategoriesChoice = ui.choiceInput('Sub Category', this.selectedSubCategory, this.qsSubCategories);
    this.subCategoriesChoice.onChanged((v) => {
      this.selectedSubCategory = this.subCategoriesChoice.value;
      this.switchButton.value ? this.updateQuestionsAcc() : this.updateQuestionsPanel();
    });
    this.subCategoriesChoice.input.style.width = '150px';
    updateDivInnerHTML(this.subCategoriesDiv, this.subCategoriesChoice.root);
  }

  private updateQuestionsAcc() {
    this.questionsAcc = ui.accordion();
    grok.shell.o = this.setPropertyPanel();
    this.updateQuestions();
    this.questions.forEach(question => {
      this.questionsAcc.addPane(`${question}`, () => this.createQuestionCharts(question));
      let questionPanel = this.questionsAcc.getPane(`${question}`);
      questionPanel.root.addEventListener('click', (event) => {
        let paneHeaderClass = 'd4-accordion-pane-header';
        if(event.composedPath()[0]['classList'].contains(paneHeaderClass)) {
          grok.shell.o = this.setPropertyPanel();
        };
      });
    });
    updateDivInnerHTML(this.questionsDiv, this.questionsAcc.root);
  }

  private updateQuestionsPanel() {
    this.updateQuestions();
    let questionsPanels = [];
    this.questions.forEach(question => {
      let questionDf = this.createQuestionDf(question);
      let questionGraph = this.selectedGraphType === 'histogram' ? this.createHistogram(questionDf) : this.createLineChart(questionDf);
      questionGraph.prepend(ui.divText(question));
      questionsPanels.push(ui.block25([questionGraph]));
    });
    updateDivInnerHTML(this.questionsDiv, ui.block(questionsPanels));
  }

  private updateQuestions() {
    this.openedQuestionsDfs = {};
    this.questionsDf = this.qsWithDm
      .groupBy([QS_CATEGORY, QS_SUB_CATEGORY, QS_TEST, QS_RES, QS_RES_N, VISIT_NUM].concat(this.splitBy))
      .where({ [QS_CATEGORY]: `${this.selectedCategory}`, [QS_SUB_CATEGORY]: `${this.selectedSubCategory}` })
      .aggregate();
    this.questions = this.questionsDf
      .col(QS_TEST)
      .categories;
  }

  private createQuestionCharts(question: string) {
    let df = this.createQuestionDf(question);
    this.openedQuestionsDfs[question] = df;
    let chart = this.isCategorical(df) ? this.createSplitByBarcharts(df) : this.createLineChart(df);
    return ui.splitH([chart], { style: { width: '100%' } });
  }

  private createQuestionDf(question: string) {
    return this.questionsDf
      .groupBy(this.questionsDf.columns.names())
      .where({ [QS_TEST]: `${question}` })
      .aggregate();
  }

  private createHistogram(questionDf: DG.DataFrame) {
    let histogram = DG.Viewer.fromType(DG.VIEWER.BAR_CHART, questionDf, {
      split: VISIT_NUM,
      barSortType: 'by category',
      barSortOrder: 'asc',
      stack: QS_RES,
    }).root;
    return histogram;
  }

  private createLineChart(questionDf: DG.DataFrame) {
    let linechart = DG.Viewer.fromType(DG.VIEWER.LINE_CHART, questionDf, {
      splitColumnName: this.selectedSplitBy,
      xColumnName: VISIT_NUM,
      yColumnNames: [QS_RES_N],
      whiskersType: 'Med | Q1, Q3',
    }).root;
    return linechart;
  }


  private isCategorical(df: DG.DataFrame){
    return df.col(QS_RES_N).isNone(0);
  }

  private createSplitByBarcharts(df: DG.DataFrame){
    const splittedDiv = ui.splitH([]);
    const splittedDf = df.groupBy([ this.selectedSplitBy ]).getGroups();
    Object.keys(splittedDf).forEach(cat => {
      const divWithCategory = ui.divV([
        this.createHistogram(splittedDf[cat]),
        ui.divText(`${cat}`)
      ])
      splittedDiv.append(divWithCategory);
    });
    return splittedDiv;
  }

  private setPropertyPanel() {
    let propPanel = async () => {
        grok.shell.o = await this.propertyPanel();
    }
    propPanel();
  }


  override async propertyPanel() {

    const acc = this.createAccWithTitle(this.name);

    const expandedPanes = this.questionsAcc.panes.filter(it => it.expanded);

    let getPaneContent = (df, rowNum) => {
      if (!rowNum) {
        return ui.divText('No records found');
      } else {
        let grid = df.plot.grid();
        grid.root.style.width = '250px';
        return ui.div(grid.root);
      }
  }

    expandedPanes.forEach(pane => {
      const df = this.openedQuestionsDfs[pane.name];
      const rowNum = df.rowCount === 1 && df.col(QS_CATEGORY).isNone(0) ? 0 : df.rowCount;
      acc.addCountPane(`${pane.name}`,
      () => getPaneContent(df, rowNum),
      () => rowNum,
      true);

    let dfPanel = acc.getPane(`${pane.name}`);
    //@ts-ignore
    $(dfPanel.root).css('display', 'flex');
    //@ts-ignore
    $(dfPanel.root).css('opacity', '1');
    });

    return acc.root;
  }

}