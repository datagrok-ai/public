/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {getMonomerPairs, getRules, LinkRuleRowArgs, Rules} from './pt-rules';
import {_package, PackageFunctions} from '../../package';
import {getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {doPolyToolConvert} from './pt-conversion';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {RuleCards} from './pt-rule-cards';
import {getOverriddenLibrary} from './pt-synthetic';
import {MmcrTemps} from '@datagrok-libraries/bio/src/utils/cell-renderer-consts';
import {helmToMol} from './pt-atomic';
import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';
import {getSeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {getReactionEditor, ReactionEditorProps} from './rule-reaction-editor';


const TAB_LINKS = 'Links';
const TAB_REACTIONS = 'Reactions';
const TAB_DIMERS = 'Dimers';

export class RulesManager {
  rules: Rules;
  linkRuleDataFrame: DG.DataFrame;
  synthRuleDataFrame: DG.DataFrame;
  fileName: string;
  private v: DG.ViewBase | null = null;

  homoDimerInput: DG.InputBase;
  heteroDimerInput: DG.InputBase;

  linkCards: RuleCards[];
  addLinkRulesFunc: (rowObj: LinkRuleRowArgs) => void;
  addSynthRulesFunc: (rowObj: ReactionEditorProps) => void;

  // every rule set will have its editor instance
  private static instances: Record<string, RulesManager> = {};

  protected constructor(rules: Rules, fileName: string) {
    this.rules = rules;
    const linksRes = this.rules.getLinkRulesDf();
    this.linkRuleDataFrame = linksRes.res;
    this.addLinkRulesFunc = linksRes.addNewRow;

    const synthRes = this.rules.getSynthesisRulesDf();
    this.addSynthRulesFunc = (r) => {
      synthRes.addNewRow(r);
      this.rules.setSynthesisRules(synthRes.df);
      this.synthRuleDataFrame = synthRes.df;
      this.substituteReactionGridDataFrame?.();
      this.save();
    };
    this.synthRuleDataFrame = synthRes.df;
    this.fileName = fileName;

    const homoValue = this.rules.homodimerCode ? this.rules.homodimerCode : '';
    const heteroValue = this.rules.heterodimerCode ? this.rules.heterodimerCode : '';
    this.homoDimerInput = ui.input.string('Homo dimer', {value: homoValue, onValueChanged: () => {}, nullable: false});
    this.heteroDimerInput = ui.input.string(
      'Hetero dimer', {value: heteroValue, onValueChanged: () => {}, nullable: false}
    );
  }

  async getAndAddView(): Promise<DG.ViewBase> {
    if (this.v) {
      try {
        //find active view; name is unique due to file names in name
        const activeView = Array.from(grok.shell.views).find((v) => v.name === this.v!.name);
        if (!activeView) {
          this.v.detach();
          this.v.close();
          throw new Error('View is closed, making it null in catch statement');
        }
        // switch to existing view
        grok.shell.v = activeView;
      } catch (_) {
        //here we only come if some error is caused due to double detaching, so no handling needed
        this.v = null;
      }
    }

    if (!this.v) {
      this.v = DG.View.create();
      this.v.name = `Manage Polytool Rules - ${this.fileName}`;
      this.v.append(await this.getForm());
      this.v = grok.shell.addView(this.v);
      grok.shell.v = this.v!; // just in any case, to make sure it switches
      return this.v;
    }

    return this.v;
  }

  public static async getInstance(name: string): Promise<RulesManager> {
    // every rull will have its own instance
    if (!this.instances[name]) {
      const rules = await getRules([name]);
      this.instances[name] = new RulesManager(rules, name);
    }
    return this.instances[name]!;
  }

  save(): void {
    this.rules.homodimerCode = this.homoDimerInput.value;
    this.rules.heterodimerCode = this.heteroDimerInput.value;
    this.rules.setLinkRules(this.linkRuleDataFrame);
    this.rules.setSynthesisRules(this.synthRuleDataFrame);

    const saveable = {
      homodimerCode: this.rules.homodimerCode,
      heterodimerCode: this.rules.heterodimerCode,
      linkRules: this.rules.linkRules,
      reactionRules: this.rules.reactionRules,
    };

    const json = JSON.stringify(saveable, undefined, 2);
    _package.files.writeAsText(`polytool-rules/${this.fileName}`, json);
    grok.shell.info(`Polytool rules at ${this.fileName} was updated`);
  }

  private createGridDiv(name: string, grid: DG.Grid, tooltipMsg?: string) {
    const header = ui.h1(name, 'polytool-grid-header');
    ui.tooltip.bind(header, tooltipMsg);
    header.style.marginTop = '10px';
    header.style.marginRight = '10px';
    grid.root.style.height = '100%';

    const gridDiv = ui.splitV([
      ui.box(
        header,
        {style: {maxHeight: '60px'}},
      ),
      grid.root,
    ]);

    gridDiv.style.height = '100%';

    return gridDiv;
  };

  async getLinkExamplesGrid() {
    const seqs: string[] = [];

    for (let i = 0; i < this.rules.linkRules.length; i++) {
      const code = this.rules.linkRules[i].code;
      const [firstMonomers, secondMonomers] = getMonomerPairs(this.rules.linkRules[i]);
      for (let j = 0; j < firstMonomers.length; j++) {
        const seq = `${firstMonomers[j]}(${code})-A-A-A-A-${secondMonomers[j]}(${code})`;
        seqs.push(seq);
      }
    }
    const helmHelper = await getHelmHelper();

    const [helms, isLinear, positionMaps] = doPolyToolConvert(seqs, this.rules, helmHelper);

    const initCol = DG.Column.fromStrings('Monomers', seqs);
    const helmCol = DG.Column.fromStrings('Helm', helms);

    PackageFunctions.applyNotationProviderForCyclized(initCol, '-');
    initCol.semType = DG.SEMTYPE.MACROMOLECULE;

    helmCol.semType = DG.SEMTYPE.MACROMOLECULE;
    helmCol.meta.units = NOTATION.HELM;
    helmCol.setTag(DG.TAGS.CELL_RENDERER, 'helm');
    return DG.DataFrame.fromColumns([
      initCol, helmCol
    ]).plot.grid();
  }

  async getReactionExamplesGrid() {
    const seqs: string[] = [];

    for (let i = 0; i < this.rules.reactionRules.length; i++) {
      const code = this.rules.reactionRules[i].code;
      const [firstMonomers, secondMonomers] = getMonomerPairs(this.rules.reactionRules[i]);
      for (let j = 0; j < firstMonomers.length; j++) {
        const seq = `${firstMonomers[j]}(${code})-A-A-A-A-${secondMonomers[j]}(${code})`;
        seqs.push(seq);
      }
    }
    const helmHelper = await getHelmHelper();

    const [helms, isLinear, positionMaps] = doPolyToolConvert(seqs, this.rules, helmHelper);

    const initCol = DG.Column.fromStrings('Monomers', seqs);
    const helmCol = DG.Column.fromStrings('Helm', helms);
    const df = DG.DataFrame.fromColumns([
      initCol, helmCol
    ]);
    initCol.semType = DG.SEMTYPE.MACROMOLECULE;
    const rdKitModule = await getRdKitModule();
    const seqHelper = await getSeqHelper();
    // initialize seqHandler for this column
    PackageFunctions.applyNotationProviderForCyclized(initCol, '-');
    initCol.tags[DG.TAGS.CELL_RENDERER] = 'Sequence';

    helmCol.semType = DG.SEMTYPE.MACROMOLECULE;
    helmCol.meta.units = NOTATION.HELM;
    helmCol.setTag(DG.TAGS.CELL_RENDERER, 'helm');

    const lib = await getOverriddenLibrary(this.rules);
    const resHelmColTemp = helmCol.temp;
    resHelmColTemp[MmcrTemps.overriddenLibrary] = lib;
    helmCol.temp = resHelmColTemp;

    const resMolCol = await helmToMol(helmCol, seqs,
      isLinear, true, false, false, lib, rdKitModule, seqHelper);
    resMolCol.name = `molfile(sequence)`;
    resMolCol.semType = DG.SEMTYPE.MOLECULE;
    df.columns.add(resMolCol);

    return df.plot.grid();
  }

  private substituteReactionGridDataFrame: (() => void) | null = null;

  async getForm() {
    const gridOptions: Partial<DG.IGridSettings> = {showAddNewRowIcon: false, allowEdit: false, rowHeight: 60};

    const linksGrid = this.linkRuleDataFrame.plot.grid(gridOptions);
    linksGrid.onCellDoubleClick.subscribe(() => {
      if (!linksGrid.dataFrame || linksGrid.dataFrame.currentRowIdx == -1 || linksGrid.dataFrame.currentRowIdx == undefined)
        return;
      const idx = linksGrid.dataFrame.currentRowIdx;
      const editArgs: LinkRuleRowArgs = {
        row: idx,
        code: this.linkRuleDataFrame.get('code', idx),
        firstMonomers: this.linkRuleDataFrame.get('firstMonomers', idx),
        secondMonomers: this.linkRuleDataFrame.get('secondMonomers', idx),
        firstLinkingGroup: this.linkRuleDataFrame.get('firstLinkingGroup', idx),
        secondLinkingGroup: this.linkRuleDataFrame.get('secondLinkingGroup', idx),
      };
      this.getAddNewLinkRuleDialog(editArgs);
    });

    const linksGridDiv = this.createGridDiv('Rules',
      linksGrid,
      'specification for monomers to link and linking positions');
    const linkExamples = this.createGridDiv('Examples', await this.getLinkExamplesGrid(),
      'specification for monomers to link and linking positions');
    linksGridDiv.style.width = '50%';
    linkExamples.style.width = '50%';
    const header = ui.h1('Monomers', 'polytool-grid-header');
    ui.tooltip.bind(header, 'Click different cobination to see how monomers will link');
    this.linkCards = await this.rules.getLinkCards();
    const gridDiv: HTMLElement = ui.splitV([
      ui.box(
        header,
        {style: {maxHeight: '30px'}},
      ),
      this.linkCards[0].root,
    ]);

    this.linkCards[0].render();
    await this.linkCards[0].reset();

    this.linkRuleDataFrame.currentRowIdx = 0;
    this.linkRuleDataFrame.onCurrentRowChanged.subscribe(async () => {
      const idx = this.linkRuleDataFrame.currentRowIdx;
      if (idx !== -1 && idx != undefined) {
        ui.empty(gridDiv);
        gridDiv.append(ui.splitV([
          ui.box(
            header,
            {style: {maxHeight: '30px'}},
          ),
          this.linkCards[idx].root,
        ]));
        this.linkCards[idx].render();
        await this.linkCards[idx].reset();
      }
    });

    const links = ui.splitH([linksGridDiv, gridDiv], null, true);
    const synthesisGrid = this.synthRuleDataFrame.plot.grid({showAddNewRowIcon: false, allowEdit: false, rowHeight: 130});
    synthesisGrid.onCellDoubleClick.subscribe(() => {
      if (!synthesisGrid.dataFrame || synthesisGrid.dataFrame.currentRowIdx == -1 || synthesisGrid.dataFrame.currentRowIdx == undefined)
        return;
      const idx = synthesisGrid.dataFrame.currentRowIdx;
      const editArgs: ReactionEditorProps = {
        rowIndex: idx,
        code: this.synthRuleDataFrame.get('code', idx),
        firstMonomers: this.synthRuleDataFrame.get('firstMonomers', idx).split(',').map((s: string) => s.trim()).filter((s: string) => s),
        secondMonomers: this.synthRuleDataFrame.get('secondMonomers', idx).split(',').map((s: string) => s.trim()).filter((s: string) => s),
        resultMonomerName: this.synthRuleDataFrame.get('name', idx),
        firstReactantSmiles: this.synthRuleDataFrame.get('firstReactant', idx),
        secondReactantSmiles: this.synthRuleDataFrame.get('secondReactant', idx),
        productSmiles: this.synthRuleDataFrame.get('product', idx),
      };
      getReactionEditor((props: ReactionEditorProps) => this.addSynthRulesFunc(props), editArgs);
    });
    const reactionsGridDiv = this.createGridDiv('Rules', synthesisGrid);

    const reactionExamplesGrid = await this.getReactionExamplesGrid();
    const reactionExamples = this.createGridDiv('Examples', reactionExamplesGrid);
    this.substituteReactionGridDataFrame = async () => {
      const newGrid = await this.getReactionExamplesGrid();
      reactionExamplesGrid.dataFrame = newGrid.dataFrame;
    };
    reactionsGridDiv.style.width = '50%';
    reactionExamples.style.width = '50%';
    const reactions = ui.divH([reactionsGridDiv, reactionExamples]);

    const dimerInputsDiv = ui.divV([
      this.homoDimerInput,
      this.heteroDimerInput,
    ]);

    const inputsTabControl = ui.tabControl({
      'Links': links,
      'Reactions': reactions,
      'Dimers': dimerInputsDiv,
    }, false);

    ui.tooltip.bind(inputsTabControl.getPane(TAB_LINKS).header,
      'Specify rules to link monomers based on HELM notation');
    ui.tooltip.bind(inputsTabControl.getPane(TAB_REACTIONS).header,
      'Specify rules to perform reactions within monomers');
    ui.tooltip.bind(inputsTabControl.getPane(TAB_DIMERS).header,
      'Specify symbols for homodimeric and heterodimeric codes');

    inputsTabControl.root.style.height = '100%';
    inputsTabControl.root.style.width = '100%';
    inputsTabControl.root.classList.add('rules-manager-form-tab-control');
    inputsTabControl.header.style.marginBottom = '10px';

    const panel = ui.divV([
      inputsTabControl.root,
    ]);

    const saveButton = ui.bigButton('Save', () => { this.save(); });
    const addButton = ui.button('Add rule', () => {
      const currentTab = inputsTabControl.currentPane.name;
      if (currentTab == TAB_LINKS)
        this.getAddNewLinkRuleDialog();
      else if (currentTab == TAB_REACTIONS)
        getReactionEditor((props: ReactionEditorProps) => this.addSynthRulesFunc(props));
    });
    const removeButton = ui.button('Remove rule', () => {
      const currentTab = inputsTabControl.currentPane.name;
      if (currentTab == TAB_LINKS) {
        if (this.linkRuleDataFrame == null || this.linkRuleDataFrame.currentRowIdx == -1 || this.linkRuleDataFrame.currentRowIdx == undefined)
          return;
        const idx = linksGrid.dataFrame.currentRowIdx;
        ui.dialog('Are you sure you want to remove the rule?')
          .add(ui.divText('This action is irreversible!'))
          .onOK(() => {
            this.linkRuleDataFrame.rows.removeAt(idx);
            this.rules.setLinkRules(this.linkRuleDataFrame);
            this.save();
          }).show();
      } else if (currentTab == TAB_REACTIONS) {
        if (this.synthRuleDataFrame == null || this.synthRuleDataFrame.currentRowIdx == -1 || this.synthRuleDataFrame.currentRowIdx == undefined)
          return;
        const idx = synthesisGrid.dataFrame.currentRowIdx;
        ui.dialog('Are you sure you want to remove the rule?')
          .add(ui.divText('This action is irreversible!'))
          .onOK(() => {
            this.synthRuleDataFrame.rows.removeAt(idx);
            this.rules.setSynthesisRules(this.synthRuleDataFrame);
            this.substituteReactionGridDataFrame?.();
            this.save();
          }).show();
      }
    });

    const topPanel = [saveButton, addButton, removeButton];
    this.v!.setRibbonPanels([topPanel]);

    panel.style.height = '100%';
    panel.style.alignItems = 'center';

    return inputsTabControl.root;
  }

  getAddNewLinkRuleDialog(preset?: Partial<LinkRuleRowArgs>): void {
    const codeInput = ui.input.int('Code', {nullable: false, value: preset?.code});
    const firstMonomersInput = ui.input.string('First monomers', {placeholder: 'E.g. C,D,E', value: preset?.firstMonomers,
      tooltipText: 'Comma separated list of first monomers applicable for the rule. If left empty, all monomers will be considered', nullable: true});
    const secondMonomersInput = ui.input.string('Second monomers', {placeholder: 'E.g. C,D,E', value: preset?.secondMonomers,
      tooltipText: 'Comma separated list of second monomers applicable for the rule. If left empty, all monomers will be considered', nullable: true});
    const firstLinkingGroup = preset?.firstLinkingGroup ? `R${preset.firstLinkingGroup}` : 'R3';
    const secondLinkingGroup = preset?.secondLinkingGroup ? `R${preset.secondLinkingGroup}` : 'R3';
    const firstLinkingGroupInput = ui.input.choice('First linking group', {value: firstLinkingGroup, items: ['R1', 'R2', 'R3', 'R4'],
      tooltipText: 'Specifies which R-group of the first monomer will be used for linking', nullable: false});
    const secondLinkingGroupInput = ui.input.choice('Second linking group', {value: secondLinkingGroup, items: ['R1', 'R2', 'R3', 'R4'],
      tooltipText: 'Specifies which R-group of the second monomer will be used for linking', nullable: false});
    ui.dialog('Add new link rule')
      .add(codeInput)
      .add(firstMonomersInput)
      .add(secondMonomersInput)
      .add(firstLinkingGroupInput)
      .add(secondLinkingGroupInput)
      .onOK(async () => {
        // we rely on validation of inputs by DG inputs
        const code = codeInput.value!;
        const firstMonomers = (firstMonomersInput.value ?? '').split(',').map((s) => s.trim()).filter((s) => s).join(',');
        const secondMonomers = (secondMonomersInput.value ?? '').split(',').map((s) => s.trim()).filter((s) => s).join(',');
        const firstLinkingGroup = parseInt(firstLinkingGroupInput.value!.substring(1));
        const secondLinkingGroup = parseInt(secondLinkingGroupInput.value!.substring(1));
        this.addLinkRulesFunc({code, firstMonomers: firstMonomers ?? '', secondMonomers: secondMonomers ?? '',
          firstLinkingGroup, secondLinkingGroup, row: preset?.row});

        this.rules.setLinkRules(this.linkRuleDataFrame);
        this.linkCards = await this.rules.getLinkCards();
        this.save();
      }).show();
  }
}

