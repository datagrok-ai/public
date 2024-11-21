import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {getMonomerPairs, getRules, Rules} from './pt-rules';
import {_package, applyNotationProviderForCyclized} from '../../package';
import {getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {doPolyToolConvert} from './pt-conversion';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';

export class RulesManager {
  rules: Rules;
  linkRuleDataFrame: DG.DataFrame;
  synthRuleDataFrame: DG.DataFrame;
  fileName: string;
  private v: DG.ViewBase | null = null;

  homoDimerInput: DG.InputBase;
  heteroDimerInput: DG.InputBase;

  currentTab = '';

  // every rule set will have its editor instance
  private static instances: Record<string, RulesManager> = {};

  protected constructor(rules: Rules, fileName: string) {
    this.rules = rules;
    this.linkRuleDataFrame = this.rules.getLinkRulesDf();
    this.synthRuleDataFrame = this.rules.getSynthesisRulesDf();
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

  private createGridDiv(name: string, grid: DG.Grid, buttons: HTMLButtonElement[] = []) {
    const header = ui.h1(name, 'polytool-grid-header');
    header.style.marginTop = '10px';
    header.style.marginRight = '10px';
    grid.root.style.height = '100%';

    buttons.forEach((b) => {
      b.style.marginLeft = '5px';
      b.style.marginRight = '10px';
    });

    const gridDiv = ui.splitV([
      ui.box(
        ui.divH([header, ...buttons]),
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

    const helms = doPolyToolConvert(seqs, this.rules, helmHelper);

    const initCol = DG.Column.fromStrings('monomers', seqs);
    const helmCol = DG.Column.fromStrings('helm', helms);

    applyNotationProviderForCyclized(initCol, '-');
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

    const helms = doPolyToolConvert(seqs, this.rules, helmHelper);

    const initCol = DG.Column.fromStrings('monomers', seqs);
    const helmCol = DG.Column.fromStrings('helm', helms);

    initCol.semType = DG.SEMTYPE.MACROMOLECULE;
    applyNotationProviderForCyclized(initCol, '-');

    helmCol.semType = DG.SEMTYPE.MACROMOLECULE;
    helmCol.meta.units = NOTATION.HELM;
    helmCol.setTag(DG.TAGS.CELL_RENDERER, 'helm');
    return DG.DataFrame.fromColumns([
      initCol, helmCol
    ]).plot.grid();
  }

  async getForm() {
    inputsTabControl: DG.TabControl;

    const saveLinksButton = ui.bigButton('Save', () => { this.save(); });
    const addLinkButton = ui.button('Add rule', () => { this.linkRuleDataFrame.rows.addNew(); });
    const linksGridDiv = this.createGridDiv('Rules',
      this.linkRuleDataFrame.plot.grid({showAddNewRowIcon: true}), [saveLinksButton, addLinkButton]);
    const linkExamples = this.createGridDiv('Examples', await this.getLinkExamplesGrid());
    linksGridDiv.style.width = '50%';
    linkExamples.style.width = '50%';
    const links = ui.divH([linksGridDiv, linkExamples]);

    const saveReactionsButton = ui.bigButton('Save', () => { this.save(); });
    const addReactionButton = ui.button('Add rule', () => { this.synthRuleDataFrame.rows.addNew(); });
    const reactionsGridDiv = this.createGridDiv('Rules',
      this.synthRuleDataFrame.plot.grid({showAddNewRowIcon: true}), [saveReactionsButton, addReactionButton]);
    const reactionExamples = this.createGridDiv('Examples', await this.getReactionExamplesGrid());
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

    inputsTabControl.root.style.height = '100%';
    inputsTabControl.root.style.width = '100%';
    inputsTabControl.root.classList.add('rules-manager-form-tab-control');
    inputsTabControl.header.style.marginBottom = '10px';

    const panel = ui.divV([
      inputsTabControl.root,
    ]);

    panel.style.height = '100%';
    panel.style.alignItems = 'center';

    return inputsTabControl.root;
  }
}

