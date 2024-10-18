import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {dfFromSynthesisRules, getRules, Rules, synthesisRulesFromDf} from './pt-rules';
import {_package} from '../../package';

export class RulesManager {
  rules: Rules;
  ruleDataFrame: DG.DataFrame;
  fileName: string;
  tv: DG.TableView | null = null;

  private static instance: RulesManager;

  protected constructor(rules: Rules, fileName: string) {
    this.rules = rules;
    this.ruleDataFrame = dfFromSynthesisRules(this.rules.reactionRules);
    this.fileName = fileName;
  }

  show(): void {
    if (!this.tv)
      this.tv = DG.TableView.create(this.ruleDataFrame, false);

    grok.shell.v = this.tv;
    const form = this.getForm();
    this.tv.dockManager.dock(form, DG.DOCK_TYPE.LEFT, null, '', 0.4);
  }

  public static async getInstance(name: string): Promise<RulesManager> {
    if (!this.instance) {
      const rules = await getRules([name]);
      this.instance = new RulesManager(rules, name);
    }
    return this.instance;
  }

  getForm() {
    inputsTabControl: DG.TabControl;
    const homoValue = this.rules.homodimerCode ? this.rules.homodimerCode : '';
    const heteroValue = this.rules.heterodimerCode ? this.rules.heterodimerCode : '';
    const homoDimerInput = ui.input.string('Homo dimer', {value: homoValue, onValueChanged: () => {}, nullable: false});
    const heteroDimerInput = ui.input.string('Hetero dimer', {value: heteroValue, onValueChanged: () => {}, nullable: false});

    const dimerInputsDiv = ui.divV([
      homoDimerInput,
      heteroDimerInput,
    ]);

    const dimerInputsDiv2 = ui.divV([]);

    const inputsTabControl = ui.tabControl({
      'Dimers': dimerInputsDiv,
      'Reactions': dimerInputsDiv2
    }, false);

    inputsTabControl.root.classList.add('rules-manager-form-tab-control');
    inputsTabControl.header.style.marginBottom = '10px';


    const saveButton = ui.button('Save changes', () => {
      this.rules.homodimerCode = homoDimerInput.value;
      this.rules.heterodimerCode = heteroDimerInput.value;
      this.rules.reactionRules = synthesisRulesFromDf(this.ruleDataFrame);

      const rrrr = JSON.stringify(this.rules, undefined, 2);
      _package.files.writeAsText(`polytool-rules/${this.fileName}`, rrrr);
      grok.shell.info(`Polytool rules at ${this.fileName} was updated`);
    });

    return ui.divV([
      inputsTabControl.root,
      saveButton
    ], {classes: 'ui-form', style: {paddingLeft: '10px', overflow: 'scroll'}});
  }
}

