import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {ActiveFiles} from '@datagrok-libraries/utils/src/settings/active-files-base';
import {RulesManager} from './rule-manager';

export const RULES_PATH = 'System:AppData/SequenceTranslator/polytool-rules/';
export const RULES_STORAGE_NAME = 'Polytool';
export const RULES_TYPE_LINK = 'link';
export const RULES_TYPE_REACTION = 'reaction';
export const RULES_TYPE_HOMODIMER = 'fragmentDuplication';
export const RULES_TYPE_HETERODIMER = 'differentFragments';

export class RuleInputs extends ActiveFiles {
  constructor(
    path: string, userStorageName: string, ext: string,
    options?: { onValueChanged: (value: string[]) => void }
  ) {
    super(path, userStorageName, ext, options);
  }

  override createInput(available: string, isChecked: boolean): DG.InputBase<boolean> {
    const res = super.createInput(available, isChecked);

    const editIcon = ui.icons.edit(async () => {
      const rulesManager = await RulesManager.getInstance(available);
      await rulesManager.show();
    }, 'Edit rules');

    res.addOptions(editIcon);

    return res;
  }
}

export type Rules = {
  homodimerCode: string | null,
  heterodimerCode: string | null,
  linkRules: RuleLink[],
  reactionRules: RuleReaction[]
}

export type RuleLink = {
  code: number,
  firstMonomer: string,
  secondMonomer: string,
  firstSubstitution: string,
  secondSubstitution: string,
  firstLinkingGroup: number,
  secondLinkingGroup: number
}

export type RuleReaction = {
  code: number,
  firstMonomer: string,
  secondMonomer: string,
  reaction: string,
  name: string
}

export function dfFromSynthesisRules(rules: RuleReaction []) : DG.DataFrame {
  const length = rules.length;
  const codeCol = DG.Column.int('code', length);
  const firstMonomerCol = DG.Column.string('firstMonomer', length);
  const secondMonomerCol = DG.Column.string('secondMonomer', length);
  const name = DG.Column.string('name', length);
  const firstReactant = DG.Column.string('firstReactant', length);
  const secondReactant = DG.Column.string('secondReactant', length);
  const product = DG.Column.string('product', length);

  for (let i = 0; i < length; i++) {
    codeCol.set(i, rules[i].code);
    firstMonomerCol.set(i, rules[i].firstMonomer);
    secondMonomerCol.set(i, rules[i].secondMonomer);
    name.set(i, rules[i].name);

    const reaction = rules[i].reaction.split('>>');
    const reactants = reaction[0].split('.');

    firstReactant.set(i, reactants[0]);
    secondReactant.set(i, reactants[1]);
    product.set(i, reaction[1]);
  }
  firstReactant.semType = DG.SEMTYPE.MOLECULE;
  secondReactant.semType = DG.SEMTYPE.MOLECULE;
  product.semType = DG.SEMTYPE.MOLECULE;


  return DG.DataFrame.fromColumns([name, firstReactant, secondReactant, product, codeCol, firstMonomerCol, secondMonomerCol]);
}

export function synthesisRulesFromDf(df: DG.DataFrame) : RuleReaction [] {
  const length = df.rowCount;
  const rules: RuleReaction [] = new Array<RuleReaction>(length);
  const codeCol = df.columns.byName('code');
  const firstMonomerCol = df.columns.byName('firstMonomer');
  const secondMonomerCol = df.columns.byName('secondMonomer');
  const name = df.columns.byName('name');
  const firstReactant = df.columns.byName('firstReactant');
  const secondReactant = df.columns.byName('secondReactant');
  const product = df.columns.byName('product');

  for (let i = 0; i < length; i++) {
    const smartsReaction = `${firstReactant.get(i)}.${secondReactant.get(i)}>>${product.get(i)}`;

    const rule = {
      code: codeCol.get(i),
      firstMonomer: firstMonomerCol.get(i),
      secondMonomer: secondMonomerCol.get(i),
      reaction: smartsReaction,
      name: name.get(i)
    };

    rules[i] = rule;
  }

  return rules;
}

export async function getRules(ruleFiles: string[]): Promise<Rules> {
  const fileSource = new DG.FileSource(RULES_PATH);
  const linkRules: RuleLink[] = [];
  const reactionRules: RuleReaction[] = [];
  const rules: Rules = {homodimerCode: null, heterodimerCode: null, linkRules: linkRules, reactionRules: reactionRules};

  for (let i = 0; i < ruleFiles.length; i++) {
    const rulesRaw = await fileSource.readAsText(ruleFiles[i].replace(RULES_PATH, ''));
    const ruleSingle : Rules = JSON.parse(rulesRaw);

    rules.homodimerCode = ruleSingle.homodimerCode;
    rules.heterodimerCode = ruleSingle.heterodimerCode;

    for (let j = 0; j < ruleSingle.linkRules.length; j++)
      linkRules.push(ruleSingle.linkRules[j]);

    for (let j = 0; j < ruleSingle.reactionRules.length; j++)
      reactionRules.push(ruleSingle.reactionRules[j]);
  }

  return rules;
}
