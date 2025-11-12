/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {ActiveFiles} from '@datagrok-libraries/utils/src/settings/active-files-base';
import {RulesManager} from './rule-manager';
import {RuleCards} from './pt-rule-cards';
import {getMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {
  _package,
  PackageFunctions} from '../../package';
import {MmcrTemps} from '@datagrok-libraries/bio/src/utils/cell-renderer-consts';
import {ReactionEditorProps} from './rule-reaction-editor';

export const RULES_PATH = 'System:AppData/SequenceTranslator/polytool-rules/';
export const RULES_STORAGE_NAME = 'Polytool';
export const RULES_TYPE_LINK = 'link';
export const RULES_TYPE_REACTION = 'reaction';
export const RULES_TYPE_HOMODIMER = 'fragmentDuplication';
export const RULES_TYPE_HETERODIMER = 'differentFragments';

const NAME_CODE = 'code';
const NAME_FIRST_MONOMERS = 'firstMonomers';
const NAME_SECOND_MONOMERS = 'secondMonomers';
const NAME_REACTION_NAME = 'name';
const NAME_FIRST_LINK = 'firstLinkingGroup';
const NAME_SECOND_LINK = 'secondLinkingGroup';

export class RuleInputs extends ActiveFiles {
  static ruleDescriptions: Record<string, Promise<HTMLElement>> = {};
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
      //await rulesManager.show();
      // close the dialogs, its easier if we just close it from active dialogs with filtering
      DG.Dialog.getOpenDialogs()?.filter((d) => d.root.contains(editIcon)).forEach((d) => d.close());
      await rulesManager.getAndAddView();
    }, 'Edit rules');

    res.addOptions(editIcon);

    return res;
  }

  override async getForm(): Promise<HTMLDivElement> {
    const form = await super.getForm();
    this.processRulesForm(form);
    return form;
  }

  private async processRulesForm(rulesForm: HTMLElement) {
    const ruleNameLabels = Array.from(rulesForm.getElementsByTagName('label') ?? [])
      .filter((l) => l.textContent?.endsWith('.json') && l.parentElement != null);
    for (const label of ruleNameLabels) {
      const ruleName = label.textContent!.trim();
      const inputRoot = label.parentElement!.getElementsByTagName('input')[0]!;
      ui.tooltip.bind(inputRoot, ui.wait(async () => {
        RuleInputs.ruleDescriptions[ruleName] ??= (async () => {
          try {
            const rule = JSON.parse(await grok.dapi.files.readAsText(`${RULES_PATH}/${ruleName}`));
            const descDiv = ui.divV([ui.h1(ruleName)]);
            const linkRules: RuleLink[] = rule.linkRules ?? [];
            // link rules
            if (linkRules.length > 0) {
              descDiv.appendChild(ui.h3('Linkage Rules'));
              const table = ui.table(linkRules, (item) => [
                item.code?.toString() ?? '',
                shortenToMaxLen((item.firstMonomers ?? []).join(', ')),
                shortenToMaxLen((item.secondMonomers ?? []).join(', ')),
                item.firstLinkingGroup?.toString() ?? '',
                item.secondLinkingGroup?.toString() ?? ''
              ], ['Code', 'M1', 'M2', 'L1', 'L2']);
              descDiv.appendChild(table);
            }
            // reaction rules
            const reactionRules: RuleReaction[] = rule.reactionRules ?? [];
            if (reactionRules.length > 0) {
              descDiv.appendChild(ui.h3('Synthesis Rules'));
              const table = ui.table(reactionRules, (item) => [
                item.code?.toString() ?? '',
                shortenToMaxLen((item.firstMonomers ?? []).join(', ')),
                shortenToMaxLen((item.secondMonomers ?? []).join(', ')),
                shortenToMaxLen(item.name ?? ''),
                grok.chem.drawMolecule(item.reaction?.split('>>')[1] ?? '', 70, 70)
              ], ['Code', 'M1', 'M2', 'Name', 'Product']);
              descDiv.appendChild(table);
            }
            return descDiv;
          } catch (e: any) {
            _package.logger.error(`Failed to load rule ${ruleName}: ${e?.toString?.()}`);
            console.error(e);
            return ui.divText(`Failed to load rule ${ruleName}`);
          }
        })();
        return await RuleInputs.ruleDescriptions[ruleName];
      }));
    }
  }
}


function shortenToMaxLen(s: string, maxLen: number = 30): string {
  if (s.length <= maxLen) return s;
  return `${s.substring(0, maxLen)}...`;
}


export type RuleLink = {
  code: number,
  firstMonomers: string[],
  secondMonomers: string[],
  firstLinkingGroup: number,
  secondLinkingGroup: number
}

export type RuleReaction = {
  code: number,
  firstMonomers: string[],
  secondMonomers: string[],
  reaction: string,
  name: string
}

export type LinkRuleRowArgs = {
  row?: number,
  code: number,
  firstMonomers: string,
  secondMonomers: string,
  firstLinkingGroup: number,
  secondLinkingGroup: number
}

export class Rules {
  homodimerCode: string | null;
  heterodimerCode: string | null;
  linkRules: RuleLink[];
  reactionRules: RuleReaction[];

  constructor(homodimerCode: string | null, heterodimerCode: string | null,
    linkRules: RuleLink[], reactionRules: RuleReaction[]) {
    this.homodimerCode = homodimerCode;
    this.heterodimerCode = heterodimerCode;
    this.linkRules = linkRules;
    this.reactionRules = reactionRules;
  }

  set homodimer(code: string) {
    this.homodimerCode = code;
  }

  set heterodimer(code: string) {
    this.heterodimerCode = code;
  }

  addLinkRules(rules: RuleLink[]): void {
    for (let j = 0; j < rules.length; j++)
      this.linkRules.push(rules[j]);
  }

  addSynthesisRules(rules: RuleReaction[]): void {
    for (let j = 0; j < rules.length; j++)
      this.reactionRules.push(rules[j]);
  }

  getLinkRulesDf() {
    const length = this.linkRules.length;
    const codeCol = DG.Column.int(NAME_CODE, length);
    codeCol.setTag('friendlyName', 'Code');
    //ui.tooltip.bind(codeCol.root, 'Click to zoom');
    const firstMonomerCol = DG.Column.string(NAME_FIRST_MONOMERS, length);
    firstMonomerCol.temp[MmcrTemps.fontSize] = 16;
    firstMonomerCol.setTag('friendlyName', 'First monomers');
    firstMonomerCol.semType = DG.SEMTYPE.MACROMOLECULE;
    PackageFunctions.applyNotationProviderForCyclized(firstMonomerCol, ',');

    const secondMonomerCol = DG.Column.string(NAME_SECOND_MONOMERS, length);
    secondMonomerCol.setTag('friendlyName', 'Second monomers');
    secondMonomerCol.temp[MmcrTemps.fontSize] = 16;
    secondMonomerCol.semType = DG.SEMTYPE.MACROMOLECULE;
    PackageFunctions.applyNotationProviderForCyclized(secondMonomerCol, ',');

    const firstLinkingGroupCol = DG.Column.int(NAME_FIRST_LINK, length);
    firstLinkingGroupCol.setTag('friendlyName', 'First group');
    const secondLinkingGroupCol = DG.Column.int(NAME_SECOND_LINK, length);
    secondLinkingGroupCol.setTag('friendlyName', 'Second group');

    for (let i = 0; i < length; i++) {
      codeCol.set(i, this.linkRules[i].code);
      firstMonomerCol.set(i, this.linkRules[i].firstMonomers.toString());
      secondMonomerCol.set(i, this.linkRules[i].secondMonomers.toString());
      firstLinkingGroupCol.set(i, this.linkRules[i].firstLinkingGroup);
      secondLinkingGroupCol.set(i, this.linkRules[i].secondLinkingGroup);
    }

    const res = DG.DataFrame.fromColumns([
      codeCol, firstMonomerCol, secondMonomerCol, firstLinkingGroupCol, secondLinkingGroupCol
    ]);

    const addNewRow = (rowObj: LinkRuleRowArgs) => {
      if (rowObj.row != undefined && rowObj.row < res.rowCount && rowObj.row >= 0) {
        codeCol.set(rowObj.row, rowObj.code);
        firstMonomerCol.set(rowObj.row, rowObj.firstMonomers);
        secondMonomerCol.set(rowObj.row, rowObj.secondMonomers);
        firstLinkingGroupCol.set(rowObj.row, rowObj.firstLinkingGroup);
        secondLinkingGroupCol.set(rowObj.row, rowObj.secondLinkingGroup);
      } else {
        const {code, firstMonomers, secondMonomers, firstLinkingGroup, secondLinkingGroup} = rowObj;
        res.rows.addNew([code, firstMonomers ?? '', secondMonomers ?? '', firstLinkingGroup ?? 1, secondLinkingGroup ?? 2]);
      }
    };

    return {res, addNewRow};
  }

  async getLinkCards(): Promise<RuleCards[]> {
    const length = this.linkRules.length;
    const cards: RuleCards[] = new Array<RuleCards>(length);
    const monomerLibHelper = await getMonomerLibHelper();
    const systemMonomerLib = monomerLibHelper.getMonomerLib();

    for (let i = 0; i < length; i++) {
      cards[i] = new RuleCards(this.linkRules[i].firstMonomers, this.linkRules[i].secondMonomers,
        systemMonomerLib, this.linkRules[i].code, this);
    }

    return cards;
  }

  getSynthesisRulesDf() {
    const length = this.reactionRules.length;
    const codeCol = DG.Column.int(NAME_CODE, length);
    codeCol.setTag('friendlyName', 'Code');
    const firstMonomerCol = DG.Column.string(NAME_FIRST_MONOMERS, length);
    firstMonomerCol.setTag('friendlyName', 'First monomers');
    const secondMonomerCol = DG.Column.string(NAME_SECOND_MONOMERS, length);
    secondMonomerCol.setTag('friendlyName', 'Second monomers');
    // set these columns as macromolecule for better visualization
    firstMonomerCol.semType = DG.SEMTYPE.MACROMOLECULE;
    firstMonomerCol.temp[MmcrTemps.fontSize] = 16;
    PackageFunctions.applyNotationProviderForCyclized(firstMonomerCol, ',');
    secondMonomerCol.semType = DG.SEMTYPE.MACROMOLECULE;
    secondMonomerCol.temp[MmcrTemps.fontSize] = 16;
    PackageFunctions.applyNotationProviderForCyclized(secondMonomerCol, ',');
    secondMonomerCol.setTag(DG.TAGS.CELL_RENDERER, 'Sequence');
    firstMonomerCol.setTag(DG.TAGS.CELL_RENDERER, 'Sequence');
    const name = DG.Column.string(NAME_REACTION_NAME, length);
    name.setTag('friendlyName', 'Name');
    const firstReactant = DG.Column.string('firstReactant', length);
    firstReactant.setTag('friendlyName', 'First reactant');
    const secondReactant = DG.Column.string('secondReactant', length);
    secondReactant.setTag('friendlyName', 'Second reactant');
    const product = DG.Column.string('product', length);
    product.setTag('friendlyName', 'Product');

    for (let i = 0; i < length; i++) {
      codeCol.set(i, this.reactionRules[i].code);
      firstMonomerCol.set(i, this.reactionRules[i].firstMonomers.toString());
      secondMonomerCol.set(i, this.reactionRules[i].secondMonomers.toString());
      name.set(i, this.reactionRules[i].name);

      const reaction = this.reactionRules[i].reaction.split('>>');
      const reactants = reaction[0].split('.');

      firstReactant.set(i, reactants[0]);
      secondReactant.set(i, reactants[1]);
      product.set(i, reaction[1]);
    }
    firstReactant.semType = DG.SEMTYPE.MOLECULE;
    secondReactant.semType = DG.SEMTYPE.MOLECULE;
    product.semType = DG.SEMTYPE.MOLECULE;


    const df = DG.DataFrame.fromColumns([
      name, firstReactant, secondReactant, product, codeCol, firstMonomerCol, secondMonomerCol
    ]);

    const addNewRow = (rowObj: ReactionEditorProps) => {
      if (rowObj.rowIndex != undefined && rowObj.rowIndex < df.rowCount && rowObj.rowIndex >= 0) {
        name.set(rowObj.rowIndex, rowObj.resultMonomerName);
        codeCol.set(rowObj.rowIndex, rowObj.code);
        firstMonomerCol.set(rowObj.rowIndex, rowObj.firstMonomers.join(','));
        secondMonomerCol.set(rowObj.rowIndex, rowObj.secondMonomers.join(','));
        firstReactant.set(rowObj.rowIndex, rowObj.firstReactantSmiles);
        secondReactant.set(rowObj.rowIndex, rowObj.secondReactantSmiles);
        product.set(rowObj.rowIndex, rowObj.productSmiles);
      } else {
        const {resultMonomerName, code, firstMonomers, secondMonomers,
          firstReactantSmiles, secondReactantSmiles, productSmiles} = rowObj;
        df.rows.addNew([resultMonomerName, firstReactantSmiles, secondReactantSmiles,
          productSmiles, code, firstMonomers.join(','), secondMonomers.join(',')]);
      }
    };
    return {df, addNewRow};
  }

  setLinkRules(df: DG.DataFrame) : void {
    const length = df.rowCount;
    const rules: RuleLink [] = new Array<RuleLink>(length);
    const codeCol = df.columns.byName(NAME_CODE);
    const firstMonomerCol = df.columns.byName(NAME_FIRST_MONOMERS);
    const secondMonomerCol = df.columns.byName(NAME_SECOND_MONOMERS);
    const firstLink = df.columns.byName(NAME_FIRST_LINK);
    const secondLink = df.columns.byName(NAME_SECOND_LINK);


    for (let i = 0; i < length; i++) {
      const fSplit = firstMonomerCol.get(i).split(',');
      const sSplit = secondMonomerCol.get(i).split(',');

      const rule = {
        code: codeCol.get(i),
        firstMonomers: fSplit[0] !== '' ? fSplit : [],
        secondMonomers: sSplit[0] !== '' ? sSplit : [],
        firstLinkingGroup: firstLink.get(i),
        secondLinkingGroup: secondLink.get(i)
      };

      rules[i] = rule;
    }

    this.linkRules = rules;
  }

  setSynthesisRules(df: DG.DataFrame) : void {
    const length = df.rowCount;
    const rules: RuleReaction [] = new Array<RuleReaction>(length);
    const codeCol = df.columns.byName(NAME_CODE);
    const firstMonomerCol = df.columns.byName(NAME_FIRST_MONOMERS);
    const secondMonomerCol = df.columns.byName(NAME_SECOND_MONOMERS);
    const name = df.columns.byName(NAME_REACTION_NAME);
    const firstReactant = df.columns.byName('firstReactant');
    const secondReactant = df.columns.byName('secondReactant');
    const product = df.columns.byName('product');

    for (let i = 0; i < length; i++) {
      const smartsReaction = `${firstReactant.get(i)}.${secondReactant.get(i)}>>${product.get(i)}`;
      const fSplit = firstMonomerCol.get(i).split(',');
      const sSplit = secondMonomerCol.get(i).split(',');

      const rule = {
        code: codeCol.get(i),
        firstMonomers: fSplit[0] !== '' ? fSplit : [],
        secondMonomers: sSplit[0] !== '' ? sSplit : [],
        reaction: smartsReaction,
        name: name.get(i)
      };

      rules[i] = rule;
    }

    this.reactionRules = rules;
  }
}

export async function getRules(ruleFiles: string[]): Promise<Rules> {
  const fileSource = new DG.FileSource(RULES_PATH);
  const rules: Rules = new Rules(null, null, [], []);

  for (let i = 0; i < ruleFiles.length; i++) {
    const rulesRaw = await fileSource.readAsText(ruleFiles[i].replace(RULES_PATH, ''));
    const ruleSingle : Rules = JSON.parse(rulesRaw);

    rules.homodimer = ruleSingle.homodimerCode!;
    rules.heterodimer = ruleSingle.heterodimerCode!;
    rules.addLinkRules(ruleSingle.linkRules);
    rules.addSynthesisRules(ruleSingle.reactionRules);
  }

  return rules;
}

export function getMonomerPairs(rule: RuleLink | RuleReaction) : [string[], string[]] {
  const allPairsNum = rule.firstMonomers.length*rule.secondMonomers.length;
  const firstMonomers = new Array<string>(allPairsNum);
  const secondMonomers = new Array<string>(allPairsNum);

  let counter = 0;
  for (let i = 0; i < rule.firstMonomers.length; i++) {
    for (let j = 0; j < rule.secondMonomers.length; j++) {
      firstMonomers[counter] = rule.firstMonomers[i];
      secondMonomers[counter] = rule.secondMonomers[j];
      counter++;
    }
  }

  return [firstMonomers, secondMonomers];
}
