/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import './style.css';

import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types';
import {Rules} from './pt-rules';
import {getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {doPolyToolConvert} from './pt-conversion';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';
import {getSeqHelper, ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getOverriddenLibrary} from './pt-synthetic';
import {helmToMol} from './pt-atomic';
import {execMonomerHoverLinks} from '@datagrok-libraries/bio/src/monomer-works/monomer-hover';

class MonomerCard {
  root: HTMLElement = ui.divV([], {classes: 'monomer-card-rule-root'});

  private _selected: boolean = false;
  get selected(): boolean { return this._selected; }
  set selected(value: boolean) {
    this._selected = value;
    this.root.style.border = value ? '2px solid var(--green-2)' : '2px solid var(--grey-2)';
  }

  constructor(public monomer: Monomer | null, public monomerSymbol: string) {}

  render() {
    ui.empty(this.root);
    this.root.appendChild(ui.h2(this.monomerSymbol, {style: {textAlign: 'center'}}));
    if (this.monomer) {
      const monomerMolSvg = this.monomer.smiles && grok.chem.checkSmiles(this.monomer.smiles) ?
        grok.chem.drawMolecule(this.monomer.smiles, 150, 120) :
        grok.chem.drawMolecule(this.monomer.molfile ?? '', 150, 120);
      this.root.appendChild(monomerMolSvg);
      const monomerName =
        ui.divH([ui.divText('Monomer Name: '), ui.divText(this.monomer.name)], {classes: 'monomer-card-info-rules'});

      this.root.appendChild(monomerName);
      ui.tooltip.bind(monomerName, this.monomer.name);
      if (this.monomer.lib?.source) {
        const monomerSource =
          ui.divH([ui.divText('Source: '), ui.divText(this.monomer.lib.source)], {classes: 'monomer-card-info-rules'});
        this.root.appendChild(monomerSource);
        ui.tooltip.bind(monomerSource, this.monomer.lib.source);
      }
      const monomerType = ui.divH([
        ui.divText('Polymer Type: '),
        ui.divText(this.monomer.polymerType)
      ], {classes: 'monomer-card-info-rules'});
      this.root.appendChild(monomerType);
      ui.tooltip.bind(monomerType, this.monomer.polymerType);

      ui.tooltip.bind(this.root, 'Select Monomer');
    } else {
      this.root.appendChild(ui.divV([`Monomer ${this.monomerSymbol} not found in libraries`],
        {style: {textAlign: 'center', height: '100%'}}));
    }
  }
}

export class RuleCards {
  root: HTMLElement = ui.divH([],
    {style: {
      width: '100%',
      overflow: 'hidden',
      visibility: 'visible',
    }, classes: 'monomer-cards'}
  );
  cardsFirst: MonomerCard[];
  cardsSecond: MonomerCard[];
  firstCard: MonomerCard;
  secondCard: MonomerCard;
  resulting: HTMLElement;
  actionable: boolean;

  constructor(public firstMonomers: string[], public secondMonomers: string[],
    lib: IMonomerLib, public code: number, public rules: Rules) {
    const emptyFirst = !firstMonomers?.length;
    const emptySecond = !secondMonomers?.length;
    // if one of the lists is empty, it means that the rule is active for all monomers
    if (emptyFirst)
      this.firstMonomers = firstMonomers = ['C', 'D', 'E', 'F', 'G'];
    if (emptySecond)
      this.secondMonomers = secondMonomers = ['D', 'C', 'E', 'F', 'G'];
    this.actionable = true;
    const monomerGalleryFirst = ui.divH([], {style: {overflowX: 'auto', width: '100%'}});
    const monomerGallerySecond = ui.divH([], {style: {overflowX: 'auto', width: '100%'}});

    this.resulting = ui.divH([], {style: {overflowX: 'auto', width: '100%', minHeight: '150px'}});
    const headerFirst = ui.h2('First Monomer' + (!emptyFirst ? '' : ' (Applied to all monomers, example list shown)'),
      {style: {width: '100%'}});
    const headerSecond = ui.h2('Second Monomer' + (!emptySecond ? '' : ' (Applied to all monomers, example list shown)'),
      {style: {width: '100%'}});
    const galleries = ui.divV([headerFirst, monomerGalleryFirst, headerSecond, monomerGallerySecond, this.resulting]);
    this.cardsFirst = firstMonomers.map((symbol) => new MonomerCard(lib.getMonomer('PEPTIDE', symbol), symbol));
    this.cardsSecond = secondMonomers.map((symbol) => new MonomerCard(lib.getMonomer('PEPTIDE', symbol), symbol));

    this.cardsFirst.forEach((card) => {
      card.root.onclick = () => {
        this.cardsFirst.forEach((c) => c.selected = false);
        card.selected = true;
        this.firstCard = card;
        this.reset();
        //onMonomerSelected(card.monomer);
      };
      monomerGalleryFirst.appendChild(card.root);
    });

    this.cardsSecond.forEach((card) => {
      card.root.onclick = () => {
        this.cardsSecond.forEach((c) => c.selected = false);
        card.selected = true;
        this.secondCard = card;
        this.reset();
        //onMonomerSelected(card.monomer);
      };
      monomerGallerySecond.appendChild(card.root);
    });

    this.cardsFirst[0].selected = true;
    this.cardsSecond[0].selected = true;
    this.firstCard = this.cardsFirst[0];
    this.secondCard = this.cardsSecond[0];
    this.root.appendChild(galleries);
  }

  async reset() {
    if (this.actionable) {
      const seqs: string[] = [
        `${this.firstCard.monomer?.symbol}(${this.code})-A-A-A-A-${this.secondCard.monomer?.symbol}(${this.code})`
      ];
      'PEPTIDE1{[NH2].A.A.A.A.D}$PEPTIDE1,PEPTIDE1,1:R2-6:R3$$$V2.0';

      const helmHelper = await getHelmHelper();
      const [helms, isLinear, positionMaps] = doPolyToolConvert(seqs, this.rules, helmHelper);

      const resHelmCol = DG.Column.fromType(DG.COLUMN_TYPE.STRING, 'helm', helms.length)
        .init((rowIdx: number) => { return helms[rowIdx]; });
      resHelmCol.semType = DG.SEMTYPE.MACROMOLECULE;
      resHelmCol.meta.units = NOTATION.HELM;
      resHelmCol.setTag(DG.TAGS.CELL_RENDERER, 'helm');

      const rdKitModule: RDModule = await getRdKitModule();
      const seqHelper: ISeqHelper = await getSeqHelper();

      const lib = await getOverriddenLibrary(this.rules);
      const resHelmColTemp = resHelmCol.temp;
      resHelmCol.temp = resHelmColTemp;
      const resMolCol = await helmToMol(resHelmCol, helms,
        isLinear, true, false, false, lib, rdKitModule, seqHelper);

      const mol = resMolCol.get(0);
      const monomerMolSvg = mol && grok.chem.checkSmiles(mol) ?
        grok.chem.drawMolecule(mol, 200, 200) :
        grok.chem.drawMolecule(mol ?? '', 200, 200);

      if (mol) {
        monomerMolSvg.style.cursor = 'pointer';
        monomerMolSvg.onclick = () => {
          const width = window.innerWidth - 200;
          const height = window.innerHeight - 200;
          const bigMol = grok.chem.checkSmiles(mol) ?
            grok.chem.drawMolecule(mol, width, height) :
            grok.chem.drawMolecule(mol ?? '', width, height);
          ui.dialog({title: 'Molecule'}).add(bigMol).showModal(true);
        };
        ui.tooltip.bind(monomerMolSvg, 'Click to expand');
      }

      const exampleDiv = ui.divV([
        ui.h2('Example Result:'), ui.h2(seqs[0])
      ], {style: {width: '200px'}});

      ui.empty(this.resulting);
      this.resulting.append(ui.divH([exampleDiv, monomerMolSvg], {style: {overflowX: 'auto', width: '100%', minHeight: '150px'}}));
    }
  }

  render() {
    if (this.actionable) {
      this.cardsFirst.forEach((card) => card.render());
      this.cardsSecond.forEach((card) => card.render());
    }
  }
}
