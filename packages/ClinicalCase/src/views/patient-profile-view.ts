import { ClinicalCaseView } from "../clinical-case-view";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import {InputBase} from "datagrok-api/dg";

export class PatientProfileView extends DG.ViewBase {

  id: HTMLHeadingElement = ui.h2('');
  sex: InputBase = ui.stringInput('Sex', '');
  age: InputBase = ui.stringInput('Age', '');
  race: InputBase = ui.stringInput('Race', '');
  arm: InputBase = ui.stringInput('ARM', '');
  // height: InputBase = ui.stringInput('Height', '');
  // weight: InputBase = ui.stringInput('Weight', '');

  constructor() {
    super();

    this.root.appendChild(ui.divH([
      ui.button('<', () => this.refresh(study.domains.dm.currentRowIdx - 1)),
      ui.button('>', () => this.refresh(study.domains.dm.currentRowIdx + 1))
    ]))

    this.root.appendChild(this.id);
    this.root.appendChild(ui.inputs([this.sex, this.age, this.race, this.arm]));
    this.refresh(0);
  }

  refresh(idx: number) {
    study.domains.dm.currentRowIdx = Math.max(0, Math.min(idx, study.domains.dm.rowCount - 1));
    let row: DG.Row = study.domains.dm.currentRow;
    this.id.textContent = row['usubjid'];
    this.sex.value = row['sex'];
    this.race.value = row['race'];
    this.age.value = row['age'].toString();
    this.arm.value = row['arm'];
  }
}