import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { ValidationHelper } from "../helpers/validation-helper";
import { createValidationErrorsDiv, getRequiredColumnsByView, updateDivInnerHTML } from "../views/utils";

export class ClinicalCaseViewBase extends DG.ViewBase {

    constructor(name) {
        super({});
        this.name = name;
      }

      loaded = false;
      domainsAndColumnsForValidation: {};
      validator: ValidationHelper;
      optDomainsWithMissingCols: string[];

      load(): void {
        this.validator = new ValidationHelper(getRequiredColumnsByView()[this.name]);
        if(this.validator.validate()) {
          this.optDomainsWithMissingCols = Object.keys(this.validator.missingColumnsInOptDomains)
            .filter(it => this.validator.missingColumnsInOptDomains[it].length);
          this.createView();
          this.loaded = true;
        } else {
          updateDivInnerHTML(this.root, createValidationErrorsDiv(this.validator.missingDomains, this.validator.missingColumnsInReqDomains, this.validator.missingColumnsInOptDomains));
        }
      }

      createView() {}

      async propertyPanel() {
        return await this.createAccWithTitle(this.name).root;
      }

      createAccWithTitle(panelName: string, title = null){
        let acc = ui.accordion(`${panelName} panel`);
        let accIcon = ui.element('i');
        accIcon.className = 'grok-icon svg-icon svg-view-layout';
        acc.addTitle(ui.span([accIcon, ui.label(`${title ?? panelName}`)]));
        return acc;
      }

}