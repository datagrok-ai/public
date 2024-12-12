
import { delay, DockManager, GridCell } from "datagrok-api/dg";
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import { loadIssueData } from "./api/data";
import { _package } from "./package-test";
import { AuthCreds } from "./api/objects";
import { getJiraCreds } from "./app/credentials";
import { issueData } from "./package";

export class JiraGridCellHandler extends DG.ObjectHandler {

    get type(): string {
        return 'JIRA Ticket';
    }


    isApplicable(x: any): boolean {
        return x instanceof DG.SemanticValue && this.type === x.semType;
    }

    renderProperties(semValue: DG.SemanticValue, context: any = null): HTMLElement {
        let resultElement = ui.div();
        issueData(semValue.cell.value).then((issue) => {
            console.log(issue);
            resultElement.innerText = 'hello';
        });
        return resultElement;
    }

    renderTooltip(semValue: DG.SemanticValue, context?: any): HTMLElement {
        let resultElement = ui.div();
        issueData(semValue.cell.value).then((issue) => {
            console.log(issue);
            resultElement.innerText = 'hello';
        });
        return resultElement;
    }
}