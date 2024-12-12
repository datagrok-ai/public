
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
        const acc = ui.panels.infoPanel(semValue);// this retrieves default properties
        const firstPane = acc.panes[0];
        acc.addPane('Jira Ticket', () => this.renderCard(semValue, context), true, firstPane, false);
        return acc.root;
    }

    renderCard(x: DG.SemanticValue, context?: any): HTMLElement {
        return ui.wait(async () => {
            const issue = await issueData(x.value);
            if (!issue)
                return ui.divText('Issue not found');
            
            const object: Record<string, WithCalcButtonType> = {};
            const column = x && x.cell && x.cell.dart && x.cell.column ? x.cell.column : null; // need to check dart, sometimes its null
            // TODO: hardcoded fields
            object.Key = withCalcForWholeColumn('', ui.link(issue.key, `https://reddata.atlassian.net/browse/${issue.key}`, 'Open in Jira'), 'Key', null);
            object.Summary = withCalcForWholeColumn('summary', textRenderer(issue.fields?.summary ?? ''), 'Summary', column);
            object.Updated = withCalcForWholeColumn('updated',textRenderer(issue.fields?.updated ?? ''), 'Updated', column);
            object.Status = withCalcForWholeColumn('status:name', textRenderer(issue.fields?.status?.name ?? ''), 'Status', column);
            object.Assignee = withCalcForWholeColumn('assignee:displayName',textRenderer(issue.fields?.assignee?.displayName ?? ''), 'Assignee', column);
            object.Creator = withCalcForWholeColumn('creator:displayName', textRenderer(issue.fields?.creator?.displayName ?? ''), 'Creator', column);
            const card = ui.card(ui.table(Object.entries(object), (item) => {
                return [item[0], item[1].element];
            }));
            ui.tools.setHoverVisibility(card , Object.values(object).filter((e: any) => e.calcButton).map((e: any) => e.calcButton));
            card.style.width = 'unset';
            card.style.margin = '0';
            card.addEventListener('click', () => grok.shell.o = x);
            return card;
        });
        
    }

    renderTooltip(semValue: DG.SemanticValue, context?: any): HTMLElement {
        return this.renderCard(semValue, context);
    }
}

type WithCalcButtonType = {
    element: HTMLElement,
    calcButton: HTMLElement | null
}

function withCalcForWholeColumn(keys: string, element: HTMLElement, fieldName: string, column?: DG.Column | null): WithCalcButtonType {
    if (!column || !column.dataFrame) {
        return {element, calcButton: null};
    }
    const func = DG.Func.find({name: 'getJiraField'})[0];
    if (!func) {
        return {element, calcButton: null};
    }
    const calcButton = ui.icons.add(async () => {
        const pg = DG.TaskBarProgressIndicator.create('Retrieving Data...');
        try {
            const resCol: DG.Column | null = await func.apply({ticketColumn: column, field: keys})
            if (resCol && resCol.length === column.length) {
                const newName = column.dataFrame.columns.getUnusedName(fieldName);
    
                column.dataFrame.columns.addNewString(newName).init((i) => resCol.get(i)?.toString());
            } else {
                grok.shell.warning(`Failed to retrieve '${fieldName}' for whole table`);
            }
        } catch (e) {
            grok.shell.error('Error while retrieving data');
            console.error(e);
        } finally {
            pg.close();
        }
        
    }, `Retrieve '${fieldName}' for whole table`);
    calcButton.style.marginLeft = '8px';

    const container = ui.divH([element, calcButton]);
    container.style.alignItems = 'center';
    return {element: container, calcButton};
}


function textRenderer(value: string | number, withTooltip = true) {
    const nameHost = ui.div(value.toString());
    nameHost.style.maxWidth = '200px';
    nameHost.style.overflow = 'hidden';
    nameHost.style.whiteSpace = 'nowrap';
    nameHost.style.textOverflow = 'ellipsis';
  
    // ability to copy
    const menu = DG.Menu.popup();
    menu.item('Copy', () => {
      navigator?.clipboard?.writeText(value.toString());
    });
    nameHost.addEventListener('contextmenu', (e) => {
      e.preventDefault();
      e.stopImmediatePropagation();
      setTimeout(() => menu.show());
    });
  
    // approximately what will fit in 150 px
    if (value.toString().length > 20 && withTooltip) ui.tooltip.bind(nameHost, value.toString());
    return nameHost;
  }