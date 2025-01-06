
import { delay, DockManager, GridCell, GridCellRenderer, GridCellStyle, LruCache, x } from "datagrok-api/dg";
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import { loadIssueData } from "./api/data";
import { _package } from "./package-test";
import { AuthCreds, JiraIssue } from "./api/objects";
import { getJiraCreds } from "./app/credentials";
import { cache, issueData, queried } from "./package";
import imgBug from './images/bug.png';
import imgEpic from './images/epic.png';
import imgFeature from './images/feature.png';
import imgImprovement from './images/improvement.png';
import imgTask from './images/task.png';

function image(url: string) {
    const img = new Image();
    img.src = url;
    return img;
}

const images = {
    Bug: image(imgBug),
    'New Feature': image(imgFeature),
    Improvement: image(imgImprovement),
    Epic: image(imgEpic),
    Task: image(imgTask)
}

export class JiraGridCellHandler extends DG.ObjectHandler {

    get type(): string {
        return 'JIRA Ticket';
    }

    isApplicable(x: any): boolean {
        return x instanceof DG.SemanticValue && this.type === x.semType;
    }

    getGridCellRenderer(): GridCellRenderer | null {
        return new JiraTicketGridCellRenderer();
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
            const creds = await getJiraCreds();
            if (!issue)
                return ui.divText('Issue not found');

            const object: Record<string, WithCalcButtonType> = {};
            const column = x && x.cell && x.cell.dart && x.cell.column ? x.cell.column : null; // need to check dart, sometimes its null
            // TODO: hardcoded fields
            object.Key = withCalcForWholeColumn('', ui.link(issue.key, `https://${creds.host}/browse/${issue.key}`, 'Open in Jira'), 'Key', null);
            object.Summary = withCalcForWholeColumn('summary', textRenderer(issue.fields?.summary ?? ''), 'Summary', column);
            object.Created = withCalcForWholeColumn('created', textRenderer(issue.fields?.created ?? ''), 'Created', column);
            object.Updated = withCalcForWholeColumn('updated', textRenderer(issue.fields?.updated ?? ''), 'Updated', column);
            object.Labels = withCalcForWholeColumn('labels', listRenderer(issue.fields?.labels ?? ['']), 'Labels', column);

            object.Status = withCalcForWholeColumn('status:name', textRenderer(issue.fields?.status?.name ?? ''), 'Status', column);
            object.Assignee = withCalcForWholeColumn('assignee:displayName', textRenderer(issue.fields?.assignee?.displayName ?? ''), 'Assignee', column);
            object.Creator = withCalcForWholeColumn('creator:displayName', textRenderer(issue.fields?.creator?.displayName ?? ''), 'Creator', column);
            const card = ui.card(ui.table(Object.entries(object), (item) => {
                return [item[0], item[1].element];
            }));
            ui.tools.setHoverVisibility(card, Object.values(object).filter((e: any) => e.calcButton).map((e: any) => e.calcButton));
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

function getIssueFields(issue: JiraIssue): Record<string, string> {
    return {
        'Summary': issue.fields.summary ?? '',
        'Type': issue.fields.issuetype?.name ?? '',
        'Priority': issue.fields.priority?.name ?? '',
        'Created': issue.fields.created ?? '',
        'Updated': issue.fields.updated ?? '',
        'Status': issue.fields.status?.name ?? '',
        'Assignee': issue.fields.assignee?.displayName ?? '',
        'Creator': issue.fields.creator?.displayName ?? '',
        'Labels': issue.fields.labels?.join(', ') ?? '',
    }
}

type WithCalcButtonType = {
    element: HTMLElement,
    calcButton: HTMLElement | null
}

function withCalcForWholeColumn(keys: string, element: HTMLElement, fieldName: string, column?: DG.Column | null): WithCalcButtonType {
    if (!column || !column.dataFrame) {
        return { element, calcButton: null };
    }
    const func = DG.Func.find({ name: 'getJiraField' })[0];
    if (!func) {
        return { element, calcButton: null };
    }
    const calcButton = ui.icons.add(async () => {
        const pg = DG.TaskBarProgressIndicator.create('Retrieving Data...');
        try {
            const resCol: DG.Column | null = await func.apply({ ticketColumn: column, field: keys })
            if (resCol && resCol.length === column.length) {
                const newName = column.dataFrame.columns.getUnusedName(fieldName);

                column.dataFrame.columns.addNew(newName, resCol.type).init((i) => resCol.get(i));
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
    return { element: container, calcButton };
}


function textRenderer(value: string | number, withTooltip = true) {
    const nameHost = ui.div(value.toString());
    nameHost.style.maxWidth = '200px';
    nameHost.style.overflow = 'hidden';
    nameHost.style.whiteSpace = 'nowrap';
    nameHost.style.textOverflow = 'ellipsis';

    const menu = DG.Menu.popup();
    menu.item('Copy', () => {
        navigator?.clipboard?.writeText(value.toString());
    });
    nameHost.addEventListener('contextmenu', (e) => {
        e.preventDefault();
        e.stopImmediatePropagation();
        setTimeout(() => menu.show());
    });

    if (value.toString().length > 20 && withTooltip)
        ui.tooltip.bind(nameHost, value.toString());
    return nameHost;
}

function listRenderer(list: string[], withTooltip = true) {
    const nameHost = ui.list(list);

    nameHost.style.maxWidth = '200px';
    nameHost.style.overflow = 'hidden';
    nameHost.style.whiteSpace = 'nowrap';
    nameHost.style.textOverflow = 'ellipsis';

    const menu = DG.Menu.popup();
    menu.item('Copy', () => {
        navigator?.clipboard?.writeText(list.join(' ').toString());
    });
    nameHost.addEventListener('contextmenu', (e) => {
        e.preventDefault();
        e.stopImmediatePropagation();
        setTimeout(() => menu.show());
    });

    if (list.join(' ').toString().length > 20 && withTooltip)
        ui.tooltip.bind(nameHost, list.join(' ').toString());
    return nameHost;
}


class JiraTicketGridCellRenderer extends DG.GridCellRenderer {

    get name(): string { return 'JIRA Ticket'; }
    get cellType(): string { return 'JIRA Ticket'; }
    get defaultWidth(): number | null { return 100; }

    render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
        const key = gridCell.cell.valueString;
        const stRen = DG.GridCellRenderer.byName('string')!;
        const keyHeight = Math.min(25, h);
        const ticket = cache.get(key);
        const props = ticket ? getIssueFields(ticket) : null;

        const renderKey = () => {
            stRen.render(g, x + 17, y, w, keyHeight, gridCell, cellStyle);
            if (props && (images as any)[props.Type]) {
                g.drawImage((images as any)[props.Type], x + 4, y + 4, 16, 16);
            }
        };

        if (!cache.has(key)) {
            issueData(key).then((_) => gridCell.render());
            renderKey();
        }
        else {
            cellStyle.textColor = ticket ? cellStyle.textColor : DG.Color.fromHtml('red');
            renderKey();

            if (ticket && props) {
                g.fillStyle = '#4d5261';
                g.textAlign =  "left";

                if (w > 150) {
                    g.fillText(props.Summary, x + 20, y + 30);
                    delete props.Summary;
                }

                if (h > 40) {
                    const fields = Object.getOwnPropertyNames(props);
                    for (let i = 0; i < fields.length; i++) {
                        g.fillText(fields[i], x + 20, y + 50 + 20 * i, 100);
                        g.fillText('' + props[fields[i]], x + 80, y + 50 + 20 * i);
                    }
                }
            }
        }
    }

    // onMouseMove(gridCell: GridCell, e: MouseEvent) {
    //     DG.GridCellRenderer.byName('string')?.onMouseMove(gridCell, e);
    // }
}