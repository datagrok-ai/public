import { category } from "@datagrok-libraries/utils/src/test";
import { delay, DockManager, GridCell } from "datagrok-api/dg";
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

export class StackTraceHandler extends DG.ObjectHandler {
    get type(): string {
        return 'stackTrace';
    }

    isApplicable(x: any): boolean {
        return x instanceof DG.SemanticValue && x.semType === 'stackTrace';
    }

    renderTooltip(semValue: DG.SemanticValue<string>, context: any = null): HTMLElement {
        let lines: string[] = semValue.value.split('\\n');
        lines = lines.slice(0, Math.min(lines.length, 2));
        return ui.span([lines.join('\\n')]);
    }

    renderProperties(semValue: DG.SemanticValue, context: any = null): HTMLElement {
        let d = ui.div([], {classes: 'div-textarea,ui-report-textarea,ui-log-textarea'});
        let lines: string[] = semValue.value.split('\\n');
        for (var i = 0; i < lines.length; i++)
            d.appendChild(ui.span([lines[i], "br"]));
        return d;
    }
}