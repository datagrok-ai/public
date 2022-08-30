import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

import {IC50SemType} from "../../IC50SemType";
import {IC50Entity} from "../../IC50Entity";
import {DGIC50GridCellRenderer} from "./DGIC50GridCellRenderer";

export class DGIC50Handler extends DG.ObjectHandler
{
    get type()
    {
        return IC50SemType.Instance;
    }

    isApplicable(x)
    {
        let b =  x instanceof IC50Entity;
        return b;
    }
    getGridCellRenderer(x)
    {
        return new DGIC50GridCellRenderer();
    }

    renderIcon(x)
    {
        return ui.iconFA('ic50-alt');
    }
    renderMarkup(x)
    {
        let m = ui.span([this.renderIcon(), ' ', x.getDefaultValue()]); $(m).css('color', "white"); return m;
    }

        renderProperties(x)
    {
        return ui.divText(`Properties for ${x.getDefaultValue()}`);
    }
    renderTooltip(x)
    {
        return ui.divText(`${x.getDefaultValue()} is in the air!`);
    }
    renderCard(x, context) {
        return ui.bind(x, ui.divV([
            this.renderMarkup(x),
            ui.divText(`Context: ${context}`)
        ]), 'd4-gallery-item');
    }

    init()
    {
     this.registerParamFunc('Eat', (entity) => {
            grok.shell.info(`Activity ${entity.getDefaultValue()}`);
        });
    }

}