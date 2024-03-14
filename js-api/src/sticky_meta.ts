import {toDart, toJs} from "./wrappers";
import {Schema} from "./entities";
import {IDartApi} from "./api/grok_api.g";
import { Column, DataFrame } from "./dataframe";

const api: IDartApi = <any>window;

export class StickyMeta {

    async getSchemas(): Promise<Schema[]> {
        return toJs(await api.grok_Sticky_GetSchemas());
    }

    setAllValues(schema: Schema, keys: Column, values: DataFrame): Promise<void> {
        return api.grok_Sticky_SetAllValues(toDart(schema), toDart(keys), toDart(values));
    }

    getAllValues(schema: Schema, keys: Column): Promise<DataFrame> {
        return toJs(api.grok_Sticky_GetAllValues(toDart(schema), toDart(keys)));
    }

    saveSchema(schema: Schema): Promise<void> {
        return api.grok_Sticky_SaveSchema(toDart(schema));
    }

    deleteSchema(id: string): Promise<void> {
        return api.grok_Sticky_DeleteSchema(id);
    }
}