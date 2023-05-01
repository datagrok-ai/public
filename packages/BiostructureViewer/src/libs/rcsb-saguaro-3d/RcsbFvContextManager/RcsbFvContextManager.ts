import {Subject, Subscription} from 'rxjs';
import {RcsbFvStructureInterface} from "../RcsbFvStructure/RcsbFvStructure";
import {RcsbFvSequenceInterface} from "../RcsbFvSequence/RcsbFvSequence";
import {PluginContext} from "molstar/lib/mol-plugin/context";

/**Main Event Data Object Interface*/
export interface RcsbFvContextManagerInterface {
    eventType: EventType;
    eventData: string | UpdateConfigInterface | ((plugin: PluginContext) => void);
}

/**Event types*/
export enum EventType {
    UPDATE_CONFIG = "updateBoardConfig",
    PLUGIN_CALL = "pluginCall"
}

export interface UpdateConfigInterface {
    structurePanelConfig?:RcsbFvStructureInterface;
    sequencePanelConfig?:RcsbFvSequenceInterface;
}

/**rxjs Event Handler Object. It allows objects to subscribe methods and then, get(send) events to(from) other objects*/
export class RcsbFvContextManager {
    private readonly subject: Subject<RcsbFvContextManagerInterface> = new Subject<RcsbFvContextManagerInterface>();
    /**Call other subscribed methods
     * @param obj Event Data Structure Interface
     * */
    public next( obj: RcsbFvContextManagerInterface ):void {
        this.subject.next(obj);
    }
    /**Subscribe loadMethod
     * @return Subscription
     * */
    public subscribe(f:(x:RcsbFvContextManagerInterface)=>void):Subscription {
        return this.subject.asObservable().subscribe(f);
    }
    /**Unsubscribe all methods*/
    public unsubscribeAll():void {
        this.subject.unsubscribe();
    }
}
