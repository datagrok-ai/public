import {ProgressPctRequestMessage} from "./ProgressPctRequestMessage";
import {ProgressPctResponseMessage} from "./ProgressPctResponseMessage";
import {MessageSerializer} from "../../../concurrent/messaging/MessageSerializer";

export class ProgressPctSerializer extends MessageSerializer
{
    constructor() {
        super("ProgressPct", ProgressPctRequestMessage, ProgressPctResponseMessage);
    }

    toRequestArgs(message)
    {
        let args = super.toRequestArgs(message);
        args[ProgressPctSerializer.PARAM_PERCENT] = message.getPct();
        return args;
    }

    toResponseArgs(message) {
        let args = super.toResponseArgs(message);
        return args;
    }


    toRequestMessage(args)
    {
        let clazz = this.getRequestMessageClass();
        let nPct = args[ProgressPctSerializer.PARAM_PERCENT];
        let message = new clazz(nPct);
        return message;
    }

    toResponseMessage(args) {

        if(super.toResponseMessage(args) === null)
            return null;

        let clazz = this.getResponseMessageClass();
        return new clazz();
    }
}

ProgressPctSerializer.Instance = new ProgressPctSerializer();
ProgressPctSerializer.PARAM_PERCENT = "percent";
