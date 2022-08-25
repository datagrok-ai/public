import {ProgressRequestMessage} from "./ProgressRequestMessage";
import {ProgressResponseMessage} from "./ProgressResponseMessage";
import {MessageSerializer} from "./MessageSerializer";

export class ProgressSerializer extends MessageSerializer
{
    constructor() {
        super("Progress", ProgressRequestMessage, ProgressResponseMessage);
    }

    toRequestArgs(message)
    {
        const args = super.toRequestArgs(message);
        args[ProgressSerializer.PARAM_JOB_COUNT] = message.getJobCount();

        return args;
    }


    toResponseArgs(message) {
        const args = super.toResponseArgs(message);
        const bCancel = message.isCancel();
        args[ProgressSerializer.PARAM_CANCEL] = bCancel;
        return args;
    }


    toRequestMessage(args) {

        if(super.toRequestMessage(args) === null)
            return null;

        const clazz = this.getRequestMessageClass();
        const nJobCount = args[ProgressSerializer.PARAM_JOB_COUNT];
        const message = new clazz(nJobCount);
        return message;
    }

    toResponseMessage(args) {

        if(super.toResponseMessage(args) === null)
            return null;

        const bCancel = args[ProgressSerializer.PARAM_CANCEL];
        const clazz = this.getResponseMessageClass();
        return new clazz(bCancel);
    }
}

ProgressSerializer.Instance = new ProgressSerializer();
ProgressSerializer.PARAM_JOB_COUNT = "job_count";
ProgressSerializer.PARAM_CANCEL = "cancel";
