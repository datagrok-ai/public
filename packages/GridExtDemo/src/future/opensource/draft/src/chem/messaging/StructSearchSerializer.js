import {MessageSerializer} from "../../concurrent/messaging/MessageSerializer";
import {StructSearchRequestMessage} from "./StructSearchRequestMessage";
import {StructSearchResponseMessage} from "./StructSearchResponseMessage";
import {TextUtils} from "../../utils/TextUtils";

export class StructSearchSerializer extends MessageSerializer
{
    constructor()
    {
        super("StructSearch", StructSearchRequestMessage, StructSearchResponseMessage);
    }

    toRequestArgs(message)
    {
        const args = super.toRequestArgs(message);

        const bufSmilesFrag = StructSearchSerializer.Encoder.encode(message.getFragmentSmiles()).buffer;
        args[StructSearchSerializer.PARAM_FRAGMENT] = bufSmilesFrag;

        return [args, [bufSmilesFrag]];
    }

    toResponseArgs(message)
    {
        const args = super.toResponseArgs(message);

        const arFlags = message.getFlags();
        const nMolCount = arFlags.length;
        const viewFlags = new Uint8Array(nMolCount);
        for(var n=0; n<nMolCount; ++n)
        {
            viewFlags[n] = arFlags[n];
        }

        args[StructSearchSerializer.PARAM_FLAGS] = viewFlags.buffer;
        return [args, viewFlags.buffer];
    }

    toRequestMessage(args)
    {
        let clazz = this.getRequestMessageClass();

        let bufSmilesFrag = args[StructSearchSerializer.PARAM_FRAGMENT];
        let strSmilesFrag = TextUtils.buf2Str(bufSmilesFrag);

        let message = new clazz(strSmilesFrag);
        return message;
    }

    toResponseMessage(args)
    {
        const clazz = this.getResponseMessageClass();
        args = args[0];
        const buFlags = args[StructSearchSerializer.PARAM_FLAGS];
        const viewFlags = new Uint8Array(buFlags);
        const nMolCount =  viewFlags.length;
        const arFlags = new Array(nMolCount);

        for(var n=0; n<nMolCount; ++n)
        {

            if(viewFlags[n] !== 0 && viewFlags[n] !== 1)
                throw new Error("Flag value " + viewFlags[n] + " must be either 0 or 1, but it is " + viewFlags[n]);

            arFlags[n] = viewFlags[n] === 1;
        }

        const message = new clazz(arFlags);
        return message;
    }
}

StructSearchSerializer.Instance = new StructSearchSerializer();
StructSearchSerializer.PARAM_FRAGMENT = "fragment";
StructSearchSerializer.PARAM_FLAGS = "flags";