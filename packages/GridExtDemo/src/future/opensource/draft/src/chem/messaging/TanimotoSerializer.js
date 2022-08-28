import {MessageSerializer} from "../../concurrent/messaging/MessageSerializer";
import {TanimotoRequestMessage} from "./TanimotoRequestMessage";
import {TanimotoResponseMessage} from "./TanimotoResponseMessage";
import {TextUtils} from "../../utils/TextUtils";

export class TanimotoSerializer extends MessageSerializer {
    constructor()
    {
        super("Tanimoto", TanimotoRequestMessage, TanimotoResponseMessage);
    }

    toRequestArgs(message)
    {
        const args = super.toRequestArgs(message);

        const bufSmilesFrag = MessageSerializer.Encoder.encode(message.getFragmentSmiles()).buffer;
        args[TanimotoSerializer.PARAM_FRAGMENT] = bufSmilesFrag;
        return [args, [bufSmilesFrag]];
    }

    toResponseArgs(message)
    {
        const args = super.toResponseArgs(message);
        const arScores = message.getScores();
        const nMolCount = arScores.length;
        const viewScores = new Float32Array(nMolCount);
        for(var n=0; n<nMolCount; ++n)
        {
            viewScores[n] = arScores[n];
        }

        args[TanimotoSerializer.PARAM_SCORES] = viewScores.buffer;
        return [args, [viewScores.buffer]];
    }

    toRequestMessage(args)
    {
        const clazz = this.getRequestMessageClass();

        const bufSmilesFrag = args[TanimotoSerializer.PARAM_FRAGMENT]
        const strSmilesFrag = TextUtils.buf2Str(bufSmilesFrag);

        const message = new clazz(strSmilesFrag);
        return message;
    }

    toResponseMessage(args)
    {
        const clazz = this.getResponseMessageClass();

        args = args[0];
        const bufScores = args[TanimotoSerializer.PARAM_SCORES];
        const viewScores = new Float32Array(bufScores);
        const nMolCount =  viewScores.length;
        const arScores = new Array(nMolCount);

        for(var n=0; n<nMolCount; ++n)
        {
            arScores[n] = viewScores[n];

            if(arScores[n] === null || 0.0 > arScores[n] || arScores[n] > 1.0)
                throw new Error("Tanimoto score value " + arScores[n] + " is out of range [0,1] for record index " + n);
        }

        const message = new clazz(arScores);
        return message;
    }


}

TanimotoSerializer.Instance = new TanimotoSerializer();
TanimotoSerializer.PARAM_FRAGMENT = "fragment";
TanimotoSerializer.PARAM_SCORES = "scores";