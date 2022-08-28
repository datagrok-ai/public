import {MessageSerializer} from "../../concurrent/messaging/MessageSerializer";
import {InitRequestMessage} from "./InitRequestMessage";
import {InitResponseMessage} from "./InitResponseMessage";
import {TextUtils} from "../../utils/TextUtils";

export class InitSerializer extends MessageSerializer {
    constructor() {
        super("Init", InitRequestMessage, InitResponseMessage);
    }

    toRequestArgs(message)
    {
        const arArgs = super.toRequestArgs(message);
        //let obArgOne = arArgs[0];
        const bufWebRoot = MessageSerializer.Encoder.encode(message.getWebRoot()).buffer;
        arArgs[InitSerializer.PARAM_WEBROOT] = bufWebRoot;

        const nLength = message.m_nTo - message.m_nFr + 2;
        const arSmilesBufs = new Array(nLength);
        for (var n = message.m_nFr; n <= message.m_nTo; ++n)
        {
            arSmilesBufs[n - message.m_nFr] = MessageSerializer.Encoder.encode(message.m_arSmilesOrMolFiles[n]).buffer;
        }

        arSmilesBufs[nLength -1] = bufWebRoot;

        arArgs[InitSerializer.PARAM_SMILES_OR_MOLFILE] = arSmilesBufs;
        arArgs[InitSerializer.PARAM_IS_MOLFILE] = message.isMolFile();
        arArgs[InitSerializer.PARAM_CPU_ID] = message.getCPUId();

        return [arArgs, arSmilesBufs];
    }

    toResponseArgs(message) {
        const args = super.toResponseArgs(message);
        args[InitSerializer.PARAM_MOL_COUNT] = message.getMolCount();
        args[InitSerializer.CANCELLED] = message.isCancelled();

        return args;
    }

    toRequestMessage(args) {
        const clazz = this.getRequestMessageClass();

        const bufWebRoot = args[InitSerializer.PARAM_WEBROOT];
        const strWebRoot = TextUtils.buf2Str(bufWebRoot);

        const arSmilesBufs = args[InitSerializer.PARAM_SMILES_OR_MOLFILE];
        const bMolFile = args[InitSerializer.PARAM_IS_MOLFILE];
        const nCPUId = args[InitSerializer.PARAM_CPU_ID];

        arSmilesBufs.pop();//last was the webRoot

        const nMolCount = arSmilesBufs.length;
        //this.m_nMolCount = nMolCount;

        for (let nMol = 0; nMol < nMolCount; ++nMol)
        {
          arSmilesBufs[nMol] = TextUtils.buf2Str(arSmilesBufs[nMol]);
        }

        const message = new clazz(strWebRoot, arSmilesBufs, 0, nMolCount - 1, bMolFile, nCPUId);
        return message;
    }


    toResponseMessage(args) {
        const clazz = this.getResponseMessageClass();
        const nMolCount = args[InitSerializer.PARAM_MOL_COUNT];
        const bCancelled = args[InitSerializer.CANCELLED];
        const message = new clazz(nMolCount, bCancelled);
        return message;
    }

}
InitSerializer.Instance = new InitSerializer();
InitSerializer.PARAM_WEBROOT = "webRoot";
InitSerializer.PARAM_SMILES_OR_MOLFILE = "smiles_or_molfiles";
InitSerializer.PARAM_IS_MOLFILE = "is_molfile";
InitSerializer.PARAM_MOL_COUNT = "mol_count";
InitSerializer.PARAM_CPU_ID = "cpu_idx";
InitSerializer.CANCELLED = "cancelled";
