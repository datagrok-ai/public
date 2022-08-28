import {MessageSerializer} from "../../../concurrent/messaging/MessageSerializer";
import {InitRequestMessage} from "./InitRequestMessage";
import {InitResponseMessage} from "./InitResponseMessage";

export class InitSerializer extends MessageSerializer
{
    constructor()
    {
        super("Init", InitRequestMessage, InitResponseMessage);
    }

    toRequestArgs(message)
    {
        let args = super.toRequestArgs(message);
        let eCanvas = message.getCanvas();
        eCanvas = eCanvas.transferControlToOffscreen();

        args[InitSerializer.PARAM_CANVAS] = eCanvas;
        return [args, [eCanvas]];
    }

    toResponseArgs(message) {
        let args = super.toResponseArgs(message);
        return args;
    }

    toRequestMessage(args)
    {
        let clazz = this.getRequestMessageClass();

        let eCanvas = args[InitSerializer.PARAM_CANVAS]

        let message = new clazz(eCanvas);
        return message;
    }

    toResponseMessage(args) {

        if(super.toResponseMessage(args) === null)
            return null;

        let clazz = this.getResponseMessageClass();
        return new clazz();
    }
}

InitSerializer.Instance = new InitSerializer();
InitSerializer.PARAM_CANVAS = "offscreen_canvas";
