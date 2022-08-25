import {MessageSerializer} from "../../../concurrent/messaging/MessageSerializer";
import {TextRequestMessage} from "./TextRequestMessage";
import {TextResponseMessage} from "./TextResponseMessage";
import {TextUtils} from "../../../utils/TextUtils";

export class TextSerializer extends MessageSerializer
{
    constructor()
    {
        super("Text", TextRequestMessage, TextResponseMessage);
    }

    toRequestArgs(message)
    {
        let args = super.toRequestArgs(message);

        let bufText = TextSerializer.Encoder.encode(message.getText()).buffer;
        args[TextSerializer.PARAM_TEXT] = bufText;
        let bUp = message.isUp();
        args[TextSerializer.PARAM_FLAG_UP] = bUp;

        return [args, [bufText]];
    }

    toResponseArgs(message)
    {
        let args = super.toResponseArgs(message);
        return args;
    }

    toRequestMessage(args)
    {
        let clazz = this.getRequestMessageClass();

        let bufText = args[TextSerializer.PARAM_TEXT];
        let strText = TextUtils.buf2Str(bufText);
        let bUp = args[TextSerializer.PARAM_FLAG_UP];

        let message = new clazz(strText, bUp);
        return message;
    }

    toResponseMessage(args)
    {
        let message = super.toResponseMessage(args);
        return message;
    }
}

TextSerializer.Instance = new TextSerializer();
TextSerializer.PARAM_TEXT = "text";
TextSerializer.PARAM_FLAG_UP = "is_up";