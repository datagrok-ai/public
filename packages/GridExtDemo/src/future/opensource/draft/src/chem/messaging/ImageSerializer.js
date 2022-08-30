import {MessageSerializer} from "../../concurrent/messaging/MessageSerializer";
import {ImageRequestMessage} from "./ImageRequestMessage";
import {ImageResponseMessage} from "./ImageResponseMessage";

export class ImageSerializer extends MessageSerializer
{
    constructor() {
        super("Image", ImageRequestMessage, ImageResponseMessage);
    }

    toRequestArgs(message)
    {
        const args = super.toRequestArgs(message);

        args[ImageSerializer.PARAM_WIDTH] = message.getW();
        args[ImageSerializer.PARAM_HEIGHT] = message.getH();
        args[ImageSerializer.PARAM_MOL_INDEX] = message.getMolIndex();

        return args;
    }

    toResponseArgs(message)
    {
        const args = super.toResponseArgs(message);

        args[ImageSerializer.PARAM_IMAGE] = message.getBitmapImage();
        return [args, [message.getBitmapImage()]];
    }


    toRequestMessage(args)
    {
        const clazz = this.getRequestMessageClass();

        const nW = args[ImageSerializer.PARAM_WIDTH];
        const nH = args[ImageSerializer.PARAM_HEIGHT];
        const nMolIndex = args[ImageSerializer.PARAM_MOL_INDEX];

        const message = new clazz(nMolIndex, nW, nH);
        return message;
    }

    toResponseMessage(args)
    {
        args = args[0];

        const clazz = this.getResponseMessageClass();
        const bitmap = args[ImageSerializer.PARAM_IMAGE];
        const message = new clazz(bitmap);
        return message;
    }
}

ImageSerializer.Instance = new ImageSerializer();
ImageSerializer.PARAM_WIDTH = "width";
ImageSerializer.PARAM_HEIGHT = "height";
ImageSerializer.PARAM_MOL_INDEX = "mol_index";
ImageSerializer.PARAM_IMAGE = "image";