import {MessageSerializer} from "./MessageSerializer";
import {DisposeRequestMessage} from "./DisposeRequestMessage";
import {DisposeResponseMessage} from "./DisposeResponseMessage";

export class DisposeSerializer extends MessageSerializer {
    constructor() {
        super("Dispose", DisposeRequestMessage, DisposeResponseMessage);
    }
}

DisposeSerializer.Instance = new DisposeSerializer();