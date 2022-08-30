export class MessageSerializer
{
    constructor(strName, classRequestMsg, classResponseMsg)
    {
        this.m_strName = strName;
        this.m_classRequestMsg = classRequestMsg;
        this.m_classResponseMsg = classResponseMsg;
    }

    getName() {return this.m_strName;}

    getRequestMessageClass()
    {
        return this.m_classRequestMsg;
    }

    getResponseMessageClass()
    {
        return this.m_classResponseMsg;
    }

    toRequestArgs(message)
    {
        if((!message instanceof this.m_classRequestMsg))
            throw new Error("Message object is not an instanceof '" + this.m_classRequestMsg.constructor.name + "'");

        let args = {msg : this.m_strName};//[{msg : this.m_strName}, []];
        return args;
    }

    toResponseArgs(message)
    {
        if((!message instanceof this.m_classResponseMsg))
            throw new Error("Message object is not an instanceof '" + this.m_classResponseMsg.constructor.name + "'");

        let args = {msg : this.m_strName};
        return args;//[{msg : this.m_strName}, []];
    }

    toRequestMessage(args)
    {
        if(args.msg !== this.m_strName)
            return null;

       // let clazz =  this.getRequestMessageClass();
        //let message = new clazz();
        return undefined;
    }

    toResponseMessage(args)
    {
        if(args.msg !== this.m_strName)
            return null;
         /*
        let clazz =  this.getResponseMessageClass();
        let message = new clazz();
        return message;*/
        return undefined;
    }
}

MessageSerializer.Encoder = new TextEncoder();