export class ErrorUtils {
    constructor() {
        throw new Error("Never create instances of this class");
    }
}

ErrorUtils.getCallerMethodName = function(nLevel)
{
    const e = new Error("dummy");
    let stack = e.stack .split('\n');
    stack = stack[4];
        // " at functionName ( ..." => "functionName"
     stack = stack.replace(/^\s+at\s+(.+?)\s.+/g, '$1');

    return stack;
}

ErrorUtils.verifyNotUndefinedNotNull = function (ob)
{
    if(ob === undefined || ob === null)
        throw new Error("Object cannot be undefined or null.");
}

ErrorUtils.verifyClass = function (ob, clazz)
{
    if(!(ob instanceof clazz))
        throw new Error("Object must be an instance of the '" + clazz.name + "' class, but it is " + ob);
}

ErrorUtils.verifyType = function (ob, obType)
{
    if(typeof ob !== obType.name.toLowerCase())
        throw new Error("Object " + ob + " must be of type '" + obType + "', but it is of type " + (typeof ob));
}
