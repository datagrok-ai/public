export class AsyncUtils
{
    constructor() {
        throw new Error("Never create instances of this class");
    }
}

AsyncUtils.sleep = function(ms) {
    return new Promise(resolve => setTimeout(resolve, ms));
}
