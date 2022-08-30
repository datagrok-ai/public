import {ErrorUtils} from "../utils/ErrorUtils";

export class NotImplementedError extends Error
{
    constructor() {
        super(NotImplementedError.message());
    }
}

NotImplementedError.message = function ()
{
  const strMathodName = ErrorUtils.getCallerMethodName();
  return "Method '" + strMathodName + "()' is not implemented."
}