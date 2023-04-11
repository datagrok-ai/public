import {MolfileHandlerBase} from './molfile-handler-base';
import {MolfileV2KHandler} from './molfile-v2k-handler';
import {MolfileV3KHandler} from './molfile-v3k-handler';

/** Defines the proper parser handler based on the molfile type  */
export class MolfileHandler {
  private constructor() {}

  static getInstance(molfile: string): MolfileHandlerBase {
    if (MolfileV2KHandler.validate(molfile))
      return new MolfileV2KHandler(molfile);
    else if (MolfileV3KHandler.validate(molfile))
      return new MolfileV3KHandler(molfile);
    else
      throw new Error('Malformed molfile');
  }
}
