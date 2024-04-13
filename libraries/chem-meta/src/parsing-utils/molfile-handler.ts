import {MolfileHandlerBase} from './molfile-handler-base';
import {MolfileV2KHandler} from './molfile-v2k-handler';
import {MolfileV3KHandler} from './molfile-v3k-handler';

/** Defines the proper handler based on the molfile type  */
export class MolfileHandler {
  private constructor() {}

  static getInstance(molfile: string): MolfileHandlerBase {
    if (MolfileHandler.isMolfileV2K(molfile))
      return new MolfileV2KHandler(molfile);
    else if (MolfileHandler.isMolfileV3K(molfile))
      return new MolfileV3KHandler(molfile);
    else
      throw new Error('Malformed molfile');
  }

  static isMolfileV2K(molfile: string): boolean {
    return MolfileV2KHandler.isValidMolfile(molfile);
  }

  static isMolfileV3K(molfile: string): boolean {
    return MolfileV3KHandler.isValidMolfile(molfile);
  }
}
