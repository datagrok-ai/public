import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package';

/* eslint-disable camelcase */
export class DummyWebEditorMonomer {
  public readonly backgroundcolor: string = '#FF4444';
  public readonly linecolor: string = '#800000';
  public readonly textcolor: string = '#FFFFFF';

  // no effect
  // public readonly na: string = 'NA';
  // public readonly nature: string = 'NA';

  /** R-Group index os single digit only is allowed in Pistoia code */
  public readonly at: { [rg: string]: string } = {
    R1: 'H', R2: 'H', R3: 'H', R4: 'H', R5: 'H', R6: 'H', R7: 'H', R8: 'H', R9: 'H'
  };

  constructor(
    public readonly biotype: string,
    /** symbol, elem */ public readonly id: string,
    /** name */ public readonly n: string = 'missing',
    /** molfile*/ public readonly m?: string,
  ) {
    if (!this.id)
      throw new Error('Not supported undefined [id].');

    return new Proxy(this, {
      get: (target: any, prop: string | symbol, _receiver: any) => {
        const logPrefix = `${this.toLog()}.proxy.get( '${String(prop)}' )`;
        switch (prop) {
        case 'id':
        case 'at':
        case 'na':
        case 'nature':
        case 'backgroundcolor':
        case 'linecolor':
        case 'textcolor': {
          // these attributes are known and handled
          return target[prop];
        }
        case 'rs': {
          /* R-Group count (?), 2 is default for monomer "?" */
          throw new Error('Not supported get "rs" attribute.');
        }
        default:
          // Unexpected attribute requested
          _package.logger.warning(`${logPrefix}`);
          return target[prop];
        }
      },
      set: (target, prop, value: any, _receiver: any): boolean => {
        throw new Error(`Not supported set '${String(prop)}' attribute, any.`);
      },
      apply: (target: any, thisArg: any, argArray: any[]): any | undefined => {
        throw new Error(`Not supported call (apply) '${thisArg.toString()}' method, any.`);
      }
    });
  }

  private static objCounter: number = -1;
  private readonly objId: number = ++DummyWebEditorMonomer.objCounter;

  protected toLog(): string {
    return `Helm: DummyWebEditorMonomer<${this.objId}>`;
  }
}
