import * as DG from 'datagrok-api/dg';
import { WebLogo, SplitterFunc } from '@datagrok-libraries/bio/src/viewers/web-logo';
import * as grok from 'datagrok-api/grok';

export const HELM_CORE_LIB_FILENAME = '/samples/HELMCoreLibrary.json';
export const HELM_CORE_LIB_MONOMER_SYMBOL = 'symbol';
export const HELM_CORE_LIB_MOLFILE = 'molfile';
export const HELM_CORE_FIELDS = ['symbol', 'molfile', 'rgroups', 'name'];

export function getMolfilesFromSeq(col: DG.Column, monomersLibObject: any[]): any[][] | null {
    const units = col.tags[DG.TAGS.UNITS];
    const sep = col.getTag('separator');
    const splitterFunc: SplitterFunc = WebLogo.getSplitter(units, sep);
    const monomersDict = createMomomersMolDict(monomersLibObject);
    const molFiles = [];
    for (let i = 0; i < col.length; ++i) {
        const monomers = splitterFunc(col.get(i));
        const molFilesForSeq = [];
        for (let j = 0; j < monomers.length; ++j) {
            if (monomers[j]) {
                if (!monomersDict[monomers[j]]) {
                    grok.shell.warning(`Monomer ${monomers[j]} is missing in HELM library. Structure cannot be created`);
                    return null;
                }
                molFilesForSeq.push(monomersDict[monomers[j]])
            }
        }
        molFiles.push(molFilesForSeq);
    }
    return molFiles;
  }

export function createMomomersMolDict(lib: any[]): {[key: string]: string | any} {
    const dict: {[key: string]: string | any} = {};
    lib.forEach(it => {
        const monomerObject: {[key: string]: any} = {};
        HELM_CORE_FIELDS.forEach(field => {
            monomerObject[field] = it[field];
        })
        dict[it[HELM_CORE_LIB_MONOMER_SYMBOL]] = monomerObject;
    })
    return dict;
}