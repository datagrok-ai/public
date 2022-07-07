import * as DG from 'datagrok-api/dg';
import { WebLogo, SplitterFunc } from '@datagrok-libraries/bio/src/viewers/web-logo';
import * as grok from 'datagrok-api/grok';

export const HELM_CORE_LIB_MONOMER_COL = 'symbol';
export const HELM_CORE_LIB_MOLFILE_COL = 'molfile';
export const HELM_CORE_LIB_FILENAME = '/samples/HELMCoreLibrary.json';

export function getMolfilesFromSeq(col: DG.Column, monomersLib: DG.DataFrame): string[][] | null {
    const units = col.tags[DG.TAGS.UNITS];
    const sep = col.getTag('separator');
    const splitterFunc: SplitterFunc = WebLogo.getSplitter(units, sep);
    const monomersDict = createMomomersMolDict(monomersLib);
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

export function createMomomersMolDict(lib: DG.DataFrame): {[key: string]: string} {
    const dict: {[key: string]: string} = {};
    const monmersCol = lib.col(HELM_CORE_LIB_MONOMER_COL);
    const molCol = lib.col(HELM_CORE_LIB_MOLFILE_COL);
    for (let i = 0; i < lib.rowCount; ++i) {
        dict[monmersCol!.get(i)] = molCol!.get(i);
    }
    return dict;
}