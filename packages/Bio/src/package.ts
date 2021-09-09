/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { SequenceAlignment, Aligned } from './seq_align';
export let _package = new DG.Package();


//name: sequenceAlignment
//input: string alignType {choices: ['Local alignment', 'Global alignment']}
//input: string alignTable {choices: ['AUTO', 'NUCLEOTIDES', 'BLOSUM45', 'BLOSUM50','BLOSUM62','BLOSUM80','BLOSUM90','PAM30','PAM70','PAM250','SCHNEIDER','TRANS']}
//input: double gap
//input: string seq1
//input: string seq2
//output: object res 
export function sequenceAlignment(alignType:string, alignTable:string, gap:number, seq1:string, seq2:string){
    let res = alignType == 'Local alignment' ? toAlign.smithWaterman() : toAlign.needlemanWunch();
    return res
}
