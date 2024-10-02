export type mmDistanceFunctionType = (seq1: string, seq2: string) => number;
export type mmDistanceFunctionArgs = {
    scoringMatrix: number[][];
    alphabetIndexes: {[id:string]:number};
    threshold?: number;
    maxLength?: number;
};
