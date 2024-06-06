export async function swapTwoWays(stringToChange: string, characterToSwap1: string, characterToSwap2: string): Promise<string> {
    let result = stringToChange.replaceAll(characterToSwap1, '\t' + characterToSwap1 + '\t')
    result = result.replaceAll(characterToSwap2, '\t' + characterToSwap2 + '\t')
    result = result.replaceAll('\t' + characterToSwap2 + '\t', characterToSwap1)
    result = result.replaceAll('\t' + characterToSwap1 + '\t', characterToSwap2)

    return result;
}