export function endsWithAny(optionsArray: string[], value: string){
    return optionsArray.some(char => value.endsWith(char))
}

export function normalizeStrings(values: any[]){
    if (values.every(item => typeof item === 'string')){
        return values.map(it => it.toLowerCase());
    }
    return values;
}