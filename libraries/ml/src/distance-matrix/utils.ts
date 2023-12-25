export const isNil = (x: any) => x === null || x === undefined;

export function insertSmaller(distancesAr: number[], indexes: number[], num: number, index: number) {
    if (num > distancesAr[distancesAr.length-1]) {
        return;
    }
    const newPosition = distancesAr.findIndex((v) => num < v);
    distancesAr.pop();
    distancesAr.splice(newPosition, 0, num);
    indexes.pop();
    indexes.splice(newPosition, 0, index);
}

export function insertLarger(distancesAr: number[], indexes: number[], num: number, index: number) {
    if (num < distancesAr[distancesAr.length-1]) {
        return;
    }
    const newPosition = distancesAr.findIndex((v) => num > v);
    distancesAr.pop();
    distancesAr.splice(newPosition, 0, num);
    indexes.pop();
    indexes.splice(newPosition, 0, index);
}
