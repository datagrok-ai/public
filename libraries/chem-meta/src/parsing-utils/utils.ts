export function isAlpha(charCode: number): boolean {
    return (charCode > 64 && charCode < 91) ||
        (charCode > 96 && charCode < 123);
}
