export function isAlpha(sym: string): boolean {
    const charCode = sym.charCodeAt(0);
    return (charCode > 64 && charCode < 91) ||
        (charCode > 96 && charCode < 123);
}
