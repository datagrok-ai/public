export function tryParseJson(s: string) : any {
    try {
        return JSON.parse(s);
    }
    catch (_) {
        return null;
    }
}