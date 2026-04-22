export function convertIdentifierFormatToRegexp(formats: string[]): string {
    if (!formats || formats.length === 0) {
        return '';
    }

    const regexps = formats.map(format => {

        let regexp = format.replace(/[.*+?^$()|[\]\\]/g, '\\$&');

        regexp = regexp.replace(/\{#+\}/g, (match) => {
            // Count the number of # symbols
            const hashCount = (match.match(/#/g) || []).length;
            // Replace with \d{count} pattern
            return `\\d{${hashCount}}`;
        });

        return `^${regexp}$`;
    });

    // Join all regexps with OR condition
    return regexps.join('|');
}