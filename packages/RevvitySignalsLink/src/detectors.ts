import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

export async function addRevvityDetector(semTypeName: string, regexp: string) {
    DG.SemanticValue.registerRegExpDetector(semTypeName, regexp);
}

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