import {describe, expect, it} from 'vitest';
import {isPackageU2Usage} from '../utils/test-utils';

const bundle = (pkg: string) =>
  `http://localhost:8889/api/packages/published/files/${pkg}/1.2.0.X-624dc27e/OgYP0H7y1Wr4/1/dist/package.js`;
const usage = (...pkgs: string[]) => 'u2 usage: owner-missing: lC was built with no owner and now owns something ' +
  `nobody will dispose. (at ${pkgs.map((p, i) => `lC._mount (${bundle(p)}:2:${1000 + i})`).join(' < ')})`;

describe('grok test: the u2 usage lines that fail a package\'s run', () => {
  it('counts a line with a frame in the tested package\'s bundle, whatever the case of its name', () => {
    expect(isPackageU2Usage(usage('U2demo'), 'U2demo')).toBe(true);
    expect(isPackageU2Usage(usage('Bio', 'HitTriage'), 'hittriage')).toBe(true);
  });

  it('leaves out a line whose frames are all in other packages\' bundles', () => {
    expect(isPackageU2Usage(usage('Bio'), 'HitTriage')).toBe(false);
    expect(isPackageU2Usage(usage('Bio', 'PowerPack'), 'Chem')).toBe(false);
    expect(isPackageU2Usage(usage('HitTriageX'), 'HitTriage')).toBe(false);
  });

  it('counts a line that carries no frame, and every line when the tested package is unknown', () => {
    expect(isPackageU2Usage('u2: U2UsageError: u2 usage: owner-disposed: Form was built under a scope already ' +
      'disposed (campaign)', 'Chem')).toBe(true);
    expect(isPackageU2Usage(usage('Bio'), undefined)).toBe(true);
  });
});
