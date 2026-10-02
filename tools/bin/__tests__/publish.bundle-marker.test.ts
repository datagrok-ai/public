import {describe, it, expect} from 'vitest';
import {BUNDLE_MARKER, BUNDLE_MARKER_NAME, predatesBundleDetection} from '../commands/publish';

describe('bundle marker for servers before 1.28.0', () => {
  it('marks older servers', () => {
    for (const v of ['1.27.5', '1.27.11', '1.0.0', '0.9.9'])
      expect(predatesBundleDetection(v), v).toBe(true);
  });

  it('leaves 1.28.0 and newer alone', () => {
    for (const v of ['1.28.0', '1.28.1', '1.29.0', '2.0.0'])
      expect(predatesBundleDetection(v), v).toBe(false);
  });

  it('compares numerically, not lexically', () => {
    expect(predatesBundleDetection('1.9.0')).toBe(true);
    expect(predatesBundleDetection('1.100.0')).toBe(false);
  });

  it('treats an unusable version as older, so the marker is never skipped by mistake', () => {
    for (const v of ['', '1.28', 'bleeding-edge', 'x.y.z', undefined as any, null as any])
      expect(predatesBundleDetection(v), String(v)).toBe(true);
  });

  it('names the marker so that an older server recognises a bundled package', () => {
    expect(BUNDLE_MARKER_NAME).toBe('webpack.config.js');
    expect(BUNDLE_MARKER).toContain('grok publish');
    expect(() => new Function(BUNDLE_MARKER)).not.toThrow();
  });
});
