import {describe, it, test, expect} from 'vitest';
import {Ids, fileId} from '../utils/ids';

describe('ids', () => {
  it('prefixes a file', () => {
    expect(fileId('a.ts')).toBe('file:a.ts');
  });

  describe('parse', () => {
    it('reads the scheme', () => {
      expect(new Ids().parse('file:a.ts')).toBe('file');
    });

    it.skip('reads a prefix', () => {});
    // it('reads a bare id', () => {});
  });

  test.each([['a', 'file:a'], ['b', 'file:b']])('fileId(%s) is %s', (file, id) => {
    expect(fileId(file)).toBe(id);
  });
});

test(`stands alone ${process.platform}`, () => {});
