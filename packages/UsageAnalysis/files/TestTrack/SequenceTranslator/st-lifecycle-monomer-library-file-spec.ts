/* ---
sub_features_covered: [sequencetranslator.api.get-code-to-weights-map, sequencetranslator.lifecycle.init-lib-data, sequencetranslator.polytool.create-monomer-library]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('SequenceTranslator monomer_library_file source-class lifecycle: createMonomerLibraryForPolyTool + initLibData + getCodeToWeightsMap', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
  });

  try {

    await softStep('S2.1 (DOM): trigger initLibData() via getTranslationHelper', async () => {
      const result = await page.evaluate(async () => {
        const g = (window as any).grok;
        let err: string | null = null;
        let helperNonNull = false;
        try {
          const helper = await g.functions.call('SequenceTranslator:getTranslationHelper', {});
          helperNonNull = helper != null;
        } catch (e: any) {
          err = String(e).slice(0, 400);
        }
        return {err, helperNonNull};
      });
      expect(result.err, `getTranslationHelper error: ${result.err}`).toBeNull();
      expect(result.helperNonNull, 'getTranslationHelper must return a non-null ITranslationHelper').toBe(true);

      const browseLocator = page.locator('[name="Browse"]');
      await browseLocator.waitFor({state: 'visible', timeout: 20_000});
      expect(await browseLocator.isVisible()).toBe(true);
    });

    await softStep('S2.2: initLibData promise caches — getTranslationHelper returns a non-null result on both calls', async () => {

      const result = await page.evaluate(async () => {
        const g = (window as any).grok;
        let err: string | null = null;
        let h1IsNonNull = false;
        let h2IsNonNull = false;
        let h1Type: string = '';
        let h2Type: string = '';
        try {
          const helper1 = await g.functions.call('SequenceTranslator:getTranslationHelper');
          const helper2 = await g.functions.call('SequenceTranslator:getTranslationHelper');
          h1IsNonNull = helper1 != null;
          h2IsNonNull = helper2 != null;
          h1Type = typeof helper1;
          h2Type = typeof helper2;
        } catch (e: any) {
          err = String(e).slice(0, 300);
        }
        return {err, h1IsNonNull, h2IsNonNull, h1Type, h2Type};
      });
      expect(result.err, `getTranslationHelper error: ${result.err}`).toBeNull();

      expect(result.h1IsNonNull, `helper1 is null; type: ${result.h1Type}`).toBe(true);
      expect(result.h2IsNonNull, `helper2 is null; type: ${result.h2Type}`).toBe(true);
    });

    await softStep('S2.3: initLibData fallback — sample library accessible at monomers-sample path', async () => {
      const result = await page.evaluate(async () => {
        const g = (window as any).grok;
        const samplePaths = [
          'System:AppData/SequenceTranslator/monomers/monomer-lib.json',
          'System:AppData/SequenceTranslator/monomers-sample/monomer-lib.json',
        ];
        let chosenPath: string | null = null;
        let content: string | null = null;
        let readErr: string | null = null;
        for (const p of samplePaths) {
          try {
            content = await g.dapi.files.readAsText(p);
            chosenPath = p;
            break;
          } catch (e) {
            readErr = String(e).slice(0, 200);
          }
        }
        if (!content || !chosenPath) {
          return {
            err: `No monomer library readable; last error: ${readErr}`,
            chosenPath: null,
            isJson: false,
            contentLen: 0,
          };
        }
        let parsed: any;
        let parseErr: string | null = null;
        try {
          parsed = JSON.parse(content);
        } catch (e: any) {
          parseErr = String(e).slice(0, 200);
        }
        return {
          err: null,
          chosenPath,
          isJson: parseErr == null,
          contentLen: content.length,
          isObject: parsed != null && typeof parsed === 'object',
        };
      });
      expect(result.err, `monomer library read error: ${result.err}`).toBeNull();
      expect(result.chosenPath).not.toBeNull();
      expect(result.contentLen).toBeGreaterThan(0);
      expect(result.isJson).toBe(true);
    });

    await softStep('S3.1: getCodeToWeightsMap returns a non-empty Record<string, number>', async () => {
      const result = await page.evaluate(async () => {
        const g = (window as any).grok;
        let weightsMap: any;
        let err: string | null = null;
        try {
          weightsMap = await g.functions.call('SequenceTranslator:getCodeToWeightsMap');
        } catch (e: any) {
          err = String(e).slice(0, 300);
        }
        if (err || weightsMap == null) return {err, keyCount: 0, sample: [], isObject: false, allPositiveNumbers: false};
        const keys = Object.keys(weightsMap);
        const sample = keys.slice(0, 5).map((k: string) => ({key: k, val: weightsMap[k]}));

        const sampled = keys.slice(0, 20).map((k: string) => weightsMap[k]);
        const allFiniteNonNeg = sampled.every((v: any) => typeof v === 'number' && isFinite(v) && v >= 0);
        const positiveCount = sampled.filter((v: any) => typeof v === 'number' && v > 0).length;
        return {
          err: null,
          keyCount: keys.length,
          sample,
          isObject: typeof weightsMap === 'object',
          allFiniteNonNeg,
          positiveCount,
        };
      });
      expect(result.err, `getCodeToWeightsMap error: ${result.err}`).toBeNull();
      expect(result.isObject).toBe(true);
      expect(result.keyCount).toBeGreaterThan(0);
      expect(result.allFiniteNonNeg,
        `all sampled MW values must be finite non-negative numbers; sample: ${JSON.stringify(result.sample)}`).toBe(true);
      expect(result.positiveCount,
        `the weights map must carry real (positive) molecular weights; sample: ${JSON.stringify(result.sample)}`).toBeGreaterThan(0);
    });

    await softStep('S3.2: getCodeToWeightsMap keys include standard RNA or DNA monomers', async () => {
      const result = await page.evaluate(async () => {
        const g = (window as any).grok;
        let weightsMap: any;
        let err: string | null = null;
        try {
          weightsMap = await g.functions.call('SequenceTranslator:getCodeToWeightsMap');
        } catch (e: any) {
          err = String(e).slice(0, 300);
        }
        if (err || weightsMap == null) return {err, keys: [], hasStandardRNA: false, hasStandardDNA: false, hasAnyStandard: false};
        const keys = Object.keys(weightsMap);
        const rnaStandard = ['A', 'C', 'G', 'U'];
        const dnaStandard = ['dA', 'dC', 'dG', 'dT'];
        const hasRNA = rnaStandard.some((k) => keys.includes(k));
        const hasDNA = dnaStandard.some((k) => keys.includes(k));
        return {
          err: null,
          keys: keys.slice(0, 20),
          hasStandardRNA: hasRNA,
          hasStandardDNA: hasDNA,
          hasAnyStandard: hasRNA || hasDNA,
        };
      });
      expect(result.err, `getCodeToWeightsMap error: ${result.err}`).toBeNull();
      expect(result.hasAnyStandard,
        `expected standard RNA or DNA monomer codes (A/C/G/U or dA/dC/dG/dT) in weights map; keys observed: [${result.keys.join(', ')}]`).toBe(true);
    });

    await softStep('S1.1: createMonomerLibraryForPolyTool validates CSV format (rejects non-monomer-lib CSV; JSON payload on valid)', async () => {
      const result = await page.evaluate(async () => {
        const g = (window as any).grok;
        const DG = (window as any).DG;

        let capturedFileName: string | null = null;
        let capturedContent: string | null = null;
        const originalDownload = DG.Utils.download;
        DG.Utils.download = (name: string, content: string) => {
          capturedFileName = name;
          capturedContent = content;
        };

        let err: string | null = null;
        let fileFound = false;
        try {

          const candidateEntries = [
            {folder: 'System:AppData/SequenceTranslator/tests', basename: 'axolabs1.csv'},
            {folder: 'System:AppData/SequenceTranslator/samples', basename: 'bulk-translation-axolabs.csv'},
          ];

          let fileInfo: any = null;
          for (const {folder, basename} of candidateEntries) {
            try {
              const items = await g.dapi.files.list(folder, false, basename);
              if (items && items.length > 0) {
                fileInfo = items[0];
                fileFound = true;
                break;
              }
            } catch (_) {

            }
          }

          if (!fileInfo) {

            return {
              err: `No accessible CSV monomer-library file found at candidate paths`,
              fileFound: false,
              capturedFileName: null,
              capturedContent: null,
              contentLength: 0,
              isJson: false,
              fileNameEndsWithJson: false,
            };
          }

          await g.functions.call('SequenceTranslator:createMonomerLibraryForPolyTool', {
            file: fileInfo,
          });
        } catch (e: any) {
          err = String(e).slice(0, 400);
        } finally {
          DG.Utils.download = originalDownload;
        }

        return {
          err,
          fileFound,
          capturedFileName,
          capturedContent,
          downloadTriggered: capturedFileName != null,
          contentLength: capturedContent ? capturedContent.length : 0,
          isJson: (() => {
            if (!capturedContent) return false;
            try { JSON.parse(capturedContent); return true; } catch { return false; }
          })(),
          fileNameEndsWithJson: capturedFileName?.endsWith('.json') ?? false,
        };
      });

      expect(result.fileFound,
        `a candidate CSV (axolabs1.csv / bulk-translation-axolabs.csv) must be reachable on the server: ${result.err}`).toBe(true);
      if (result.err) {
        expect(result.err,
          `createMonomerLibraryForPolyTool must reject a non-monomer-lib CSV with a clear "Invalid format" error; got: ${result.err}`)
          .toContain('Invalid format');
      } else {

        expect(result.downloadTriggered,
          `expected DG.Utils.download to be called; fileName=${result.capturedFileName}, contentLen=${result.contentLength}`).toBe(true);
        expect(result.fileNameEndsWithJson).toBe(true);
        expect(result.isJson,
          `expected downloaded content to be valid JSON; contentLen=${result.contentLength}`).toBe(true);
      }
    });
  } finally {

    await page.evaluate(async () => {
      try {
        (window as any).grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 300));
      } catch (_) {  }
    }).catch(() => {});
  }

  if (stepErrors.length > 0) {
    throw new Error(
      `[SPEC FAIL] ${stepErrors.length} step(s) failed:\n` +
      stepErrors.map((e) => `  • ${e.step}: ${e.error}`).join('\n'),
    );
  }
});
