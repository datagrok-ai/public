/* eslint-disable */
/**
 * Runtime APIs that current browsers still lack in some contexts. Imported as the first import from
 * `package.ts` / `package-test.ts`, so these are in place before any other module runs.
 */

(function installPolyfills() {
  // crypto.randomUUID exists only in a secure context: a stand served over plain HTTP on a hostname has none.
  const cryptoObj: any = (typeof crypto !== 'undefined') ? crypto : undefined;
  if (cryptoObj && typeof cryptoObj.randomUUID !== 'function') {
    try {
      cryptoObj.randomUUID = function randomUUID(): string {
        return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, (c) => {
          const r = Math.random() * 16 | 0;
          return (c === 'x' ? r : (r & 0x3 | 0x8)).toString(16);
        });
      };
    } catch (e) {
      console.warn('Grokky: could not polyfill crypto.randomUUID', e);
    }
  }

  // SpeechRecognition: Chrome ships only webkitSpeechRecognition
  const W: any = (typeof window !== 'undefined') ? window : undefined;
  if (W && !W.SpeechRecognition && W.webkitSpeechRecognition)
    W.SpeechRecognition = W.webkitSpeechRecognition;
})();
