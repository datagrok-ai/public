/* eslint-disable */
/**
 * Runtime APIs that current browsers still lack in some contexts. Imported as the first import from
 * `package.ts` / `package-test.ts`, so these are in place before any other module runs.
 */

(function installPolyfills() {
  // SpeechRecognition: Chrome ships only webkitSpeechRecognition
  const W: any = (typeof window !== 'undefined') ? window : undefined;
  if (W && !W.SpeechRecognition && W.webkitSpeechRecognition)
    W.SpeechRecognition = W.webkitSpeechRecognition;
})();
