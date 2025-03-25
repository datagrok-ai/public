class NlpPackageDetectors extends DG.Package {
  // `isTextFile` appears in JavaScript panel conditions
  // `supportedExt` is used in Python scripts

  //input: file file
  //output: bool result
  isTextFile(file) {
    let extensions = ["doc", "docx", "odt", "pdf", "rtf", "txt"];
    return extensions.includes(file.extension);
  }

  //input: string filename
  //output: bool result
  supportedExt(filename) {
    let extensions = ["doc", "docx", "odt", "pdf", "rtf", "txt"];
    return extensions.some(ext => filename.endsWith(ext));
  }

  //input: string filename
  //output: bool result
  supportedExtExtract(filename) {
    let extensions = ["csv", "doc", "docx", "eml", "epub", "gif", "htm", "html",
      "jpeg", "jpg", "json", "log", "mp3", "msg", "odt", "ogg",
      "pdf", "png", "pptx", "ps", "psv", "rtf", "tff", "tif",
      "tiff", "tsv", "txt", "wav", "xls", "xlsx"];
    return extensions.some(ext => filename.endsWith(ext));
  }

  //name: detectText
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  //meta.skipTest: GROK-17630
  detectText(col) {
    if (col.type !== 'string')
      return null;
    if (DG.Detector.sampleCategories(col, (s) => s.includes('M  END'), 1, undefined, 0.5))
      return null;
    const likelyTextColumnNames = ['question', 'comment', 'description', 'messige', 'text'];
    const isLikelyColumn = likelyTextColumnNames.some((name) => col.name.toLowerCase().includes(name));
    const ratio = isLikelyColumn ? 0.7 : 0.85;
    const minWordCount = isLikelyColumn ? 8 : 12;
    const delimiters = ` ,`;
    const isText = (s) => {
      let wordStartIdx = -1;
      let maxLength = 0;
      let wordCount = 0;
      let delimCount = 0;

      for (let i = 0; i <= s.length; i++) {
        const c = i == s.length ? ' ' : s[i];
        const isDelim = delimiters.includes(c);
        delimCount += isDelim ? 1 : 0;

        if (wordStartIdx == -1 && /[a-zA-Z]/.test(c))
          wordStartIdx = i;
        else if (isDelim && wordStartIdx > -1) {
          wordCount++;
          maxLength = Math.max(maxLength, i - wordStartIdx);
          wordStartIdx = -1;

          if (maxLength > 25)
            return false;
        }
      }

      return wordCount > minWordCount && delimCount > 3;
    };
    return DG.Detector.sampleCategories(col, isText, 1, undefined, ratio, 1) ? DG.SEMTYPE.TEXT : null;
  }
}
