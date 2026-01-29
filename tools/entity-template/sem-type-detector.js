
  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detect#{PACKAGE_DETECTORS_NAME}(col) {
    if (col.name.startsWith('#{NAME_PREFIX}')) {
      col.semType = '#{NAME}';
      return col.semType;
    }
    return null;
  }
