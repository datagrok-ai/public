class UsageAnalysisPackageDetectors extends DG.Package {

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectUserIds(col) {
    if (col.type === DG.TYPE.STRING && col.name === 'user_id')
      col.semType = 'user_id';
    return col.semType;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectErrorLog(col) {
    if (col.type === DG.TYPE.STRING) {
      if (col.name === 'error_message')
        col.semType = 'ErrorMessage';
      else if (col.name === 'error_stack_trace')
        col.semType = 'ErrorStackTrace';
    }
    return col.semType;
  }
}
