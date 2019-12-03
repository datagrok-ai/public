class UsageAnalysisPackageDetectors extends GrokPackage {

    //tags: semTypeDetector
    //input: column col
    //output: string semType
    detectErrorLog(col) {
        if (col.type === TYPE_STRING) {
            if (col.name === 'error_message')
                col.semType = 'ErrorMessage';
            else if (col.name === 'error_stack_trace')
                col.semType = 'ErrorStackTrace';
        }
        return col.semType;
    }
}
