/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */
class StanislavSequencePackageDetectors extends DG.Package {
    //tags: semTypeDetector
    //input: column col
    //output: string semType
    detectNucleotides(col) {
        console.log(col);
        if (col.type == "string" && DG.Detector.sampleCategories(col, (s) => {
            console.log(s.length === 7);
            console.log(new RegExp("[A-Z]*").test(s));
            return s.length === 7 && new RegExp("[A-Z]*").test(s);
        }, 1)) {
            col.semType = 'dna_nucleotide';
            return col.semType;
        }
        return null;
    }
}