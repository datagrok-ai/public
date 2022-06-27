/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */
class ArtemSequensePackageDetectors extends DG.Package {
    //tags: semTypeDetector
    //input: column col
    //output: string semType
    detectNucleotides(col) {
        const alphabet = new Set(['A', 'C', 'G', 'T', 'U', '-', 'a', 'c', 'g', 't', 'u', ' ']);
        return DG.Detector.sampleCategories(col, (s) => {
            return s.split('').every((symbol) => alphabet.has(symbol));
        }, 3) ? 'dna_nucleotide' : null;
    }

    //input: string str
    //output: bool result
    isPotentialENAId(str) {
        // returns true, if name is of the form [A-Z]{2}[0-9]{6}
        return /^[A-Z]{2}[0-9]{6}$/.test(str);
    }
}
