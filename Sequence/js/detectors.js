class SequencePackageDetectors extends GrokPackage {

    //tags: semTypeDetector
    //input: column col
    //output: string semType
    detectNucleotides(col) {
        if (col.name.startsWith('nuc')) {
            col.semType = 'nucleotides';
            return 'nucleotides';
        }

        return null;
    }

    //input: string s
    //output: string ss
    dupXX(s) {
        return s + '_' + s + ' (sequence)';
    }


    //input: string s
    //output: string ss
    dup(s) {
        return s + '_' + s + ' (sequence)';
    }

}