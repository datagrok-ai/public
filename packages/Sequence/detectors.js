class SequencePackageDetectors extends grok.Package {

    //tags: autostart
    autostartTest() {
        console.log('sequence autostarted.');
    }

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

    //tags: semTypeDetector
    //input: column col
    //output: string semType
    detectPdb(col) {
        if (col.name === "pdb") {
            col.semType = 'pdb_id';
            return 'pdb_id';
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