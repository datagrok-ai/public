import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Aminoacids} from './src/aminoacids';
import {Nucleotides} from './src/nucleotides';

/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */
class BioPackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectNucleotidesSeq(col) {
    const alphabet = Nucleotides.Names + {'-': 'gap'};
    return DG.Detector.sampleCategories(col, (s) => {
      return s.split('').every((n) => n in alphabet);
    }) ? Nucleotides.SemType : null;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectAminoacidsSeq(col) {
    const alphabet = Aminoacids.Names + {'-': 'gap'};
    return DG.Detector.sampleCategories(col, (s) => {
      return s.split('').every((aa) => aa in alphabet);
    }) ? Aminoacids.SemType : null;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectMultipleAlignmentNucleotidesSeq(col) {
    const alphabet = Nucleotides.Names + {'-': 'gap'};
    const len = col.get(0).length;
    return DG.Detector.sampleCategories(col, (s) => {
      return s.length == len && s.split('').every((n) => n in alphabet);
    }) ? Nucleotides.SemTypeMultipleAlignment : null;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectMultipleAlignmentAminoacidsSeq(col) {
    const alphabet = Aminoacids.Names + {'-': 'gap'};
    const len = col.get(0).length;
    return DG.Detector.sampleCategories(col, (s) => {
      return s.length == len && s.split('').every((aa) => aa in alphabet);
    }) ? Aminoacids.SemTypeMultipleAlignment : null;
  }
}
