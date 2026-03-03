/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {ALPHABET, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

import {_package} from '../package';
import {AnnotationRenderer} from '@datagrok-libraries/bio/src/utils/cell-renderer-annotations';
import {getColumnAnnotations} from './annotations/annotation-manager';
import {AnnotationCategory} from '@datagrok-libraries/bio/src/utils/macromolecule/annotations';


export function addCopyMenuUI(cell: DG.Cell, menu: DG.Menu, seqHelper: ISeqHelper): void {
  const tgtNotationList: string[] = Object.values(NOTATION).filter((v) => v !== NOTATION.CUSTOM);

  menu.group('Copy')
    .items(tgtNotationList, (tgtNotation) => {
      const srcCol = cell.column;
      const srcRowIdx = cell.rowIndex;
      const srcSh = seqHelper.getSeqHandler(srcCol);
      const separator = tgtNotation === NOTATION.SEPARATOR ? (_package.properties.defaultSeparator ?? '-') : undefined;
      const joiner = srcSh.getJoiner({notation: tgtNotation as NOTATION, separator});
      const srcSS = srcSh.getSplitted(srcRowIdx);
      const tgtSeq = joiner(srcSS);

      if (!navigator.clipboard)
        grok.shell.warning('The clipboard functionality requires a secure origin — either HTTPS or localhost');
      else {
        navigator.clipboard.writeText(tgtSeq);
        grok.shell.info(`Value of notation '${tgtNotation}' copied to clipboard`);
      }
    });

  // Annotation context menu items
  const srcCol = cell.column;
  const annotations = getColumnAnnotations(srcCol);
  if (annotations.length > 0) {
    const annotRenderer = new AnnotationRenderer(srcCol);
    if (annotRenderer.hasAnnotations()) {
      // Add annotation info submenu
      const structAnnots = annotations.filter((a) => a.category === AnnotationCategory.Structure);
      if (structAnnots.length > 0) {
        const annotMenu = menu.group('Annotations');

        // Extract region actions (uses per-row data for unaligned sequences)
        for (const annot of structAnnots) {
          if (annot.start && annot.end) {
            annotMenu.item(`Extract ${annot.name} as Column`, () => {
              import('./annotations/annotation-actions').then((m) => {
                try {
                  m.extractAnnotatedRegion(srcCol.dataFrame, srcCol, annot.name, seqHelper);
                } catch (err: any) {
                  grok.shell.error(`Failed to extract region: ${err.message ?? err}`);
                }
              });
            });
          }
        }

        // Filter by liability hits
        const liabAnnots = annotations.filter((a) => a.category === AnnotationCategory.Liability);
        if (liabAnnots.length > 0) {
          annotMenu.item('Filter: Rows with Liability Hits', () => {
            import('./annotations/annotation-actions').then((m) => {
              m.filterByLiabilityHits(srcCol.dataFrame, srcCol);
            });
          });
        }
      }
    }
  }
}
