import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {
  SeqAnnotation, AnnotationCategory, AnnotationVisualType,
} from '@datagrok-libraries/bio/src/utils/macromolecule/annotations';
import {NumberingScheme} from '@datagrok-libraries/bio/src/utils/macromolecule/numbering-schemes';
import {setColumnAnnotations, getColumnAnnotations} from './annotation-manager';
import {_package} from '../../package';

export function showNumberingSchemeDialog(): void {
  const df = grok.shell.tv?.dataFrame;
  if (!df) {
    grok.shell.warning('No table open');
    return;
  }

  const seqCols = df.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE);
  if (seqCols.length === 0) {
    grok.shell.warning('No macromolecule columns found');
    return;
  }

  const schemeChoices = Object.values(NumberingScheme);

  const tableInput = ui.input.table('Table', {value: df});
  const seqInput = ui.input.column('Sequence', {
    table: df, value: seqCols[0],
    filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE,
  });
  const schemeInput = ui.input.choice('Scheme', {value: NumberingScheme.IMGT, items: schemeChoices});
  const populateRegions = ui.input.bool('Populate FR/CDR regions', {value: true});
  const openVdRegions = ui.input.bool('Open VD Regions viewer', {value: true});

  const dialog = ui.dialog({title: 'Apply Antibody Numbering'})
    .add(ui.inputs([tableInput, seqInput, schemeInput, populateRegions, openVdRegions]))
    .add(ui.div([
      ui.divText('Uses AntPack (Python) for numbering. No external dependencies required.', {style: {fontSize: '11px', color: '#888', marginTop: '8px'}}),
    ]))
    .onOK(async () => {
      try {
        const seqCol = seqInput.value!;
        const schemeName = schemeInput.value!;

        grok.shell.info(`Running AntPack numbering (${schemeName})...`);

        // Call Python script
        const result: DG.DataFrame = await grok.functions.call(
          'Bio:NumberAntibodySequences', {
            df: df,
            seqCol: seqCol,
            scheme: schemeName.toLowerCase(),
          },
        );

        if (!result || result.rowCount === 0) {
          grok.shell.warning('No numbering results returned');
          return;
        }

        // Apply results to the sequence column
        const posNamesCol = result.getCol('position_names');
        const chainTypeCol = result.getCol('chain_type');
        const annotJsonCol = result.getCol('annotations_json');

        // Use the first non-empty result for column-level position names
        let positionNames = '';
        let chainType = '';
        let annotationsJson = '[]';

        for (let i = 0; i < result.rowCount; i++) {
          const pn = posNamesCol.get(i);
          if (pn && pn.length > 0) {
            positionNames = pn;
            chainType = chainTypeCol.get(i) ?? '';
            annotationsJson = annotJsonCol.get(i) ?? '[]';
            break;
          }
        }

        if (!positionNames) {
          grok.shell.warning('AntPack could not number the sequences. Check that they are valid antibody variable region sequences.');
          return;
        }

        // Set position names
        seqCol.setTag(bioTAGS.positionNames, positionNames);
        seqCol.setTag(bioTAGS.numberingScheme, schemeName);

        // Set region annotations
        if (populateRegions.value) {
          try {
            const regionAnnotations: SeqAnnotation[] = JSON.parse(annotationsJson);
            const existing = getColumnAnnotations(seqCol).filter((a) => a.category !== AnnotationCategory.Structure);
            setColumnAnnotations(seqCol, [...existing, ...regionAnnotations]);
          } catch (err) {
            console.warn('Failed to parse region annotations from AntPack result:', err);
          }
        }

        df.fireValuesChanged();
        grok.shell.info(`Numbering applied: ${schemeName}, chain type: ${chainType}`);

        // Open VD Regions viewer
        if (openVdRegions.value && grok.shell.tv) {
          try {
            await grok.shell.tv.dataFrame.plot.fromType('VdRegions', {});
          } catch (err) {
            console.warn('Could not open VD Regions viewer:', err);
          }
        }
      } catch (err: any) {
        grok.shell.error(`Numbering failed: ${err.message ?? err}`);
        console.error(err);
      }
    });

  dialog.show();
}
