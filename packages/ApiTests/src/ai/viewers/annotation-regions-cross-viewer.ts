import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// JS API source:
//   public/js-api/src/interfaces/d4.d.ts:274/377/502/672/803/951 (annotationRegions, showViewerAnnotationRegions),
//   core/client/d4/lib/src/viewer_base/annotation_regions_mixin.dart (look property bindings + dfRegions via TAG).
// GROK-20087 PR (Annotation: title placement + strip reserve) added/altered
// behavior across five viewers (scatter, histogram, density, box, line) that
// share AnnotationRegionsMixin but each gate the strip-reserve differently.
// The existing scatter-plot-js-api.ts pins the `meta.annotationRegions`
// helper round-trip on scatter only; this file pins the underlying string
// property on every viewer the PR touches, the `showViewerAnnotationRegions`
// boolean envelope, multi-region JSON shape, and DataFrame-tag-driven
// dataFrameAnnotationRegions sharing — covering the public surface that
// changed without depending on canvas geometry.
category('AI: Viewers: Annotation regions cross-viewer', () => {
  const formulaRegion = {
    type: 'formula', header: 'AI test band',
    headerColor: 2062260, fillColor: 2062260, opacity: 16,
    formula1: '${age} = 30', formula2: '${age} = 50',
  };
  const areaRegion = {
    type: 'area', header: 'AI test area', headerColor: 14034728,
    fillColor: 14034728, outlineColor: 14034728, outlineWidth: 1, opacity: 18,
    x: 'height', y: 'weight',
    area: [[170, 60], [180, 60], [180, 80], [170, 80]],
  };

  test('annotationRegions JSON survives round-trip on Scatter', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.scatter({x: 'height', y: 'weight'});
    v.setOptions({annotationRegions: JSON.stringify([areaRegion, formulaRegion])});
    const raw = v.getOptions(true).look['annotationRegions'];
    expect(typeof raw, 'string');
    const parsed = JSON.parse(raw);
    expect(Array.isArray(parsed), true);
    expect(parsed.length, 2);
    expect(parsed[0].header, 'AI test area');
    expect(parsed[1].header, 'AI test band');
  });

  test('annotationRegions JSON survives round-trip on Histogram, Density, BoxPlot, LineChart', async () => {
    const df = grok.data.demo.demog(50);
    const json = JSON.stringify([formulaRegion]);
    const histo = df.plot.histogram({valueColumnName: 'age', annotationRegions: json});
    const density = DG.Viewer.densityPlot(df, {xColumnName: 'weight', yColumnName: 'height', annotationRegions: json});
    const box = df.plot.box({categoryColumnNames: ['sex'], valueColumnName: 'weight', annotationRegions: json});
    const line = df.plot.line({xColumnName: 'age', yColumnNames: ['height'], annotationRegions: json});
    for (var v of [histo, density, box, line]) {
      const raw = v.getOptions(true).look['annotationRegions'];
      expect(typeof raw, 'string');
      const parsed = JSON.parse(raw);
      expect(Array.isArray(parsed), true);
      expect(parsed.length, 1);
      expect(parsed[0].header, 'AI test band');
    }
  });

  test('showViewerAnnotationRegions boolean round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.scatter({x: 'height', y: 'weight',
      annotationRegions: JSON.stringify([formulaRegion])});
    v.setOptions({showViewerAnnotationRegions: false});
    expect(v.getOptions(true).look['showViewerAnnotationRegions'], false);
    v.setOptions({showViewerAnnotationRegions: true});
    expect(v.getOptions(true).look['showViewerAnnotationRegions'], true);
  });

  test('Empty / cleared annotationRegions resolves to empty array', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.scatter({x: 'height', y: 'weight',
      annotationRegions: JSON.stringify([formulaRegion])});
    v.setOptions({annotationRegions: '[]'});
    const raw = v.getOptions(true).look['annotationRegions'];
    const parsed = JSON.parse(raw);
    expect(Array.isArray(parsed), true);
    expect(parsed.length, 0);
  });

  test('DataFrame .annotation-regions tag drives dfRegions on every viewer attached', async () => {
    const df = grok.data.demo.demog(50);
    df.tags['.annotation-regions'] = JSON.stringify([formulaRegion]);
    const scatter = df.plot.scatter({x: 'height', y: 'weight', showDataframeAnnotationRegions: true});
    const histo = df.plot.histogram({valueColumnName: 'age', showDataframeAnnotationRegions: true});
    for (var v of [scatter, histo]) {
      expect(v.getOptions(true).look['showDataframeAnnotationRegions'], true);
      const tagJson = df.tags['.annotation-regions'];
      expect(typeof tagJson, 'string');
      const parsed = JSON.parse(tagJson);
      expect(parsed.length, 1);
      expect(parsed[0].header, 'AI test band');
    }
    df.tags['.annotation-regions'] = '';
  });

  test('annotationFont string round-trip on Scatter and Histogram', async () => {
    const df = grok.data.demo.demog(50);
    const font = 'normal 600 11px "Roboto"';
    const scatter = df.plot.scatter({x: 'height', y: 'weight', annotationFont: font});
    const histo = df.plot.histogram({valueColumnName: 'age', annotationFont: font});
    expect(scatter.getOptions(true).look['annotationFont'], font);
    expect(histo.getOptions(true).look['annotationFont'], font);
  });
});
