/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */
//TODO: add unit test for coordDetector
class GisPackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectGisCoord(col) {
    let colSemType = null;
    const colname = col.name.toLowerCase();
    //TODO: add pattern like (51° 28′ 38″ N)
    if (col.type != DG.COLUMN_TYPE.FLOAT) return null;
    if (col.stats.min < -180) return null;
    if (col.stats.max > 180) return null;

    //TODO: change to map or array or pattern search
    if ((colname.includes('lon')) || (colname.includes('lng'))) colSemType = 'gis_longitude';
    if ((colname.includes('lat')) || (colname.includes('ltt')))
      if ((col.stats.min > -90) && (col.stats.max < 90)) colSemType = 'gis_latitude';
    return colSemType;
  }
}
