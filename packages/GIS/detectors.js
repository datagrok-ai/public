/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */
// import * as GisTypes from '../src/gis-semtypes';
// import {SEMTYPEGIS, countriesList} from '../src/gis-constants';

const SEMTYPEGIS = {
  LONGITUDE: 'gis-longitude',
  LATIITUDE: 'gis-latitude',
  ALTITUDE: 'gis-altitude',
  GISPOINT: 'gis-point',
  GISAREA: 'gis-gis_area',
  GISCOUNTRY: 'gis-country',
  GISSTATE: 'gis-state',
  GISADDRESS: 'gis-address',
  GISZIPCODE: 'gis-zipcode',
};

//TODO: add unit test for Detectors
class GisPackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectGisCoord(col) {
    let colSemType = null;
    const colName = col.name.toLowerCase();
    //TODO: add pattern like (51° 28′ 38″ N)
    if (col.type != DG.COLUMN_TYPE.FLOAT) return null;
    if (col.stats.min < -180) return null;
    if (col.stats.max > 180) return null;

    //TODO: change to map or array or pattern search
    if ((colName.includes('lon')) || (colName.includes('lng'))) colSemType = SEMTYPEGIS.LONGITUDE;
    if ((colName.includes('lat')) || (colName.includes('ltt')))
      if ((col.stats.min > -90) && (col.stats.max < 90)) colSemType = SEMTYPEGIS.LATIITUDE;
    return colSemType;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectGisCountry(col) {
    const colName = col.name.toLowerCase();
    if ((colName.includes('country') || colName.includes('countri')) && col.type === DG.TYPE.STRING) {
      col.semType = SEMTYPEGIS.GISCOUNTRY;
      return col.semType;
    }
    //TODO: add checking for country name or id
    return null;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectGisZipcode(col) {
    let estCoeff = 0; //coefficient of estimation [0-100] the more value - the more probability

    const colName = col.name.toLowerCase();
    //TODO: check incorrect length of column values (<3? of >9?)
    for (let i = 0; i < col.categories.length; i++) {
      if (((col.categories[i].length > 0) && (col.categories[i].length < 3)) ||
        (col.categories[i].length > 10)) return null;
    }
    // if (DG.Detector.sampleCategories(col, (s) => s.includes('M  END'), 1)) {}
    //US Zipcode format checking
    // const zipUS1 = /\b\d{5}\b/i;
    const zipUS2 = /\b[0-9]{5}-[0-9]{4}\b/i;
    const zipEU1 = /\b\d{4,6}\b/i;
    const zipBRAZ = /\b\d{5}-\d{3}\b/i;
    const zipCAN = /\b[a-z]\d[a-z]\s\d[a-z]\d\b/i;
    const zipJPN = /\b\d{3}-\d{4}\b/i;
    //TODO: add checking for zip codes of other countries Great Britain, AZ, AG, GR, SW, Livan, Islands, NL, PL, PT
    let catNumber = col.categories.length;
    if (catNumber > 100) catNumber = 100; //shorten amount of checking categories to 100
    const caseWeight = 80 / (catNumber+1);
    for (let i = 0; i < catNumber; i++) {
      if ((col.categories[i].length > 0) &&
        (col.categories[i].match(zipUS2) != null)||
        (col.categories[i].match(zipEU1) != null)||
        (col.categories[i].match(zipBRAZ) != null)||
        (col.categories[i].match(zipJPN) != null)||
        (col.categories[i].match(zipCAN) != null)) estCoeff += caseWeight;
      else estCoeff -= caseWeight*2;
    }

    //TODO: should we add checking for "Почтовый индекс"
    if (colName.includes('zip') || colName.includes('code') || colName.includes('post'))
      estCoeff += 40;

    if (estCoeff > 75) {
      col.semType = SEMTYPEGIS.GISZIPCODE;
      return col.semType;
    }

    return null;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectGisAddress(col) {
    let estCoeff = 0; //coefficient of estimation [0-100] the more value - the more probability
    const colName = col.name.toLowerCase();
    if (col.type !== DG.TYPE.STRING) return null;

    if ((colName.includes('address') || colName.includes('street')) && (col.type === DG.TYPE.STRING))
      estCoeff += 40;
    // console.log('col names estCoeff=' + estCoeff);

    //Address patterns
    const Addr1 = /ave/i;
    const Addr2 = /blvd/i;
    const Addr3 = /str./i;
    const Addr4 = /street/i;
    const Addr5 = /square/i;
    const Addr6 = /rd./i;
    const Addr7 = /road/i;
    const Addr8 = /boulevad/i;

    let catNumber = col.categories.length;
    if (catNumber > 100) catNumber = 100;
    const caseWeight = 80 / (catNumber+1);
    for (let i = 0; i < catNumber; i++) {
      if (col.categories[i].length > 6) {
        if (col.categories[i].match(Addr1) != null) estCoeff += caseWeight;
        else if (col.categories[i].match(Addr2) != null) estCoeff += caseWeight;
        else if (col.categories[i].match(Addr3) != null) estCoeff += caseWeight;
        else if (col.categories[i].match(Addr4) != null) estCoeff += caseWeight;
        else if (col.categories[i].match(Addr5) != null) estCoeff += caseWeight;
        else if (col.categories[i].match(Addr6) != null) estCoeff += caseWeight;
        else if (col.categories[i].match(Addr7) != null) estCoeff += caseWeight;
        else if (col.categories[i].match(Addr8) != null) estCoeff += caseWeight;
      } else estCoeff -= caseWeight*2;
    }
    // console.log('worlds match estCoeff=' + estCoeff);

    if (estCoeff > 75) {
      col.semType = SEMTYPEGIS.GISADDRESS;
      return col.semType;
    }
    return null;
  }
//end of GisPackageDetectors class
}
