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
  LONGITUDE: 'Longitude', //'gis-longitude',
  LATIITUDE: 'Latitude', //'gis-latitude',
  ALTITUDE: 'gis-altitude',
  GISPOINT: 'gis-point',
  GISAREA: 'gis-area',
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
  //description: detect longitude/latitude coordinates
  detectGisCoord(col) {
    //TODO: add pattern like (51° 28′ 38″ N) (? if we need it)
    if (col.type !== DG.COLUMN_TYPE.FLOAT)
      return null;

    let estCoeff = 0;
    let colSemType = null;
    const colName = col.name.toLowerCase();
    //the check below leads to miss defining of coordinates type in case of at least one error in data>>
    // if ((col.stats.min < -180) || (col.stats.max > 180)) return null;

    // const begTime = Date.now();

    if ((colName.includes('lon')) || (colName.includes('lng')) || (colName === 'x')) {
      colSemType = SEMTYPEGIS.LONGITUDE;
      estCoeff += 40;
    } else if ((colName.includes('lat')) || (colName.includes('ltt')) || (colName === 'y')) {
      colSemType = SEMTYPEGIS.LATIITUDE;
      estCoeff += 40;
    }

    if (colSemType === null)
      return null;

    const samplesNum = Math.min(col.length, 50);
    const step = Math.round(col.length / (samplesNum + 1));
    const caseWeight = 80 / (samplesNum + 1);
    const columnArr = col.getRawData();
    for (let i = 0; i < columnArr.length; i += step) {
      if (colSemType === SEMTYPEGIS.LATIITUDE) {
        if ((columnArr[i] > -90) && (columnArr[i] < 90))
          estCoeff += caseWeight;
        else
          estCoeff -= caseWeight * 2;
      }
      if (colSemType === SEMTYPEGIS.LONGITUDE) {
        if ((columnArr[i] > -180) && (columnArr[i] < 180))
          estCoeff += caseWeight;
        else
          estCoeff -= caseWeight * 2;
      }
    }

    // const endTime = Date.now();
    // console.log('coord detector: ' + (endTime - begTime));
    if (estCoeff > 75)
      return colSemType;

    return null;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  //description: detector of countries
  detectGisCountry(col) {
    if (col.type !== DG.COLUMN_TYPE.STRING)
      return null;
    const colName = col.name.toLowerCase();
    if ((colName.includes('country') || colName.includes('countri'))) {
      col.semType = SEMTYPEGIS.GISCOUNTRY;
      return col.semType;
    }
    //TODO: add checking for country name or id (?)
    return null;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  //description: detector of ZIP codes for different countries
  detectGisZipcode(col) {
    if ((col.type !== DG.COLUMN_TYPE.STRING)) // || (col.type !== DG.COLUMN_TYPE.INT))
      return null;
    let estCoeff = 0; //coefficient of estimation [0 - >100] the more value - the more probability

    const zipcodePatterns = new Map([
      ['zipUS1EU1', /[\b\d{4,6}\b|\b\d{9}\b]/i],
      ['zipUS1nodash', /\b\d{9}\b/i],
      ['zipUS2BRZ', /\b\d{5}-\d{3,4}\b/i],
      ['zipJPN', /\b\d{3}-\d{3,4}\b/i],
      ['zipCAN', /\b[a-z]\d[a-z]\s\d[a-z]\d\b/i],
    ]);

    const zipcodeLengthMap = new Map([
      ['zipUS1EU1', 7],
      ['zipUS1nodash', 10],
      ['zipUS2BRZ', 11],
      ['zipJPN', 9],
      ['zipCAN', 8],
    ]);

    //TODO: add checking for zip codes of other countries Great Britain, AZ, AG, GR, SW, Livan, Islands, NL, PL, PT ?

    const samplesNum = Math.min(col.categories.length, 30);
    const caseWeight = 75 / (samplesNum + 1);
    const step = Math.round(col.categories.length / (samplesNum + 1));
    for (let i = 0; i < col.categories.length; i += step) {
      //checking for incorrect length of column values (<3 or >10)
      if ((col.categories[i].length < 3) || (col.categories[i].length > 10)) {
        estCoeff -= caseWeight * 2;
        continue;
      }

      const oldValue = estCoeff;
      for (const [key, value] of zipcodePatterns) {
        if ((col.categories[i].match(value) !== null) && (col.categories[i].length < zipcodeLengthMap.get(key))) {
          estCoeff += caseWeight;
          break;
        }
      }
      if (estCoeff === oldValue) estCoeff -= caseWeight * 2;
    }

    //TODO: should we add checking for "Почтовый индекс"?
    const colName = col.name.toLowerCase();
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
  //description: detector address
  detectGisAddress(col) {
    if (col.type !== DG.TYPE.STRING)
      return null;

    let estCoeff = 0; //coefficient of estimation [0-100] the more value - the more probability
    const colName = col.name.toLowerCase();

    if (colName.includes('address') || colName.includes('street') || colName.includes('location'))
      estCoeff += 40;

    //Address patterns
    const addressPatterns = new Map([
      ['street', /[str.|street]/i],
      ['boulevard', /[blvd|boulevad]/i],
      ['square', /square/i],
      ['road', /[rd.|road]/i],
    ]);

    const samplesNum = Math.min(col.categories.length, 50);
    const caseWeight = 75 / (samplesNum + 1);
    for (let i = 0; i < samplesNum; i++) {
      if (col.categories[i].length < 6) {
        estCoeff -= caseWeight * 2;
        continue;
      }

      const oldValue = estCoeff;
      for (const value of addressPatterns.values()) {
        if ((col.categories[i].match(value) !== null)) {
          estCoeff += caseWeight;
          break;
        }
      }
      if (estCoeff === oldValue) estCoeff -= caseWeight * 2;
    }

    if (estCoeff > 75) {
      col.semType = SEMTYPEGIS.GISADDRESS;
      return col.semType;
    }
    return null;
  }
//end of GisPackageDetectors class
}
