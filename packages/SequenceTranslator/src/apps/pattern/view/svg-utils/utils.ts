import {NUCLEOTIDES} from '../../../common/model/const';
import {AXOLABS_STYLE_MAP as styleMap} from '../../../common/data-loader/json-loader';
import {SVGElementFactory} from './svg-element-factory';
import {isOverhangNucleotide} from '../../model/helpers';
import { STRAND, STRANDS, TERMINUS, TERMINI } from '../../model/const';
import {PatternConfiguration, StrandType, TerminalType} from '../../model/types';
import {STRAND_END, STRAND_ENDS, LUMINANCE_COEFFICIENTS, TEXT_COLOR, SVG_CIRCLE_SIZES, SVG_TEXT_FONT_SIZES, SVG_ELEMENT_COLORS, STRAND_END_LABEL_TEXT, NUMERIC_LABEL_POSITION_OFFSET, DEFAULT_FONT_FAMILY, Y_POSITIONS_FOR_STRAND_ELEMENTS} from './const';

export function computeCommentYPosition(isAntisenseStrandActive: boolean): number {
  return (isAntisenseStrandActive ? 11 : 8.5) * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS;
}

export function computeLegendCircleYPosition(isAntisenseStrandActive: boolean): number {
  return (isAntisenseStrandActive ? 9.5 : 6) * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS;
}

export function computeLegendTextYPosition(isAntisenseStrandActive: boolean): number {
  const position = isAntisenseStrandActive ? 10 * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS : Y_POSITIONS_FOR_STRAND_ELEMENTS[STRAND.ANTISENSE].NUCLEOBASE_CIRCLE;
  return position - 3;
}

export function computeTotalSVGHeight(isAntisenseStrandActive: boolean): number {
  return (isAntisenseStrandActive ? 11 : 9) * SVG_CIRCLE_SIZES.NUCLEOBASE_RADIUS;
}

export function isSingleDigitNumber(n: number): boolean {
  return n >= 0 && n < 10;
}

export function countOverhangNucleotidesAtStrandEnd(modifications: string[]): number {
  const lastIdx = modifications.length - 1;
  let count = 0;
  while (count <= lastIdx && isOverhangNucleotide(modifications[count])) {
    count++;
  }
  return count === lastIdx + 1 ? 0 : count;
}

export function computeTextWidthInPixels(text: string, fontSize: number, fontFamily: string): number {
  const canvas = document.createElement('canvas');
  const context = canvas.getContext('2d');
  if (context) {
    context.font = `${fontSize}px ${fontFamily}`;
    const metrics = context.measureText(text);
    return 2 * metrics.width;
  }
  return 0;
}

export function getNucleobaseLabelForCircle(nucleobase: string): string {
  const criterion = !isOverhangNucleotide(nucleobase) && NUCLEOTIDES.includes(nucleobase);

  return criterion ? nucleobase : '';
}

export function computeTextColorForNucleobaseLabel(nucleobase: string): string {
  const nucleobaseColor = styleMap[nucleobase]?.color || '';

  const rgbValues = nucleobaseColor.match(/\d+/g)?.map(Number);
  if (!rgbValues || rgbValues.length < 3) {
    return TEXT_COLOR.LIGHT;
  }

  const [r, g, b] = rgbValues;
  const luminance = r * LUMINANCE_COEFFICIENTS.RED + g * LUMINANCE_COEFFICIENTS.GREEN + b * LUMINANCE_COEFFICIENTS.BLUE;
  return luminance > LUMINANCE_COEFFICIENTS.THRESHOLD ? TEXT_COLOR.DARK : TEXT_COLOR.LIGHT;
}

export function getNucleobaseColorFromStyleMap(nucleobase: string): string {
  return styleMap[nucleobase].color;
}

export function computeMaxWidthForStrandEnd(strandEnd: STRAND_END): number {
  return Math.max(
    ...STRANDS.map(strand =>
      computeTextWidthInPixels(
        STRAND_END_LABEL_TEXT[strandEnd][strand],
        SVG_TEXT_FONT_SIZES.NUCLEOBASE,
        DEFAULT_FONT_FAMILY
      )
    )
  );
}
