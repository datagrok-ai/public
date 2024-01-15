import {PatternConfiguration} from '../../model/types';
import {SVGRenderer} from './svg-renderer';

export function renderNucleotidePattern(patternConfiguration: PatternConfiguration): Element {
  const renderer = new SVGRenderer(patternConfiguration);
  const svg = renderer.renderNucleotidePattern();
  return svg;
}
