import {PatternConfiguration} from '../../model/types';
import {NucleotidePatternSVGRenderer} from './svg-renderer';

export function renderNucleotidePattern(patternConfiguration: PatternConfiguration): Element {
  const renderer = new NucleotidePatternSVGRenderer(patternConfiguration);
  const svg = renderer.renderPattern();
  return svg;
}
