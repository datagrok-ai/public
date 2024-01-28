import {PatternConfiguration} from '../../model/types';
import {NucleotidePatternSVGRenderer} from './svg-renderer';

// todo: to be deleted after removing legacy code
export function renderNucleotidePattern(patternConfiguration: PatternConfiguration): Element {
  const renderer = new NucleotidePatternSVGRenderer(patternConfiguration);
  const svg = renderer.renderPattern();
  return svg;
}
