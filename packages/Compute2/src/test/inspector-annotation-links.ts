import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {isLinkShown} from '../components/Inspector/annotation-links';

category('Inspector: annotation links', () => {
  test('Missing-value checks are never shown, other annotation links only with the toggle', async () => {
    expect(isLinkShown(undefined, false), true);
    expect(isLinkShown('rule', false), false);
    expect(isLinkShown('check', false), false);
    expect(isLinkShown('required', false), false);
    expect(isLinkShown('rule', true), true);
    expect(isLinkShown('check', true), true);
    expect(isLinkShown('required', true), false);
  });
});
