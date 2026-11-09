// Enter and Space press a focused element like a button; the key stops there, so a surrounding widget
// (e.g. a tree row) doesn't also act on it
export const onPress = (action: (e: KeyboardEvent) => void) => (e: KeyboardEvent) => {
  if (e.key !== 'Enter' && e.key !== ' ')
    return;
  e.preventDefault();
  e.stopPropagation();
  action(e);
};
