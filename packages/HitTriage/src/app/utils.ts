export const toFormatedDateString = (d: Date): string => {
  return `${d.getFullYear()}/${d.getMonth() + 1}/${d.getDate()}`;
};

export function modifyUrl(key: string, value: string) {
  const title = document.title;
  const url = window.location.href.split('?')[0] + '?' + key + '=' + value;
  if (history.replaceState) {
    const obj = {
      Title: title,
      Url: url,
    };
    history.replaceState(obj, obj.Title, obj.Url);
  }
}

export function hideComponents(styles: CSSStyleDeclaration[]) {
  styles.forEach((style) => {
    style.display = 'none';
  });
}

export function showComponents(styles: CSSStyleDeclaration[]) {
  styles.forEach((style) => {
    style.display = 'flex';
  });
}
