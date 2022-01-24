export namespace ArrayUtils {

  export function indexesOf<T>(values: T[], filter: (item: T) => boolean): number[] {
    const indexes = [];
    for (let i = 0; i < values.length; i++)
      if (filter(values[i]))
        indexes.push(i);
    return indexes;
  }
}