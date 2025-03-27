export enum StringMetricsNames {
    Levenshtein = 'Levenshtein',
    JaroWinkler = 'Jaro-Winkler',
    Manhattan = 'Manhattan String',
    Onehot = 'One-Hot',
  }

export enum VectorMetricsNames {
    Euclidean = 'Euclidean',
    Manhattan = 'Manhattan',
    Cosine = 'Vector Cosine',
  }

export enum BitArrayMetricsNames {
    Tanimoto = 'Tanimoto',
    Dice = 'Dice',
    Asymmetric = 'Asymmetric',
    BraunBlanquet = 'Braun-Blanquet',
    Cosine = 'Cosine',
    Kulczynski = 'Kulczynski',
    McConnaughey = 'Mc-Connaughey',
    RogotGoldberg = 'Rogot-Goldberg',
    Russel = 'Russel',
    Sokal = 'Sokal',
    Hamming = 'Hamming',
    Euclidean = 'Euclidean BitArray',
  }

export enum IntArrayMetricsNames {
  TanimotoIntArray = 'TanimotoIntArray',
}

export enum DistanceMetricsSubjects {
    Vector = 'Vector',
    String = 'String',
    BitArray = 'BitArray',
    MacroMolecule = 'MacroMolecule',
    Number = 'Number',
    IntArray = 'IntArray',
    NumberArray = 'NumberArray',
  }

export enum NumberMetricsNames {
  Difference = 'Difference',
}

export enum NumberArrayMetricsNames {
  CommonItems = 'Common Items',
}
