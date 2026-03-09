#name: PythonEnvCondaExtras
#description: Environment with extra conda-only dependencies beyond base
#language: python
#environment: env_conda_extras
#input: string text
#output: string result

import scipy
import sklearn
import numpy as np
from sklearn.feature_extraction.text import CountVectorizer
from scipy.spatial.distance import cosine

# Verify scipy works
vec = np.array([1.0, 2.0, 3.0])
dist = cosine(vec, vec * 2)

# Verify sklearn works
vectorizer = CountVectorizer()
X = vectorizer.fit_transform([text, 'test document'])

result = f"scipy={scipy.__version__},sklearn={sklearn.__version__},dist={dist:.2f},features={X.shape[1]}"
