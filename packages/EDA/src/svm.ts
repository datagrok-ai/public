/* Support vector machine (SVM) tools.

   Training & predicting are provided by wasm-computations.

   Least square support vector machine (LS-SVM) is implemented:
     [1] Suykens, J., Vandewalle, J. "Least Squares Support Vector Machine Classifiers",
	       Neural Processing Letters 9, 293-300 (1999). https://doi.org/10.1023/A:1018628609742 
*/

// 1. CONSTANTS

// kernel types
export const LINEAR = 0;
export const POLYNOMIAL = 1;
export const RBF = 2;
export const SIGMOID = 3;

// limitations CURRENTLY, IT'S NOT IMPLEMENTED!
//const SAMPLES_MAX = 10000;
//const FEATURES_MAX = 1000000;

// output-related
const CONFUSION_MATR_SIZE = 4;
const NORMALIZED_DATA_INDEX = 0;
const MEANS_INDEX = 1;
const STD_DEVS_INDEX = 2;
const MODEL_PARAMS_INDEX = 3;
const MODEL_WEIGHTS_INDEX = 4;
const PREDICTED_LABELS_INDEX = 5;
const CORRECTNESS_INDEX = 6;
const CONFUSION_MATRIX_INDEX = 7;
const TRUE_POSITIVE_INDEX = 0;
const FALSE_NEGATIVE_INDEX = 1;
const FALSE_POSITIVE_INDEX = 2;
const TRUE_NEGATIVE_INDEX = 3;

// kernel parameters indeces
const RBF_SIGMA_INDEX = 0;
const POLYNOMIAL_C_INDEX = 0;
const POLYNOMIAL_D_INDEX = 1;
const SIGMOID_KAPPA_INDEX = 0;
const SIGMOID_THETA_INDEX = 1;

// hyperparameters limits
const GAMMA_INFIMUM_LIMIT = 0;
const RBF_SIGMA_INFIMUM_LIMIT = 0;
const POLYNOMIAL_C_INFIMUM_LIMIT = 0;
const POLYNOMIAL_D_INFIMUM_LIMIT = 0;

// error messages
const WRONG_GAMMA_MESSAGE = 'gamma must be strictly positive.';
const WRONG_RBF_SIGMA_MESSAGE = 'sigma must be strictly positive.';
const WRONG_POLYNOMIAL_C_MESSAGE = 'c must be strictly positive.';
const WRONG_POLYNOMIAL_D_MESSAGE = 'd must be strictly positive.';
const WRONG_KERNEL_MESSAGE = 'incorrect kernel.';
const TOO_BIG_DATA_MESSAGE = 'Training data is too big.'

// names
const NORMALIZED = 'Normalized';
const CHARACTERISTICS = 'mean, deviation';
const PARAMETERS = 'Parameters';
const WEIGHTS = 'Weights';
const LABELS = 'Labels';
const PREDICTED = 'predicted';
const CORRECTNESS = 'correctness';
const CONFUSION_MATRIX_NAME = 'Confusion matrix';
const MEAN = 'mean';
const STD_DEV = 'std dev';
const MODEL_PARAMS_NAME = 'alpha';
const MODEL_WEIGHTS_NAME = 'weight';
const GAMMA = 'gamma';
const KERNEL = 'kernel';
const KERNEL_PARAMS = 'kernel params'; 
const KERNEL_PARAM_1 = 'kernel param 1';
const KERNEL_PARAM_2 = 'kernel param 2';
const FEATURES_COUNT_NAME = 'features count';
const TRAIN_SAMPLES_COUNT_NAME = 'train samples count';
const TRAIN_ERROR = 'Train error,%';
const MODEL_INFO = 'Model info';
const MODEL_INFO_FULL = 'Model full info';
const KERNEL_TYPE_TO_NAME_MAP = ['linear', 'polynomial', 'RBF', 'sigmoid'];
const POSITIVE_NAME = 'positive (P)';
const NEGATIVE_NAME = 'negative (N)';
const PREDICTED_POSITIVE_NAME = 'predicted positive (PP)';
const PREDICTED_NEGATIVE_NAME = 'predicted negative (PN)';
const SENSITIVITY = 'Sensitivity';
const SPECIFICITY = 'Specificity';
const BALANCED_ACCURACY = 'Balanced accuracy';
const POSITIVE_PREDICTIVE_VALUE = 'Positive predicitve value';
const NEGATIVE_PREDICTIVE_VALUE = 'Negative predicitve value';
const ML_REPORT = 'Model report';
const ML_REPORT_PREDICTED_LABELS = 'Predicted labels';
const ML_REPORT_TRAIN_LABELS = 'Train labels';
const ML_REPORT_CORRECTNESS = 'Prediction correctness';
const PREDICTION = 'prediction';

// Pack/unpack constants
const BYTES = 4;
const INTS_COUNT = 3;
const KER_PARAMS_COUNT = 2;
const MODEL_KERNEL_INDEX = 0;
const SAMPLES_COUNT_INDEX = 1;
const FEATURES_COUNT_INDEX = 2;

// misc
const INIT_VALUE = 0; // any number can be used
const LS_SVM_ADD_CONST = 1; // see [1] for more details
