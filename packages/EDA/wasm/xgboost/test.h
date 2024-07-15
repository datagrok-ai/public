#pragma once
#include <iostream>
#include <cstdlib>
#include <cmath>
using namespace std;

//#include <xgboost/c_api.h>
#include "../../../XGBoost/xgboost/include/xgboost/c_api.h"

#include "xgboost-api.h"

int firstTest() {

	// 1. DATA

	const int ROWS = 5, COLS = 3;
	float* features = new float[ROWS * COLS];
	float* labels = new float[ROWS];
	const float MISSING_VALUE = -1;

	for (int i = 0; i < ROWS; ++i) {
		cout << "sample " << i << ": (";

		float sum = 0;

		for (int j = 0; j < COLS; ++j) {
			features[COLS * i + j] = static_cast<float>(rand() % 10);
			cout << features[COLS * i + j] << (j < COLS - 1 ? ", " : "");
			sum += features[COLS * i + j];
		}

		labels[i] = static_cast<float>(static_cast<int>(sum) % 3);

		cout << ")  <-->  " << labels[i] << endl;
	}

	// 2. XGBooster data matrix
	DMatrixHandle dmatrix;
	XGDMatrixCreateFromMat(features, ROWS, COLS, MISSING_VALUE, &dmatrix);

	// 3. XGBooster labels
	XGDMatrixSetFloatInfo(dmatrix, "label", labels, ROWS);

	// 4. Create model
	BoosterHandle booster;
	const int eval_dmats_size = 2;	

	// We assume that training and test data have been loaded into 'train' and 'test'
	DMatrixHandle eval_dmats[eval_dmats_size] = { dmatrix, dmatrix };
	XGBoosterCreate(eval_dmats, eval_dmats_size, &booster);

	// Settings
	XGBoosterSetParam(booster, "booster", "gbtree");
	XGBoosterSetParam(booster, "eta", "0.3");
	XGBoosterSetParam(booster, "max_depth", "3");
	XGBoosterSetParam(booster, "lambda", "1");
	XGBoosterSetParam(booster, "alpha", "0");

	// 5. Train & evaluate
	int num_of_iterations = 20;

	const char* eval_names[eval_dmats_size] = { "train", "test" };
	const char* eval_result = NULL;

	for (int i = 0; i < num_of_iterations; ++i) {
		// Update the model performance for each iteration
		XGBoosterUpdateOneIter(booster, i, dmatrix);

		// Give the statistics for the learner for training & testing dataset in terms of error after each iteration
		XGBoosterEvalOneIter(booster, i, eval_dmats, eval_names, eval_dmats_size, &eval_result);
		printf("%s\n", eval_result);
	}

	// 6. Save model
	const char* model_path = "model.json";
	XGBoosterSaveModel(booster, model_path);

	// 6*. Pack model
	char const packConfig[] = "{\"format\": \"json\"}";
	uint64_t out_len;
	char const * out_dptr = NULL;
	XGBoosterSaveModelToBuffer(booster, packConfig, &out_len, &out_dptr);

	cout << "Packed size: " << out_len << endl;

	char* copiedBytes = new char[out_len];
	for (uint64_t i = 0; i < out_len; ++i)
		copiedBytes[i] = out_dptr[i];

	// 7*. Unpack model
	BoosterHandle unpackedModel;
	XGBoosterCreate(NULL, 0, &unpackedModel);
	cout << "Unpacking: " << XGBoosterLoadModelFromBuffer(unpackedModel, reinterpret_cast<void *>(copiedBytes), out_len) << endl;

	// 7. Load predictor
	BoosterHandle loadedModel;
	XGBoosterCreate(NULL, 0, &loadedModel);
	// load model
	XGBoosterLoadModel(loadedModel, model_path);

	// 8. Predict by LOADED
	cout << "\nBy LOADED:\n";
	char const config[] =
		"{\"training\": false, \"type\": 0, "
		"\"iteration_begin\": 0, \"iteration_end\": 0, \"strict_shape\": false}";
	/* Shape of output prediction */
	uint64_t const* out_shape;
	/* Dimension of output prediction */
	uint64_t out_dim;
	/* Pointer to a thread local contiguous array, assigned in prediction function. */
	float const* out_result = NULL;
	XGBoosterPredictFromDMatrix(loadedModel, dmatrix, config, &out_shape, &out_dim, &out_result);

	for (unsigned int i = 0; i < ROWS; i++) {
		printf("prediction[%i] = %f \n", i, out_result[i]);
	}

	cout << "Predictions: ";
	for (unsigned int i = 0; i < ROWS; i++)
		cout << out_result[i] << (i == ROWS - 1 ? "\n" : ", ");

	// 8. Predict by UNPACKED

	cout << "\nBy UNPACKED:\n";
	/* Shape of output prediction */
	uint64_t const* out_shape1;
	/* Dimension of output prediction */
	uint64_t out_dim1;
	/* Pointer to a thread local contiguous array, assigned in prediction function. */
	float const* out_result1 = NULL;
	XGBoosterPredictFromDMatrix(unpackedModel, dmatrix, config, &out_shape1, &out_dim1, &out_result1);

	for (unsigned int i = 0; i < ROWS; i++) {
		printf("prediction[%i] = %f \n", i, round(out_result1[i]));
	}

	cout << "Predictions: ";
	for (unsigned int i = 0; i < ROWS; i++)
		cout << out_result1[i] << (i == ROWS - 1 ? "\n" : ", ");

	XGDMatrixFree(dmatrix);
	XGBoosterFree(booster);
	XGBoosterFree(loadedModel);
	XGBoosterFree(unpackedModel);

	delete[] features;
	delete[] labels;
	delete[] copiedBytes;

	return 0;
}

int mainTest() {
	int iterations = 20;
	float eta = 0.3f;
	int maxDepth = 6;
	float lambda = 1.0f;
	float alpha = 0.0f;

	const int ROWS = 5, COLS = 3;
	float* features = new float[ROWS * COLS];
	float* labels = new float[ROWS];
	int* predictions = new int[ROWS];
	const float MISSING_VALUE = -1;

	for (int i = 0; i < ROWS; ++i) {
		cout << "sample " << i;

		float sum = 0;

		for (int j = 0; j < COLS; ++j) {
			features[COLS * i + j] = static_cast<float>(rand() % 10);
			//cout << features[COLS * i + j] << (j < COLS - 1 ? ", " : "");
			sum += features[COLS * i + j];
		}

		labels[i] = static_cast<float>(static_cast<int>(sum) % 3);

		cout << "  <-->  " << labels[i] << endl;
	}

	int* modelSizePtr = new int[1];

	int reserved = 1000000;
	int* modelBytes = new int[reserved];

	cout << "\nTraining model: "
		<< train(features, ROWS, COLS, MISSING_VALUE, labels, ROWS,
			iterations, eta, maxDepth, lambda, alpha,
			modelSizePtr, 1,
			modelBytes, reserved)
		<< endl;

	cout << "\nPredicting: "
		<< predict(features, ROWS, COLS, MISSING_VALUE,
			modelBytes, modelSizePtr[0],
			predictions, ROWS) << endl;

	cout << "Predictions: ";
	for (unsigned int i = 0; i < ROWS; i++)
		cout << predictions[i] << (i == ROWS - 1 ? "\n" : ", ");

	delete[] features;
	delete[] labels;
	delete[] predictions;
	delete[] modelSizePtr;
	delete[] modelBytes;

	return 0;
}