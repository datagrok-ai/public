// newFAE.cpp

#include <cmath>
#include <vector>
#include <string>
using std::vector;
using std::pow;

#include <emscripten.h>

extern "C" {
    int solveFAEnew(float t0, float t1, float h, 
                    int timesCount, int varsCount,
                    float * resultCols, int resultRowCount, int resultColCount) noexcept;
}

double sqrt(double x)
{
    if(x >= 0.0)
      return std::sqrt(x);
    return 0.0;
}

void Solve(
  std::vector<std::vector<double>> *result,
  std::vector<double> *times,
  std::vector<double> *_parameters,
  double _tolerance)
{
  std::vector<double>::iterator _itB_times = times->begin();
  std::vector<double>::iterator _itE_times = times->end();
  double _time = *_itB_times;
  double _askedTime = *_itB_times;
  double step = 0;

  //parameters
  double qin = (*_parameters)[0];
  double percentO2saturation = (*_parameters)[1];
  double yO2in = (*_parameters)[2];
  double pKa2MEA = (*_parameters)[3];
  double H = (*_parameters)[4];
  double T = (*_parameters)[5];
  double R = (*_parameters)[6];
  double P = (*_parameters)[7];
  //constants
  double VLinitial = 6.6;
  double Vtotalvessel = 10;
  double AgitatorSpeed = 200;
  double AgitatorDiameter = 2.36;
  double AgitatorPowerNumber = 2.1;
  double pH = 7.4;
  double k1red = 0.000934 * 60;
  double k1ox = 0.00018 * 60;
  double k2Fd = 0.0225 * 60;
  double k2Fa = 1840000 * 60;
  double k2Kd = 0.000673 * 60;
  double k2Ka = 2000000 * 60;
  double k3FKa = 3020000 * 60;
  double k3FKd = 0.000198 * 60;
  double k4ox = 0.00018 * 60;
  double k4red = 0.000934 * 60;
  double kthiolox = 0.04;
  double krcyst = 0;
  //allocation for selected variables
  std::vector<double>* MEAthiol = new std::vector<double>();
  std::vector<double>* VL = new std::vector<double>();
  std::vector<double>* klasurface = new std::vector<double>();
  std::vector<double>* CO2 = new std::vector<double>();
  std::vector<double>* yO2P = new std::vector<double>();
  std::vector<double>* FFox = new std::vector<double>();
  std::vector<double>* KKox = new std::vector<double>();
  std::vector<double>* FFred = new std::vector<double>();
  std::vector<double>* KKred = new std::vector<double>();
  std::vector<double>* Ffree = new std::vector<double>();
  std::vector<double>* Kfree = new std::vector<double>();
  std::vector<double>* FKred = new std::vector<double>();
  std::vector<double>* FKox = new std::vector<double>();
  std::vector<double>* Cystamine = new std::vector<double>();
  std::vector<double>* tioBalance = new std::vector<double>();
  //initial expressions and starting values
  double MEAthiol_t = 34;
  double MEAthiolate_t = MEAthiol_t * pow(10, (pH - pKa2MEA));
  double VL_t = 6.6;
  double klasurface_t = 0.33 * pow(VL_t, (-0.281)) * pow((AgitatorPowerNumber * pow((AgitatorSpeed / 60), 3) * pow(AgitatorDiameter, 5) * pow(2.54, 3) / (pow(39.37, 2) * 1000) / (VL_t / 1000)), 0.36);
  double CO2_t = 0.22;
  double yO2P_t = 0.209;
  double qout_t = qin - klasurface_t * (yO2P_t * H - CO2_t) * VL_t * R * T / (P * 1000);
  double OTR_t = klasurface_t * (yO2P_t * H - CO2_t);
  double Vg_t = Vtotalvessel - VL_t;
  double CO2in_t = percentO2saturation * 7.17 / (32 * 100);
  double Fin_t = _time < 120 ? (0) : (0.025);
  double Fpermeate_t = _time < 120 ? (0.025) : (Fin_t);
  double Vres_t = VLinitial / VL_t;
  double FFox_t = 0.268;
  double KKox_t = 0.268;
  double FFred_t = 0;
  double KKred_t = 0;
  double Ffree_t = 0;
  double Kfree_t = 0;
  double FKred_t = 0;
  double FKox_t = 0;
  double Cystamine_t = 0;
  double FFox_to_FFred_t = k1red * FFox_t * Vres_t * pow(MEAthiolate_t, 2) * pow(Vres_t, 2);
  double FFred_to_FFox_t = k1ox * FFred_t * Vres_t;
  double KKox_to_KKred_t = k1red * KKox_t * Vres_t * pow(MEAthiolate_t, 2) * pow(Vres_t, 2);
  double KKred_to_KKox_t = k1ox * KKred_t * Vres_t;
  double FFred_to_Ffree_t = k2Fd * FFred_t * Vres_t;
  double Ffree_to_FFred_t = k2Fa * pow(Ffree_t, 2) * pow(Vres_t, 2) * pow(MEAthiolate_t, 2) * pow(Vres_t, 2);
  double KKred_to_Kfree_t = k2Kd * KKred_t * Vres_t;
  double Kfree_to_KKred_t = k2Ka * pow(Kfree_t, 2) * pow(Vres_t, 2) * pow(MEAthiolate_t, 2) * pow(Vres_t, 2);
  double free_to_FKred_t = k3FKa * Ffree_t * Vres_t * Kfree_t * Vres_t;
  double FKred_to_free_t = k3FKd * FKred_t * Vres_t;
  double FKred_to_FKox_t = k4ox * FKred_t * Vres_t * pow(Cystamine_t, 2) * pow(Vres_t, 2);
  double FKox_to_FKred_t = k4red * FKox_t * Vres_t * pow(MEAthiolate_t, 2) * pow(Vres_t, 2);
  double o2balance_t = klasurface_t * (yO2P_t * H - CO2_t);
  double tioBalance_t = kthiolox * pow(MEAthiol_t, 2) * pow(Vres_t, 2) * sqrt(CO2_t) * sqrt(Vres_t);
  //selected variables to memory
  MEAthiol->push_back(MEAthiol_t);
  VL->push_back(VL_t);
  klasurface->push_back(klasurface_t);
  CO2->push_back(CO2_t);
  yO2P->push_back(yO2P_t);
  FFox->push_back(FFox_t);
  KKox->push_back(KKox_t);
  FFred->push_back(FFred_t);
  KKred->push_back(KKred_t);
  Ffree->push_back(Ffree_t);
  Kfree->push_back(Kfree_t);
  FKred->push_back(FKred_t);
  FKox->push_back(FKox_t);
  Cystamine->push_back(Cystamine_t);
  tioBalance->push_back(tioBalance_t);
  ++_itB_times;

  //ODE Solver: explicit adaptive Cash-Karp 4 order method
  step = 0.000001;
  bool _increase = true;
  _time = 0 + step;
  double _timePrev = 0;
  _askedTime = *_itB_times;
  double _infinitezimal = 1e-20;
  double FFox_scale = 0;
  double KKox_scale = 0;
  double FFred_scale = 0;
  double KKred_scale = 0;
  double Ffree_scale = 0;
  double Kfree_scale = 0;
  double FKred_scale = 0;
  double FKox_scale = 0;
  double MEAthiol_scale = 0;
  double CO2_scale = 0;
  double yO2P_scale = 0;
  double Cystamine_scale = 0;
  double VL_scale = 0;
  double FFox_4it = 0;
  double KKox_4it = 0;
  double FFred_4it = 0;
  double KKred_4it = 0;
  double Ffree_4it = 0;
  double Kfree_4it = 0;
  double FKred_4it = 0;
  double FKox_4it = 0;
  double MEAthiol_4it = 0;
  double CO2_4it = 0;
  double yO2P_4it = 0;
  double Cystamine_4it = 0;
  double VL_4it = 0;
  double FFox_5it = 0;
  double KKox_5it = 0;
  double FFred_5it = 0;
  double KKred_5it = 0;
  double Ffree_5it = 0;
  double Kfree_5it = 0;
  double FKred_5it = 0;
  double FKox_5it = 0;
  double MEAthiol_5it = 0;
  double CO2_5it = 0;
  double yO2P_5it = 0;
  double Cystamine_5it = 0;
  double VL_5it = 0;
  double FFox_k1 = 0;
  double KKox_k1 = 0;
  double FFred_k1 = 0;
  double KKred_k1 = 0;
  double Ffree_k1 = 0;
  double Kfree_k1 = 0;
  double FKred_k1 = 0;
  double FKox_k1 = 0;
  double MEAthiol_k1 = 0;
  double CO2_k1 = 0;
  double yO2P_k1 = 0;
  double Cystamine_k1 = 0;
  double VL_k1 = 0;
  double FFox_k2 = 0;
  double KKox_k2 = 0;
  double FFred_k2 = 0;
  double KKred_k2 = 0;
  double Ffree_k2 = 0;
  double Kfree_k2 = 0;
  double FKred_k2 = 0;
  double FKox_k2 = 0;
  double MEAthiol_k2 = 0;
  double CO2_k2 = 0;
  double yO2P_k2 = 0;
  double Cystamine_k2 = 0;
  double VL_k2 = 0;
  double FFox_k3 = 0;
  double KKox_k3 = 0;
  double FFred_k3 = 0;
  double KKred_k3 = 0;
  double Ffree_k3 = 0;
  double Kfree_k3 = 0;
  double FKred_k3 = 0;
  double FKox_k3 = 0;
  double MEAthiol_k3 = 0;
  double CO2_k3 = 0;
  double yO2P_k3 = 0;
  double Cystamine_k3 = 0;
  double VL_k3 = 0;
  double FFox_k4 = 0;
  double KKox_k4 = 0;
  double FFred_k4 = 0;
  double KKred_k4 = 0;
  double Ffree_k4 = 0;
  double Kfree_k4 = 0;
  double FKred_k4 = 0;
  double FKox_k4 = 0;
  double MEAthiol_k4 = 0;
  double CO2_k4 = 0;
  double yO2P_k4 = 0;
  double Cystamine_k4 = 0;
  double VL_k4 = 0;
  double FFox_k5 = 0;
  double KKox_k5 = 0;
  double FFred_k5 = 0;
  double KKred_k5 = 0;
  double Ffree_k5 = 0;
  double Kfree_k5 = 0;
  double FKred_k5 = 0;
  double FKox_k5 = 0;
  double MEAthiol_k5 = 0;
  double CO2_k5 = 0;
  double yO2P_k5 = 0;
  double Cystamine_k5 = 0;
  double VL_k5 = 0;
  double FFox_k6 = 0;
  double KKox_k6 = 0;
  double FFred_k6 = 0;
  double KKred_k6 = 0;
  double Ffree_k6 = 0;
  double Kfree_k6 = 0;
  double FKred_k6 = 0;
  double FKox_k6 = 0;
  double MEAthiol_k6 = 0;
  double CO2_k6 = 0;
  double yO2P_k6 = 0;
  double Cystamine_k6 = 0;
  double VL_k6 = 0;
  double FFox_add = 0;
  double KKox_add = 0;
  double FFred_add = 0;
  double KKred_add = 0;
  double Ffree_add = 0;
  double Kfree_add = 0;
  double FKred_add = 0;
  double FKox_add = 0;
  double MEAthiol_add = 0;
  double CO2_add = 0;
  double yO2P_add = 0;
  double Cystamine_add = 0;
  double VL_add = 0;
  double MEAthiolate_add = 0;
  double klasurface_add = 0;
  double qout_add = 0;
  double OTR_add = 0;
  double Vg_add = 0;
  double CO2in_add = 0;
  double Fin_add = 0;
  double Fpermeate_add = 0;
  double Vres_add = 0;
  double FFox_to_FFred_add = 0;
  double FFred_to_FFox_add = 0;
  double KKox_to_KKred_add = 0;
  double KKred_to_KKox_add = 0;
  double FFred_to_Ffree_add = 0;
  double Ffree_to_FFred_add = 0;
  double KKred_to_Kfree_add = 0;
  double Kfree_to_KKred_add = 0;
  double free_to_FKred_add = 0;
  double FKred_to_free_add = 0;
  double FKred_to_FKox_add = 0;
  double FKox_to_FKred_add = 0;
  double o2balance_add = 0;
  double tioBalance_add = 0;
  while (_itB_times != _itE_times)
  {
    step = _time - _timePrev;
    FFox_k1 = -FFox_to_FFred_t + FFred_to_FFox_t;
    KKox_k1 = -KKox_to_KKred_t + KKred_to_KKox_t;
    FFred_k1 = FFox_to_FFred_t - FFred_to_FFox_t - FFred_to_Ffree_t + Ffree_to_FFred_t;
    KKred_k1 = KKox_to_KKred_t - KKred_to_KKox_t - KKred_to_Kfree_t + Kfree_to_KKred_t;
    Ffree_k1 = 2 * FFred_to_Ffree_t - 2 * Ffree_to_FFred_t - free_to_FKred_t + FKred_to_free_t;
    Kfree_k1 = 2 * KKred_to_Kfree_t - 2 * Kfree_to_KKred_t - free_to_FKred_t + FKred_to_free_t;
    FKred_k1 = free_to_FKred_t - FKred_to_free_t - FKred_to_FKox_t + FKox_to_FKred_t;
    FKox_k1 = FKred_to_FKox_t - FKox_to_FKred_t;
    MEAthiol_k1 = -2 * FFox_to_FFred_t + 2 * FFred_to_FFox_t - 2 * KKox_to_KKred_t + 2 * KKred_to_KKox_t + 2 * FFred_to_Ffree_t + 2 * KKred_to_Kfree_t - 2 * Ffree_to_FFred_t - 2 * Kfree_to_KKred_t;
    CO2_k1 = Fin_t * CO2in_t / VL_t - Fpermeate_t * CO2_t / VL_t + o2balance_t - Fpermeate_t * CO2_t / VL_t - 0.5 * tioBalance_t;
    yO2P_k1 = yO2in * qin - yO2P_t * qout_t - o2balance_t * (VL_t / Vg_t) * R * T * P;
    Cystamine_k1 = tioBalance_t - krcyst * Cystamine_t * Vres_t - 2 * k4ox * FKred_t * Vres_t * pow(Cystamine_t, 2) * pow(Vres_t, 2) + 2 * k4ox * FKred_t * Vres_t * pow(Cystamine_t, 2) * pow(Vres_t, 2) - Cystamine_t * Fin_t / VL_t - Cystamine_t * Fpermeate_t / VL_t;
    VL_k1 = Fin_t - Fpermeate_t;
    FFox_add = FFox_t + 0.2 * step * FFox_k1;
    KKox_add = KKox_t + 0.2 * step * KKox_k1;
    FFred_add = FFred_t + 0.2 * step * FFred_k1;
    KKred_add = KKred_t + 0.2 * step * KKred_k1;
    Ffree_add = Ffree_t + 0.2 * step * Ffree_k1;
    Kfree_add = Kfree_t + 0.2 * step * Kfree_k1;
    FKred_add = FKred_t + 0.2 * step * FKred_k1;
    FKox_add = FKox_t + 0.2 * step * FKox_k1;
    MEAthiol_add = MEAthiol_t + 0.2 * step * MEAthiol_k1;
    CO2_add = CO2_t + 0.2 * step * CO2_k1;
    yO2P_add = yO2P_t + 0.2 * step * yO2P_k1;
    Cystamine_add = Cystamine_t + 0.2 * step * Cystamine_k1;
    VL_add = VL_t + 0.2 * step * VL_k1;
    MEAthiolate_add = MEAthiol_add * pow(10, (pH - pKa2MEA));
    klasurface_add = 0.33 * pow(VL_add, (-0.281)) * pow((AgitatorPowerNumber * pow((AgitatorSpeed / 60), 3) * pow(AgitatorDiameter, 5) * pow(2.54, 3) / (pow(39.37, 2) * 1000) / (VL_add / 1000)), 0.36);
    qout_add = qin - klasurface_add * (yO2P_add * H - CO2_add) * VL_add * R * T / (P * 1000);
    OTR_add = klasurface_add * (yO2P_add * H - CO2_add);
    Vg_add = Vtotalvessel - VL_add;
    CO2in_add = percentO2saturation * 7.17 / (32 * 100);
    Fin_add = _time < 120 ? (0) : (0.025);
    Fpermeate_add = _time < 120 ? (0.025) : (Fin_add);
    Vres_add = VLinitial / VL_add;
    FFox_to_FFred_add = k1red * FFox_add * Vres_add * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    FFred_to_FFox_add = k1ox * FFred_add * Vres_add;
    KKox_to_KKred_add = k1red * KKox_add * Vres_add * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    KKred_to_KKox_add = k1ox * KKred_add * Vres_add;
    FFred_to_Ffree_add = k2Fd * FFred_add * Vres_add;
    Ffree_to_FFred_add = k2Fa * pow(Ffree_add, 2) * pow(Vres_add, 2) * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    KKred_to_Kfree_add = k2Kd * KKred_add * Vres_add;
    Kfree_to_KKred_add = k2Ka * pow(Kfree_add, 2) * pow(Vres_add, 2) * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    free_to_FKred_add = k3FKa * Ffree_add * Vres_add * Kfree_add * Vres_add;
    FKred_to_free_add = k3FKd * FKred_add * Vres_add;
    FKred_to_FKox_add = k4ox * FKred_add * Vres_add * pow(Cystamine_add, 2) * pow(Vres_add, 2);
    FKox_to_FKred_add = k4red * FKox_add * Vres_add * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    o2balance_add = klasurface_add * (yO2P_add * H - CO2_add);
    tioBalance_add = kthiolox * pow(MEAthiol_add, 2) * pow(Vres_add, 2) * sqrt(CO2_add) * sqrt(Vres_add);

    FFox_k2 = -FFox_to_FFred_add + FFred_to_FFox_add;
    KKox_k2 = -KKox_to_KKred_add + KKred_to_KKox_add;
    FFred_k2 = FFox_to_FFred_add - FFred_to_FFox_add - FFred_to_Ffree_add + Ffree_to_FFred_add;
    KKred_k2 = KKox_to_KKred_add - KKred_to_KKox_add - KKred_to_Kfree_add + Kfree_to_KKred_add;
    Ffree_k2 = 2 * FFred_to_Ffree_add - 2 * Ffree_to_FFred_add - free_to_FKred_add + FKred_to_free_add;
    Kfree_k2 = 2 * KKred_to_Kfree_add - 2 * Kfree_to_KKred_add - free_to_FKred_add + FKred_to_free_add;
    FKred_k2 = free_to_FKred_add - FKred_to_free_add - FKred_to_FKox_add + FKox_to_FKred_add;
    FKox_k2 = FKred_to_FKox_add - FKox_to_FKred_add;
    MEAthiol_k2 = -2 * FFox_to_FFred_add + 2 * FFred_to_FFox_add - 2 * KKox_to_KKred_add + 2 * KKred_to_KKox_add + 2 * FFred_to_Ffree_add + 2 * KKred_to_Kfree_add - 2 * Ffree_to_FFred_add - 2 * Kfree_to_KKred_add;
    CO2_k2 = Fin_add * CO2in_add / VL_add - Fpermeate_add * CO2_add / VL_add + o2balance_add - Fpermeate_add * CO2_add / VL_add - 0.5 * tioBalance_add;
    yO2P_k2 = yO2in * qin - yO2P_add * qout_add - o2balance_add * (VL_add / Vg_add) * R * T * P;
    Cystamine_k2 = tioBalance_add - krcyst * Cystamine_add * Vres_add - 2 * k4ox * FKred_add * Vres_add * pow(Cystamine_add, 2) * pow(Vres_add, 2) + 2 * k4ox * FKred_add * Vres_add * pow(Cystamine_add, 2) * pow(Vres_add, 2) - Cystamine_add * Fin_add / VL_add - Cystamine_add * Fpermeate_add / VL_add;
    VL_k2 = Fin_add - Fpermeate_add;
    FFox_add = FFox_t + step * (0.075 * FFox_k1 + 0.225 * FFox_k2);
    KKox_add = KKox_t + step * (0.075 * KKox_k1 + 0.225 * KKox_k2);
    FFred_add = FFred_t + step * (0.075 * FFred_k1 + 0.225 * FFred_k2);
    KKred_add = KKred_t + step * (0.075 * KKred_k1 + 0.225 * KKred_k2);
    Ffree_add = Ffree_t + step * (0.075 * Ffree_k1 + 0.225 * Ffree_k2);
    Kfree_add = Kfree_t + step * (0.075 * Kfree_k1 + 0.225 * Kfree_k2);
    FKred_add = FKred_t + step * (0.075 * FKred_k1 + 0.225 * FKred_k2);
    FKox_add = FKox_t + step * (0.075 * FKox_k1 + 0.225 * FKox_k2);
    MEAthiol_add = MEAthiol_t + step * (0.075 * MEAthiol_k1 + 0.225 * MEAthiol_k2);
    CO2_add = CO2_t + step * (0.075 * CO2_k1 + 0.225 * CO2_k2);
    yO2P_add = yO2P_t + step * (0.075 * yO2P_k1 + 0.225 * yO2P_k2);
    Cystamine_add = Cystamine_t + step * (0.075 * Cystamine_k1 + 0.225 * Cystamine_k2);
    VL_add = VL_t + step * (0.075 * VL_k1 + 0.225 * VL_k2);
    MEAthiolate_add = MEAthiol_add * pow(10, (pH - pKa2MEA));
    klasurface_add = 0.33 * pow(VL_add, (-0.281)) * pow((AgitatorPowerNumber * pow((AgitatorSpeed / 60), 3) * pow(AgitatorDiameter, 5) * pow(2.54, 3) / (pow(39.37, 2) * 1000) / (VL_add / 1000)), 0.36);
    qout_add = qin - klasurface_add * (yO2P_add * H - CO2_add) * VL_add * R * T / (P * 1000);
    OTR_add = klasurface_add * (yO2P_add * H - CO2_add);
    Vg_add = Vtotalvessel - VL_add;
    CO2in_add = percentO2saturation * 7.17 / (32 * 100);
    Fin_add = _time < 120 ? (0) : (0.025);
    Fpermeate_add = _time < 120 ? (0.025) : (Fin_add);
    Vres_add = VLinitial / VL_add;
    FFox_to_FFred_add = k1red * FFox_add * Vres_add * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    FFred_to_FFox_add = k1ox * FFred_add * Vres_add;
    KKox_to_KKred_add = k1red * KKox_add * Vres_add * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    KKred_to_KKox_add = k1ox * KKred_add * Vres_add;
    FFred_to_Ffree_add = k2Fd * FFred_add * Vres_add;
    Ffree_to_FFred_add = k2Fa * pow(Ffree_add, 2) * pow(Vres_add, 2) * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    KKred_to_Kfree_add = k2Kd * KKred_add * Vres_add;
    Kfree_to_KKred_add = k2Ka * pow(Kfree_add, 2) * pow(Vres_add, 2) * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    free_to_FKred_add = k3FKa * Ffree_add * Vres_add * Kfree_add * Vres_add;
    FKred_to_free_add = k3FKd * FKred_add * Vres_add;
    FKred_to_FKox_add = k4ox * FKred_add * Vres_add * pow(Cystamine_add, 2) * pow(Vres_add, 2);
    FKox_to_FKred_add = k4red * FKox_add * Vres_add * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    o2balance_add = klasurface_add * (yO2P_add * H - CO2_add);
    tioBalance_add = kthiolox * pow(MEAthiol_add, 2) * pow(Vres_add, 2) * sqrt(CO2_add) * sqrt(Vres_add);

    FFox_k3 = -FFox_to_FFred_add + FFred_to_FFox_add;
    KKox_k3 = -KKox_to_KKred_add + KKred_to_KKox_add;
    FFred_k3 = FFox_to_FFred_add - FFred_to_FFox_add - FFred_to_Ffree_add + Ffree_to_FFred_add;
    KKred_k3 = KKox_to_KKred_add - KKred_to_KKox_add - KKred_to_Kfree_add + Kfree_to_KKred_add;
    Ffree_k3 = 2 * FFred_to_Ffree_add - 2 * Ffree_to_FFred_add - free_to_FKred_add + FKred_to_free_add;
    Kfree_k3 = 2 * KKred_to_Kfree_add - 2 * Kfree_to_KKred_add - free_to_FKred_add + FKred_to_free_add;
    FKred_k3 = free_to_FKred_add - FKred_to_free_add - FKred_to_FKox_add + FKox_to_FKred_add;
    FKox_k3 = FKred_to_FKox_add - FKox_to_FKred_add;
    MEAthiol_k3 = -2 * FFox_to_FFred_add + 2 * FFred_to_FFox_add - 2 * KKox_to_KKred_add + 2 * KKred_to_KKox_add + 2 * FFred_to_Ffree_add + 2 * KKred_to_Kfree_add - 2 * Ffree_to_FFred_add - 2 * Kfree_to_KKred_add;
    CO2_k3 = Fin_add * CO2in_add / VL_add - Fpermeate_add * CO2_add / VL_add + o2balance_add - Fpermeate_add * CO2_add / VL_add - 0.5 * tioBalance_add;
    yO2P_k3 = yO2in * qin - yO2P_add * qout_add - o2balance_add * (VL_add / Vg_add) * R * T * P;
    Cystamine_k3 = tioBalance_add - krcyst * Cystamine_add * Vres_add - 2 * k4ox * FKred_add * Vres_add * pow(Cystamine_add, 2) * pow(Vres_add, 2) + 2 * k4ox * FKred_add * Vres_add * pow(Cystamine_add, 2) * pow(Vres_add, 2) - Cystamine_add * Fin_add / VL_add - Cystamine_add * Fpermeate_add / VL_add;
    VL_k3 = Fin_add - Fpermeate_add;
    FFox_add = FFox_t + step * (0.3 * FFox_k1 - 0.9 * FFox_k2 + 1.2 * FFox_k3);
    KKox_add = KKox_t + step * (0.3 * KKox_k1 - 0.9 * KKox_k2 + 1.2 * KKox_k3);
    FFred_add = FFred_t + step * (0.3 * FFred_k1 - 0.9 * FFred_k2 + 1.2 * FFred_k3);
    KKred_add = KKred_t + step * (0.3 * KKred_k1 - 0.9 * KKred_k2 + 1.2 * KKred_k3);
    Ffree_add = Ffree_t + step * (0.3 * Ffree_k1 - 0.9 * Ffree_k2 + 1.2 * Ffree_k3);
    Kfree_add = Kfree_t + step * (0.3 * Kfree_k1 - 0.9 * Kfree_k2 + 1.2 * Kfree_k3);
    FKred_add = FKred_t + step * (0.3 * FKred_k1 - 0.9 * FKred_k2 + 1.2 * FKred_k3);
    FKox_add = FKox_t + step * (0.3 * FKox_k1 - 0.9 * FKox_k2 + 1.2 * FKox_k3);
    MEAthiol_add = MEAthiol_t + step * (0.3 * MEAthiol_k1 - 0.9 * MEAthiol_k2 + 1.2 * MEAthiol_k3);
    CO2_add = CO2_t + step * (0.3 * CO2_k1 - 0.9 * CO2_k2 + 1.2 * CO2_k3);
    yO2P_add = yO2P_t + step * (0.3 * yO2P_k1 - 0.9 * yO2P_k2 + 1.2 * yO2P_k3);
    Cystamine_add = Cystamine_t + step * (0.3 * Cystamine_k1 - 0.9 * Cystamine_k2 + 1.2 * Cystamine_k3);
    VL_add = VL_t + step * (0.3 * VL_k1 - 0.9 * VL_k2 + 1.2 * VL_k3);
    MEAthiolate_add = MEAthiol_add * pow(10, (pH - pKa2MEA));
    klasurface_add = 0.33 * pow(VL_add, (-0.281)) * pow((AgitatorPowerNumber * pow((AgitatorSpeed / 60), 3) * pow(AgitatorDiameter, 5) * pow(2.54, 3) / (pow(39.37, 2) * 1000) / (VL_add / 1000)), 0.36);
    qout_add = qin - klasurface_add * (yO2P_add * H - CO2_add) * VL_add * R * T / (P * 1000);
    OTR_add = klasurface_add * (yO2P_add * H - CO2_add);
    Vg_add = Vtotalvessel - VL_add;
    CO2in_add = percentO2saturation * 7.17 / (32 * 100);
    Fin_add = _time < 120 ? (0) : (0.025);
    Fpermeate_add = _time < 120 ? (0.025) : (Fin_add);
    Vres_add = VLinitial / VL_add;
    FFox_to_FFred_add = k1red * FFox_add * Vres_add * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    FFred_to_FFox_add = k1ox * FFred_add * Vres_add;
    KKox_to_KKred_add = k1red * KKox_add * Vres_add * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    KKred_to_KKox_add = k1ox * KKred_add * Vres_add;
    FFred_to_Ffree_add = k2Fd * FFred_add * Vres_add;
    Ffree_to_FFred_add = k2Fa * pow(Ffree_add, 2) * pow(Vres_add, 2) * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    KKred_to_Kfree_add = k2Kd * KKred_add * Vres_add;
    Kfree_to_KKred_add = k2Ka * pow(Kfree_add, 2) * pow(Vres_add, 2) * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    free_to_FKred_add = k3FKa * Ffree_add * Vres_add * Kfree_add * Vres_add;
    FKred_to_free_add = k3FKd * FKred_add * Vres_add;
    FKred_to_FKox_add = k4ox * FKred_add * Vres_add * pow(Cystamine_add, 2) * pow(Vres_add, 2);
    FKox_to_FKred_add = k4red * FKox_add * Vres_add * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    o2balance_add = klasurface_add * (yO2P_add * H - CO2_add);
    tioBalance_add = kthiolox * pow(MEAthiol_add, 2) * pow(Vres_add, 2) * sqrt(CO2_add) * sqrt(Vres_add);

    FFox_k4 = -FFox_to_FFred_add + FFred_to_FFox_add;
    KKox_k4 = -KKox_to_KKred_add + KKred_to_KKox_add;
    FFred_k4 = FFox_to_FFred_add - FFred_to_FFox_add - FFred_to_Ffree_add + Ffree_to_FFred_add;
    KKred_k4 = KKox_to_KKred_add - KKred_to_KKox_add - KKred_to_Kfree_add + Kfree_to_KKred_add;
    Ffree_k4 = 2 * FFred_to_Ffree_add - 2 * Ffree_to_FFred_add - free_to_FKred_add + FKred_to_free_add;
    Kfree_k4 = 2 * KKred_to_Kfree_add - 2 * Kfree_to_KKred_add - free_to_FKred_add + FKred_to_free_add;
    FKred_k4 = free_to_FKred_add - FKred_to_free_add - FKred_to_FKox_add + FKox_to_FKred_add;
    FKox_k4 = FKred_to_FKox_add - FKox_to_FKred_add;
    MEAthiol_k4 = -2 * FFox_to_FFred_add + 2 * FFred_to_FFox_add - 2 * KKox_to_KKred_add + 2 * KKred_to_KKox_add + 2 * FFred_to_Ffree_add + 2 * KKred_to_Kfree_add - 2 * Ffree_to_FFred_add - 2 * Kfree_to_KKred_add;
    CO2_k4 = Fin_add * CO2in_add / VL_add - Fpermeate_add * CO2_add / VL_add + o2balance_add - Fpermeate_add * CO2_add / VL_add - 0.5 * tioBalance_add;
    yO2P_k4 = yO2in * qin - yO2P_add * qout_add - o2balance_add * (VL_add / Vg_add) * R * T * P;
    Cystamine_k4 = tioBalance_add - krcyst * Cystamine_add * Vres_add - 2 * k4ox * FKred_add * Vres_add * pow(Cystamine_add, 2) * pow(Vres_add, 2) + 2 * k4ox * FKred_add * Vres_add * pow(Cystamine_add, 2) * pow(Vres_add, 2) - Cystamine_add * Fin_add / VL_add - Cystamine_add * Fpermeate_add / VL_add;
    VL_k4 = Fin_add - Fpermeate_add;
    FFox_add = FFox_t + step * (-0.2037037037 * FFox_k1 + 2.5 * FFox_k2 - 2.592592593 * FFox_k3 + 1.296296296 * FFox_k4);
    KKox_add = KKox_t + step * (-0.2037037037 * KKox_k1 + 2.5 * KKox_k2 - 2.592592593 * KKox_k3 + 1.296296296 * KKox_k4);
    FFred_add = FFred_t + step * (-0.2037037037 * FFred_k1 + 2.5 * FFred_k2 - 2.592592593 * FFred_k3 + 1.296296296 * FFred_k4);
    KKred_add = KKred_t + step * (-0.2037037037 * KKred_k1 + 2.5 * KKred_k2 - 2.592592593 * KKred_k3 + 1.296296296 * KKred_k4);
    Ffree_add = Ffree_t + step * (-0.2037037037 * Ffree_k1 + 2.5 * Ffree_k2 - 2.592592593 * Ffree_k3 + 1.296296296 * Ffree_k4);
    Kfree_add = Kfree_t + step * (-0.2037037037 * Kfree_k1 + 2.5 * Kfree_k2 - 2.592592593 * Kfree_k3 + 1.296296296 * Kfree_k4);
    FKred_add = FKred_t + step * (-0.2037037037 * FKred_k1 + 2.5 * FKred_k2 - 2.592592593 * FKred_k3 + 1.296296296 * FKred_k4);
    FKox_add = FKox_t + step * (-0.2037037037 * FKox_k1 + 2.5 * FKox_k2 - 2.592592593 * FKox_k3 + 1.296296296 * FKox_k4);
    MEAthiol_add = MEAthiol_t + step * (-0.2037037037 * MEAthiol_k1 + 2.5 * MEAthiol_k2 - 2.592592593 * MEAthiol_k3 + 1.296296296 * MEAthiol_k4);
    CO2_add = CO2_t + step * (-0.2037037037 * CO2_k1 + 2.5 * CO2_k2 - 2.592592593 * CO2_k3 + 1.296296296 * CO2_k4);
    yO2P_add = yO2P_t + step * (-0.2037037037 * yO2P_k1 + 2.5 * yO2P_k2 - 2.592592593 * yO2P_k3 + 1.296296296 * yO2P_k4);
    Cystamine_add = Cystamine_t + step * (-0.2037037037 * Cystamine_k1 + 2.5 * Cystamine_k2 - 2.592592593 * Cystamine_k3 + 1.296296296 * Cystamine_k4);
    VL_add = VL_t + step * (-0.2037037037 * VL_k1 + 2.5 * VL_k2 - 2.592592593 * VL_k3 + 1.296296296 * VL_k4);
    MEAthiolate_add = MEAthiol_add * pow(10, (pH - pKa2MEA));
    klasurface_add = 0.33 * pow(VL_add, (-0.281)) * pow((AgitatorPowerNumber * pow((AgitatorSpeed / 60), 3) * pow(AgitatorDiameter, 5) * pow(2.54, 3) / (pow(39.37, 2) * 1000) / (VL_add / 1000)), 0.36);
    qout_add = qin - klasurface_add * (yO2P_add * H - CO2_add) * VL_add * R * T / (P * 1000);
    OTR_add = klasurface_add * (yO2P_add * H - CO2_add);
    Vg_add = Vtotalvessel - VL_add;
    CO2in_add = percentO2saturation * 7.17 / (32 * 100);
    Fin_add = _time < 120 ? (0) : (0.025);
    Fpermeate_add = _time < 120 ? (0.025) : (Fin_add);
    Vres_add = VLinitial / VL_add;
    FFox_to_FFred_add = k1red * FFox_add * Vres_add * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    FFred_to_FFox_add = k1ox * FFred_add * Vres_add;
    KKox_to_KKred_add = k1red * KKox_add * Vres_add * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    KKred_to_KKox_add = k1ox * KKred_add * Vres_add;
    FFred_to_Ffree_add = k2Fd * FFred_add * Vres_add;
    Ffree_to_FFred_add = k2Fa * pow(Ffree_add, 2) * pow(Vres_add, 2) * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    KKred_to_Kfree_add = k2Kd * KKred_add * Vres_add;
    Kfree_to_KKred_add = k2Ka * pow(Kfree_add, 2) * pow(Vres_add, 2) * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    free_to_FKred_add = k3FKa * Ffree_add * Vres_add * Kfree_add * Vres_add;
    FKred_to_free_add = k3FKd * FKred_add * Vres_add;
    FKred_to_FKox_add = k4ox * FKred_add * Vres_add * pow(Cystamine_add, 2) * pow(Vres_add, 2);
    FKox_to_FKred_add = k4red * FKox_add * Vres_add * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    o2balance_add = klasurface_add * (yO2P_add * H - CO2_add);
    tioBalance_add = kthiolox * pow(MEAthiol_add, 2) * pow(Vres_add, 2) * sqrt(CO2_add) * sqrt(Vres_add);

    FFox_k5 = -FFox_to_FFred_add + FFred_to_FFox_add;
    KKox_k5 = -KKox_to_KKred_add + KKred_to_KKox_add;
    FFred_k5 = FFox_to_FFred_add - FFred_to_FFox_add - FFred_to_Ffree_add + Ffree_to_FFred_add;
    KKred_k5 = KKox_to_KKred_add - KKred_to_KKox_add - KKred_to_Kfree_add + Kfree_to_KKred_add;
    Ffree_k5 = 2 * FFred_to_Ffree_add - 2 * Ffree_to_FFred_add - free_to_FKred_add + FKred_to_free_add;
    Kfree_k5 = 2 * KKred_to_Kfree_add - 2 * Kfree_to_KKred_add - free_to_FKred_add + FKred_to_free_add;
    FKred_k5 = free_to_FKred_add - FKred_to_free_add - FKred_to_FKox_add + FKox_to_FKred_add;
    FKox_k5 = FKred_to_FKox_add - FKox_to_FKred_add;
    MEAthiol_k5 = -2 * FFox_to_FFred_add + 2 * FFred_to_FFox_add - 2 * KKox_to_KKred_add + 2 * KKred_to_KKox_add + 2 * FFred_to_Ffree_add + 2 * KKred_to_Kfree_add - 2 * Ffree_to_FFred_add - 2 * Kfree_to_KKred_add;
    CO2_k5 = Fin_add * CO2in_add / VL_add - Fpermeate_add * CO2_add / VL_add + o2balance_add - Fpermeate_add * CO2_add / VL_add - 0.5 * tioBalance_add;
    yO2P_k5 = yO2in * qin - yO2P_add * qout_add - o2balance_add * (VL_add / Vg_add) * R * T * P;
    Cystamine_k5 = tioBalance_add - krcyst * Cystamine_add * Vres_add - 2 * k4ox * FKred_add * Vres_add * pow(Cystamine_add, 2) * pow(Vres_add, 2) + 2 * k4ox * FKred_add * Vres_add * pow(Cystamine_add, 2) * pow(Vres_add, 2) - Cystamine_add * Fin_add / VL_add - Cystamine_add * Fpermeate_add / VL_add;
    VL_k5 = Fin_add - Fpermeate_add;
    FFox_add = FFox_t + step * (0.0294958044 * FFox_k1 + 0.341796875 * FFox_k2 + 0.04158229679 * FFox_k3 + 0.4003454138 * FFox_k4);
    KKox_add = KKox_t + step * (0.0294958044 * KKox_k1 + 0.341796875 * KKox_k2 + 0.04158229679 * KKox_k3 + 0.4003454138 * KKox_k4);
    FFred_add = FFred_t + step * (0.0294958044 * FFred_k1 + 0.341796875 * FFred_k2 + 0.04158229679 * FFred_k3 + 0.4003454138 * FFred_k4);
    KKred_add = KKred_t + step * (0.0294958044 * KKred_k1 + 0.341796875 * KKred_k2 + 0.04158229679 * KKred_k3 + 0.4003454138 * KKred_k4);
    Ffree_add = Ffree_t + step * (0.0294958044 * Ffree_k1 + 0.341796875 * Ffree_k2 + 0.04158229679 * Ffree_k3 + 0.4003454138 * Ffree_k4);
    Kfree_add = Kfree_t + step * (0.0294958044 * Kfree_k1 + 0.341796875 * Kfree_k2 + 0.04158229679 * Kfree_k3 + 0.4003454138 * Kfree_k4);
    FKred_add = FKred_t + step * (0.0294958044 * FKred_k1 + 0.341796875 * FKred_k2 + 0.04158229679 * FKred_k3 + 0.4003454138 * FKred_k4);
    FKox_add = FKox_t + step * (0.0294958044 * FKox_k1 + 0.341796875 * FKox_k2 + 0.04158229679 * FKox_k3 + 0.4003454138 * FKox_k4);
    MEAthiol_add = MEAthiol_t + step * (0.0294958044 * MEAthiol_k1 + 0.341796875 * MEAthiol_k2 + 0.04158229679 * MEAthiol_k3 + 0.4003454138 * MEAthiol_k4);
    CO2_add = CO2_t + step * (0.0294958044 * CO2_k1 + 0.341796875 * CO2_k2 + 0.04158229679 * CO2_k3 + 0.4003454138 * CO2_k4);
    yO2P_add = yO2P_t + step * (0.0294958044 * yO2P_k1 + 0.341796875 * yO2P_k2 + 0.04158229679 * yO2P_k3 + 0.4003454138 * yO2P_k4);
    Cystamine_add = Cystamine_t + step * (0.0294958044 * Cystamine_k1 + 0.341796875 * Cystamine_k2 + 0.04158229679 * Cystamine_k3 + 0.4003454138 * Cystamine_k4);
    VL_add = VL_t + step * (0.0294958044 * VL_k1 + 0.341796875 * VL_k2 + 0.04158229679 * VL_k3 + 0.4003454138 * VL_k4);
    MEAthiolate_add = MEAthiol_add * pow(10, (pH - pKa2MEA));
    klasurface_add = 0.33 * pow(VL_add, (-0.281)) * pow((AgitatorPowerNumber * pow((AgitatorSpeed / 60), 3) * pow(AgitatorDiameter, 5) * pow(2.54, 3) / (pow(39.37, 2) * 1000) / (VL_add / 1000)), 0.36);
    qout_add = qin - klasurface_add * (yO2P_add * H - CO2_add) * VL_add * R * T / (P * 1000);
    OTR_add = klasurface_add * (yO2P_add * H - CO2_add);
    Vg_add = Vtotalvessel - VL_add;
    CO2in_add = percentO2saturation * 7.17 / (32 * 100);
    Fin_add = _time < 120 ? (0) : (0.025);
    Fpermeate_add = _time < 120 ? (0.025) : (Fin_add);
    Vres_add = VLinitial / VL_add;
    FFox_to_FFred_add = k1red * FFox_add * Vres_add * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    FFred_to_FFox_add = k1ox * FFred_add * Vres_add;
    KKox_to_KKred_add = k1red * KKox_add * Vres_add * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    KKred_to_KKox_add = k1ox * KKred_add * Vres_add;
    FFred_to_Ffree_add = k2Fd * FFred_add * Vres_add;
    Ffree_to_FFred_add = k2Fa * pow(Ffree_add, 2) * pow(Vres_add, 2) * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    KKred_to_Kfree_add = k2Kd * KKred_add * Vres_add;
    Kfree_to_KKred_add = k2Ka * pow(Kfree_add, 2) * pow(Vres_add, 2) * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    free_to_FKred_add = k3FKa * Ffree_add * Vres_add * Kfree_add * Vres_add;
    FKred_to_free_add = k3FKd * FKred_add * Vres_add;
    FKred_to_FKox_add = k4ox * FKred_add * Vres_add * pow(Cystamine_add, 2) * pow(Vres_add, 2);
    FKox_to_FKred_add = k4red * FKox_add * Vres_add * pow(MEAthiolate_add, 2) * pow(Vres_add, 2);
    o2balance_add = klasurface_add * (yO2P_add * H - CO2_add);
    tioBalance_add = kthiolox * pow(MEAthiol_add, 2) * pow(Vres_add, 2) * sqrt(CO2_add) * sqrt(Vres_add);

    FFox_k6 = -FFox_to_FFred_add + FFred_to_FFox_add;
    KKox_k6 = -KKox_to_KKred_add + KKred_to_KKox_add;
    FFred_k6 = FFox_to_FFred_add - FFred_to_FFox_add - FFred_to_Ffree_add + Ffree_to_FFred_add;
    KKred_k6 = KKox_to_KKred_add - KKred_to_KKox_add - KKred_to_Kfree_add + Kfree_to_KKred_add;
    Ffree_k6 = 2 * FFred_to_Ffree_add - 2 * Ffree_to_FFred_add - free_to_FKred_add + FKred_to_free_add;
    Kfree_k6 = 2 * KKred_to_Kfree_add - 2 * Kfree_to_KKred_add - free_to_FKred_add + FKred_to_free_add;
    FKred_k6 = free_to_FKred_add - FKred_to_free_add - FKred_to_FKox_add + FKox_to_FKred_add;
    FKox_k6 = FKred_to_FKox_add - FKox_to_FKred_add;
    MEAthiol_k6 = -2 * FFox_to_FFred_add + 2 * FFred_to_FFox_add - 2 * KKox_to_KKred_add + 2 * KKred_to_KKox_add + 2 * FFred_to_Ffree_add + 2 * KKred_to_Kfree_add - 2 * Ffree_to_FFred_add - 2 * Kfree_to_KKred_add;
    CO2_k6 = Fin_add * CO2in_add / VL_add - Fpermeate_add * CO2_add / VL_add + o2balance_add - Fpermeate_add * CO2_add / VL_add - 0.5 * tioBalance_add;
    yO2P_k6 = yO2in * qin - yO2P_add * qout_add - o2balance_add * (VL_add / Vg_add) * R * T * P;
    Cystamine_k6 = tioBalance_add - krcyst * Cystamine_add * Vres_add - 2 * k4ox * FKred_add * Vres_add * pow(Cystamine_add, 2) * pow(Vres_add, 2) + 2 * k4ox * FKred_add * Vres_add * pow(Cystamine_add, 2) * pow(Vres_add, 2) - Cystamine_add * Fin_add / VL_add - Cystamine_add * Fpermeate_add / VL_add;
    VL_k6 = Fin_add - Fpermeate_add;

    FFox_scale = abs(FFox_t) + step * abs(FFox_k1) + _infinitezimal;
    KKox_scale = abs(KKox_t) + step * abs(KKox_k1) + _infinitezimal;
    FFred_scale = abs(FFred_t) + step * abs(FFred_k1) + _infinitezimal;
    KKred_scale = abs(KKred_t) + step * abs(KKred_k1) + _infinitezimal;
    Ffree_scale = abs(Ffree_t) + step * abs(Ffree_k1) + _infinitezimal;
    Kfree_scale = abs(Kfree_t) + step * abs(Kfree_k1) + _infinitezimal;
    FKred_scale = abs(FKred_t) + step * abs(FKred_k1) + _infinitezimal;
    FKox_scale = abs(FKox_t) + step * abs(FKox_k1) + _infinitezimal;
    MEAthiol_scale = abs(MEAthiol_t) + step * abs(MEAthiol_k1) + _infinitezimal;
    CO2_scale = abs(CO2_t) + step * abs(CO2_k1) + _infinitezimal;
    yO2P_scale = abs(yO2P_t) + step * abs(yO2P_k1) + _infinitezimal;
    Cystamine_scale = abs(Cystamine_t) + step * abs(Cystamine_k1) + _infinitezimal;
    VL_scale = abs(VL_t) + step * abs(VL_k1) + _infinitezimal;
    FFox_4it = FFox_t + step * (0.09788359788 * FFox_k1 + 0.4025764895 * FFox_k3 + 0.2104377104 * FFox_k4 + 0.2891022021 * FFox_k6);
    KKox_4it = KKox_t + step * (0.09788359788 * KKox_k1 + 0.4025764895 * KKox_k3 + 0.2104377104 * KKox_k4 + 0.2891022021 * KKox_k6);
    FFred_4it = FFred_t + step * (0.09788359788 * FFred_k1 + 0.4025764895 * FFred_k3 + 0.2104377104 * FFred_k4 + 0.2891022021 * FFred_k6);
    KKred_4it = KKred_t + step * (0.09788359788 * KKred_k1 + 0.4025764895 * KKred_k3 + 0.2104377104 * KKred_k4 + 0.2891022021 * KKred_k6);
    Ffree_4it = Ffree_t + step * (0.09788359788 * Ffree_k1 + 0.4025764895 * Ffree_k3 + 0.2104377104 * Ffree_k4 + 0.2891022021 * Ffree_k6);
    Kfree_4it = Kfree_t + step * (0.09788359788 * Kfree_k1 + 0.4025764895 * Kfree_k3 + 0.2104377104 * Kfree_k4 + 0.2891022021 * Kfree_k6);
    FKred_4it = FKred_t + step * (0.09788359788 * FKred_k1 + 0.4025764895 * FKred_k3 + 0.2104377104 * FKred_k4 + 0.2891022021 * FKred_k6);
    FKox_4it = FKox_t + step * (0.09788359788 * FKox_k1 + 0.4025764895 * FKox_k3 + 0.2104377104 * FKox_k4 + 0.2891022021 * FKox_k6);
    MEAthiol_4it = MEAthiol_t + step * (0.09788359788 * MEAthiol_k1 + 0.4025764895 * MEAthiol_k3 + 0.2104377104 * MEAthiol_k4 + 0.2891022021 * MEAthiol_k6);
    CO2_4it = CO2_t + step * (0.09788359788 * CO2_k1 + 0.4025764895 * CO2_k3 + 0.2104377104 * CO2_k4 + 0.2891022021 * CO2_k6);
    yO2P_4it = yO2P_t + step * (0.09788359788 * yO2P_k1 + 0.4025764895 * yO2P_k3 + 0.2104377104 * yO2P_k4 + 0.2891022021 * yO2P_k6);
    Cystamine_4it = Cystamine_t + step * (0.09788359788 * Cystamine_k1 + 0.4025764895 * Cystamine_k3 + 0.2104377104 * Cystamine_k4 + 0.2891022021 * Cystamine_k6);
    VL_4it = VL_t + step * (0.09788359788 * VL_k1 + 0.4025764895 * VL_k3 + 0.2104377104 * VL_k4 + 0.2891022021 * VL_k6);
    MEAthiolate_t = MEAthiol_t * pow(10, (pH - pKa2MEA));
    klasurface_t = 0.33 * pow(VL_t, (-0.281)) * pow((AgitatorPowerNumber * pow((AgitatorSpeed / 60), 3) * pow(AgitatorDiameter, 5) * pow(2.54, 3) / (pow(39.37, 2) * 1000) / (VL_t / 1000)), 0.36);
    qout_t = qin - klasurface_t * (yO2P_t * H - CO2_t) * VL_t * R * T / (P * 1000);
    OTR_t = klasurface_t * (yO2P_t * H - CO2_t);
    Vg_t = Vtotalvessel - VL_t;
    CO2in_t = percentO2saturation * 7.17 / (32 * 100);
    Fin_t = _time < 120 ? (0) : (0.025);
    Fpermeate_t = _time < 120 ? (0.025) : (Fin_t);
    Vres_t = VLinitial / VL_t;
    FFox_to_FFred_t = k1red * FFox_t * Vres_t * pow(MEAthiolate_t, 2) * pow(Vres_t, 2);
    FFred_to_FFox_t = k1ox * FFred_t * Vres_t;
    KKox_to_KKred_t = k1red * KKox_t * Vres_t * pow(MEAthiolate_t, 2) * pow(Vres_t, 2);
    KKred_to_KKox_t = k1ox * KKred_t * Vres_t;
    FFred_to_Ffree_t = k2Fd * FFred_t * Vres_t;
    Ffree_to_FFred_t = k2Fa * pow(Ffree_t, 2) * pow(Vres_t, 2) * pow(MEAthiolate_t, 2) * pow(Vres_t, 2);
    KKred_to_Kfree_t = k2Kd * KKred_t * Vres_t;
    Kfree_to_KKred_t = k2Ka * pow(Kfree_t, 2) * pow(Vres_t, 2) * pow(MEAthiolate_t, 2) * pow(Vres_t, 2);
    free_to_FKred_t = k3FKa * Ffree_t * Vres_t * Kfree_t * Vres_t;
    FKred_to_free_t = k3FKd * FKred_t * Vres_t;
    FKred_to_FKox_t = k4ox * FKred_t * Vres_t * pow(Cystamine_t, 2) * pow(Vres_t, 2);
    FKox_to_FKred_t = k4red * FKox_t * Vres_t * pow(MEAthiolate_t, 2) * pow(Vres_t, 2);
    o2balance_t = klasurface_t * (yO2P_t * H - CO2_t);
    tioBalance_t = kthiolox * pow(MEAthiol_t, 2) * pow(Vres_t, 2) * sqrt(CO2_t) * sqrt(Vres_t);

    FFox_5it = FFox_t + step * (0.1021773727 * FFox_k1 + 0.3839079034 * FFox_k3 + 0.2445927373 * FFox_k4 + 0.01932198661 * FFox_k5 + 0.25 * FFox_k6);
    KKox_5it = KKox_t + step * (0.1021773727 * KKox_k1 + 0.3839079034 * KKox_k3 + 0.2445927373 * KKox_k4 + 0.01932198661 * KKox_k5 + 0.25 * KKox_k6);
    FFred_5it = FFred_t + step * (0.1021773727 * FFred_k1 + 0.3839079034 * FFred_k3 + 0.2445927373 * FFred_k4 + 0.01932198661 * FFred_k5 + 0.25 * FFred_k6);
    KKred_5it = KKred_t + step * (0.1021773727 * KKred_k1 + 0.3839079034 * KKred_k3 + 0.2445927373 * KKred_k4 + 0.01932198661 * KKred_k5 + 0.25 * KKred_k6);
    Ffree_5it = Ffree_t + step * (0.1021773727 * Ffree_k1 + 0.3839079034 * Ffree_k3 + 0.2445927373 * Ffree_k4 + 0.01932198661 * Ffree_k5 + 0.25 * Ffree_k6);
    Kfree_5it = Kfree_t + step * (0.1021773727 * Kfree_k1 + 0.3839079034 * Kfree_k3 + 0.2445927373 * Kfree_k4 + 0.01932198661 * Kfree_k5 + 0.25 * Kfree_k6);
    FKred_5it = FKred_t + step * (0.1021773727 * FKred_k1 + 0.3839079034 * FKred_k3 + 0.2445927373 * FKred_k4 + 0.01932198661 * FKred_k5 + 0.25 * FKred_k6);
    FKox_5it = FKox_t + step * (0.1021773727 * FKox_k1 + 0.3839079034 * FKox_k3 + 0.2445927373 * FKox_k4 + 0.01932198661 * FKox_k5 + 0.25 * FKox_k6);
    MEAthiol_5it = MEAthiol_t + step * (0.1021773727 * MEAthiol_k1 + 0.3839079034 * MEAthiol_k3 + 0.2445927373 * MEAthiol_k4 + 0.01932198661 * MEAthiol_k5 + 0.25 * MEAthiol_k6);
    CO2_5it = CO2_t + step * (0.1021773727 * CO2_k1 + 0.3839079034 * CO2_k3 + 0.2445927373 * CO2_k4 + 0.01932198661 * CO2_k5 + 0.25 * CO2_k6);
    yO2P_5it = yO2P_t + step * (0.1021773727 * yO2P_k1 + 0.3839079034 * yO2P_k3 + 0.2445927373 * yO2P_k4 + 0.01932198661 * yO2P_k5 + 0.25 * yO2P_k6);
    Cystamine_5it = Cystamine_t + step * (0.1021773727 * Cystamine_k1 + 0.3839079034 * Cystamine_k3 + 0.2445927373 * Cystamine_k4 + 0.01932198661 * Cystamine_k5 + 0.25 * Cystamine_k6);
    VL_5it = VL_t + step * (0.1021773727 * VL_k1 + 0.3839079034 * VL_k3 + 0.2445927373 * VL_k4 + 0.01932198661 * VL_k5 + 0.25 * VL_k6);
    MEAthiolate_t = MEAthiol_t * pow(10, (pH - pKa2MEA));
    klasurface_t = 0.33 * pow(VL_t, (-0.281)) * pow((AgitatorPowerNumber * pow((AgitatorSpeed / 60), 3) * pow(AgitatorDiameter, 5) * pow(2.54, 3) / (pow(39.37, 2) * 1000) / (VL_t / 1000)), 0.36);
    qout_t = qin - klasurface_t * (yO2P_t * H - CO2_t) * VL_t * R * T / (P * 1000);
    OTR_t = klasurface_t * (yO2P_t * H - CO2_t);
    Vg_t = Vtotalvessel - VL_t;
    CO2in_t = percentO2saturation * 7.17 / (32 * 100);
    Fin_t = _time < 120 ? (0) : (0.025);
    Fpermeate_t = _time < 120 ? (0.025) : (Fin_t);
    Vres_t = VLinitial / VL_t;
    FFox_to_FFred_t = k1red * FFox_t * Vres_t * pow(MEAthiolate_t, 2) * pow(Vres_t, 2);
    FFred_to_FFox_t = k1ox * FFred_t * Vres_t;
    KKox_to_KKred_t = k1red * KKox_t * Vres_t * pow(MEAthiolate_t, 2) * pow(Vres_t, 2);
    KKred_to_KKox_t = k1ox * KKred_t * Vres_t;
    FFred_to_Ffree_t = k2Fd * FFred_t * Vres_t;
    Ffree_to_FFred_t = k2Fa * pow(Ffree_t, 2) * pow(Vres_t, 2) * pow(MEAthiolate_t, 2) * pow(Vres_t, 2);
    KKred_to_Kfree_t = k2Kd * KKred_t * Vres_t;
    Kfree_to_KKred_t = k2Ka * pow(Kfree_t, 2) * pow(Vres_t, 2) * pow(MEAthiolate_t, 2) * pow(Vres_t, 2);
    free_to_FKred_t = k3FKa * Ffree_t * Vres_t * Kfree_t * Vres_t;
    FKred_to_free_t = k3FKd * FKred_t * Vres_t;
    FKred_to_FKox_t = k4ox * FKred_t * Vres_t * pow(Cystamine_t, 2) * pow(Vres_t, 2);
    FKox_to_FKred_t = k4red * FKox_t * Vres_t * pow(MEAthiolate_t, 2) * pow(Vres_t, 2);
    o2balance_t = klasurface_t * (yO2P_t * H - CO2_t);
    tioBalance_t = kthiolox * pow(MEAthiol_t, 2) * pow(Vres_t, 2) * sqrt(CO2_t) * sqrt(Vres_t);

    FFox_t = FFox_5it;
    KKox_t = KKox_5it;
    FFred_t = FFred_5it;
    KKred_t = KKred_5it;
    Ffree_t = Ffree_5it;
    Kfree_t = Kfree_5it;
    FKred_t = FKred_5it;
    FKox_t = FKox_5it;
    MEAthiol_t = MEAthiol_5it;
    CO2_t = CO2_5it;
    yO2P_t = yO2P_5it;
    Cystamine_t = Cystamine_5it;
    VL_t = VL_5it;

    MEAthiolate_t = MEAthiol_t * pow(10, (pH - pKa2MEA));
    klasurface_t = 0.33 * pow(VL_t, (-0.281)) * pow((AgitatorPowerNumber * pow((AgitatorSpeed / 60), 3) * pow(AgitatorDiameter, 5) * pow(2.54, 3) / (pow(39.37, 2) * 1000) / (VL_t / 1000)), 0.36);
    qout_t = qin - klasurface_t * (yO2P_t * H - CO2_t) * VL_t * R * T / (P * 1000);
    OTR_t = klasurface_t * (yO2P_t * H - CO2_t);
    Vg_t = Vtotalvessel - VL_t;
    CO2in_t = percentO2saturation * 7.17 / (32 * 100);
    Fin_t = _time < 120 ? (0) : (0.025);
    Fpermeate_t = _time < 120 ? (0.025) : (Fin_t);
    Vres_t = VLinitial / VL_t;
    FFox_to_FFred_t = k1red * FFox_t * Vres_t * pow(MEAthiolate_t, 2) * pow(Vres_t, 2);
    FFred_to_FFox_t = k1ox * FFred_t * Vres_t;
    KKox_to_KKred_t = k1red * KKox_t * Vres_t * pow(MEAthiolate_t, 2) * pow(Vres_t, 2);
    KKred_to_KKox_t = k1ox * KKred_t * Vres_t;
    FFred_to_Ffree_t = k2Fd * FFred_t * Vres_t;
    Ffree_to_FFred_t = k2Fa * pow(Ffree_t, 2) * pow(Vres_t, 2) * pow(MEAthiolate_t, 2) * pow(Vres_t, 2);
    KKred_to_Kfree_t = k2Kd * KKred_t * Vres_t;
    Kfree_to_KKred_t = k2Ka * pow(Kfree_t, 2) * pow(Vres_t, 2) * pow(MEAthiolate_t, 2) * pow(Vres_t, 2);
    free_to_FKred_t = k3FKa * Ffree_t * Vres_t * Kfree_t * Vres_t;
    FKred_to_free_t = k3FKd * FKred_t * Vres_t;
    FKred_to_FKox_t = k4ox * FKred_t * Vres_t * pow(Cystamine_t, 2) * pow(Vres_t, 2);
    FKox_to_FKred_t = k4red * FKox_t * Vres_t * pow(MEAthiolate_t, 2) * pow(Vres_t, 2);
    o2balance_t = klasurface_t * (yO2P_t * H - CO2_t);
    tioBalance_t = kthiolox * pow(MEAthiol_t, 2) * pow(Vres_t, 2) * sqrt(CO2_t) * sqrt(Vres_t);

    _timePrev = _time;
    if (_time == *_itB_times) {
      MEAthiol->push_back(MEAthiol_t);
      VL->push_back(VL_t);
      klasurface->push_back(klasurface_t);
      CO2->push_back(CO2_t);
      yO2P->push_back(yO2P_t);
      FFox->push_back(FFox_t);
      KKox->push_back(KKox_t);
      FFred->push_back(FFred_t);
      KKred->push_back(KKred_t);
      Ffree->push_back(Ffree_t);
      Kfree->push_back(Kfree_t);
      FKred->push_back(FKred_t);
      FKox->push_back(FKox_t);
      Cystamine->push_back(Cystamine_t);
      tioBalance->push_back(tioBalance_t);
      ++_itB_times;
    }
    double _delta = 0;
    double _add = 0;
    double _power = _increase ? 0.2 : 0.25;
    _add = abs(FFox_5it - FFox_4it) / FFox_scale;
    _delta = _add > _delta ? _add : _delta;
    _add = abs(KKox_5it - KKox_4it) / KKox_scale;
    _delta = _add > _delta ? _add : _delta;
    _add = abs(FFred_5it - FFred_4it) / FFred_scale;
    _delta = _add > _delta ? _add : _delta;
    _add = abs(KKred_5it - KKred_4it) / KKred_scale;
    _delta = _add > _delta ? _add : _delta;
    _add = abs(Ffree_5it - Ffree_4it) / Ffree_scale;
    _delta = _add > _delta ? _add : _delta;
    _add = abs(Kfree_5it - Kfree_4it) / Kfree_scale;
    _delta = _add > _delta ? _add : _delta;
    _add = abs(FKred_5it - FKred_4it) / FKred_scale;
    _delta = _add > _delta ? _add : _delta;
    _add = abs(FKox_5it - FKox_4it) / FKox_scale;
    _delta = _add > _delta ? _add : _delta;
    _add = abs(MEAthiol_5it - MEAthiol_4it) / MEAthiol_scale;
    _delta = _add > _delta ? _add : _delta;
    _add = abs(CO2_5it - CO2_4it) / CO2_scale;
    _delta = _add > _delta ? _add : _delta;
    _add = abs(yO2P_5it - yO2P_4it) / yO2P_scale;
    _delta = _add > _delta ? _add : _delta;
    _add = abs(Cystamine_5it - Cystamine_4it) / Cystamine_scale;
    _delta = _add > _delta ? _add : _delta;
    _add = abs(VL_5it - VL_4it) / VL_scale;
    _delta = _add > _delta ? _add : _delta;
    double _newStep = 0.8 * step * pow(_tolerance / _delta, _power);
    _increase = _newStep > step ? true : false;
    if (_itB_times != _itE_times && _time + _newStep > *_itB_times) {
      _time = *_itB_times;
    }
    else {
      _time += _newStep;
    }
    if (_time != _time)
    {
      break;
    }
  } // Solve

  //selected variables to memory
  result->push_back(*MEAthiol);
  delete MEAthiol;
  result->push_back(*VL);
  delete VL;
  result->push_back(*klasurface);
  delete klasurface;
  result->push_back(*CO2);
  delete CO2;
  result->push_back(*yO2P);
  delete yO2P;
  result->push_back(*FFox);
  delete FFox;
  result->push_back(*KKox);
  delete KKox;
  result->push_back(*FFred);
  delete FFred;
  result->push_back(*KKred);
  delete KKred;
  result->push_back(*Ffree);
  delete Ffree;
  result->push_back(*Kfree);
  delete Kfree;
  result->push_back(*FKred);
  delete FKred;
  result->push_back(*FKox);
  delete FKox;
  result->push_back(*Cystamine);
  delete Cystamine;
  result->push_back(*tioBalance);
  delete tioBalance;
}

//name: solveFAEnew
//input: double t0 = 0 {units: minutes; caption: initial; category: Time}
//input: double t1 = 10 {units: minutes; caption: final; category: Time}
//input: double h = 0.1 {units: minutes; caption: step; category: Time}
//input: int timesCount
//input: int varsCount
//output: column_list resultCols [new(timesCount, varsCount); 't, time (minutes)'; 'MEAthiol'; 'VL'; 'klasurface'; 'CO2'; 'yO2P'; 'FFox'; 'KKox'; 'FFred'; 'KKred'; 'Ffree'; 'Kfree'; 'FKred'; 'FKox'; 'Cystamine'; 'tioBalance']
//output: dataframe solution [resultCols] {caption: Solution; viewer: Line chart(x: "t, time (minutes)", sharex: "true", multiAxis: "true", yGlobalScale: "true", multiAxisLegendPosition: "RightCenter") | Grid(block: 100) }
//editor: Compute:RichFunctionViewEditor
EMSCRIPTEN_KEEPALIVE
int solveFAEnew(float t0, float t1, float h, 
  int timesCount, int varsCount,
  float * resultCols, int resultRowCount, int resultColCount) noexcept
{
    // times
    vector<double> times(timesCount);

	// initialization of times
    times[0] = t0;
	for (int i = 1; i < timesCount; i++)
		times[i] = times[i - 1] + h;

    // parameters
    int numOfParams = 8;
	vector<double> _parameters(numOfParams);
	_parameters[0] = 1.0;
	_parameters[1] = 100.0;
	_parameters[2] = 0.209;
	_parameters[3] = 8.19;
	_parameters[4] = 7.17 / 32 / 0.209;
	_parameters[5] = 300;
	_parameters[6] = 0.082;
	_parameters[7] = 1.0;

    // numerical solution
	vector<vector<double>> result;

    double _tolerance = 0.000001;

    // solving the system
    Solve(&result, &times, &_parameters, _tolerance);

    for(int i = 0; i < timesCount; i++)
      resultCols[i] = static_cast<float>(times[i]);

    for(int j = 1; j < varsCount; j++)
      for(int i = 0; i < timesCount; i++)
        resultCols[i + j * timesCount] = result[j - 1][i];

    return 0;
} // solveFAE

