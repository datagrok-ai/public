# PopPK Studio — Glossary / Глоссарий

| Term (EN) | Definition (EN) | Значение (RU) |
|-----------|-----------------|---------------|
| ADVAN | NONMEM subroutine defining the PK model structure (e.g., ADVAN2 = one-compartment oral) | Подпрограмма NONMEM, определяющая структуру PK-модели (напр., ADVAN2 = однокомпартментная пероральная) |
| AIC (Akaike Information Criterion) | Model selection metric balancing goodness of fit and complexity; lower is better | Критерий выбора модели, балансирующий качество подгонки и сложность; меньше — лучше |
| AMT | Amount of dose administered (NMTRAN data column) | Количество введённой дозы (колонка данных NMTRAN) |
| AUC (Area Under the Curve) | Total drug exposure over time, calculated as the integral of the concentration-time curve | Общая экспозиция препарата за время, интеграл кривой концентрация-время |
| BID | Twice daily dosing (bis in die) | Двукратный приём в сутки (bis in die) |
| BLQ / BQL (Below Quantification Limit) | Observed concentration below the assay's limit of detection | Наблюдаемая концентрация ниже предела определения аналитического метода |
| Bootstrap | Resampling method: repeatedly re-estimate the model on resampled datasets to obtain parameter confidence intervals | Метод ресэмплинга: повторная оценка модели на пересобранных датасетах для получения доверительных интервалов параметров |
| BPV (Between-Patient Variability) | See IIV | См. IIV |
| CDISC (Clinical Data Interchange Standards Consortium) | Standards for clinical trial data formatting (SDTM, ADaM) | Стандарты форматирования данных клинических исследований (SDTM, ADaM) |
| CL (Clearance) | Volume of plasma completely cleared of drug per unit time | Объём плазмы, полностью очищенный от препарата за единицу времени |
| Cmax | Maximum observed or predicted drug concentration | Максимальная наблюдаемая или предсказанная концентрация препарата |
| CMT | Compartment number in NMTRAN data (identifies which compartment a dose or observation applies to) | Номер компартмента в данных NMTRAN (указывает, к какому компартменту относится доза или наблюдение) |
| Compartmental model | Mathematical model describing the body as interconnected compartments (e.g., central, peripheral) | Математическая модель, описывающая организм как взаимосвязанные компартменты (напр., центральный, периферический) |
| Control stream | NONMEM input file specifying model, data, and estimation settings in NM-TRAN syntax | Входной файл NONMEM, задающий модель, данные и настройки оценки в синтаксисе NM-TRAN |
| Covariate | Patient characteristic (weight, age, GFR, genotype) that may explain variability in PK/PD parameters | Характеристика пациента (вес, возраст, СКФ, генотип), которая может объяснять вариабельность PK/PD параметров |
| CRO (Contract Research Organization) | Company providing outsourced pharmaceutical research services | Контрактная исследовательская организация, оказывающая аутсорсинговые услуги по разработке лекарств |
| CSR (Clinical Study Report) | Comprehensive regulatory document summarizing a clinical trial's design, conduct, and results | Регуляторный документ, резюмирующий дизайн, проведение и результаты клинического исследования |
| Ctrough | Drug concentration at the end of a dosing interval (just before next dose) | Концентрация препарата в конце дозового интервала (непосредственно перед следующей дозой) |
| CWRES (Conditional Weighted Residuals) | Diagnostic metric: standardized residuals conditional on individual parameter estimates; should follow N(0,1) | Диагностическая метрика: стандартизированные остатки, условные на индивидуальные оценки параметров; должны следовать N(0,1) |
| DSL (Domain-Specific Language) | Specialized programming language designed for a specific problem domain | Специализированный язык программирования для конкретной предметной области |
| DV (Dependent Variable) | Observed data value in NMTRAN format (typically drug concentration) | Наблюдаемое значение в формате NMTRAN (обычно концентрация препарата) |
| E-R (Exposure-Response) | Analysis linking drug exposure metrics to efficacy or safety endpoints | Анализ, связывающий метрики экспозиции препарата с конечными точками эффективности или безопасности |
| EBE (Empirical Bayes Estimate) | Individual parameter estimate obtained by combining population information with individual data | Индивидуальная оценка параметра, полученная путём объединения популяционной информации с индивидуальными данными |
| eCTD (Electronic Common Technical Document) | Standardized format for regulatory submissions to FDA/EMA/PMDA | Стандартизированный формат регуляторных подач в FDA/EMA/PMDA |
| eGFR (Estimated Glomerular Filtration Rate) | Kidney function metric used as a covariate in PK models | Показатель функции почек, используемый как ковариата в PK-моделях |
| Emax | Maximum drug effect in the Emax pharmacodynamic model | Максимальный эффект препарата в фармакодинамической модели Emax |
| EOP2 (End of Phase 2) | FDA meeting to discuss Phase 2 results and Phase 3 design | Встреча с FDA для обсуждения результатов Phase 2 и дизайна Phase 3 |
| Eta (η) | Random effect representing individual deviation from the population parameter | Случайный эффект, представляющий индивидуальное отклонение от популяционного параметра |
| EVID | Event ID in NMTRAN data (0=observation, 1=dose, 2=other, 3=reset, 4=reset+dose) | ID события в данных NMTRAN (0=наблюдение, 1=доза, 2=другое, 3=сброс, 4=сброс+доза) |
| Fixed effects (θ) | Population-average model parameters estimated from data | Среднепопуляционные параметры модели, оцениваемые из данных |
| FO (First-Order) | Simplest NLME approximation method; linearizes around η=0 | Простейший метод аппроксимации NLME; линеаризация вокруг η=0 |
| FOCE (First-Order Conditional Estimation) | NLME method linearizing around individual conditional η estimates | Метод NLME с линеаризацией вокруг индивидуальных условных оценок η |
| FOCE-I (FOCE with Interaction) | FOCE accounting for η–ε interaction in the residual variance model | FOCE с учётом взаимодействия η–ε в модели остаточной дисперсии |
| Forest plot | Visualization showing covariate effects on PK parameters with confidence intervals | Визуализация ковариатных эффектов на PK-параметры с доверительными интервалами |
| GOF (Goodness of Fit) | Diagnostic plots assessing how well the model describes observed data | Диагностические графики, оценивающие, насколько хорошо модель описывает наблюдаемые данные |
| GxP | Umbrella term for regulatory quality practices (GLP, GCP, GMP, etc.) | Общий термин для регуляторных практик качества (GLP, GCP, GMP и др.) |
| ICH M15 | International guideline on general principles for Model-Informed Drug Development | Международное руководство по общим принципам разработки лекарств на основе моделей (MIDD) |
| IIV (Inter-Individual Variability) | Between-subject variability in PK/PD parameters, described by the Ω matrix | Межиндивидуальная вариабельность PK/PD параметров, описываемая матрицей Ω |
| IOV (Inter-Occasion Variability) | Within-subject, between-occasion variability in PK/PD parameters | Внутрисубъектная вариабельность PK/PD параметров между визитами/случаями |
| IPRED (Individual Prediction) | Model prediction using individual (post-hoc) parameter estimates | Предсказание модели с использованием индивидуальных (post-hoc) оценок параметров |
| ITT (Intent-To-Treat) | Analysis population including all randomized subjects regardless of compliance | Анализируемая популяция, включающая всех рандомизированных субъектов независимо от комплаентности |
| IWRES (Individual Weighted Residuals) | Residuals scaled by individual-level residual variance estimate | Остатки, масштабированные индивидуальной оценкой остаточной дисперсии |
| Ka (Absorption rate constant) | First-order rate constant for drug absorption from the GI tract | Константа скорости абсорбции первого порядка из ЖКТ |
| MAP (Maximum A Posteriori) | Bayesian estimation of individual parameters using population priors | Байесовская оценка индивидуальных параметров с использованием популяционных априорных распределений |
| MCMC (Markov Chain Monte Carlo) | Sampling-based Bayesian inference method | Метод байесовского вывода, основанный на сэмплировании |
| MDV (Missing Dependent Variable) | Flag in NMTRAN data indicating the observation should not contribute to the likelihood | Флаг в данных NMTRAN, указывающий, что наблюдение не должно вносить вклад в функцию правдоподобия |
| MIDD (Model-Informed Drug Development) | FDA framework for using quantitative models to inform drug development decisions | Фреймворк FDA по использованию количественных моделей для обоснования решений в разработке лекарств |
| MMRM (Mixed Model for Repeated Measures) | Linear mixed-effects model for longitudinal continuous endpoints in clinical trials | Линейная смешанная модель для продольных непрерывных конечных точек в клинических исследованиях |
| M&S (Modeling and Simulation) | Umbrella term for quantitative model-based analysis in drug development | Общий термин для количественного модельного анализа в разработке лекарств |
| NCA (Non-Compartmental Analysis) | Model-independent PK analysis calculating AUC, Cmax, t½ directly from data | Модельно-независимый PK-анализ, рассчитывающий AUC, Cmax, t½ непосредственно из данных |
| NDA (New Drug Application) | Regulatory submission to FDA for marketing approval of a new drug | Регуляторная подача в FDA для получения разрешения на маркетинг нового лекарства |
| NLME (Nonlinear Mixed-Effects Model) | Statistical model combining fixed and random effects with a nonlinear structural model | Статистическая модель, объединяющая фиксированные и случайные эффекты с нелинейной структурной моделью |
| NMTRAN | NONMEM translator: converts user-friendly model code to NONMEM's internal format; also refers to the standard data format | Транслятор NONMEM: преобразует пользовательский код модели; также обозначает стандартный формат данных |
| NONMEM | NONlinear Mixed-Effects Modeling — the original and most widely used pharmacometrics software | NONlinear Mixed-Effects Modeling — оригинальное и наиболее используемое ПО фармакометрики |
| NPDE (Normalized Prediction Distribution Errors) | Simulation-based diagnostic; if model is correct, NPDE ~ N(0,1) | Диагностика на основе симуляции; при корректной модели NPDE ~ N(0,1) |
| ODE (Ordinary Differential Equation) | Mathematical equation describing rate of change; used to define compartmental PK/PD models | Обыкновенное дифференциальное уравнение, описывающее скорость изменения; используется для определения компартментных PK/PD моделей |
| OFV (Objective Function Value) | -2 × log-likelihood; minimized during NLME estimation; used for model comparison (lower = better fit) | -2 × логарифм правдоподобия; минимизируется при NLME-оценке; используется для сравнения моделей (меньше = лучше подгонка) |
| Omega (Ω) | Variance-covariance matrix of inter-individual random effects | Ковариационная матрица межиндивидуальных случайных эффектов |
| pcVPC (Prediction-Corrected VPC) | VPC variant normalizing data by population prediction to remove covariate-driven variability | Вариант VPC, нормализующий данные по популяционному предсказанию для устранения ковариатной вариабельности |
| PD (Pharmacodynamics) | Study of drug effects on the body (what the drug does to the body) | Изучение эффектов препарата на организм (что препарат делает с организмом) |
| PK (Pharmacokinetics) | Study of drug absorption, distribution, metabolism, and excretion (what the body does to the drug) | Изучение абсорбции, распределения, метаболизма и экскреции препарата (что организм делает с препаратом) |
| PopPK (Population Pharmacokinetics) | Application of NLME models to characterize PK variability across a patient population | Применение NLME-моделей для характеристики вариабельности PK в популяции пациентов |
| PRED (Population Prediction) | Model prediction using population-average (typical) parameters (η=0) | Предсказание модели с использованием среднепопуляционных (типичных) параметров (η=0) |
| PREDPP | NONMEM library of pre-programmed PK model subroutines | Библиотека NONMEM с предварительно запрограммированными подпрограммами PK-моделей |
| PsN (Perl-speaks-NONMEM) | Automation toolkit for NONMEM providing bootstrap, VPC, SCM, and cluster computing | Инструментарий автоматизации NONMEM для bootstrap, VPC, SCM и кластерных вычислений |
| QD | Once daily dosing (quaque die) | Однократный приём в сутки (quaque die) |
| QSP (Quantitative Systems Pharmacology) | Mechanistic modeling approach integrating systems biology with pharmacology | Механистический подход к моделированию, объединяющий системную биологию с фармакологией |
| Random effects (η, ε) | Model components describing unexplained variability (between-subject and residual) | Компоненты модели, описывающие необъяснённую вариабельность (межсубъектную и остаточную) |
| REML (Restricted Maximum Likelihood) | Estimation method for variance parameters that accounts for degrees of freedom lost in estimating fixed effects | Метод оценки параметров дисперсии с учётом степеней свободы, потраченных на оценку фиксированных эффектов |
| SAEM (Stochastic Approximation EM) | NLME estimation via iterating between MCMC sampling of individual parameters and closed-form population updates | Оценка NLME через итерации MCMC-сэмплирования индивидуальных параметров и замкнутых обновлений популяционных оценок |
| SAP (Statistical Analysis Plan) | Pre-specified document describing all planned analyses for a clinical trial | Предварительно определённый документ, описывающий все запланированные анализы клинического исследования |
| SCM (Stepwise Covariate Modeling) | Automated covariate selection via forward inclusion then backward elimination based on OFV changes | Автоматический отбор ковариат прямым включением и обратным исключением на основе изменений OFV |
| SDTM (Study Data Tabulation Model) | CDISC standard for organizing clinical trial raw data | Стандарт CDISC для организации сырых данных клинических исследований |
| Shrinkage | Degree to which individual parameter estimates are pulled toward the population mean; high shrinkage (>30%) reduces reliability of individual-level diagnostics | Степень «притяжения» индивидуальных оценок параметров к популяционному среднему; высокий shrinkage (>30%) снижает надёжность диагностики на индивидуальном уровне |
| Sigma (Σ) | Parameter(s) describing residual unexplained variability | Параметр(ы), описывающие остаточную необъяснённую вариабельность |
| Stiff ODE | ODE system with widely separated time scales requiring specialized solvers (BDF, Rosenbrock) | Система ОДУ с сильно различающимися временными масштабами, требующая специализированных солверов (BDF, Rosenbrock) |
| TDM (Therapeutic Drug Monitoring) | Clinical practice of measuring drug concentrations to optimize individual dosing | Клиническая практика измерения концентраций препарата для оптимизации индивидуальной дозировки |
| Theta (θ) | Fixed-effect population parameter (e.g., typical CL, typical V) | Параметр фиксированного эффекта популяции (напр., типичный CL, типичный V) |
| TMDD (Target-Mediated Drug Disposition) | PK model for drugs whose disposition is significantly affected by binding to their pharmacological target | PK-модель для препаратов, чья диспозиция существенно зависит от связывания с фармакологической мишенью |
| TRANS | NONMEM parameterization translator (e.g., TRANS2: CL,V; TRANS4: CL,V1,Q,V2) | Транслятор параметризации NONMEM (напр., TRANS2: CL,V; TRANS4: CL,V1,Q,V2) |
| V (Volume of distribution) | Apparent volume relating drug amount in the body to plasma concentration | Кажущийся объём, связывающий количество препарата в организме с концентрацией в плазме |
| VPC (Visual Predictive Check) | Model qualification plot comparing observed data percentiles against simulated prediction intervals | График квалификации модели, сравнивающий перцентили наблюдаемых данных с интервалами симулированных предсказаний |
| WASM (WebAssembly) | Binary instruction format enabling near-native performance in web browsers | Бинарный формат инструкций, обеспечивающий производительность, близкую к нативной, в веб-браузерах |
