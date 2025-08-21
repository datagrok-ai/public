#include "updated-sampler.cpp"
#include <iostream>
#include <fstream>
#include <chrono>

const double TOL = 1e-6;

const double DEFAULT_LOW_BOUND = -10;
const double DEFAULT_UPPER_BOUND = 10;

const unsigned SEED = 42;
const unsigned THINNING = 50;
const int NPROJ = 8;
const int WARM_STEPS = 2000;
const double FEAS_TOL = 1e-8;
const int FEAS_MAX_ITER = 4000;

double samplerRandomTest(const int metabolitesCount = 2, const int reactionsCount = 4, const int samplesCount = 1000) {
	Eigen::VectorXd beq(metabolitesCount);
	beq.setConstant(0.0);

	Eigen::VectorXd lb(reactionsCount);
	lb.setConstant(DEFAULT_LOW_BOUND);

	Eigen::VectorXd ub(reactionsCount);
	ub.setConstant(DEFAULT_UPPER_BOUND);

	//Eigen::MatrixXd S(metabolitesCount, reactionsCount);
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> S(metabolitesCount, reactionsCount);
	//S << -2, 1, 0, 1, 0, 1, -2, 1;
	S.setRandom(); // [-1, 1]

	//std::cout << "S=\n" << S << "\n";

	LinearPolytope P;
	P.A = Eigen::MatrixXd();
	P.b = Eigen::VectorXd();
	P.Aeq = S; P.beq = beq;
	P.lb = lb; P.ub = ub;

	OptGPOptions opt;
	opt.seed = SEED;
	opt.thinning = THINNING;
	opt.nproj = NPROJ;
	opt.warmupSteps = WARM_STEPS;
	opt.feasTol = FEAS_TOL;
	opt.feasMaxIter = FEAS_MAX_ITER;

	OptGPSampler sampler(P, opt);

	auto start = std::chrono::high_resolution_clock::now();

	auto samples = sampler.sample(samplesCount);
	
	auto end = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double, std::milli> elapsed = end - start; // In milliseconds
	
	/*for (int i = 0; i < samplesCount; ++i)
		std::cout << i << "-th: " << samples[i].transpose() << "\n";*/

	double mad = 0;

	for (int i = 0; i < samplesCount; ++i) {
		mad = (S * samples[i]).norm();
		if (mad > TOL) {
			throw std::runtime_error("Invalid sample: " + std::to_string(i) + ", too big deviation: " + std::to_string(mad));
		}

		for (int k = 0; k < reactionsCount; ++k) {
			if ( (samples[i](k) > ub(k)) || (samples[i](k) < lb(k)) ) {
				throw std::runtime_error("Invalid sample: value is out of range");
			}
		}
	}	

	return elapsed.count();
}


double samplerChompact(const int samplesCount = 1000) {
	const int metabolitesCount = 101;
	const int reactionsCount = 144;

	Eigen::VectorXd beq(metabolitesCount);
	beq.setConstant(0.0);

	Eigen::VectorXd lb(reactionsCount);
	//lb.setConstant(DEFAULT_LOW_BOUND);
	//lb.setConstant(-0.0000001);

	double val = 0;

	std::ifstream lbFile("reaction_lower_bounds.txt", std::ios::in);
	for (int i = 0; i < reactionsCount; ++i) {
		lbFile >> val;		
		lb(i) = val < 0.0 ? val : -1e-10;//-0.0000001;
	//	std::cout << " " << lb(i);
	}

	Eigen::VectorXd ub(reactionsCount);
	ub.setConstant(DEFAULT_UPPER_BOUND);
	//ub.setConstant(2.0);

	Eigen::MatrixXd S(metabolitesCount, reactionsCount);
	
	std::ifstream file("chompact.txt", std::ios::in);	

	for (int i = 0; i < metabolitesCount; ++i) {
		for (int j = 0; j < reactionsCount; ++j) {
			file >> val;

			S(i, j) = val;

			//std::cout << " " << val;
		}

		//std::cout << "\n";
	}

	int initOptRowCount = 288;
	int initOptColCount = 144;

	Eigen::MatrixXd InitOpt(initOptRowCount, initOptColCount);

	std::ifstream initOptFile("init-opt.txt", std::ios::in);

	for (int i = 0; i < initOptRowCount; ++i) {
		for (int j = 0; j < initOptColCount; ++j) {
			initOptFile >> val;
			InitOpt(i, j) = val;
		}
	}

	/*std::cout << S * InitOpt.row(0).transpose() << "\n";
	std::cout << InitOpt.transpose().rowwise().mean().size() << "\n";
	std::cout << (S * InitOpt.transpose().rowwise().mean()).norm() << "\n";
	std::cout << "|S * InitOpt_tr| = " << (S * InitOpt.transpose()).norm() << "\n";
	return 0;*/



	//std::cout << "S=\n" << S << "\n";

	//return 0;

	//Eigen::VectorXd x(reactionsCount);

	//// mean
	//// x << 2.88696837413413, 2.82100848365597, 2.79721931959462, 2.51076098626129, 5.30798030585591, 1.65593085470754, 8.29781909375042, 0.070493213472858, 0.399393690745588, 0.0868055555555555, 9.14106825347296, 8.60800609468012, 1.36517010719548, 1.45142386336462, -7.4870821177433, 8.30331575017093, 0.0759898698933729, 2.45193276260344, 2.74764865732274, 0.0694444444444445, -8.86251611121454, 0.416839941682045, 0.225093354979181, -8.41808157122955, -9.10049099589946, -3.09363120287671, 0.0184739731457025, 0.054780627725467, -1.23808106787509, 0.12230214225842, -9.65616263910463, -9.56820511753623, 1.32025712825345, -0.379311803799166, 0.114039081163163, 0.337093444163699, 0.119235592018175, 0.20125358471021, 0.0618872846693813, 0.0346926195800632, 0.173611111111111, 0.0347222222222222, 0.138888888888889, 0.069394276472952, 0.607777230499511, 0.377978384242748, 0.229798846256763, 0.528222214291086, 0.460453857955481, 0.0820179926920349, -9.54894282114226, -9.39626262299088, 0.043402777777774, 0.297095696236722, 3.6504811263257, -8.1396295149746, 0.185917832481748, 0.260589941682053, 0.347395497237601, -0.0845145940387369, 0.00984248483843235, 0.00164041413973873, 0.00765863224965369, 0.00765863224965369, 0.00601936048491425, 0.0247148014191648, 0.0129038196130459, 0.0118109818061188, 0.0300658263063023, 0.0300658263063023, 0.0144394682350605, 0.00313163823774526, 0.002202841844792, 0.00534306319800613, 0.00314022135321413, 0.194406666502026, 0.194374778477652, 0.000641291711640513, 0.104265018948946, 0.0184103131338416, 0.00621346483144762, 0.00467160513281025, -9.55327383997151, 0.00287933879139529, 0.00287933879139529, 0.000759427211332416, 0, 0.00200211655244578, 0.00814808795018311, 0.234344877105532, 2.72025666641734, 0.286458333333334, 9.81692366695694, 0.542780487800658, 0.00120461017076561, 0.0108750982876006, 0.0129140657737514, 0.001568324822668, 0.00226487791947547, 0.00102020276635806, 6.42040759194236E-06, 0.00117197683351213, 9.04810852842805, 0.720486111111111, 2.88696837413413, 9.72067882985945, -9.10049099589946, 0.342096119890277, -4.65124633602191, -3.28275713748243, 0.0312967019755127, -9.4261147726933, 0.0212385248716719, -9.8323985042735, 8.65848990724926, 9.37899953131425, 0.356071728881824, 0.167436616133692, 0.83984278812672, 0.291582147183778, 0.127361240012397, 0.630418271200858, -9.13345183792803, -3.25107807341801, 9.78369465834856, 0.0732629346501976, 0.0437787017516424, 0.407483946975114, 0.312588135190033, 0.759214220056152, 8.08126822968287, 0.444434539984996, 0.138888888888889, 0.229798846256763, 2.72025666641734, 0.656819568963821, -8.91043968150091, 9.81692366695694, -9.55327383997151, 0.0858547058151044, 0.185917832481748, 0.234344877105532, 0.00117197683351213, 1.38310185185185;

	//x << 10, 10, 10, 0, 10, 10, 10, 0, 2E-15, 0, 10, 10, 0, 1.2E-15, -10, 10, 0, -2.8E-14, -5.68E-14, 0, -10, 0, 0, -10, -10, 1.4E-15, 4E-16, -1.2E-15, -4E-16, 0, -10, -10, -1.3E-14, 2.84E-14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.2E-15, 1.2E-15, 0, 0, 0, 0, -10, -10, 0, 0, 0, -10, 0, 0, 0, 2E-16, -2E-16, 0, 0, 0, 0, -4E-16, -4E-16, -2E-16, -6E-16, -6E-16, -4E-16, 0, 0, 0, 0, -3.2E-15, -3.2E-15, 2E-15, 0, -4E-16, 0, 0, -8E-16, 0, 0, 0, 0, 0, -2E-16, -4.6E-15, 0, 10, 10, -2.92E-14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10, -10, 2.4E-15, -10, 6E-15, -1.4E-15, -10, -2.8E-15, -10, 10, 10, -6E-16, -1.2E-15, -1.2E-15, -2.4E-15, -6E-16, -1.2E-15, -10, 2.52E-14, 10, 1.8E-15, -8E-16, -2E-15, 2E-15, -2E-15, 10, 0, 0, 0, 0, -2.92E-14, -3E-14, 10, -8E-16, 4E-16, 0, -4.6E-15, 0, 0;

	//std::cout << "|S * x| = " << (S * x).norm() << "\n";

	//for (int i = 0; i < reactionsCount; ++i)
	//	std::cout << i << "-th condition: " << (x(i) >= lb(i)) << ", " << (x(i) <= ub(i)) << "\n";

	//return 0;

	LinearPolytope P;
	P.A = Eigen::MatrixXd();
	P.b = Eigen::VectorXd();
	P.Aeq = S; P.beq = beq;
	P.lb = lb; P.ub = ub;
	P.InitOpt = InitOpt;

	OptGPOptions opt;
	//opt.seed = 1983;
	opt.seed = SEED;
	opt.thinning = THINNING;
	opt.nproj = NPROJ;
	opt.warmupSteps = WARM_STEPS;
	opt.feasTol = FEAS_TOL;
	opt.feasMaxIter = FEAS_MAX_ITER;
	opt.spreadRatio = 0.5;

	OptGPSampler sampler(P, opt);

	auto start = std::chrono::high_resolution_clock::now();

	auto samples = sampler.sample(samplesCount);

	auto end = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double, std::milli> elapsed = end - start; // In milliseconds	

	//std::ofstream resFile("samples (FullPivLU, seed " + std::to_string(opt.seed) + ", spread " + std::to_string(opt.spreadRatio) + ").csv", std::ios::out);
	std::ofstream resFile("samples.csv", std::ios::out);

	for (int k = 0; k < reactionsCount; ++k)
		resFile << "reaction " << k << ",";
	resFile << "\n";

	for (int i = 0; i < samplesCount; ++i) {
		//std::cout << i << "-th: " << samples[i].transpose() << "\n\n";

		for (int k = 0; k < reactionsCount; ++k)
			resFile << samples[i](k) << ",";

		resFile << "\n";
	}

	double mad = 0;

	for (int i = 0; i < samplesCount; ++i) {
		mad = (S * samples[i]).norm();
		if (mad > TOL) {
			throw std::runtime_error("Invalid sample: " + std::to_string(i) + ", too big deviation: " + std::to_string(mad));
		}

		for (int k = 0; k < reactionsCount; ++k) {
			if ((samples[i](k) > ub(k)) || (samples[i](k) < lb(k))) {
				throw std::runtime_error("Invalid sample: value is out of range");
			}
		}
	}

	for (int i = 0; i < samplesCount; ++i)
		std::cout << "Min: " << samples[i].minCoeff() << "   max: " << samples[i].maxCoeff() << "\n";

	return elapsed.count();
}

/*
   NOTES:
   
     1) in sampler.cpp, the path to the Eigne lib is manually specified (see line 32)
   
     2) command for gcc: 
	     g++ -O3 -std=c++17 sampler.cpp main.cpp -o check-performance
*/

int main() {
	std::cout << samplerChompact()/* samplerRandomTest(10, 20, 1) */ << " ms" << " --- " << "SUCCESS!\n";
	return 0;

	int reactMin = 0;
	int reactMax = 0;

	int metabMin = 0;
	int metabMax = 0;

	int step = 0;

	int samples = 0;

	std::cout << "Reactions min: ";
	std::cin >> reactMin;

	std::cout << "Reactions max: ";
	std::cin >> reactMax;

	std::cout << "Metabolites min: ";
	std::cin >> metabMin;

	std::cout << "Metabolites max: ";
	std::cin >> metabMax;

	std::cout << "Step: ";
	std::cin >> step;

	std::cout << "Samples: ";
	std::cin >> samples;

	std::ofstream file("eigen-" + std::to_string(samples) + "-samples-ms.csv", std::ios::out);

	file << "metabolites\\reactions";

	for (int i = reactMin; i <= reactMax; i += step)
		file << ',' << i;

	file << "\n";

	for (int metabolitesCount = metabMin; metabolitesCount <= metabMax; metabolitesCount += step) {
		file << metabolitesCount;

		for (int reactionsCount = reactMin; reactionsCount <= reactMax; reactionsCount += step) {
			std::cout << metabolitesCount << " x " << reactionsCount << ": ";
			file << ",";

			try {
				double execTime = samplerRandomTest(metabolitesCount, reactionsCount, samples);
				std::cout << execTime << " ms" << " --- " << "SUCCESS!\n";
				file << std::ceil(execTime);
			} catch (const std::runtime_error& e) {
				std::cerr << "Error: " << e.what() << std::endl;				
			}
			catch (...) { // Catch-all for any other unhandled exceptions
				std::cerr << "An unknown error occurred." << std::endl;
			}
		}

		file << "\n";
	}

	std::cout << "\nPress any key...";

	//
	//// Simple chain Sv=0 with bounds 0<=v<=10
	//Eigen::MatrixXd S(2, 4); S << -2, 1, 0, 1, 0, 1, -2, 1; // Aeq
	//Eigen::VectorXd beq(2); beq << 0, 0;
	//Eigen::VectorXd lb(4); lb << -10, -10, -10, -10;
	//Eigen::VectorXd ub(4); ub << 10, 10, 10, 10;

	//LinearPolytope P;
	//P.A = Eigen::MatrixXd();
	//P.b = Eigen::VectorXd();
	//P.Aeq = S; P.beq = beq;
	//P.lb = lb; P.ub = ub;

	//OptGPOptions opt; 
	//opt.seed = 42;
	//opt.thinning = 50;
	//opt.nproj = 8;
	//opt.warmupSteps = 2000;
	//opt.feasTol = 1e-8;
	//opt.feasMaxIter = 4000;

	//OptGPSampler sampler(P, opt);
	//auto samples = sampler.sample(1000);
	//std::cout << "N=" << samples.size() << " first sample: " << samples[0].transpose() << "\n";

	//std::cout << S * samples[0];

	std::cin.get();
}
