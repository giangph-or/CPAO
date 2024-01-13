#include "SolutionReader.h"
#include <iomanip>

void Report::create_report(vector<string> instance, vector<string> noPay, vector<string> budget, string nProducts) {
	report_file = nProducts + ".csv";
	ofstream reportFile(report_file, ofstream::out);
	reportFile << "noPay,budget,set,obj,runtime\n";
	for (int n = 0; n < noPay.size(); ++n) {
		for (int b = 0; b < budget.size(); ++b) {
			for (int i = 0; i < instance.size(); ++i) {
				ifstream resultFile;
				string name = instance[i] + "_" + noPay[n] + "_" + budget[b] + ".txt";
				cout << name << endl;
				resultFile.open(name, ofstream::in);
				double obj, time;
				resultFile >> obj >> time;
				if (!resultFile) {
					cout << "Cannot open file " << name << " ..." << endl;
				}
				resultFile.close();
				reportFile << setprecision(5) << fixed << noPay[n] << "," << budget[b] << "," << instance[i] << "," << obj << "," << time << endl;
			}
		}
	}
	reportFile.close();
}