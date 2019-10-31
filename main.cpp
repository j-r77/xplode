#include <fstream>
#include <iostream>
#include "util/output.h"
#include "data/databasereader.h"
#include "algorithms/xplode.h"

int main(int argc, char *argv[]) {
    if (argc != 5) {
        std::cout << "Usage: ./XPlode dirtyfile cleanfile epsilon delta" << std::endl;
		std::cout << "\t where dirtyfile is a dirty dataset in csv format," << std::endl;
		std::cout << "\t cleanfile is a (partially) cleaned version of this dataset," << std::endl;
		std::cout << "\t and epsilon and delta are confidence and support parameters for the underlying rule mining algorithm"
        std::cout << "\t For Example: ./XPlode abalone-dirty.csv abalone-cleaned.csv 0.1 100" << std::endl;
    }
    else {
        std::ifstream dirtyFile(argv[1]);
        Database dirtyDb = DatabaseReader::fromTable(dirtyFile, ',');
        std::ifstream cleanFile(argv[2]);
        Database cleanDb = DatabaseReader::fromTable(cleanFile, dirtyDb, ',');
        double epsilon = atof(argv[3]);
        int delta = atoi(argv[4]);
        XPlode xplode(dirtyDb, cleanDb);
        CFD expl = xplode.explain(delta, 1.0-epsilon);
        std::cout << "Best Explanation: " << std::endl;
        Output::printCFD(expl, cleanDb);
    }
    return 0;
}
