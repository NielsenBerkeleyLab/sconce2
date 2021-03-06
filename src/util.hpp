#ifndef UTIL_HPP
#define UTIL_HPP

#include <iostream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <regex>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

typedef boost::minstd_rand base_generator_type;

void printMatrix(gsl_matrix* mat);
void printMatrix(FILE* stream, gsl_matrix* mat);
void printRowVector(gsl_vector* vec);
void printRowVector(FILE* stream, gsl_vector* vec);
void printColVector(gsl_vector* vec);
void printColVector(FILE* stream, gsl_vector* vec);
void createDiagFromCol(gsl_matrix* dest, int colNum, gsl_matrix* src);
void raiseMatrixToPower(gsl_matrix* dest, int power, gsl_matrix* src);
double matrixExponential(gsl_matrix* destP, gsl_matrix* srcQ, double t, gsl_eigen_nonsymmv_workspace* eigenWorkspace, gsl_vector_complex* eval, gsl_matrix_complex* evec, gsl_matrix* realEvecMat, gsl_matrix* diagMat, gsl_matrix* LUdecompMat, gsl_permutation* perm, gsl_matrix* inverseMat);
void vectorAbsoluteValue(gsl_vector* vec);
std::vector<std::string>* createStateSpace(int numCells, int maxPloidy);
std::vector<int>* parseStateToInt(std::string state);
bool compareDoubles(double depth1, double depth2, double epsilon = 1e-6);
std::string parseSampleName(std::string filename);
double calcChiSqOfMatrix(gsl_matrix* expected, gsl_matrix* observed);
double calcSqDiffOfMatrix(gsl_matrix* expected, gsl_matrix* observed);
bool checkUniformMatrix(gsl_matrix* mat, double epsilon = 1e-8);
bool checkNearDouble(double val1, double val2, double epsilon = 0.1);
void conflicting_options(const boost::program_options::variables_map& vm, const char* opt1, const char* opt2);
double calcEuclideanDistOfVectors(std::vector<int>* p, std::vector<int>* q);
std::vector<int>* genUniqueRandomSeqInRange(int lower, int upper, int exclude, base_generator_type* generator);
void deleteChrWindowIdxLineNumMaps(std::pair<std::unordered_map<std::string, std::vector<long long int>*>*, std::unordered_map<std::string, std::vector<int>*>*>* chrWindowIdxLineNumMaps);

#endif

