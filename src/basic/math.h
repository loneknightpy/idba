/**
 * @file math.h
 * @brief Basic math functions.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.10
 * @date 2012-04-15
 */

#ifndef __BASIC_MATH_H_

#define __BASIC_MATH_H_

#include <vector>

double NormalRand();
double NormalCDF(double x);
void GenerateNormalCDFTable(int times, int unit, std::vector<double> &table);

#endif

