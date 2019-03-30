#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
std::string arond_st(double x, int n){
	std::ostringstream convert;
	convert << std::fixed << std::setprecision(n) << x;
	return convert.str();
}
