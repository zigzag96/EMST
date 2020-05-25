#ifndef MORTON
#define MORTON
#include <bitset>
#include "head.h"

typedef union {
	float f;
	struct {
		unsigned int mantisa : 23;
		unsigned int exponent : 8;
		unsigned int sign : 1;
	} binary;
} float_cast;

bool MortonCompare(std::vector<float> p, std::vector<float> q);

bool operator<=(const std::vector<float>& v1,const std::vector<float>& v2);
bool operator>=(const std::vector<float>& v1,const std::vector<float>& v2);
bool operator<(const std::vector<float>& v1,const std::vector<float>& v2);
bool operator>(const std::vector<float>& v1,const std::vector<float>& v2);
bool operator==(const std::vector<float>& v1,const std::vector<float>& v2);

float Vol(std::vector<float> x1, std::vector<float> x2);
std::vector<float> Vol(std::vector<std::vector<float>> points);



#endif // !MORTON

