#ifndef MORTON
#define MORTON
#include "head.h"
#include <bitset>

typedef union {
	float f;
	struct {
		unsigned int mantisa : 23;
		unsigned int exponent : 8;
		unsigned int sign : 1;
	} binary;
} float_cast;

bool MortonCompare(std::vector<float> p, std::vector<float> q);

void QuickSort(std::vector<std::vector<float>>& points, int left, int right);
void QuickSort(std::vector<int>& points, int left, int right);

float Vol(std::vector<float> x1, std::vector<float> x2);
std::vector<float> Vol(std::vector<std::vector<float>> points);



#endif // !MORTON

