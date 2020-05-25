#include "Morton.h"

//return true if p<=q in Morton order
bool MortonCompare(std::vector<float> p, std::vector<float> q)
{
	int diffBit = -150;
	bool result = true;
	for (int i = 0; i < 3; i++) {
		float_cast pi, qi;
		int diffbit = -150;
		pi.f = p[i];
		qi.f = q[i];
		if (pi.binary.exponent == qi.binary.exponent) {
			auto xorMantisa = pi.binary.mantisa ^ qi.binary.mantisa;
			for (int j = 22; j > 0; j--) {
				if (((xorMantisa >> j) & 1) == 1) {
					diffbit = pi.binary.exponent - 127 - (23 - j);
					break;
				}
			}
		}
		else if (pi.binary.exponent < qi.binary.exponent) {
			diffbit = qi.binary.exponent - 127;
		}
		else {
			diffbit = pi.binary.exponent - 127;
		}
		if (diffbit > diffBit) {
			diffBit = diffbit;
			if (pi.f < qi.f) {
				result = true;
			}
			else {
				result = false;
			}
		}
	}
	return result;
}

void QuickSort(std::vector<std::vector<float>>& points, int left, int right) {
	auto mid = points[right];
	int i = left;
	int j = right;
	while (i != j) {
		while (MortonCompare(points[i], mid) && i < j) {
			i++;
		}
		while (MortonCompare(mid, points[j]) && i < j) {
			j--;
		}
		if (i < j) {
			swap(points[i], points[j]);
		}
	}
	swap(points[j], points[right]);
	if (j  - left > 1) {
		QuickSort(points, left, j - 1);
	}
	if (right - j  > 1) {
		QuickSort(points, j + 1, right);
	}
}

void QuickSort(std::vector<int>& points, int left, int right) {
	auto mid = points[right];
	int i = left;
	int j = right;
	while (i != j) {
		while (points[i]<=mid && i < j) {
			i++;
		}
		while (mid<=points[j] && i < j) {
			j--;
		}
		if (i < j) {
			std::swap(points[i], points[j]);
		}
	}
	std::swap(points[j], points[right]);
	if (j  - left > 1) {
		QuickSort(points, left, j - 1);
	}
	if (right - j  > 1) {
		QuickSort(points, j + 1, right);
	}
}

float Vol(std::vector<float> x1, std::vector<float> x2) {
	float volume;
	int diffBit = -150;
	for (int i = 0; i < 3; i++) {
		float_cast pi, qi;
		int diffbit = -150;
		pi.f = x1[i];
		qi.f = x2[i];
		if (pi.binary.exponent == qi.binary.exponent) {
			auto xorMantisa = pi.binary.mantisa ^ qi.binary.mantisa;
			for (int j = 22; j > 0; j--) {
				if (((xorMantisa >> j) & 1) == 1) {
					diffbit = pi.binary.exponent - 127 - (23 - j);
					break;
				}
			}
		}
		else if (pi.binary.exponent < qi.binary.exponent) {
			diffbit = qi.binary.exponent - 127;
		}
		else {
			diffbit = pi.binary.exponent - 127;
		}
		if (diffbit > diffBit) {
			diffBit = diffbit;
		}
	}
	volume = pow(2.0f, diffBit + 1);
	return volume;
}

std::vector<float> Vol(std::vector<std::vector<float>> points) {
	std::vector<float> volumes;
	for (auto i = points.begin(); i != points.end()-1; i++) {
		volumes.push_back(Vol(*i, *(i + 1)));
	}
	return volumes;
}
