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

bool operator<=(const std::vector<float>& v1,const std::vector<float>& v2){
	if(v1.size()==3&&v2.size()==3){
		bool result=MortonCompare(v1,v2);
		return result;
	}
	else{
		std::cerr<<"point size wrong!\n"<<"p1:";
		for(int i =0;i<v1.size();i++){
			std::cerr<<v1[i];
		}
		std::cerr<<"\np2:";
		for(int i =0;i<v2.size();i++){
			std::cerr<<v2[i];
		}
		std::cerr<<std::endl;
		exit(1);
	}
}

bool operator>=(const std::vector<float>& v1,const std::vector<float>& v2){
	if(v1.size()==3&&v2.size()==3){
		bool result=MortonCompare(v2,v1);
		return result;
	}
	else{
		std::cerr<<"point size wrong!\n"<<"p1:";
		for(int i =0;i<v1.size();i++){
			std::cerr<<v1[i];
		}
		std::cerr<<"\np2:";
		for(int i =0;i<v2.size();i++){
			std::cerr<<v2[i];
		}
		std::cerr<<std::endl;
		exit(1);
	}
}

bool operator<(const std::vector<float>& v1,const std::vector<float>& v2){
	if(v1.size()==3&&v2.size()==3){
		bool result=(!MortonCompare(v2,v1));
		return result;
	}
	else{
		std::cerr<<"point size wrong!\n"<<"p1:";
		for(int i =0;i<v1.size();i++){
			std::cerr<<v1[i];
		}
		std::cerr<<"\np2:";
		for(int i =0;i<v2.size();i++){
			std::cerr<<v2[i];
		}
		std::cerr<<std::endl;
		exit(1);
	}
}

bool operator>(const std::vector<float>& v1,const std::vector<float>& v2){
	if(v1.size()==3&&v2.size()==3){
		bool result=(!MortonCompare(v1,v2));
		return result;
	}
	else{
		std::cerr<<"point size wrong!\n"<<"p1:";
		for(int i =0;i<v1.size();i++){
			std::cerr<<v1[i];
		}
		std::cerr<<"\np2:";
		for(int i =0;i<v2.size();i++){
			std::cerr<<v2[i];
		}
		std::cerr<<std::endl;
		exit(1);
	}
}

bool operator==(const std::vector<float>& v1,const std::vector<float>& v2){
	if(v1.size()==3&&v2.size()==3){
		bool result=(MortonCompare(v1,v2) && MortonCompare(v2,v1));
		return result;
	}
	else{
		std::cerr<<"point size wrong!\n"<<"p1:";
		for(int i =0;i<v1.size();i++){
			std::cerr<<v1[i];
		}
		std::cerr<<"\np2:";
		for(int i =0;i<v2.size();i++){
			std::cerr<<v2[i];
		}
		std::cerr<<std::endl;
		exit(1);
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
