#pragma once
#include"wspd.h"

struct Point {
	Point(int id) :id(id) {};
	int id;
	Point* par=nullptr;
};

class Forests {
public:
	Forests(int PointsNum);
	void Union(int id1, int id2);
	int Find(int id1);
private:
	std::vector<Point> forests;
};

void GFKruskal(std::vector<WSP>& E, std::vector < std::vector<int>>& T, Forests& forests,int beta) {
	std::vector<WSP> El, Eu, El1, El2;
	for (auto i = E.begin(); i != E.end(); i++) {
		if()
	}
}