#ifndef WSPDH
#define WSPDH
#include "Morton.h"
#include "iofile.h"
#include "float.h"
#include "pquicksort.h"

class Box {
public:
	Box(std::vector<float> pt, float volume) :volume(volume) {
		if (volume != 0) {
			for (int i = 0; i < 3; i++) {
				min[i] = floor(pt[i] / volume)*volume;
				max[i] = ceil(pt[i] / volume)*volume;
			}
		}
		else {
			for (int i = 0; i < 3; i++) {
				min[i] = pt[i];
				max[i] = pt[i];
			}
		}
	};
	float volume;
	float min[3];
	float max[3];
};

struct Point {
	Point(int id) :id(id) {};
	int id;
	Point* par = nullptr;
};

class Forests {
public:
	Forests(int PointsNum);
	void Union(int id1, int id2);
	int Find(int id1);
private:
	std::vector<Point> forests;
};

class Node {
public:
	Node(std::vector<int> pts, float volume, std::vector<float> loc) :pts(pts), volume(volume), loc(loc) {
		left = nullptr;
		right = nullptr;
		par = nullptr;
	};
	Node* left;
	Node* right;
	Node* par;
	std::vector<int> get_pts() {
		return pts;
	};
	int get_pts(int i) {
		return pts[i];
	};
	float get_volume() {
		return volume;
	};
	std::vector<float> get_loc() {
		return loc;
	};
private:
	std::vector<float> loc;
	std::vector<int> pts;
	float volume;
};

struct BCCP {
	int np = -1;
	int nq = -1;
	float lenth = -1.0f;
};

class WSP {
public:
	WSP(Node* P, Node* Q) :P(P), Q(Q) {};
	Node* P;
	Node* Q;
	BCCP bccp;
	bool operator<=(const WSP& wsp2){
		return (this->bccp.lenth<=wsp2.bccp.lenth);
	}
	bool operator>=(const WSP& wsp2){
		return (this->bccp.lenth>=wsp2.bccp.lenth);
	}
	bool operator<(const WSP& wsp2){
		return (this->bccp.lenth<wsp2.bccp.lenth);
	}
	bool operator>(const WSP& wsp2){
		return (this->bccp.lenth>wsp2.bccp.lenth);
	}
	bool operator==(const WSP& wsp2){
		return (this->bccp.lenth==wsp2.bccp.lenth);
	}
};

class comQuadtree {
public:
	comQuadtree(std::vector<std::vector<float>>& points);
	void compute_tree();
	void compute_wspd();
	float Bccp(std::vector<WSP>::iterator wsp);
	void GFKruskal();
	void write_EMST();
	void check_node();
	void check_wspd();

private:
	std::vector<Node> tree;
	std::vector<std::vector<float>> points;
	std::vector<WSP> wspd;
	std::vector<std::vector<int>> edges;
	bool WellSeperated(Box box1, Box box2);
	void WSPD(Node* node1, Node* node2);
	void WSPD(Node* node);
	BCCP Bccp_recurrent(WSP& wsp, BCCP bccp);
	int get_nodeNum(Node* node);
	int get_nodeNum(std::vector<WSP>::iterator i);
	void GFKruskal(std::vector<std::vector<int>>& T, Forests& forests, int beta, std::vector<WSP>& E,std::vector<WSP>::iterator l);
	void Kruskal(std::vector<WSP>& E, std::vector<WSP>::iterator l,std::vector<WSP>::iterator r,std::vector<std::vector<int>>& edges, Forests& forest);
	void Filter(std::vector<WSP>& E, std::vector<WSP>::iterator l,std::vector<WSP>::iterator r,Forests& forest);

	FLAG neutral_le(std::vector<WSP>& vec, std::vector<WSP>::iterator begin, std::vector<WSP>::iterator end,std::vector<WSP>::iterator& left, std::vector<WSP>::iterator& right, int mid, FLAG f, int blocksize);
	std::vector<WSP>::iterator ppartition_le(std::vector<WSP>& vec, std::vector<WSP>::iterator vec_l, std::vector<WSP>::iterator vec_r, int mid);

	FLAG neutral_le(std::vector<WSP>& vec, std::vector<WSP>::iterator begin, std::vector<WSP>::iterator end,std::vector<WSP>::iterator& left, std::vector<WSP>::iterator& right, float mid, FLAG f, int blocksize);
	std::vector<WSP>::iterator ppartition_le(std::vector<WSP>& vec, std::vector<WSP>::iterator vec_l, std::vector<WSP>::iterator vec_r, float mid);

	FLAG neutral_eq(std::vector<WSP>& vec,std::vector<WSP>::iterator begin, std::vector<WSP>::iterator end,std::vector<WSP>::iterator& left, std::vector<WSP>::iterator& right, FLAG f, int blocksize,Forests& forest);
	std::vector<WSP>::iterator ppartition_eq(std::vector<WSP>& vec, std::vector<WSP>::iterator vec_l, std::vector<WSP>::iterator vec_r,Forests& forest);
};

void QuickSort(std::vector<WSP>& wsp, int left, int right);

float Distance(std::vector<float> p, std::vector<float> q);

float BoxDis(Box box1, Box box2);

float BoxDis(Node node1, Node node2);

#endif // !WSPDH
