#include "wspd.h"

Forests::Forests(int pointsNum)
{
	for (int i = 0; i < pointsNum; i++) {
		forests.push_back(Point(i));
	}
}

void Forests::Union(int id1, int id2)
{
	Point* p1 = &(forests[id1]);
	while (p1->par != nullptr) {
		p1 = p1->par;
	}
	Point* p2 = &(forests[id2]);
	while (p2->par != nullptr) {
		p2 = p2->par;
	}
	if (p1 != p2) {
		p1->par = p2;
	}
}

int Forests::Find(int id1)
{
	Point* p = &(forests[id1]);
	while (p->par != nullptr) {
		p = p->par;
	}
	return p->id;
}

comQuadtree::comQuadtree(std::vector<std::vector<float>>& points)
{
	QuickSort(points, 0, points.size() - 1);
	this->points = points;
	tree.push_back(Node({}, FLT_MAX, {}));
	for (int i = 0; i < points.size() - 1; i++) {
		tree.push_back(Node({ i,i+1 },Vol(points[i],points[i+1]),points[i]));
	}
	for (int i = 0; i < points.size(); i++) {
		tree.push_back(Node({ i },0.0f, { points[i] }));
	}
}

void comQuadtree::compute_tree() {
	tree[0].right = &tree[points.size()];
	Node* k = &(tree[0]);
	for (int i = 0; i < points.size() - 1; i++) {
		while (tree[i + 1].get_volume() > k->get_volume()) {
			k = k->par;
		}
		tree[i + 1].left = k->right;
		k->right->par = &(tree[i + 1]);
		tree[i + 1].right = &(tree[points.size() + i + 1]);
		tree[points.size() + i + 1].par = &(tree[i + 1]);
		k->right = &(tree[i + 1]);
		tree[i + 1].par = k;
		k = &(tree[i + 1]);
	}
	std::cout << "compute_tree done" << std::endl;
}

void comQuadtree::compute_wspd() {
	wspd = WSPD(tree[0].right);
	std::cout << "compute_wspd done" << std::endl;
}

void comQuadtree::Bccp(WSP& wsp)
{
	if (wsp.bccp.lenth != -1) {
		return ;
	}
	BCCP bccp;
	{
		auto iP = wsp.P->get_pts();
		auto iQ = wsp.Q->get_pts();
		if (iP[0] > iQ[0]) {
			std::swap(wsp.P, wsp.Q);
		}//P<=Q
		auto pright = wsp.P;
		while (pright->right != nullptr) {
			pright = pright->right;
		}
		bccp.np = pright->get_pts(0);
		auto pleft = wsp.Q;
		while (pleft->left != nullptr) {
			pleft = pleft->left;
		}
		bccp.nq = pleft->get_pts(0);
		bccp.lenth = Distance(points[bccp.np], points[bccp.nq]);
	}
	bccp=Bccp_recurrent(wsp,bccp);
	wsp.bccp = bccp;
	return ;
}

void comQuadtree::GFKruskal()
{
	Forests forests(points.size());
	GFKruskal(edges, forests, 2, wspd);
	std::cout << "EMST done" << std::endl;
}

void comQuadtree::write_EMST()
{
	draw_ply(points, edges);
}

void comQuadtree::check_node()
{
	std::ofstream ptfile;
	ptfile.open("node.txt", std::ios::out | std::ios::trunc);
	for (auto i = tree.begin(); i != tree.end(); i++) {
		for (auto j : i->get_pts()) {
			ptfile << j << " ";
		}
		ptfile << ":";
		if (i->left != nullptr) {
			for (auto j : i->left->get_pts()) {
				ptfile << j << " ";
			}
		}
		if (i->right != nullptr) {
			ptfile << "/";
			for (auto j : i->right->get_pts()) {
				ptfile << j << " ";
			}
		}
		ptfile << std::endl;
	}
	ptfile.close();
}

bool comQuadtree::WellSeperated(Box box1,Box box2) {
	float diameter = std::max(box1.volume*sqrt(3),box2.volume*sqrt(3));
	float dis = BoxDis(box1,box2);
	if (diameter <= dis) {
		return true;
	}
	return false;
}

std::vector<WSP> comQuadtree::WSPD(Node* node1, Node* node2)
{
	Box box1(node1->get_loc(),node1->get_volume());
	Box box2(node2->get_loc(), node2->get_volume());
	if (WellSeperated(box1, box2)) {
		WSP wsp(node1, node2);
		return { wsp };
	}
	else if (box1.volume >= box2.volume) {
		std::vector<WSP> wsp1 = WSPD(node1->left, node2);
		std::vector<WSP> wsp2 = WSPD(node1->right, node2);
		wsp1.insert(wsp1.end(), wsp2.begin(), wsp2.end());
		return wsp1;
	}
	else {
		std::vector<WSP> wsp1 = WSPD(node1, node2->left);
		std::vector<WSP> wsp2 = WSPD(node1, node2->right);
		wsp1.insert(wsp1.end(), wsp2.begin(), wsp2.end());
		return wsp1;
	}
}


std::vector<WSP> comQuadtree::WSPD(Node* node)
{
	std::vector<WSP> result;
	if (node->right == nullptr) {
		return result;
	}
	auto lwsp = WSPD(node->left);
	auto rwsp = WSPD(node->right);
	auto pwsp = WSPD(node->left, node->right);
	result.insert(result.end(), lwsp.begin(), lwsp.end());
	result.insert(result.end(), rwsp.begin(), rwsp.end());
	result.insert(result.end(), pwsp.begin(), pwsp.end());
	return result;
}

BCCP comQuadtree::Bccp_recurrent(WSP& wsp,BCCP bccp)
{
	if (wsp.P->left == nullptr && wsp.Q->left == nullptr) {
		float dis = Distance(points[wsp.P->get_pts(0)], points[wsp.Q->get_pts(0)]);
		if (dis <= bccp.lenth ) {
			bccp.np = wsp.P->get_pts(0);
			bccp.nq = wsp.Q->get_pts(0);
			bccp.lenth = dis;
			if (bccp.lenth == -1) {
				std::cout << "error" << std::endl;
			}
			return bccp;
		}
	}
	else {
		if (wsp.P->get_volume() < wsp.Q->get_volume()) {
			std::swap(wsp.P, wsp.Q);
		}
		Box boxpL(wsp.P->left->get_loc(),wsp.P->left->get_volume());
		Box boxpR(wsp.P->right->get_loc(), wsp.P->right->get_volume());
		Box boxq(wsp.Q->get_loc(), wsp.Q->get_volume());
		auto d1 = BoxDis(boxpL, boxq);
		auto d2 = BoxDis(boxpR, boxq);
		BCCP bccp1, bccp2;
		if (d1 <= bccp.lenth) {
			WSP wsp1(wsp.P->left,wsp.Q);
			bccp1 = Bccp_recurrent(wsp1, bccp);
		}
		if (d2 <= bccp.lenth) {
			WSP wsp2(wsp.P->right, wsp.Q);
			bccp2 = Bccp_recurrent(wsp2, bccp);
		}
		//if (d1 > bccp.lenth && d2 > bccp.lenth) {
		//	std::cout << "error" << std::endl;
		//}
		if (bccp1.lenth != -1) {
			if (bccp2.lenth != -1) {
				if (bccp1.lenth <= bccp2.lenth && bccp1.lenth<=bccp.lenth) {
					bccp = bccp1;

				}
				else {
					if (bccp2.lenth <= bccp.lenth) {
						bccp = bccp2;
					}
				}
			}
			else {
				if (bccp1.lenth <= bccp.lenth&&bccp1.lenth != -1) {
					bccp = bccp1;
				}
			}
		}
		else {
			if (bccp2.lenth <= bccp.lenth && bccp2.lenth!=-1) {
				bccp = bccp2;
			}
		}
		if (bccp.lenth == -1) {
			std::cout << "error" << std::endl;
		}
	}
	return bccp;
}

int comQuadtree::get_nodeNum(Node* node)
{
	Node* nleft=node;
	Node* nright = node;
	while (nleft->left != nullptr) {
		nleft = nleft->left;
	}
	while (nright->right != nullptr) {
		nright = nright->right;
	}
	int Num = nright->get_pts(0) - nleft->get_pts(0) + 1;
	return Num;
}

void comQuadtree::Kruskal(std::vector<WSP>& E, std::vector<std::vector<int>>& edges, Forests & forest)
{
	QuickSort(E, 0, E.size()-1);
	for (auto i = E.begin(); i != E.end(); i++) {
		if ((*i).bccp.np == -1) {
			std::cout << "error" << std::endl;
		}
		if (forest.Find((*i).bccp.np) != forest.Find((*i).bccp.nq)) {
			edges.push_back({ (*i).bccp.np ,(*i).bccp.nq });
			forest.Union((*i).bccp.np, (*i).bccp.nq);
		}
	}
}

void comQuadtree::Filter(std::vector<WSP>& E, Forests& forest)
{
	auto i = E.begin();
	while(i!=E.end()) {
		if (forest.Find(i->P->get_pts(0)) == forest.Find(i->Q->get_pts(0))){
			i=E.erase(i);
		}
		else {
			i++;
		}
	}
}

void comQuadtree::check_wspd()
{
	std::ofstream ptfile;
	ptfile.open("wspd.txt", std::ios::out | std::ios::trunc);
	for (auto i = wspd.begin(); i != wspd.end(); i++) {
		for (auto j : i->P->get_pts()) {
			ptfile << "P:" << j;
		}
		for (auto j : i->Q->get_pts()) {
			ptfile << " Q:" << j;
		}
		ptfile << std::endl;
	}
	ptfile.close();
}

void comQuadtree::GFKruskal(std::vector<std::vector<int>>& T, Forests & forests, int beta, std::vector<WSP>& E)
{
	std::vector<WSP> El,Eu,El1,El2;
	for (auto i = E.begin(); i != E.end(); i++) {
		if (get_nodeNum(i->P) + get_nodeNum(i->Q) <= beta) {
			El.push_back(*i);
		}
		else {
			Eu.push_back(*i);
		}
	}
	float ro=FLT_MAX;
	for (auto i = Eu.begin(); i != Eu.end(); i++) {
		ro = std::min(ro, BoxDis(*(i->P), *(i->Q)));
	}
	for (auto i = El.begin(); i != El.end(); i++) {
		Bccp(*i);
		if (i->bccp.lenth <= ro) {
			El1.push_back(*i);
		}
		else {
			El2.push_back(*i);
		}
	}
	Kruskal(El1, T, forests);
	El2.insert(El2.begin(), Eu.begin(), Eu.end());
	Filter(El2, forests);
	if (T.size() < points.size() - 1) {
		GFKruskal(T, forests, beta + 1, El2);
	}

}

void QuickSort(std::vector<WSP>& wsp,int left,int right)
{
	if(left>=right){
		return;
	}
	auto mid = wsp[right].bccp.lenth;
	int i = left;
	int j = right;
	while (i != j) {
		while (wsp[i].bccp.lenth <= mid && i < j) {
			i++;
		}
		while (mid <= wsp[j].bccp.lenth && i < j) {
			j--;
		}
		if (i < j) {
			std::swap(wsp[i], wsp[j]);
		}
	}
	std::swap(wsp[j], wsp[right]);
	if (j - left > 1) {
		QuickSort(wsp, left, j - 1);
	}
	if (right - j  > 1) {
		QuickSort(wsp, j + 1, right);
	}
}

float Distance(std::vector<float> p, std::vector<float> q) {
	float result = 0.0f;
	for (int i = 0; i < 3; i++) {
		result += pow(p[i] - q[i], 2);
	}
	result = sqrt(result);
	return result;
}

float BoxDis(Box box1, Box box2) {
	float d[3];
	for (int i = 0; i < 3; i++) {
		d[i] = pow(box1.min[i] - box2.min[i], 2);
		d[i] = std::min(d[i], static_cast<float>(pow(box1.min[i] - box2.max[i], 2)));
		d[i] = std::min(d[i], static_cast<float>(pow(box1.max[i] - box2.max[i], 2)));
		d[i] = std::min(d[i], static_cast<float>(pow(box1.max[i] - box2.min[i], 2)));
	}
	float dis = sqrt(d[0] + d[1] + d[2]);
	return dis;
}

float BoxDis(Node node1, Node node2) {
	Box box1(node1.get_loc(), node1.get_volume());
	Box box2(node2.get_loc(), node2.get_volume());
	float dis = BoxDis(box1, box2);
	return dis;
}
