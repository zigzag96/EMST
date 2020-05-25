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
	PQuickSort(points);
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
	#pragma omp parallel
	{
		#pragma omp single nowait
		{
			WSPD(tree[0].right);
		}
	}
	std::cout << "compute_wspd done" << std::endl;
}

float comQuadtree::Bccp(std::vector<WSP>::iterator wsp)
{
	if (wsp->bccp.lenth != -1) {
		return wsp->bccp.lenth;
	}
	BCCP bccp;
	{
		auto iP = wsp->P->get_pts();
		auto iQ = wsp->Q->get_pts();
		if (iP[0] > iQ[0]) {
			std::swap(wsp->P, wsp->Q);
		}//P<=Q
		auto pright = wsp->P;
		while (pright->right != nullptr) {
			pright = pright->right;
		}
		bccp.np = pright->get_pts(0);
		auto pleft = wsp->Q;
		while (pleft->left != nullptr) {
			pleft = pleft->left;
		}
		bccp.nq = pleft->get_pts(0);
		bccp.lenth = Distance(points[bccp.np], points[bccp.nq]);
	}
	bccp=Bccp_recurrent(*wsp,bccp);
	wsp->bccp = bccp;
	return wsp->bccp.lenth;
}

void comQuadtree::GFKruskal()
{
	Forests forests(points.size());
	GFKruskal(edges, forests, 2, wspd,wspd.begin());
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

void comQuadtree::WSPD(Node* node1, Node* node2)
{
	Box box1(node1->get_loc(),node1->get_volume());
	Box box2(node2->get_loc(), node2->get_volume());
	if (WellSeperated(box1, box2)) {
		WSP wsp(node1, node2);
		#pragma omp critical
		{
			wspd.push_back(wsp);
		}
	}
	else if (box1.volume >= box2.volume) {
		#pragma omp task
		{
			WSPD(node1->left, node2);
		}
		#pragma omp task
		{
			WSPD(node1->right, node2);
		}
	}
	else {
		#pragma omp task
		{
			WSPD(node1, node2->left);
		}
		#pragma omp task
		{
			WSPD(node1, node2->right);
		}
	}
}


void comQuadtree::WSPD(Node* node)
{
	if (node->right == nullptr) {
		return ;
	}
	#pragma omp task
	{
		WSPD(node->left);
	}
	#pragma omp task
	{
		WSPD(node->right);
	}
	#pragma omp task
	{
		WSPD(node->left, node->right);
	}
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

int comQuadtree::get_nodeNum(std::vector<WSP>::iterator i){
	return get_nodeNum(i->P)+get_nodeNum(i->Q);
}

void comQuadtree::Kruskal(std::vector<WSP>& E, std::vector<WSP>::iterator l,std::vector<WSP>::iterator r,std::vector<std::vector<int>>& edges, Forests & forest)
{
	PQuickSort(E, l, r,true);
	for (auto i = l; i <= r; i++) {
		if ((*i).bccp.np == -1) {
			std::cout << "error" << std::endl;
		}
		if (forest.Find((*i).bccp.np) != forest.Find((*i).bccp.nq)) {
			edges.push_back({ (*i).bccp.np ,(*i).bccp.nq });
			forest.Union((*i).bccp.np, (*i).bccp.nq);
		}
	}
}

void comQuadtree::Filter(std::vector<WSP>& E,std::vector<WSP>::iterator l,std::vector<WSP>::iterator r, Forests& forest)
{
	auto i = ppartition_eq(E,l,r,forest);
	E.erase(i,E.end());
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

float min_ro(std::vector<WSP>& E,std::vector<WSP>::iterator l,std::vector<WSP>::iterator r){
	auto Num_threads = omp_get_max_threads();
	int step = ((r-l+1)-((r-l+1)%Num_threads))/Num_threads;
	float result=FLT_MAX;
	if(step!=0){
		float ro[Num_threads];
		#pragma omp parallel num_threads(Num_threads)
		{
			ro[omp_get_thread_num()] = FLT_MAX;
			for(auto i=l+step*omp_get_thread_num();i<l+step*omp_get_thread_num()+step;i++){
				ro[omp_get_thread_num()] = std::min(ro[omp_get_thread_num()], BoxDis(*(i->P), *(i->Q)));
			}
		}
		for(auto i=l+step*Num_threads;i<=r;i++){
			ro[0]=std::min(ro[0], BoxDis(*(i->P), *(i->Q)));
		}
		result=ro[0];
		for (int i =1;i<Num_threads;i++){
			result=std::min(result, ro[i]);
		}
	}
	else{
		for(auto i =l;i<=r;i++){
			result= std::min(result, BoxDis(*(i->P), *(i->Q)));
		}
	}
	return result;
}

void comQuadtree::GFKruskal(std::vector<std::vector<int>>& T, Forests & forests, int beta, std::vector<WSP>& E, std::vector<WSP>::iterator l)
{
	auto Eu = ppartition_le(E,l,E.end()-1,beta);
	auto ro = min_ro(E,Eu,E.end()-1);
	auto El2= ppartition_le(E,l,Eu-1,ro);
	Kruskal(E,l,El2-1, T, forests);
	Filter(E,El2,E.end()-1,forests);
	if (T.size() < points.size() - 1) {
		GFKruskal(T, forests, beta + 1, E,El2);
	}

}

void QuickSort(std::vector<WSP>& wsp,int left,int right)
{
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

FLAG comQuadtree::neutral_le(std::vector<WSP>& vec,std::vector<WSP>::iterator begin, std::vector<WSP>::iterator end, std::vector<WSP>::iterator& left, std::vector<WSP>::iterator& right, int mid, FLAG f, int blocksize)
{
	auto pleft = left;
	auto pright = right;
	auto left_end = left;
	auto right_end = right;
	switch (f)
	{
	case leftf:
		left_end = left + blocksize;
		right_end = right - blocksize + ((end - right) % blocksize);
		break;
	case rightf:
		left_end = left + blocksize - ((left - begin) % blocksize);
		right_end = right - blocksize;
		break;
	case bothf:
		left_end = left + blocksize;
		right_end = right - blocksize;
		break;
	}
	f = nonef;
	while (f == nonef)
	{
		while (pleft != left_end && get_nodeNum(pleft) <= mid) {
			pleft++;
		}
		while (pright != right_end && get_nodeNum(pright) > mid) {
			pright--;
		}
		if (pleft == left_end) {
			if (pright == right_end) {
				f = bothf;
			}
			else {
				f = leftf;
			}
		}
		else if (pright == right_end) {
			f = rightf;
		}
		else {
			std::swap(*pleft, *pright);
		}
	}
	left = pleft;
	right = pright;
	return f;
}

std::vector<WSP>::iterator comQuadtree::ppartition_le(std::vector<WSP>& vec, std::vector<WSP>::iterator vec_l, std::vector<WSP>::iterator vec_r, int mid)
{
	int blocksize = 64;
	auto begin = vec_l;
	auto end =vec_r;
	std::vector<typename std::vector<WSP>::iterator> left_list, right_list;
	int n_threads = std::min(omp_get_num_threads(), static_cast<int>(floor(vec.size() / (2 * blocksize))));
	if (n_threads > 1) {
#pragma omp parallel shared(vec,vec_l,vec_r,mid,left_list,right_list) num_threads(n_threads)
		{
			auto _left = begin + blocksize * omp_get_thread_num();
			auto _right = end - blocksize * omp_get_thread_num();
			vec_l = vec_l + blocksize;
			vec_r = vec_r - blocksize;
			FLAG _f = bothf;
#pragma omp barrier
			while (_f != nonef)
			{
				_f = neutral_le(vec,begin,end, _left, _right, mid, _f, blocksize);
				switch (_f)
				{
				case leftf:
#pragma omp critical
				{
					if (vec_l + blocksize - 1 <= vec_r) {
						_left = vec_l;
						vec_l = vec_l + blocksize;
					}
					else {
						right_list.push_back(_right);
						_f = nonef;
					}
				}
				break;
				case rightf:
#pragma omp critical
				{
					if (vec_r - blocksize + 1 >= vec_r) {
						_right = vec_r;
						vec_r = vec_r - blocksize;
					}
					else {
						left_list.push_back(_left);
						_f = nonef;
					}
				}
				break;
				case bothf:
#pragma omp critical
				{
					if (vec_l + 2 * blocksize - 1 <= vec_r) {
						_left = vec_l;
						vec_l = vec_l + blocksize;
						_right = vec_r;
						vec_r = vec_r - blocksize;
					}
					else {
						_f = nonef;
					}
				}
				break;
				case nonef:
					break;
				}
			}

		}

		while (left_list.size() != 0 && right_list.size() != 0) {
			std::vector<typename std::vector<WSP>::iterator> _left_list, _right_list;
#pragma omp parallel for
			for (int i = 0; i < std::min(left_list.size(), right_list.size()); i++) {
				auto _left = left_list[i];
				auto _right = right_list[i];
				std::swap(*_left, *_right);
				_left++;
				_right--;
				while ((_left - begin) % blocksize != 0 && (end - _right) % blocksize != 0) {
					while ((_left - begin) % blocksize != 0 && get_nodeNum(_left) <= mid) {
						_left++;
					}
					while ((end - _right) % blocksize != 0 && get_nodeNum(_right) > mid) {
						_right--;
					}
					if ((_left - begin) % blocksize != 0 && (end - _right) % blocksize != 0) {
						std::swap(*_left, *_right);
						_left++;
						_right--;
					}
				}
				if ((_left - begin) % blocksize != 0) {
#pragma omp critical
					_left_list.push_back(_left);
				}
				if ((end - _right) % blocksize != 0) {
#pragma omp critical
					_right_list.push_back(_right);
				}
			}
			left_list = _left_list;
			right_list = _right_list;
		}


		for (int i = 0; i < left_list.size(); i++) {
			auto _left = left_list[i];
			while (get_nodeNum(vec_l) > mid && _left < vec_l) {
				vec_l--;
			}
			if (_left < vec_l) {
				std::swap(*_left, *(vec_l));
				vec_l--;
				_left++;
			}
			else {
				continue;
			}
			while ((_left - begin) % blocksize != 0) {
				while (get_nodeNum(vec_l) > mid && _left < vec_l) {
					vec_l--;
				}
				if (_left < vec_l) {
					std::swap(*_left, *(vec_l));
					vec_l--;
					_left++;
				}
				else
				{
					break;
				}
			}
		}
		for (int i = 0; i < right_list.size(); i++) {
			auto _right = right_list[i];
			while (get_nodeNum(vec_r) <= mid && vec_r < _right) {
				vec_r++;
			}
			if (vec_r < _right) {
				std::swap(*_right, *(vec_r));
				vec_r++;
				_right--;
			}
			else {
				continue;
			}
			while ((end - _right) % blocksize != 0) {
				while (get_nodeNum(vec_r) <= mid && vec_r < _right) {
					vec_r++;
				}
				if (vec_r < _right) {
					std::swap(*_right, *(vec_r));
					vec_r++;
					_right--;
				}
				else {
					break;
				}
			}
		}
	}
	if (vec_l < vec_r) {
		while (vec_l < vec_r) {
			while (vec_l < vec_r  && get_nodeNum(vec_l) <= mid) {
				vec_l++;
			}
			while (vec_l < vec_r  && get_nodeNum(vec_r) > mid) {
				vec_r--;
			}
			if (vec_l < vec_r) {
				std::swap(*vec_l, *vec_r);
				vec_l++;
				vec_r--;
			}
		}
	}
	if (get_nodeNum(vec_l) > mid) {
		return vec_l;
	}
	else {
		return vec_l + 1;
	}
}


FLAG comQuadtree::neutral_le(std::vector<WSP>& vec, std::vector<WSP>::iterator begin, std::vector<WSP>::iterator end,std::vector<WSP>::iterator& left, std::vector<WSP>::iterator& right, float mid, FLAG f, int blocksize)
{
	auto pleft = left;
	auto pright = right;
	auto left_end = left;
	auto right_end = right;
	switch (f)
	{
	case leftf:
		left_end = left + blocksize;
		right_end = right - blocksize + ((end - right) % blocksize);
		break;
	case rightf:
		left_end = left + blocksize - ((left - begin) % blocksize);
		right_end = right - blocksize;
		break;
	case bothf:
		left_end = left + blocksize;
		right_end = right - blocksize;
		break;
	}
	f = nonef;
	while (f == nonef)
	{
		while (pleft != left_end && Bccp(pleft) <= mid) {
			pleft++;
		}
		while (pright != right_end && Bccp(pright) > mid) {
			pright--;
		}
		if (pleft == left_end) {
			if (pright == right_end) {
				f = bothf;
			}
			else {
				f = leftf;
			}
		}
		else if (pright == right_end) {
			f = rightf;
		}
		else {
			std::swap(*pleft, *pright);
		}
	}
	left = pleft;
	right = pright;
	return f;
}

std::vector<WSP>::iterator comQuadtree::ppartition_le(std::vector<WSP>& vec, std::vector<WSP>::iterator vec_l, std::vector<WSP>::iterator vec_r, float mid)
{
	int blocksize = 64;
	auto begin = vec_l;
	auto end =vec_r;
	std::vector<typename std::vector<WSP>::iterator> left_list, right_list;
	int n_threads = std::min(omp_get_num_threads(), static_cast<int>(floor(vec.size() / (2 * blocksize))));
	if (n_threads > 1) {
#pragma omp parallel shared(vec,vec_l,vec_r,mid,left_list,right_list) num_threads(n_threads)
		{
			auto _left = begin + blocksize * omp_get_thread_num();
			auto _right = end - blocksize * omp_get_thread_num();
			vec_l = vec_l + blocksize;
			vec_r = vec_r - blocksize;
			FLAG _f = bothf;
#pragma omp barrier
			while (_f != nonef)
			{
				_f = neutral_le(vec,begin,end, _left, _right, mid, _f, blocksize);
				switch (_f)
				{
				case leftf:
#pragma omp critical
				{
					if (vec_l + blocksize - 1 <= vec_r) {
						_left = vec_l;
						vec_l = vec_l + blocksize;
					}
					else {
						right_list.push_back(_right);
						_f = nonef;
					}
				}
				break;
				case rightf:
#pragma omp critical
				{
					if (vec_r - blocksize + 1 >= vec_r) {
						_right = vec_r;
						vec_r = vec_r - blocksize;
					}
					else {
						left_list.push_back(_left);
						_f = nonef;
					}
				}
				break;
				case bothf:
#pragma omp critical
				{
					if (vec_l + 2 * blocksize - 1 <= vec_r) {
						_left = vec_l;
						vec_l = vec_l + blocksize;
						_right = vec_r;
						vec_r = vec_r - blocksize;
					}
					else {
						_f = nonef;
					}
				}
				break;
				case nonef:
					break;
				}
			}

		}

		while (left_list.size() != 0 && right_list.size() != 0) {
			std::vector<typename std::vector<WSP>::iterator> _left_list, _right_list;
#pragma omp parallel for
			for (int i = 0; i < std::min(left_list.size(), right_list.size()); i++) {
				auto _left = left_list[i];
				auto _right = right_list[i];
				std::swap(*_left, *_right);
				_left++;
				_right--;
				while ((_left - begin) % blocksize != 0 && (end - _right) % blocksize != 0) {
					while ((_left - begin) % blocksize != 0 && Bccp(_left) <= mid) {
						_left++;
					}
					while ((end - _right) % blocksize != 0 && Bccp(_right) > mid) {
						_right--;
					}
					if ((_left - begin) % blocksize != 0 && (end - _right) % blocksize != 0) {
						std::swap(*_left, *_right);
						_left++;
						_right--;
					}
				}
				if ((_left - begin) % blocksize != 0) {
#pragma omp critical
					_left_list.push_back(_left);
				}
				if ((end - _right) % blocksize != 0) {
#pragma omp critical
					_right_list.push_back(_right);
				}
			}
			left_list = _left_list;
			right_list = _right_list;
		}


		for (int i = 0; i < left_list.size(); i++) {
			auto _left = left_list[i];
			while (Bccp(vec_l) > mid && _left < vec_l) {
				vec_l--;
			}
			if (_left < vec_l) {
				std::swap(*_left, *(vec_l));
				vec_l--;
				_left++;
			}
			else {
				continue;
			}
			while ((_left - begin) % blocksize != 0) {
				while (Bccp(vec_l) > mid && _left < vec_l) {
					vec_l--;
				}
				if (_left < vec_l) {
					std::swap(*_left, *(vec_l));
					vec_l--;
					_left++;
				}
				else
				{
					break;
				}
			}
		}
		for (int i = 0; i < right_list.size(); i++) {
			auto _right = right_list[i];
			while (Bccp(vec_r) <= mid && vec_r < _right) {
				vec_r++;
			}
			if (vec_r < _right) {
				std::swap(*_right, *(vec_r));
				vec_r++;
				_right--;
			}
			else {
				continue;
			}
			while ((end - _right) % blocksize != 0) {
				while (Bccp(vec_r) <= mid && vec_r < _right) {
					vec_r++;
				}
				if (vec_r < _right) {
					std::swap(*_right, *(vec_r));
					vec_r++;
					_right--;
				}
				else {
					break;
				}
			}
		}
	}
	if (vec_l < vec_r) {
		while (vec_l < vec_r) {
			while (vec_l < vec_r  && Bccp(vec_l) <= mid) {
				vec_l++;
			}
			while (vec_l < vec_r  && Bccp(vec_r) > mid) {
				vec_r--;
			}
			if (vec_l < vec_r) {
				std::swap(*vec_l, *vec_r);
				vec_l++;
				vec_r--;
			}
		}
	}
	if (Bccp(vec_l) > mid) {
		return vec_l;
	}
	else {
		return vec_l + 1;
	}
}


FLAG comQuadtree::neutral_eq(std::vector<WSP>& vec, std::vector<WSP>::iterator begin, std::vector<WSP>::iterator end,std::vector<WSP>::iterator& left, std::vector<WSP>::iterator& right, FLAG f, int blocksize,Forests& forest)
{
	auto pleft = left;
	auto pright = right;
	auto left_end = left;
	auto right_end = right;
	switch (f)
	{
	case leftf:
		left_end = left + blocksize;
		right_end = right - blocksize + ((end - right) % blocksize);
		break;
	case rightf:
		left_end = left + blocksize - ((left - begin) % blocksize);
		right_end = right - blocksize;
		break;
	case bothf:
		left_end = left + blocksize;
		right_end = right - blocksize;
		break;
	}
	f = nonef;
	while (f == nonef)
	{
		while (pleft != left_end && forest.Find(pleft->P->get_pts(0)) != forest.Find(pleft->Q->get_pts(0))) {
			pleft++;
		}
		while (pright != right_end && forest.Find(pright->P->get_pts(0)) == forest.Find(pright->Q->get_pts(0))) {
			pright--;
		}
		if (pleft == left_end) {
			if (pright == right_end) {
				f = bothf;
			}
			else {
				f = leftf;
			}
		}
		else if (pright == right_end) {
			f = rightf;
		}
		else {
			std::swap(*pleft, *pright);
		}
	}
	left = pleft;
	right = pright;
	return f;
}

std::vector<WSP>::iterator comQuadtree::ppartition_eq(std::vector<WSP>& vec, std::vector<WSP>::iterator vec_l, std::vector<WSP>::iterator vec_r,Forests& forest)
{
	int blocksize = 64;
	auto begin = vec_l;
	auto end = vec_r;
	std::vector<typename std::vector<WSP>::iterator> left_list, right_list;
	int n_threads = std::min(omp_get_num_threads(), static_cast<int>(floor(vec.size() / (2 * blocksize))));
	if (n_threads > 1) {
#pragma omp parallel num_threads(n_threads)
		{
			auto _left = begin + blocksize * omp_get_thread_num();
			auto _right = end - blocksize * omp_get_thread_num();
			vec_l = vec_l + blocksize;
			vec_r = vec_r - blocksize;
			FLAG _f = bothf;
#pragma omp barrier
			while (_f != nonef)
			{
				_f = neutral_eq(vec,begin,end, _left, _right, _f, blocksize,forest);
				switch (_f)
				{
				case leftf:
#pragma omp critical
				{
					if (vec_l + blocksize - 1 <= vec_r) {
						_left = vec_l;
						vec_l = vec_l + blocksize;
					}
					else {
						right_list.push_back(_right);
						_f = nonef;
					}
				}
				break;
				case rightf:
#pragma omp critical
				{
					if (vec_r - blocksize + 1 >= vec_r) {
						_right = vec_r;
						vec_r = vec_r - blocksize;
					}
					else {
						left_list.push_back(_left);
						_f = nonef;
					}
				}
				break;
				case bothf:
#pragma omp critical
				{
					if (vec_l + 2 * blocksize - 1 <= vec_r) {
						_left = vec_l;
						vec_l = vec_l + blocksize;
						_right = vec_r;
						vec_r = vec_r - blocksize;
					}
					else {
						_f = nonef;
					}
				}
				break;
				case nonef:
					break;
				}
			}

		}

		while (left_list.size() != 0 && right_list.size() != 0) {
			std::vector<typename std::vector<WSP>::iterator> _left_list, _right_list;
#pragma omp parallel for
			for (int i = 0; i < std::min(left_list.size(), right_list.size()); i++) {
				auto _left = left_list[i];
				auto _right = right_list[i];
				std::swap(*_left, *_right);
				_left++;
				_right--;
				while ((_left - begin) % blocksize != 0 && (end - _right) % blocksize != 0) {
					while ((_left - begin) % blocksize != 0 && forest.Find(_left->P->get_pts(0)) != forest.Find(_left->Q->get_pts(0))) {
						_left++;
					}
					while ((end - _right) % blocksize != 0 && forest.Find(_right->P->get_pts(0)) == forest.Find(_right->Q->get_pts(0))) {
						_right--;
					}
					if ((_left - begin) % blocksize != 0 && (end - _right) % blocksize != 0) {
						std::swap(*_left, *_right);
						_left++;
						_right--;
					}
				}
				if ((_left - begin) % blocksize != 0) {
#pragma omp critical
					_left_list.push_back(_left);
				}
				if ((end - _right) % blocksize != 0) {
#pragma omp critical
					_right_list.push_back(_right);
				}
			}
			left_list = _left_list;
			right_list = _right_list;
		}


		for (int i = 0; i < left_list.size(); i++) {
			auto _left = left_list[i];
			while (forest.Find(vec_l->P->get_pts(0)) == forest.Find(vec_l->Q->get_pts(0)) && _left < vec_l) {
				vec_l--;
			}
			if (_left < vec_l) {
				std::swap(*_left, *(vec_l));
				vec_l--;
				_left++;
			}
			else {
				continue;
			}
			while ((_left - begin) % blocksize != 0) {
				while (forest.Find(vec_l->P->get_pts(0)) == forest.Find(vec_l->Q->get_pts(0)) && _left < vec_l) {
					vec_l--;
				}
				if (_left < vec_l) {
					std::swap(*_left, *(vec_l));
					vec_l--;
					_left++;
				}
				else
				{
					break;
				}
			}
		}
		for (int i = 0; i < right_list.size(); i++) {
			auto _right = right_list[i];
			while (forest.Find(vec_r->P->get_pts(0)) != forest.Find(vec_r->Q->get_pts(0))&& vec_r < _right) {
				vec_r++;
			}
			if (vec_r < _right) {
				std::swap(*_right, *(vec_r));
				vec_r++;
				_right--;
			}
			else {
				continue;
			}
			while ((end - _right) % blocksize != 0) {
				while (forest.Find(vec_r->P->get_pts(0)) != forest.Find(vec_r->Q->get_pts(0)) && vec_r < _right) {
					vec_r++;
				}
				if (vec_r < _right) {
					std::swap(*_right, *(vec_r));
					vec_r++;
					_right--;
				}
				else {
					break;
				}
			}
		}
	}
	if (vec_l < vec_r) {
		while (vec_l < vec_r) {
			while (vec_l < vec_r  && forest.Find(vec_l->P->get_pts(0)) != forest.Find(vec_l->Q->get_pts(0))) {
				vec_l++;
			}
			while (vec_l < vec_r  && forest.Find(vec_r->P->get_pts(0)) == forest.Find(vec_r->Q->get_pts(0))) {
				vec_r--;
			}
			if (vec_l < vec_r) {
				std::swap(*vec_l, *vec_r);
				vec_l++;
				vec_r--;
			}
		}
	}
	if (forest.Find(vec_l->P->get_pts(0)) == forest.Find(vec_l->Q->get_pts(0))) {
		return vec_l;
	}
	else {
		return vec_l + 1;
	}
}