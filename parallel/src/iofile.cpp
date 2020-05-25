#include "iofile.h"

using namespace std;
vector<vector<float>> GenPoints(int ptNum) {
	vector<vector<float>> points;
	auto seed = time(0);
	srand(seed);
	for (int i = 0; i < ptNum; i++) {
		float x = 128.0*static_cast <float> (rand()) / (static_cast <float> (RAND_MAX)+1.0f);
		float y = 128.0*static_cast <float> (rand()) / (static_cast <float> (RAND_MAX)+1.0f);
		float z = 128.0*static_cast <float> (rand()) / (static_cast <float> (RAND_MAX)+1.0f);
		points.push_back({ x,y,z });
	}
	return points;
}

void write_points(vector<vector<float>> points,string FileName)
{
	ofstream ptfile;
	ptfile.open(FileName, ios::out | ios::trunc);
	for (auto i = points.begin(); i != points.end(); i++) {
		ptfile << (*i)[0] << "," << (*i)[1] << "," << (*i)[2] << endl;
	}
	ptfile.close();
}

std::vector<std::vector<float>> read_points(string FileName) {
	std::vector<std::vector<float>> points;
	ifstream ptfile(FileName,ios::in);
	string line;
	int i = 0;
	while (getline(ptfile, line)) {
		stringstream ss(line);
		string str;
		vector<float> point;
		while (getline(ss, str, ',')) {
			point.push_back(static_cast<float>(atof(str.c_str())));
		}
		points.push_back(point);
		i++;
	}
	ptfile.close();
	return points;
}

void write_volumes(vector<float> volumes,string FileName) {
	ofstream ptfile;
	ptfile.open(FileName, ios::out | ios::trunc);
	for (auto i = volumes.begin(); i != volumes.end(); i++) {
		ptfile << *i << endl;
	}
	ptfile.close();
}

void draw_ply(vector<vector<float>> points, vector<vector<int>> edges, string FileName) {
	ofstream plyfile;
	if (edges.size() == 0) {
		plyfile.open(FileName, ios::out | ios::trunc);
		plyfile << "ply" << endl;
		plyfile << "format ascii 1.0" << endl;
		plyfile << "comment object : colored pcd" << endl;
		plyfile << "element vertex " << points.size() << endl;
		plyfile << "property float x" << endl;
		plyfile << "property float y" << endl;
		plyfile << "property float z" << endl;
		plyfile << "element edge " << points.size() - 1 << endl;
		plyfile << "property int32 vertex1" << endl;
		plyfile << "property int32 vertex2" << endl;
		plyfile << "end_header" << endl;
		for (auto i = points.begin(); i != points.end(); i++) {
			plyfile << (*i)[0] << " " << (*i)[1] << " " << (*i)[2] << endl;
		}
		for (int i = 0; i < points.size()-1; i++) {
			plyfile << i << " " << i + 1 << endl;
		}
	}
	else {
		plyfile.open(FileName, ios::out | ios::trunc);
		plyfile << "ply" << endl;
		plyfile << "format ascii 1.0" << endl;
		plyfile << "comment object : colored pcd" << endl;
		plyfile << "element vertex " << points.size() << endl;
		plyfile << "property float x" << endl;
		plyfile << "property float y" << endl;
		plyfile << "property float z" << endl;
		plyfile << "element edge " << edges.size() << endl;
		plyfile << "property int32 vertex1" << endl;
		plyfile << "property int32 vertex2" << endl;
		plyfile << "end_header" << endl;
		for (auto i = points.begin(); i != points.end(); i++) {
			plyfile << (*i)[0] << " " << (*i)[1] << " " << (*i)[2] << endl;
		}
		for (auto i = edges.begin(); i != edges.end() ; i++) {
			plyfile << (*i)[0] << " " << (*i)[1] << endl;
		}
	}
	plyfile.close();
}