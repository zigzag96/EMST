#ifndef __IOFILE__
#define __IOFILE__

#include "head.h"
#include <fstream>
#include <cstdlib> 
#include <ctime>
#include <string>
#include <sstream>

std::vector<std::vector<float>> GenPoints(int ptNum);
void write_points(std::vector<std::vector<float>> points, std::string FileName = "../out/points.csv");
std::vector<std::vector<float>> read_points(std::string FileName = "../out/points.csv");
void write_volumes(std::vector<float> volumes, std::string FileName = "../out/volumes.csv");
void draw_ply(std::vector<std::vector<float>> points, std::vector<std::vector<int>> edges = {}, std::string FileName = "../out/order.ply");

#endif // !__IOFILE__
