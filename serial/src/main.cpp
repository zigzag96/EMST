#include<sys/time.h>

#include "head.h"
#include "iofile.h"
#include "wspd.h"

int main(void) {
	float time_use = 0;
	struct timeval start;
	struct timeval end;
	// auto points=GenPoints(1000);
	// write_points(points);
	auto points=read_points();
	gettimeofday(&start,NULL);
	comQuadtree QT(points);
	QT.compute_tree();
	QT.compute_wspd();
	std::cout<<"wspd done"<<std::endl;
	QT.GFKruskal();
	gettimeofday(&end,NULL);
	time_use=(end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec);
	std::cout<<"total time:"<<time_use<<std::endl;
	QT.write_EMST();
	return 0;
}