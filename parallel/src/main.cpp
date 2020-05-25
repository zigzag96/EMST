#include<sys/time.h>

#include"iofile.h"
#include"pquicksort.h"
#include"wspd.h"


int main() {
	float time_use = 0;
	struct timeval start;
	struct timeval end;
	omp_set_num_threads(1);
	omp_set_dynamic(0);
	auto points=GenPoints(2000);
	write_points(points);
	// auto points=read_points();
	gettimeofday(&start,NULL);
	comQuadtree QT(points);
	QT.compute_tree();
	QT.compute_wspd();
	QT.GFKruskal();
	gettimeofday(&end,NULL);
	time_use=(end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec);
	std::cout<<"total time:"<<time_use<<std::endl;
	QT.write_EMST();
	return 0;
}
