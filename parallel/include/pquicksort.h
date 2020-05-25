#ifndef __P_QUICK_SORT
#define __P_QUICK_SORT
#include"head.h"

enum FLAG
{
	leftf, rightf, bothf, nonef
};

template<typename T>
FLAG neutral(std::vector<T>& vec,typename std::vector<T>::iterator begin,typename std::vector<T>::iterator end, typename std::vector<T>::iterator& left, typename std::vector<T>::iterator& right, T mid, FLAG f, int blocksize)
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
		while (pleft != left_end && *pleft < mid) {
			pleft++;
		}
		while (pright != right_end && *pright >= mid) {
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

//partition into P1 and P2: P1<mid, P2>=mid, result point to first element of P2
template<typename T>
typename std::vector<T>::iterator ppartition(std::vector<T>& vec, typename std::vector<T>::iterator vec_l, typename std::vector<T>::iterator vec_r, T mid) {
	int blocksize = 64;
	auto begin = vec_l;
	auto end = vec_r;
	std::vector<typename std::vector<T>::iterator> left_list, right_list;
	int n_threads = std::min(omp_get_num_threads(), static_cast<int>(floor((end-begin+1) / (2 * blocksize))));
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
				_f = neutral(vec,begin,end, _left, _right, mid, _f, blocksize);
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
			std::vector<typename std::vector<T>::iterator> _left_list, _right_list;
#pragma omp parallel for
			for (int i = 0; i < std::min(left_list.size(), right_list.size()); i++) {
				auto _left = left_list[i];
				auto _right = right_list[i];
				std::swap(*_left, *_right);
				_left++;
				_right--;
				while ((_left - begin) % blocksize != 0 && (end - _right) % blocksize != 0) {
					while ((_left - begin) % blocksize != 0 && *_left < mid) {
						_left++;
					}
					while ((end - _right) % blocksize != 0 && *_right >= mid) {
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
			while (*(vec_l) >= mid && _left < vec_l) {
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
				while (*(vec_l) >= mid && _left < vec_l) {
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
			while (*(vec_r) < mid && vec_r < _right) {
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
				while (*(vec_r) < mid && vec_r < _right) {
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
			while (vec_l < vec_r  && *vec_l < mid) {
				vec_l++;
			}
			while (vec_l < vec_r  && *vec_r >= mid) {
				vec_r--;
			}
			if (vec_l < vec_r) {
				std::swap(*vec_l, *vec_r);
				vec_l++;
				vec_r--;
			}
		}
	}
	if (*vec_l >= mid) {
		return vec_l;
	}
	else {
		return vec_l + 1;
	}
}

template<typename T>
void PQuickSort(std::vector<T>& vec, typename std::vector<T>::iterator left, typename std::vector<T>::iterator right) {
	if(left>=right){
		return;
	}
	T mid = *right;
	auto p = ppartition(vec, left, right-1, mid);
	std::swap(*p,*right);
	if(p>left){
#pragma omp task shared(vec) firstprivate(left,p)
		{
			PQuickSort(vec,left,p-1);
		}
	}
	if(p<right){
#pragma omp task shared(vec) firstprivate(right,p)
		{
			PQuickSort(vec,p+1,right);
		}
	}
}

template<typename T>
void PQuickSort(std::vector<T>& vec){
#pragma omp parallel shared(vec)
	{
#pragma omp single nowait
		{
			PQuickSort(vec,vec.begin(),vec.end()-1);
		}
	}
}

template<typename T> 
void PQuickSort(std::vector<T>& vec,typename std::vector<T>::iterator left, typename std::vector<T>::iterator right,bool f){
#pragma omp parallel
	{
#pragma omp single nowait
		{
			PQuickSort(vec,left,right);
		}
	}
}
#endif // !__P_QUICK_SORT

