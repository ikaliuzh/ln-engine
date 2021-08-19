#pragma once

#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>


template <size_t DIM, typename T>
struct vec{
	vec() {
		for (size_t i = 0; i < DIM; i++)
			data_[i] = T();
	}

	T& operator[](const size_t i){
		assert(i < DIM);
		return data_[i];
	}

	const T& operator[](const size_t i) const{
			assert(i < DIM);
			return data_[i];
	}
private:
	T data_[DIM];
};


template <typename T>
struct vec<2, T> {
	vec(): x(T()), y(T()) {}
	vec(T X, T Y) : x(X), y(Y) {}

	T& operator[](const size_t i){
		assert(i < 2);
		return i == 0 ? x : y;
	}

	const T& operator[](const size_t i) const{
		assert(i < 2);
		return i == 0 ? x : y;
	}
	T x, y;
};


template <typename T>
struct vec<3, T> {
	vec(): x(T()), y(T()), z(T()) {}
	vec(T X, T Y, T Z) : x(X), y(Y), z(Z) {}

	T& operator[](const size_t i){
		assert(i < 3);
		return i == 0 ? x : (i == 1? y : z);
	}

	const T& operator[](const size_t i) const{
		assert(i < 3);
		return i == 0 ? x : (i == 1? y : z);
	}

	float norm() const{
		return std::sqrt(x*x + y*y + z*z);
	}

	vec<3, T>& normalize(){
		*this = (*this) / norm();
		return *this;
	}

	T x, y, z;
};

template<size_t DIM, typename T, typename U>
vec<DIM, T> operator* (const vec<DIM, T>& lhs, U rhs){
	vec<DIM, T> res;
	for (size_t i = 0; i < DIM; ++i){
		res[i] = lhs[i] * rhs;
	}
	return res;
}


template<size_t DIM, typename T, typename U>
vec<DIM, T> operator* (U lhs, const vec<DIM, T>& rhs){
	vec<DIM, T> res;
	for (size_t i = 0; i < DIM; ++i){
		res[i] = rhs[i] * lhs;
	}
	return res;
}

template <size_t DIM, typename T>
T operator*(const vec<DIM, T>& lhs, const vec<DIM, T>& rhs){
	T res = T();
	for (size_t i = 0; i < DIM; ++i){
		res += lhs[i] * rhs[i];
	}
	return res;
}

template <typename T>
vec<3,T> cross(const vec<3,T>& v1, const vec<3,T>& v2) {
    return vec<3,T>(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

template <size_t DIM, typename T>
vec<DIM, T> operator+(const vec<DIM, T>& lhs, const vec<DIM, T>& rhs){
	vec<DIM, T> res;
	for (size_t i = 0; i < DIM; ++i){
		res[i] = lhs[i] + rhs[i];
	}
	return res;
}

template <size_t DIM, typename T>
vec<DIM, T> operator-(const vec<DIM, T>& lhs, const vec<DIM, T>& rhs){
	vec<DIM, T> res;
	for (size_t i = 0; i < DIM; ++i){
		res[i] = lhs[i] - rhs[i];
	}
	return res;
}


template<size_t DIM, typename T, typename U>
vec<DIM, T> operator/ (const vec<DIM, T>& lhs, U rhs){
	vec<DIM, T> res;
	for (size_t i = 0; i < DIM; ++i){
		res[i] = lhs[i] / rhs;
	}
	return res;
}

template<size_t DIM,typename T>
vec<DIM,T> operator-(const vec<DIM,T> &lhs) {
    return lhs*T(-1);
}

template<size_t DIM, typename T>
std::ostream& operator<<(std::ostream& out, const vec<DIM,T>& v) {
    for(size_t i = 0; i < DIM; i++) {
        out << v[i] << " " ;
    }
    return out ;
}
