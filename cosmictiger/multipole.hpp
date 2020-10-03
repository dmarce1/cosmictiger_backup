/*  
 Copyright (c) 2016 Dominic C. Marcello

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef multipole_H_
#define multipole_H_

#include <cosmictiger/position.hpp>
#include <cosmictiger/vect.hpp>
#include <array>

constexpr int MP = 17;

template<class T>
class multipole {
private:
	T data[MP+1];
public:
	CUDA_EXPORT multipole();
	CUDA_EXPORT T operator ()() const;
	CUDA_EXPORT T& operator ()();
	CUDA_EXPORT T operator ()(int i, int j) const;
	CUDA_EXPORT T& operator ()(int i, int j);
	CUDA_EXPORT T operator ()(int i, int j, int k) const;
	CUDA_EXPORT T& operator ()(int i, int j, int k);
	CUDA_EXPORT multipole<T>& operator =(const multipole<T> &other);
	CUDA_EXPORT multipole<T>& operator =(T other);
	CUDA_EXPORT multipole operator>>(const vect<T> &dX) const;
	CUDA_EXPORT multipole<T>& operator>>=(const vect<T> &Y);
	CUDA_EXPORT multipole<T> operator +(const multipole<T> &vec) const;
	CUDA_EXPORT T& operator[](int i) {
		return data[i];
	}
	CUDA_EXPORT const T operator[](int i) const {
		return data[i];
	}

	template<class A>
	void serialize( A&& arc, unsigned ) {
		for( int i = 0; i < MP; i++) {
			arc & data[i];
		}
	}

};

template<class T>
CUDA_EXPORT inline multipole<T>::multipole() {
}

template<class T>
CUDA_EXPORT inline T multipole<T>::operator ()() const {
	return data[0];
}

template<class T>
CUDA_EXPORT inline T& multipole<T>::operator ()() {
	return data[0];
}

template<class T>
CUDA_EXPORT inline T multipole<T>::operator ()(int i, int j) const {
	static constexpr size_t map2[3][3] = { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } };
	return data[1 + map2[i][j]];
}

template<class T>
CUDA_EXPORT inline T& multipole<T>::operator ()(int i, int j) {
	static constexpr size_t map2[3][3] = { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } };
	return data[1 + map2[i][j]];
}

template<class T>
CUDA_EXPORT inline T multipole<T>::operator ()(int i, int j, int k) const {
	static constexpr size_t map3[3][3][3] = { { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } }, { { 1, 3, 4 }, { 3, 6, 7 }, { 4, 7, 8 } }, { { 2, 4, 5 }, { 4, 7, 8 },
			{ 5, 8, 9 } } };

	return data[7 + map3[i][j][k]];
}
template<class T>
CUDA_EXPORT inline T& multipole<T>::operator ()(int i, int j, int k) {
	static constexpr size_t map3[3][3][3] = { { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } }, { { 1, 3, 4 }, { 3, 6, 7 }, { 4, 7, 8 } }, { { 2, 4, 5 }, { 4, 7, 8 },
			{ 5, 8, 9 } } };
	return data[7 + map3[i][j][k]];
}

#include <cstring>

template<class T>
CUDA_EXPORT inline multipole<T>& multipole<T>::operator =(const multipole<T> &other) {
	memcpy(&data[0], &other.data[0], MP * sizeof(float));
	return *this;
}

template<class T>
CUDA_EXPORT inline multipole<T>& multipole<T>::operator =(T other) {
	for (int i = 0; i < MP; i++) {
		data[i] = other;
	}
	return *this;
}

template<class T>
CUDA_EXPORT inline multipole<T> multipole<T>::operator>>(const vect<T> &dX) const {
	multipole you = *this;
	you >>= dX;
	return you;
}

template<class T>
CUDA_EXPORT inline multipole<T> multipole<T>::operator +(const multipole<T> &vec) const {
	multipole<T> C;
	for (int i = 0; i < MP; i++) {
		C[i] = data[i] + vec[i];
	}
	return C;
}

template<class T>
CUDA_EXPORT inline multipole<T>& multipole<T>::operator>>=(const vect<T> &Y) {
	multipole<T> &me = *this;
	for (int p = 0; p < 3; p++) {
		for (int q = p; q < 3; q++) {
			for (int l = q; l < 3; l++) {
				me(p, q, l) -= me() * Y[p] * Y[q] * Y[l];
				me(p, q, l) -= me(p, q) * Y[l];
				me(p, q, l) -= me(q, l) * Y[p];
				me(p, q, l) -= me(l, p) * Y[q];
			}
		}
	}
	for (int p = 0; p < 3; p++) {
		for (int q = p; q < 3; q++) {
			me(p, q) += me() * Y[p] * Y[q];
		}
	}
	return me;
}


struct multi_src {
	multipole<float> m;
	float r;
	vect<position> x;
	template<class A>
	void serialize( A&& arc, unsigned ) {
		arc & m;
		arc & r;
		arc & x;
	}
};


#endif /* multipole_H_ */
