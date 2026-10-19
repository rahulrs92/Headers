#pragma once

#include<vector>
#include<valarray>

using namespace std;

template <typename T = double> class pick
{
public:
	T &a;
	pick(T &ia) : a(ia) {}
	void operator = (T ia) { a = ia; }
	operator T& () { return a; }
};

template <typename T = double> class column
{
public:
	vector<T> &a;
	T cn;
	column(vector<T>& ia, size_t& icn) : a(ia)
	{
		if constexpr (is_same<T,string>::value)
			cn = to_string(icn);
		else
			cn = icn;
	}
	void operator = (vector<T> ia) { for(int i=0; i<(a.size()>ia.size() ? ia : a).size(); i++) a[i] = ia[i]; }
	T& operator [] (size_t ix) { return (ix==0 ? cn : a[ix-1]); }
	operator vector<T>& () { return a; }
};

template <typename T = double> class row
{
public:
	vector<pick<T>> a;
	size_t sz;
	T rn;
	row() { sz = 0; rn = T(); }
	void set_row(size_t irn)
	{
		if constexpr (is_same<T,string>::value)
			rn = to_string(irn);
		else
			rn = irn;
	}
	void push_back(T &ia) { a.push_back(ia); sz++; }
	size_t size() { return sz; }
	T& operator [] (size_t ix) { return (ix==0 ? pick<T>(rn) : a[ix-1]); }
	void operator = (T* ia) { for(int i=1; i<sz; i++) a[i] = ia[i-1]; }
	void clear() { a.clear(); }
	
	operator T* ()
	{
		T* oa = new T[a.size()];
		for(int i=1; i<a.size(); i++) oa[i-1] = a[i];
		return oa;
	}
	
	operator valarray<T> ()
	{
		valarray<T> ov = new T[a.size()];
		for(int i=1; i<a.size(); i++) ov[i-1] = a[i];
		return ov;
	}
};
 
