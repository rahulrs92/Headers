#pragma once

#include<iostream>
#include<iomanip>
#include<valarray>
#include<array>

using namespace std;

template <size_t m,size_t n,class T = double> class matrix;

template <class T> ostream& operator<< (ostream& os, const valarray<T> V) { for(int i=0; i<V.size(); i++) os<<V[i]<<" "; return os;  }

class R
{
public:
	size_t ix;
	R() { ix = 0; }
	R(size_t iix) { ix = iix; }
};

class C
{
public:
	size_t ix;
	C() {}
	C(size_t iix) { ix = iix; }
};

template <size_t m,size_t n,class T = double> class matrix_base
{
protected:
	valarray<T> V;
	
	valarray<T> m_tr(valarray<T> iV, size_t im, size_t in) { return iV[gslice(0,{in,im},{1,in})]; }
	
public:
	
	slice r(size_t iq) const { return slice(iq*n,n,1); }
	slice c(size_t iq) const { return slice(iq,m,n); }
	slice r(size_t iq) { return slice(iq*n,n,1); }
	slice c(size_t iq) { return slice(iq,m,n); }
	valarray<size_t> R, C;
	
	matrix_base() : matrix_base(T(0)) {}
	matrix_base(T iS)
	{
		V = valarray<T> (iS,m*n);
		R = valarray<size_t>(size_t(0),m*n);
		C = valarray<size_t>(size_t(0),m*n);
		for(int i=0; i<m; i++) R[r(i)] = i;
		for(int i=0; i<n; i++) C[c(i)] = i;
		//cout<<R<<endl;
		//cout<<C<<endl;
	}
	
	matrix_base(array<T,m*n> iV) : matrix_base() { for(int i=0; i<m*n; i++) V[i] = iV[i]; }
	matrix_base(valarray<T> iV) : matrix_base() { V = iV; }
	matrix_base(slice_array<T> iS) : matrix_base(valarray<T>(iS)) {}
	void operator = (matrix<m,n,T> iM) { V = iM.V; }
	
	size_t M() { return m; }
	size_t N() { return n; }

	valarray<T> operator [] (size_t q) const { return V[r(q)]; }
	slice_array<T> operator [] (size_t q) { return V[r(q)]; }
	
	valarray<T> operator () (size_t q) const { return V[c(q)]; }
	slice_array<T> operator () (size_t q) { return V[c(q)]; }
	T& operator () (size_t p, size_t q) { return V[n*p + q]; }
	
	template <class U> matrix<m,n,T> operator + (matrix_base<m,n,U> iM)
	{
		matrix<m,n,T> oM((*this)(0,0));
		for(int i=0; i<m; i++) for(int j=0; j<n; j++) oM(i,j) = (*this)(i,j) + iM(i,j); 
		return oM;
	}
	matrix<m,n,valarray<T>> operator + (matrix_base<m,n,valarray<T>> iM)
	{
		matrix<m,n,valarray<T>> oM((*this)(0,0) + iM(0,0));
		for(int i=0; i<m; i++) for(int j=0; j<n; j++) oM(i,j) = (*this)(i,j) + iM(i,j); 
		return oM;
	}
	matrix<m,n,T> operator + (T iS) { return matrix_base<m,n,T>(V+iS); }
	matrix<m,n,T> operator + (valarray<T> iV) { return matrix_base<m,n,T>(iV.size() == m*n ? V+iV : V); }
	
	template <class U> matrix<m,n,T> operator - (matrix_base<m,n,U> iM)
	{
		matrix<m,n,T> oM((*this)(0,0));
		for(int i=0; i<m; i++) for(int j=0; j<n; j++) oM(i,j) = (*this)(i,j) - iM(i,j); 
		return oM;
	}
	matrix<m,n,valarray<T>> operator - (matrix_base<m,n,valarray<T>> iM)
	{
		matrix<m,n,valarray<T>> oM((*this)(0,0) - iM(0,0));
		for(int i=0; i<m; i++) for(int j=0; j<n; j++) oM(i,j) = (*this)(i,j) - iM(i,j); 
		return oM;
	}
	matrix<m,n,T> operator - (T iS) { return matrix_base<m,n,T>(V-iS); }
	matrix<m,n,T> operator - (valarray<T> iV) { return matrix_base<m,n,T>(iV.size() == m*n ? V-iV : V); }
	
	matrix<m,n,T> operator += (matrix_base<m,n,T> iM) { return matrix_base<m,n,T>(V+iM.V); }
	matrix<m,n,T> operator += (T iS) { return matrix_base<m,n,T>(V+iS); }
	matrix<m,n,T> operator += (valarray<T> iV) { return matrix_base<m,n,T>(iV.size() == m*n ? V+iV : V); }
	
	matrix<m,n,T> operator -= (matrix_base<m,n,T> iM) { return matrix_base<m,n,T>(V-iM.V); }
	matrix<m,n,T> operator -= (T iS) { return matrix_base<m,n,T>(V-iS); }
	matrix<m,n,T> operator -= (valarray<T> iV) { return matrix_base<m,n,T>(iV.size() == m*n ? V-iV : V); }
	
	matrix<m,n,T> operator * (T iS) { return matrix_base<m,n,T>(V*iS); }
	matrix<m,n,T> operator * (valarray<T> iV) { return (iV.size() == m*n ? matrix_base<m,n,T>(V*iV) : *this); }
	
	matrix<m,n,T> operator / (valarray<T> iV) { return matrix_base<m,n,T>(iV.size() == m*n ? V/iV : V); }
	matrix<m,n,T> operator / (T iS) { return matrix_base<m,n,T>(V/iS); }
	
	valarray<T>& operator ! () { return V; }
	matrix<n,m,T> operator ~ () { return matrix_base<n,m,T>(m_tr(V,m,n)); }
	matrix<n,m,T> tp() { return matrix_base<n,m,T>(m_tr(V,m,n)); }
	
		//cout<<"oM(0,0) from matrix_base: "<<oM(0,0)<<endl;
		//cout<<"YAY!!!!!! "<<Mp<<endl;
		//cout<<"Size of iM(0,0) from matrix_base: "<<iM(0,0).size()<<endl;
		//cout<<"Mp from matrix_base: "<<Mp<<endl;
	
	template <size_t o> matrix<m,o,T> operator * (matrix_base<n,o,T> iM)
	{
		matrix_base<m,o,T> Mp(iM(0,0)*0);
		for(int i=0; i<m; i++) for(int j=0; j<o; j++) for(int k=0; k<n; k++) Mp(i,j) += (*this)(i,k) * iM(k,j);
		return Mp;
	}
	
	template <size_t o> matrix<m,o,valarray<T>> operator * (matrix_base<n,o,valarray<T>> iM)
	{
		matrix_base<m,o,valarray<T>> Mp(valarray<valarray<T>>(iM(0,0)*0,m*o));
		for(int i=0; i<m; i++) for(int j=0; j<o; j++) for(int k=0; k<n; k++) Mp(i,j) += (*this)(i,k) * iM(k,j);
		return Mp;
	}
};



template <size_t m,size_t n,class T> class matrix : public matrix_base<m,n,T>
{
private:
	valarray<T>& V = matrix_base<m,n,T>::V;
	slice r(size_t iq) { return matrix_base<m,n,T>::r(iq); }
	slice c(size_t iq) { return matrix_base<m,n,T>::c(iq); }
public:
	matrix() : matrix_base<m,n,T>() {}
	matrix(T iS) : matrix_base<m,n,T>(iS) {}
	matrix(array<T,m*n> iV) : matrix_base<m,n,T>(iV) {}
	matrix(matrix_base<m,n,T> iM) : matrix_base<m,n,T>(iM) {}
	matrix_base<m,n,T>& operator ! () { return *this; }
	void operator = (matrix<m,n,T> iM) { V = iM.V; }
	
	matrix<1,n,T> operator [] (size_t q) { return matrix_base<1,n,T>(V[r(q)]); }
	matrix<m,1,T> operator () (size_t q) { return matrix_base<m,1,T>(V[c(q)]); }
	T& operator () (size_t p, size_t q) { return V[n*p + q]; }
	
	//template <size_t o> matrix<m,o,T> operator* (matrix<n,o,T> iM) { return matrix_base<m,n,T>(*this) * matrix_base<n,o,T>(iM); }
};

template <size_t n,class T> class matrix<1,n,T> : public matrix_base<1,n,T>
{
private:
	valarray<T>& V = matrix_base<1,n,T>::V;
public:
	matrix() : matrix_base<1,n,T>() {}
	matrix(T iS) : matrix_base<1,n,T>(iS) {}
	matrix(array<T,n> iV) : matrix_base<1,n,T>(iV) {}
	matrix(matrix_base<1,n,T> iM) : matrix_base<1,n,T>(iM) {}
	matrix_base<1,n,T>& operator ! () { return *this; }
	void operator = (matrix<1,n,T> iM) { V = iM.V; }
	
	T& operator [] (size_t q) { return V[q]; }
	
	template <size_t o> matrix<1,o,T> operator * (matrix<n,o,T> iM) { return matrix_base<1,n,T>(*this) * matrix_base<n,o,T>(iM); }
	matrix<1,n,T> operator * (T iS) { return matrix_base<1,n,T>(V*iS); }
	matrix<1,n,T> operator * (valarray<T> iV) { return matrix_base<1,n,T>(iV.size() == n ? (V*iV) : V); }
	T operator* (matrix<n,1,T> iM) { return (V*(!!iM)).sum(); }
	valarray<T> operator* (matrix<n,1,valarray<T>> iM)
	{
		valarray<T> oV(0.,iM(0,0).size());
		for(int i=0; i<iM(0,0).size(); i++) for(int j=0; j<n; j++) oV[i] += V[j] * iM[j][i];
		return oV;
	}
};


template <size_t m,class T> class matrix<m,1,T> : public matrix_base<m,1,T>
{
private:
	valarray<T>& V = matrix_base<m,1,T>::V;
public:
	matrix() : matrix_base<m,1,T>() {}
	matrix(T iS) : matrix_base<m,1,T>(iS) {}
	matrix(array<T,m> iV) : matrix_base<m,1,T>(iV) {}
	matrix(matrix_base<m,1,T> iM) : matrix_base<m,1,T>(iM) {}
	matrix_base<m,1,T>& operator ! () { return *this; }
	void operator = (matrix<m,1,T> iM) { V = iM.V; }
	
	T& operator [] (size_t q) { return V[q]; }
};

template <class T> class matrix<3,1,T> : public matrix_base<3,1,T>
{
private:
	valarray<T>& V = matrix_base<3,1,T>::V;
public:
	matrix() : matrix_base<3,1,T>() {}
	matrix(T iS) : matrix_base<3,1,T>(iS) {}
	matrix(array<T,3> iV) : matrix_base<3,1,T>(iV) {}
	matrix(matrix_base<3,1,T> iM) : matrix_base<3,1,T>(iM) {}
	matrix_base<3,1,T>& operator ! () { return *this; }
	void operator = (matrix<3,1,T> iM) { V = iM.V; }
	
	T& operator [] (size_t q) { return V[q]; }
	
	template <size_t o> matrix<3,o,T> operator * (matrix<1,o,T> iM) { return matrix_base<3,1,T>(*this) * matrix_base<1,o,T>(iM); }
	template <size_t o> matrix<3,o,valarray<T>> operator * (matrix_base<1,o,valarray<T>> iM)
	{
		matrix_base<3,o,valarray<T>> Mp(valarray<valarray<T>>(iM(0,0)*0,3*o));
		for(int i=0; i<3; i++) for(int j=0; j<o; j++) for(int k=0; k<1; k++) Mp(i,j) += (*this)(i,k) * iM(k,j);
		return Mp;
	}
	matrix<3,1,T> operator * (T iS) { return matrix_base<3,1,T>(V*iS); }
	matrix<3,1,T> operator * (valarray<T> iV) { return (iV.size() == 3 ? matrix_base<3,1,T>(V*iV) : *this); }
	matrix<3,1,T> operator * (matrix<3,1,T> iM)
	{
		matrix<3,1,T> oM;
		oM[0] = (*this)[1]*iM[2] - (*this)[2]*iM[1];
		oM[1] = (*this)[2]*iM[0] - (*this)[0]*iM[2];
		oM[2] = (*this)[0]*iM[1] - (*this)[1]*iM[0];
		return oM;
	}
	
};

		//cout<<"R at minor is of "<<iV<<" is: "<<R<<endl;
		//cout<<"Size of R is: "<<R.size()<<endl;
	// cout<<"minor of "<<iV<<" at "<<im<<", "<<in<<" is: "<<minor(iV,im,in)<<endl;
	//cout<<"cofactor of "<<iV[i]<<": "<<tv<<endl;
		// cout<<"YAY!!!!!!"<<endl;
		//cout<<"calculating Determinant of: "<<iV<<" of size "<<s<<endl;
				//cout<<"cofactor of "<<iV[0]<<": "<<tv<<endl;
		//cout<<"Cofactor matrix: "<<oV<<endl;

template <size_t m,class T> class matrix<m,m,T> : public matrix_base<m,m,T>
{
private:
	valarray<T>& V = matrix_base<m,m,T>::V;
	valarray<size_t>& R = matrix_base<m,m,T>::R;
	valarray<size_t>& C = matrix_base<m,m,T>::C;
	slice r(size_t iq) { return matrix_base<m,m,T>::r(iq); }
	slice c(size_t iq) { return matrix_base<m,m,T>::c(iq); }
	
	
	size_t l(valarray<T> iV) { return sqrt(iV.size()); }
	
	T minor(valarray<T> iV,size_t im,size_t in)
	{
		return vdet(iV[valarray<size_t>(R[R<l(iV) && C<l(iV)])!=im && valarray<size_t>(C[R<l(iV) && C<l(iV)])!=in]);
	}
	T cf(valarray<T> iV,size_t im,size_t in) { return pow(-1,im+in)*minor(iV,im,in); }
	
	void circ_l(valarray<T>& iV) { for(int i=0; i<l(iV); i++) iV[slice(l(iV)*i,l(iV),1)] = valarray<T>(iV[slice(l(iV)*i,l(iV),1)]).cshift(1); }
	
	void circ_u(valarray<T>& iV) { for(int i=0; i<l(iV); i++) iV[slice(i,l(iV),l(iV))] = valarray<T>(iV[slice(i,l(iV),l(iV))]).cshift(1); }

	
	T vdet(valarray<T> iV)
	{
		T D = 0;
		size_t s = l(iV);
		for(int i=0; i<s && s>1; i++)
		{
			T tv = cf(iV,0,i);
			D = D + iV[i]*tv;
		}
		return (s>1 ? D : iV[0]);
	}
	
	valarray<T> m_cf(valarray<T> iV)
	{
		valarray<T> oV(T(0),iV.size());
		size_t s = l(iV);
		for(int i=0; i<s; i++)
		{
			for(int j=0; j<s; j++)
			{
				T tv = cf(iV,i,j);
				oV[(i*s)+j] = tv;
			}
		}
		return oV;
	}
	
	valarray<T> m_adj(valarray<T> iV) { return matrix_base<m,m,T>::m_tr(m_cf(iV),m,m); }
	
	valarray<T> m_inv(valarray<T> iV) { valarray<T> a = (m_adj(iV)/vdet(iV)); return a; }
	
	valarray<T> m_delta(size_t cn, valarray<T> iV, valarray<T> iC) { iV[slice(cn,l(iV),l(iV))] = iC; return iV; }
	
public:
	matrix() : matrix_base<m,m,T>() {}
	matrix(T iS) : matrix_base<m,m,T>(iS) {}
	matrix(array<T,m*m> iV) : matrix_base<m,m,T>(iV) {}
	matrix(matrix_base<m,m,T> iM) : matrix_base<m,m,T>(iM) {}
	matrix_base<m,m,T>& operator! () { return *this; }
	void operator = (matrix<m,m,T> iM) { V = iM.V; }
	
	T det() { return (m>1 ? vdet(V) : V[0]); }
	matrix<1,m,T> operator [] (size_t q) { return matrix_base<1,m,T>(V[r(q)]); }
	matrix<m,1,T> operator () (size_t q) { return matrix_base<m,1,T>(V[c(q)]); }
	T& operator () (size_t p, size_t q) { return V[m*p + q]; }

	
	matrix<m,m,T> operator~ () { return ~matrix_base<m,m,T>(*this); }
	
	//void operator= (matrix<m,m,T> iM) { V = !!iM; }
	//void operator= (matrix_base<m,m,T> iM) { *this = matrix<m,m,T>(iM); }
	
	matrix<m,m,T> operator / (T iS) { return matrix_base<m,m,T>(*this)/iS; }
	matrix<m,m,T> operator / (valarray<T> iV) { return matrix_base<m,m,T>(iV.size() == m*m ? V/iV : V); }
	matrix<m,m,T> operator / (matrix<m,m,T> iM) { return *this * iM^(-1); }
	
	matrix<m,m,T> operator^ (int p)
	{
		matrix<m,m,T> oM;
		(!!oM)[slice(0,m,m+1)] = 1;
		for(int i=0; i<abs(p); i++) { matrix<m,m,T> a = (p>0 ? *this : matrix<m,m,T>(m_inv(V)));  oM = oM * a;  }
		return oM;
	}
	
	T delta(size_t cn, matrix<1,m,T> iR) { return matrix<m,m,T>(m_delta(cn,V,!!iR)).det(); }
	
	slice_array<T> diag() { return V[slice(0,m,m+1)]; }
};

template <size_t m,size_t n,class T> matrix<m,n,T> sum(matrix<m,n,valarray<T>> iM)
{
	matrix<m,n,T> oM;
	for(int i=0; i<m; i++)
	{
		for(int j=0; j<n; j++)
		{
			//cout<<"iM(i,j): "<<iM(i,j).size()<<endl;
			oM(i,j) = iM(i,j).sum();
		}
	}
	return oM;
}

template <size_t m,size_t n,class T = double> ostream& operator<< (ostream& os, const matrix_base<m,n,T> iM)
{
	//os.precision(4);
	ios init(NULL);
	init.copyfmt(os);
	for(int i=0; i<m; i++)
	{
		os << '\t';
		for(int j=0; j<n; j++) os << setprecision(8) << setw(14) << iM[i][j] << " "; //<< '\t'  
		os << endl;
	}
	os.copyfmt(init);
	return os;
}
