#pragma once

#include<iostream>
#include<iomanip>
#include<valarray>
#include<array>
#include<initializer_list>
#include<tuple>

#include<pick.hh>

using namespace std;

template <class T> ostream& operator<< (ostream& os, const valarray<T> V) 
{
	for(int i=0; i<V.size(); i++) os<<V[i]<<" "; return os;
}

template <size_t m,size_t n,class T> class ref_matrix;
template <size_t m,size_t n,class T> class ref_vector;
template <size_t m,class T> class ref_select;
template <size_t m,size_t n,class T = double> class matrix;

template <size_t m,size_t n,class T = double> class matrix_base
{
protected:
	valarray<T> V;
	matrix_base(valarray<T> iV) : matrix_base() { V = iV; }
	static valarray<T> m_tr(valarray<T> iV) { return iV[gslice(0,{n,m},{1,n})]; }
	
public:
	matrix_base() : matrix_base(T(0)) {}
	matrix_base(T iS) { V = valarray<T> (iS,m*n); }
	matrix_base(initializer_list<T> iV) : matrix_base() { V = iV; }
	
	typedef T value_type;
	
	matrix<n,m,T> tp()
	{
		matrix_base<n,m,T> oM;
		valarray<T> oV = V[gslice(0,{n,m},{1,n})];
		for(int i=0; i<m; i++) for(int j=0; j<n; j++) oM[j][i] = V[i*n + j];
		//for(int i=0; i<m; i++) oM(i) = (*this)[i];
		return oM;
	}
	
	void operator = (matrix_base<m,n,T> iM) { V = iM.V; }
	slice r(size_t iq) { return slice(iq*n,n,1); }
	slice c(size_t iq) { return slice(iq,m,n); }
	
	size_t constexpr M() { return m; }
	size_t constexpr N() { return n; }
	valarray<T>& v() { return V; }
	
	ref_vector<1,n,T> operator [] (size_t q) { return ref_vector<1,n,T>(V,r(q)); }
	ref_vector<m,1,T> operator () (size_t q) { return ref_vector<m,1,T>(V,c(q)); }
	
	T& operator () (size_t p, size_t q) { return V[n*p + q]; }
	
	template <int o,int p> ref_matrix<(o<0 ? -o : o),(p<0 ? -p : p),T> sub(size_t ir=0, size_t ic=0)
	{
		size_t constexpr q = (o<0 ? -o : o), r = (p<0 ? -p : p);
		if(ir>=m || ic>=n)
		{
			cerr<<"sub error: origin out of range, assuming 0,0"<<endl;
			return ref_matrix<q,r,T>(V,gslice(0,{q,r},{n,1}));
		}
		size_t tr = (o<0 ? abs(int(o+ir)) : ir), tc = (p<0 ? abs(int(p+ic)) : ic);								// start
		size_t zr = (o<0 ? ir : m-ir), zc = (p<0 ? ic : n-ic);																		// max size
		//cout<<"tr = "<<tr<<", tc = "<<tc<<endl;
		//cout<<"zr = "<<zr<<", zc = "<<zc<<endl;
		if(q>zr || r>zc) { cerr<<"sub error: sub-matrix size out of range (max size: "<<zr<<","<<zc<<")"<<endl; }
		return ref_matrix<q,r,T>(V,gslice(tr*n + tc,{q,r},{n,1}));
	}

	ref_matrix<m-1,n-1,T> minor(size_t ir=0, size_t ic=0)
	{
		if(ir>=m || ic>=n)
		{
			cerr<<"minor error: origin out of range, assuming 0,0"<<endl;
			return ref_matrix<m-1,n-1,T>(V,gslice(1,{m-1,n-1},{n-1,1}));
		}
		matrix_base<m,n,size_t> ns;
		for(int i=0; i<m; i++) for(int j=0; j<n; j++) ns[i][j] = i*n + j;

		matrix_base<m-1,n,size_t> sel;
		for(int i=0; i<ir; i++) sel[i] = ns[i];
		for(int i=ir+1; i<m; i++) sel[i-1] = ns[i];
		matrix_base<m-1,n-1,size_t> sel2;
		for(int i=0; i<ic; i++) sel2(i) = sel(i);
		for(int i=ic+1; i<n; i++) sel2(i-1) = sel(i);

		return ref_matrix<m-1,n-1,T>(V,sel2.v());
	}
	
	matrix_base<m,n,T>& swap_rows(size_t r1,size_t r2)
	{
		matrix_base<1,n,T> tR = ref_vector<1,n,T>(V,r(r1));
		ref_vector<1,n,T>(V,r(r1)) = ref_vector<1,n,T>(V,r(r2));
		ref_vector<1,n,T>(V,r(r2)) = tR;
		return this[0];
	}
	matrix_base<m,n,T>& swap_cols(size_t c1,size_t c2)
	{
		matrix_base<1,n,T> tC = ref_vector<1,n,T>(V,c(c1));
		ref_vector<1,n,T>(V,c(c1)) = ref_vector<1,n,T>(V,c(c2));
		ref_vector<1,n,T>(V,c(c2)) = tC;
		return this[0];
	}
	
	matrix<m,n,T> operator + () { return *this; }
	matrix<m,n,T> operator + (T iS) { return matrix_base<m,n,T>(V+iS); }
	matrix<m,n,T> operator + (matrix_base<m,n,T> iM) { return matrix_base<m,n,T>(V+iM.V); }
	matrix<m,n,valarray<T>> operator + (matrix_base<m,n,valarray<T>> ivM)
	{
		matrix_base<m,n,valarray<T>> oM(ivM[0][0]);
		for(int i=0; i<m; i++) for(int j=0; j<n; j++) oM[i][j] = V[i*n + j] + ivM[i][j];
		return oM;
	}
	
	matrix<m,n,T> operator - () { return matrix_base<m,n,T>(-V); }
	matrix<m,n,T> operator - (T iS) { return matrix_base<m,n,T>(V-iS); }
	matrix<m,n,T> operator - (matrix_base<m,n,T> iM) { return matrix_base<m,n,T>(V-iM.V); }
	matrix<m,n,valarray<T>> operator - (matrix_base<m,n,valarray<T>> ivM)
	{
		matrix_base<m,n,valarray<T>> oM(ivM[0][0]);
		for(int i=0; i<m; i++) for(int j=0; j<n; j++) oM[i][j] = V[i*n + j] - ivM[i][j]; 
		return oM;
	}
	
	//matrix<m,n,T> operator * (T iS) { return matrix_base<m,n,T>(V*iS); }
	
	void operator += (matrix_base<m,n,T> iM) { V += iM.V; }
	void operator += (T iS) { V += iS; }
	
	void operator -= (matrix_base<m,n,T> iM) { V -= iM.V; }
	void operator -= (T iS) { V -= iS; }
	
	void operator *= (T iS) { V *= iS; }
	void operator /= (T iS) { V /= iS; }
	
	bool operator == (matrix_base<m,n,T> iM) { return (V == iM.v()).min(); }
	bool operator != (matrix_base<m,n,T> iM) { return !(V == iM.v()).min(); }
	
	matrix_base<m,n,bool> operator == (T it) { return matrix_base<m,n,bool>(V == it); }
	matrix_base<m,n,bool> operator != (T it) { return matrix_base<m,n,bool>(!(V == it)); }
	
	matrix<m,n,T> operator * (T iS)
	{
		matrix_base<m,n,T> Mp(V);
		for(int i=0; i<m; i++) for(int j=0; j<n; j++) Mp[i][j] *= iS;
		return Mp;
	}
	matrix<m,n,T> operator / (T iS)
	{
		matrix_base<m,n,T> Mp(V);
		for(int i=0; i<m; i++) for(int j=0; j<n; j++) Mp[i][j] /= iS;
		return Mp;
	}

	
	//valarray<T>& operator ! () { return V; }
	
	//matrix<n,m,T> operator ~ () { return matrix<n,m,T>(m_tr(V,m,n)); }
	
		//cout<<"oM(0,0) from matrix: "<<oM(0,0)<<endl;
		//cout<<"YAY!!!!!! "<<Mp<<endl;
		//cout<<"Size of iM(0,0) from matrix: "<<iM(0,0).size()<<endl;
		//cout<<"Mp from matrix: "<<Mp<<endl;
	
	
	template <size_t o> matrix<m,o,T> operator * (matrix_base<n,o,T> iM)
	{
		matrix_base<m,o,T> Mp(V[0]*0);
		//cout<<"iM[n-1][o-1]: "<<iM[n-1][o-1]<<endl;
		//cout<<"yes: "<<m<<", "<<n<<", "<<o<<endl;
		for(int i=0; i<m; i++) for(int j=0; j<o; j++) for(int k=0; k<n; k++) Mp[i][j] += V[i*n + k] * iM[k][j];
		return Mp;
	}
	template <size_t o> matrix<m,o,valarray<T>> operator * (matrix_base<n,o,valarray<T>> ivM)
	{
		matrix_base<m,o,valarray<T>> oM(ivM[0][0]*0);
		for(int i=0; i<m; i++) for(int j=0; j<o; j++) for(int k=0; k<n; k++)
					oM[i][j] += ivM[k][j] * V[i*n + k];
		return oM;
	}
	
	T norm()
	{
		T nrm(0);
		for(int i=0; i<m*n; i++) nrm += V[i]*V[i];
		return sqrt(nrm);
	}
/*
	matrix<m,1,T> operator * (matrix_base<n,1,T> iM)
	{
		matrix_base<m,1,T> Mp(V[0]*0);
		for(int i=0; i<m; i++) for(int k=0; k<n; k++) Mp(i,0) += V[i*n + k] * iM(k,0);
		return Mp;
	}
	matrix<m,1,valarray<T>> operator * (matrix_base<n,1,valarray<T>> ivM)
	{
		matrix<m,1,valarray<T>> oM(ivM(0,0)*0);
		cout<<"yes"<<endl;
		for(int i=0; i<m; i++) for(int k=0; k<n; k++) oM[i] += V[i*n + k] * valarray<T>(ivM[k]); 
		return oM;
	}
	*/
	
};

template <size_t m,size_t n,class T> matrix<m,n,T> operator + (T it,matrix_base<m,n,T> iM) { return iM + it; }
template <size_t m,size_t n,class T> matrix<m,n,T> operator - (T it,matrix_base<m,n,T> iM) { return -(iM - it); }
template <size_t m,size_t n,class T> matrix<m,n,T> operator * (T it,matrix_base<m,n,T> iM) { return iM * it; }
template <size_t m,size_t n,class T> matrix<m,n,T> operator / (T it,matrix_base<m,n,T> iM)
{
	for(int i=0; i<m; i++) for(int j=0; j<n; j++) iM[i][j] = it / iM[i][j];
	return iM;
}

//template <size_t m,size_t n,class T> matrix<m,n,T> operator * (matrix_base<m,n,T> iM,T iS) { iM.v() *= iS; return iM; }
//template <size_t m,size_t n,class T> matrix<m,n,T> operator / (matrix_base<m,n,T> iM,T iS) { iM.v() /= T(iS); return iM; }

template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator + 
															(valarray<T> iV,matrix_base<m,n,valarray<T>> iM) { return iM + iV; }
template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator - 
															(valarray<T> iV,matrix_base<m,n,valarray<T>> iM) { return -(iM - iV); }
template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator * 
															(valarray<T> iV,matrix_base<m,n,valarray<T>> iM) { return iM * iV; }

template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator + 
															(T it,matrix_base<m,n,valarray<T>> ivM)
{
	for(int i=0; i<m; i++) for(int j=0; j<n; j++) ivM[i][j] += it; 
	return ivM;
}
template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator - 
															(T it,matrix_base<m,n,valarray<T>> ivM)
{
	for(int i=0; i<m; i++) for(int j=0; j<n; j++) ivM[i][j] -= it; 
	return -ivM;
}
template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator * 
															(T it,matrix_base<m,n,valarray<T>> ivM)
{
	for(int i=0; i<m; i++) for(int j=0; j<n; j++) ivM[i][j] *= it;
	return ivM;
}
template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator / 
															(T it,matrix_base<m,n,valarray<T>> ivM)
{
	for(int i=0; i<m; i++) for(int j=0; j<n; j++) ivM[i][j] /= it;
	return ivM;
}

template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator + 
											(matrix_base<m,n,valarray<T>> ivM,matrix_base<m,n,T> iM) { return iM + ivM; }
template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator - 
											(matrix_base<m,n,valarray<T>> ivM,matrix_base<m,n,T> iM) { return -(iM - ivM); }
											
//template <size_t m,size_t n,class T> matrix_base<n,m,T> tp(matrix<m,n,T> iM) { return iM.transpose(); }

template <size_t m,size_t n,size_t o,class T> matrix<m,o,valarray<T>> operator * 
																				(matrix_base<m,n,valarray<T>> ivM,matrix_base<n,o,T> iM)
{
	matrix_base<m,o,valarray<T>> oM(ivM[0][0]*0);
	for(int i=0; i<m; i++) for(int j=0; j<o; j++) for(int k=0; k<n; k++) oM[i][j] += ivM[i][k] * iM[k][j];
	return oM;
}

template <size_t m,size_t n,class T> matrix<m,n,T> pull(matrix_base<m,n,valarray<T>> iM,size_t iq)
{
	matrix_base<m,n,T> oM;
	for(int i=0; i<m; i++) for(int j=0; j<n; j++) oM[i][j] = iM[i][j][iq];
	return oM;
}

template <size_t m,size_t n,size_t o,class T> matrix<m,n+o,T> cat(matrix_base<m,n,T> iM1,matrix_base<m,o,T> iM2)
{
	matrix_base<m,n+o,T> oM;
	for(int i=0; i<m; i++) for(int j=0; j<n; j++) oM[i][j] = iM1[i][j];
	for(int i=0; i<m; i++) for(int j=0; j<o; j++) oM[i][j+n] = iM2[i][j];
	return oM;
}

template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> pull(matrix_base<m,n,valarray<T>> iM,valarray<size_t> iv)
{
	matrix_base<m,n,valarray<T>> oM;
	for(int i=0; i<m; i++) for(int j=0; j<n; j++) oM[i][j] = valarray<T>(iM[i][j][iv]);
	return oM;
}

template <class T> matrix<3,1,T> operator*(matrix<3,1,T> a, matrix<3,1,T> b)
{
	matrix<3,1,T> c;
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
	return c;
}

template <size_t m,size_t n,class T> class matrix : public matrix_base<m,n,T>
{
public:
	matrix() : matrix_base<m,n,T>(T(0)) {}
	matrix(matrix_base<m,n,T> iM) : matrix_base<m,n,T>(iM) {}
	matrix(T iS) : matrix_base<m,n,T>(iS) {}
	matrix(initializer_list<T> iV) : matrix_base<m,n,T>(iV) {}

	void SVD(matrix<m, n, T> A, matrix<m, n, T> &U, valarray<T> &s, matrix<n, n, T> &V, T tolerance = 1e-15)
	{
		static_assert(m >= n, "SVD implementation requires rows >= columns (m >= n)");
		
		U = A;              // U starts as a copy of A
		V = T(0);           // Initialize V as zero matrix
		for(size_t i=0; i<n; ++i) V(i, i) = T(1); // Set Identity using direct access
		s = valarray<T>(T(0), n);
		
		bool changed = true;
		int max_sweeps = 100;
		int sweep = 0;
		
		while (changed && sweep < max_sweeps)
		{
			changed = false;
			sweep++;
			
			for (size_t i = 0; i < n; i++)
			{
				for (size_t j = i + 1; j < n; j++)
				{
					// Compute dot products of columns i and j of U using operator()
					T alpha = 0, beta = 0, gamma = 0;
					for (size_t k = 0; k < m; k++)
					{
						T val_i = U(k, i);
						T val_j = U(k, j);
						alpha += val_i * val_i;
						beta  += val_j * val_j;
						gamma += val_i * val_j;
					}
					
					// If columns are already orthogonal, skip rotation
					if (std::abs(gamma) <= tolerance * std::sqrt(alpha * beta)) continue;
					
					changed = true;
					
					// Jacobi rotation calculation
					T zeta = (beta - alpha) / (2.0 * gamma);
					T t = (zeta > 0 ? 1.0 : -1.0) / (std::abs(zeta) + std::sqrt(1.0 + zeta * zeta));
					T cos_theta = 1.0 / std::sqrt(1.0 + t * t);
					T sin_theta = t * cos_theta;
					
					// Update columns of U and V using direct reference access
					for (size_t k = 0; k < m; k++)
					{
						T u_ik = U(k, i);
						T u_jk = U(k, j);
						U(k, i) = cos_theta * u_ik - sin_theta * u_jk;
						U(k, j) = sin_theta * u_ik + cos_theta * u_jk;
					}
					for (size_t k = 0; k < n; k++)
					{
						T v_ik = V(k, i);
						T v_jk = V(k, j);
						V(k, i) = cos_theta * v_ik - sin_theta * v_jk;
						V(k, j) = sin_theta * v_ik + cos_theta * v_jk;
					}
				}
			}
		}
		
		// After orthogonalization, calculate singular values and normalize U
		for (size_t j = 0; j < n; j++)
		{
			T norm = 0;
			for (size_t i = 0; i < m; i++)
			{
				T val = U(i, j);
				norm += val * val;
			}
			norm = std::sqrt(norm);
			s[j] = norm;
			
			if (norm > tolerance)
			{
				for (size_t i = 0; i < m; i++) U(i, j) /= norm;
			}
		}
	}

};
/*
template <size_t m,class T> class matrix<m,m,T> : public matrix_base<m,m,T>
{
protected:
	valarray<T>& V = matrix_base<m,m,T>::V;
	T D;
	array<T,m> CE;
	bool FAD;
	
	void Fad() 
	{    
		matrix<m,m,T> I, N, tN; I.diag() = 1;
		N = I;
		for(int i=1;i<=m;i++)
		{
			tN = (*this) * N;
			CE[i-1] =  tr(tN) * (-1.0/i);
			N = tN + I*(CE[i-1]);
		}
		FAD = 1;
		D = CE[m-1]*pow(-1,m);
	}
	
	matrix<m,m,T> power(size_t ip)
	{
		matrix<m,m,T> oM; oM.diag() = 1;
		for(int i=0; i<ip; i++) oM = oM * (*this);
		return oM;
	}
	
public:
	matrix() : matrix_base<m,m,T>(T(0)), FAD(0) {}
	matrix(matrix_base<m,m,T> iM) : matrix_base<m,m,T>(iM), FAD(0) {}
	matrix(T iS) : matrix_base<m,m,T>(iS), FAD(0) {}
	matrix(initializer_list<T> iV) : matrix_base<m,m,T>(iV), FAD(0) {}
	void operator = (matrix<m,m,T> iM) { V = iM.V; FAD = iM.FAD; D = iM.D; CE = iM.CE; }
	
	ref_select<m,T> diag() { return ref_select<m,T>(V,slice(0,m,m+1)); }
	
	T determinant() { if(!FAD) Fad();return D; };
	
	matrix<m,m,T> inverse()
	{
		if(!FAD) Fad(); 
		matrix<m,m,T> oM;
		if(D==0) { cout<<"Singular Matrix!"<<endl; }
		else for(int i=0;i<m;i++) oM += power(m-i-1) * (i==0 ? 1 : CE[i-1]);
		return oM * (-1/CE[m-1]);
	}
	
	matrix<m,m,T> operator ^ (int ip)
	{
		matrix<m,m,T> oM, iM = (ip<0 ? inverse() : *this);
		oM.diag() = 1;
		for(int i=0; i<abs(ip); i++) oM = oM * iM;
		return oM;
	}
};
*/

template <size_t m,class T> class matrix<m,m,T> : public matrix_base<m,m,T>
{
protected:
	valarray<T>& V = matrix_base<m,m,T>::V;
	valarray<size_t> R, C;
	
	matrix(valarray<T> iV) : matrix_base<m,m,T>(iV) {}
	
	size_t l(valarray<T> iV) { return sqrt(iV.size()); }
	
	T minor(valarray<T> iV,size_t im,size_t in)
	{
		return vdet(iV[valarray<size_t>(R[R<l(iV) && C<l(iV)])!=im && valarray<size_t>(C[R<l(iV) && C<l(iV)])!=in]);
	}
	T cf(valarray<T> iV,size_t im,size_t in) { return pow(-1,im+in)*minor(iV,im,in); }
	
	void circ_l(valarray<T>& iV)
	{
		for(int i=0; i<l(iV); i++)
			iV[slice(l(iV)*i,l(iV),1)] = valarray<T>(iV[slice(l(iV)*i,l(iV),1)]).cshift(1);
	}
	
	void circ_u(valarray<T>& iV)
	{
		for(int i=0; i<l(iV); i++)
			iV[slice(i,l(iV),l(iV))] = valarray<T>(iV[slice(i,l(iV),l(iV))]).cshift(1);
	}
	
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
	
	valarray<T> m_adj(valarray<T> iV) { return m_cf(iV)[gslice(0,{m,m},{1,m})]; }
	
	valarray<T> m_inv(valarray<T> iV) { valarray<T> a = (m_adj(iV)/vdet(iV)); return a; }
	
	valarray<T> m_delta(size_t cn, valarray<T> iV, valarray<T> iC) { iV[slice(cn,l(iV),l(iV))] = iC; return iV; }
	
public:
	matrix() : matrix<m,m,T>(T(0)) {}
	matrix(matrix_base<m,m,T> iM) : matrix_base<m,m,T>(iM)
	{
		R = valarray<size_t>(size_t(0),m*m);
		C = valarray<size_t>(size_t(0),m*m);
		for(int i=0; i<m; i++) { R[matrix_base<m,m,T>::r(i)] = i; C[matrix_base<m,m,T>::c(i)] = i; }
		//cout<<"V = "<<V<<", iM = "<<iM<<endl;

	}
	matrix(T iS) : matrix_base<m,m,T>(iS)
	{
		R = valarray<size_t>(size_t(0),m*m);
		C = valarray<size_t>(size_t(0),m*m);
		for(int i=0; i<m; i++) { R[matrix_base<m,m,T>::r(i)] = i; C[matrix_base<m,m,T>::c(i)] = i; }
	}
	matrix(initializer_list<T> iV) : matrix<m,m,T>(T(0)) { V = iV; }
	
	matrix(const matrix<m,m,T>& iM) : matrix_base<m,m,T>(iM.V), V(matrix_base<m,m,T>::V), R(iM.R), C(iM.C) {}
	
	void operator = (matrix<m,m,T> iM) { V = iM.V; R = iM.R; C = iM.C; }
	
	//T determinant() { return (m>1 ? vdet(V) : V[0]); }
	T determinant()
	{
		matrix<m,m,T> L,D,R;
		LDR(L,D,R);
		T d = 1;
		for(int i=0; i<m; i++) d *= D.diag()[i];
		//cout<<"determinant: "<<vdet(V)<<endl;
		return d;
	}
	
	ref_select<m,T> diag() { return ref_select<m,T>(V,slice(0,m,m+1)); }
	
	//matrix<m,m,T> operator / (matrix<m,m,T> iM) { return *this * iM^(-1); }

	valarray<T>& v() { return V; } 
	
	matrix<m,m,T> operator ^ (int ip)
	{
		//matrix<m,m,T> oM, iM = (ip<0 ? matrix<m,m,T>(m_inv(V)) : *this);
		//cout<<"ip = "<<ip<<endl;
		matrix<m,m,T> oM, iM = (ip<0 ? inv() : *this);
		oM.diag() = 1;
		for(int i=0; i<abs(ip); i++) oM = oM * iM;
		return oM;
	}

	void LDR(matrix<m,m,T> &oL,matrix<m,m,T> &oD,matrix<m,m,T> &oR)
	{
		matrix<m,m,T> L[m],D[m],R[m];
		L[0].diag() = 1; D[0] = *this; R[0].diag() = 1;
		for(int i=1; i<m; i++)
		{
			L[i].diag() = 1; R[i].diag() = 1;
			for(int j=i; j<m; j++)
			{
				L[i][j][i-1] = - (D[i-1][j][i-1] / D[i-1][i-1][i-1]);
				R[i][i-1][j] = - (D[i-1][i-1][j] / D[i-1][i-1][i-1]);
			}
			D[i] = L[i] * D[i-1] * R[i];

			L[i] = - L[i]; L[i].diag() = 1;
			L[i] = L[i-1] * L[i];

			R[i] = - R[i]; R[i].diag() = 1;
			R[i] = R[i] * R[i-1];
			//cout<<"L["<<i<<"] = "<<L[i]<<endl;
			//cout<<"D["<<i<<"] = "<<D[i]<<endl;
			//cout<<"R["<<i<<"] = "<<R[i]<<endl;
		}
		oL = L[m-1]; oD = D[m-1]; oR = R[m-1];

		matrix<m,m,T> Ld, Rd, Ldi, Rdi;
		Ld.diag() = oL.diag(); Rd.diag() = oR.diag();
		for(int i=0; i<m; i++)
		{
			Ldi[i][i] = 1 / oL[i][i];
			Rdi[i][i] = 1 / oR[i][i];
		}
		oL = oL * Ldi; oD = Ld * oD * Rd; oR = Rdi * oR;
	}

	/*
	matrix<m,m,T> inv()
	{
		matrix<m,m,T> L[m],D[m],R[m],P,Q,IV;
		L[0].diag() = 1; D[0] = *this; R[0].diag() = 1;
		P.diag() = 1; Q.diag() = 1;
		for(int i=1; i<m; i++)
		{
			L[i].diag() = 1; R[i].diag() = 1;
			for(int j=i; j<m; j++)
			{
				L[i][j][i-1] = - (D[i-1][j][i-1] / D[i-1][i-1][i-1]);
				R[i][i-1][j] = - (D[i-1][i-1][j] / D[i-1][i-1][i-1]);
			}
			D[i] = L[i] * D[i-1] * R[i];
			//cout<<"D["<<i<<"] = "<<D[i]<<endl;
			P = P * R[i]; Q = L[i] * Q;
		}
		for(int i=0; i<m; i++) IV[i][i] = 1 / D[m-1][i][i];
		IV = P * IV * Q;
		
		//cout<<"I = "<<D[0]*IV<<endl;
		
		return IV;
	}
	*/

	// Optimized Inverse using Gauss-Jordan Elimination
	// Complexity: O(N^3)
	// Access Pattern: Uses (row, col) to avoid heap allocations
	matrix<m, m, T> inv()
	{
		matrix<m, m, T> A = *this;  // Copy original matrix
		matrix<m, m, T> I;
		I.diag() = 1;               // Identity matrix

		// 1. Forward Elimination
		for (size_t i = 0; i < m; i++)
		{
			T pivot = A(i, i);

			// Optional: Add pivoting here for numerical stability if pivot is close to 0

			// Normalize the pivot row
			for (size_t j = 0; j < m; j++)
			{
				A(i, j) /= pivot;
				I(i, j) /= pivot;
			}

			// Eliminate other rows
			for (size_t k = 0; k < m; k++)
			{
				if (k != i)
				{
					T factor = A(k, i);
					for (size_t j = 0; j < m; j++)
					{
						A(k, j) -= factor * A(i, j);
						I(k, j) -= factor * I(i, j);
					}
				}
			}
		}
		return I;
	}

	matrix<m,m,T> inv_adj() { return matrix<m,m,T>(m_inv(V)); }

	void Eigens(valarray<T> &evl,matrix<m,m,T> &evt)
	{
		matrix<m,m,T> Q, Qt, E, nE, L, K, A;
		Q.diag() = 1; Qt.diag() = 1;
		for(int i=0; i<m; i++) { A[i] = (*this)[i]; E[i] = (*this)[i]; }
		//for(int k=0; k<10000; k++)
		int k = 0;
		while(1)
		{
			k++;
			matrix<m,m,T> Qk, nQ;
			for(int i=0; i<m; i++)
			{
				Qk(i) = E(i);
				for(int j=0; j<i; j++) Qk(i) -= (Qk(j).tp()*E(i))/(Qk(j).tp()*Qk(j)) * Qk(j);
				Qk(i) /= pow(Qk(i).tp()*Qk(i),0.5);
			}
			//cout<<"Q*Q' = "<<Q*Q.tp()<<endl;
			nE = Qk.tp() * E * Qk;
			Qt = Qk.tp() * Qt;
			nQ = Q * Qk;

			//cout<<"E = "<<E.diag()<<endl;
			//cout<<"nE = "<<nE.diag()<<endl;

			if(E.diag() == nE.diag()) { cout<<"QR Eigenvalue Algorithm converged at loop "<<k<<endl; break; }
			//if(Q == nQ) { cout<<"QR Eigenvalue Algorithm converged at loop "<<k<<endl; break; }
			//if(k == 9999) { cout<<"QR Eigenvalue Algorithm couldn't converge in 10000 loops"<<endl; break; }
			Q = nQ;
			E = nE;
		}
		cout<<"QR Eigenvalue Algorithm exited with E = \n"<<E<<endl;
		evt = Q;
		evl = E.diag();

		matrix<m,1> v1;
		v1[0] = 1;
		for(int i=0; i<m; i++)
			for(int j=0; j<m; j++)
			{
				L[i][j] = pow(evl[j],i);
				K[i][j] = ((A^i)*v1)[j];
			}
		evt = (L^-1) * K;
		evt = evt.tp();
		for(int i=0; i<m; i++) evt(i) /= pow(evt(i).tp() * evt(i),0.5);
	}
};

template <size_t m,class T> T det(matrix<m,m,T> iM) { return iM.determinant(); }
template <size_t m,class T = double> matrix<m,m,T> eye() { matrix<m,m,T> oM; oM.diag() = 1; return oM; }
template <size_t m,class T> T tr(matrix<m,m,T> iM) { return valarray<T>(iM.diag()).sum(); }

template <size_t n,class T> class matrix<1,n,T> : public matrix_base<1,n,T>
{
protected:
	valarray<T>& V = matrix_base<1,n,T>::V;
	T zero;
public:
	matrix() : matrix_base<1,n,T>(T(0)) {}
	matrix(matrix_base<1,n,T> iM) : matrix_base<1,n,T>(iM) { /*cout<<"V at []: "<<V<<endl;*/ }
	matrix(T iS) : matrix_base<1,n,T>(iS) {}
	matrix(initializer_list<T> iV) : matrix_base<1,n,T>(iV) {}

	matrix(ref_vector<1,n,T> irV) { for(int i=0; i<n; i++) V[i] = irV[i]; }
	
	void operator = (matrix<1,n,T> iM) { V = iM.V; }
	
	T& operator [] (size_t ii) { return V[ii]; }
	
	ref_vector<1,n-1,T> minor(size_t ic=0)
	{
		if(ic>=n)
		{
			cerr<<"minor error: origin out of range, assuming 0"<<endl;
			return ref_vector<1,n-1,T>(V,slice(1,n-1,1));
		}

		valarray<size_t> ns(n);
		for(int i=0; i<n; i++) ns[i] = i;

		return ref_vector<1,n-1,T>(V,valarray<size_t>(ns[ns!=ic]));
	}

	T operator * (matrix<n,1,T> iM)
	{
		T ot = (*this)[0]*0;
		for(int k=0; k<n; k++) ot += V[k] * iM[k];
		return ot;
	}
	valarray<T> operator * (matrix<n,1,valarray<T>> ivM)
	{
		valarray<T> oV(valarray<T>(ivM[0])*0);
		for(int k=0; k<n; k++) oV += V[k] * valarray<T>(ivM[k]);
		return oV;
	}
	template <size_t o> matrix<1,o,T> operator * (matrix_base<n,o,T> iM)
	{
		matrix<1,o,T> Mp(V[0]*0);
		for(int i=0; i<o; i++) for(int k=0; k<n; k++) Mp[i] += V[k] * iM[k][i];
		return Mp;
	}
	template <size_t o> matrix<1,o,valarray<T>> operator * (matrix_base<n,o,valarray<T>> ivM)
	{
		matrix<1,o,valarray<T>> Mp(ivM[0][0]*0);
		for(int i=0; i<o; i++) for(int k=0; k<n; k++) Mp[i] += V[k] * ivM[k][i];
		return Mp;
	}
};
template <size_t n,class T> valarray<T> operator * 
											(matrix<1,n,valarray<T>> ivM,matrix<n,1,T> iM)
{
	valarray<T> oV(valarray<T>(ivM[0])*0);
	for(int k=0; k<n; k++) oV += iM[k] * valarray<T>(ivM[k]);
	return oV;
}


template <size_t m,class T> class matrix<m,1,T> : public matrix_base<m,1,T>
{
protected:
	valarray<T>& V = matrix_base<m,1,T>::V;
	T zero;
public:
	matrix() : matrix_base<m,1,T>(T(0)) {}
	matrix(matrix_base<m,1,T> iM) : matrix_base<m,1,T>(iM) {}
	matrix(T iS) : matrix_base<m,1,T>(iS) {}
	matrix(initializer_list<T> iV) : matrix_base<m,1,T>(iV) {}

	matrix(ref_matrix<m,1,T> irV) { for(int i=0; i<m; i++) V[i] = irV[i]; }
	
	void operator = (matrix<m,1,T> iM) { V = iM.V; }
	
	T& operator [] (size_t ii) { return V[ii]; }

	ref_vector<m-1,1,T> minor(size_t ir=0)
	{
		if(ir>=m)
		{
			cerr<<"minor error: origin out of range, assuming 0"<<endl;
			return ref_vector<m-1,1,T>(V,slice(1,m-1,1));
		}

		valarray<size_t> ns(m);
		for(int i=0; i<m; i++) ns[i] = i;

		return ref_vector<m-1,1,T>(V,valarray<size_t>(ns[ns!=ir]));
	}

};

template <size_t m,class T> class ref_select
{
protected:
	valarray<T> &V;
	valarray<size_t> Vi;
	void operator = (valarray<T> iV) { V[Vi] = iV; }
public:
	ref_select(valarray<T> &iV, valarray<size_t> iVi) : V(iV), Vi(iVi) {}
	ref_select(valarray<T> &iV, slice is) : V(iV), Vi(valarray<size_t>(m))
	{
		valarray<size_t> vix(V.size());
		for(int i=0; i<V.size(); i++) vix[i] = i;
		Vi = valarray<size_t>(vix[is]);
	}
	
	void operator = (T it) { V[Vi] = it; }
	void operator = (array<T,m> iA) { for(int i=0; i<m; i++) V[Vi[i]] = iA[i]; }
	void operator = (ref_select<m,T> iR) { for(int i=0; i<m; i++) V[Vi[i]] = iR[i]; }
	
	bool operator == (ref_select<m,T> iR)
	{
		bool ob = 1;
		for(int i=0; i<m; i++)
			ob = ob && (abs(V[Vi[i]] - iR[i]) < 1e-13 * abs(V[Vi[i]] + iR[i]));
		return ob;
	}

	bool operator != (ref_select<m,T> iR)
	{
		bool ob = 0;
		for(int i=0; i<m; i++) ob = ob || (V[Vi[i]] != iR[i]);
		return ob;
	}

	T& operator [] (size_t ii) { return V[Vi[ii]]; }
	
	operator matrix<m,1,T> ()
	{
		matrix<m,1,T> oM;
		for(int i=0; i<m; i++) oM[i] = V[Vi[i]];
		return oM;
	}
	operator matrix<1,m,T> ()
	{
		matrix<1,m,T> oM;
		for(int i=0; i<m; i++) oM[i] = V[Vi[i]];
		return oM;
	}
	operator valarray<T> () { return valarray<T>(V[Vi]); }
};

template <size_t m,size_t n,class T> class ref_matrix : public matrix<m,n,T>
{
protected:
	valarray<T> &rV;
	valarray<size_t> Vi;
public:
	ref_matrix(valarray<T> &iV, valarray<size_t> iVi) : matrix<m,n,T>(), rV(iV), Vi(iVi) { matrix<m,n,T>::V = rV[iVi]; }
	ref_matrix(valarray<T> &iV, gslice is) : matrix<m,n,T>(), rV(iV), Vi(valarray<size_t>(m*n))
	{
		valarray<size_t> vix(rV.size());
		for(int i=0; i<rV.size(); i++) vix[i] = i;
		Vi = valarray<size_t>(vix[is]);
		matrix<m,n,T>::V = rV[Vi];
	}
	
	void operator = (T it) { rV[Vi] = it; }
	void operator = (array<T,m*n> iA) { for(int i=0; i<m*n; i++) rV[Vi[i]] = iA[i]; }
	
	void operator = (ref_matrix<m,n,T> iM) { rV[Vi] = iM.v(); }
	
	void operator = (matrix<m,n,T> iM) { rV[Vi] = iM.v(); }
	
	ref_vector<1,n,T> operator [] (size_t ir) { return ref_vector<1,n,T>(rV,Vi[slice(n*ir,n,1)]); }
	ref_vector<m,1,T> operator () (size_t ic) { return ref_vector<m,1,T>(rV,Vi[slice(ic,m,n)]); }
	
	ref_select<(m>n ? n : m),T> diag() { return ref_select<(m>n ? n : m),T>(rV,Vi[slice(0,(m>n ? n : m),n+1)]); }
	
	void operator += (T it) { for(int i=0; i<m*n; i++) rV[Vi[i]] += it; }
	void operator -= (T it) { for(int i=0; i<m*n; i++) rV[Vi[i]] -= it; }
	void operator *= (T it) { for(int i=0; i<m*n; i++) rV[Vi[i]] *= it; }
	void operator /= (T it) { for(int i=0; i<m*n; i++) rV[Vi[i]] /= it; }
	void operator += (ref_matrix<m,n,T> iM) { rV[Vi] += iM.v(); }
	void operator -= (ref_matrix<m,n,T> iM) { rV[Vi] -= iM.v(); }
	void operator *= (ref_matrix<n,n,T> iM) { rV[Vi] = ((*this) * matrix<m,n,T>(iM)).v(); }
	void operator += (matrix<m,n,T> iM) { rV[Vi] += iM.v(); }
	void operator -= (matrix<m,n,T> iM) { rV[Vi] -= iM.v(); }
	void operator *= (matrix<n,n,T> iM) { rV[Vi] = ((*this) * iM).v(); }
	/*
	operator matrix<m,n,T> ()
	{
		array<T,m*n> oA;
		cout<<"V at the base class = "<<V<<endl;
		for(int i=0; i<m*n; i++) {  oA[i] = V[Vi[i]];}
		return matrix<m,n,T>(oA);
	}
	*/
};

template <class T> class ref_matrix<1,1,T>
{
protected:
	valarray<T> &rV;
	size_t x;
public:
	ref_matrix(valarray<T> &iV, valarray<size_t> iVi) : rV(iV) { x = iVi[0]; }
	ref_matrix(valarray<T> &iV, gslice is) : rV(iV), x(0)
	{
		valarray<size_t> vix(rV.size());
		for(int i=0; i<rV.size(); i++) vix[i] = i;
		x = valarray<size_t>(vix[is])[0];
	}
	ref_matrix(valarray<T> &iV, size_t ix) : rV(iV), x(ix) {}
	
	void operator = (T it) { rV[x] = it; }
	
	void operator = (ref_matrix<1,1,T>& iM) { rV[x] = iM; }
	
	void operator += (T it) { rV[x] += it; }
	void operator -= (T it) { rV[x] -= it; }
	void operator *= (T it) { rV[x] *= it; }
	void operator /= (T it) { rV[x] /= it; }
	
	ref_vector<1,1,T> operator [] (size_t ir) { return ref_vector<1,1,T>(rV,x); }
	ref_vector<1,1,T> operator () (size_t ic) { return ref_vector<1,1,T>(rV,x); }
	
	operator T& () { return rV[x]; }
};

template <size_t m,size_t n,class T> class ref_vector : public ref_matrix<m,n,T>
{
protected:
	valarray<T> &V = ref_matrix<m,n,T>::rV;
	valarray<size_t> &Vi = ref_matrix<m,n,T>::Vi;
public:
	ref_vector(valarray<T> &iV, valarray<size_t> iVi) : ref_matrix<m,n,T>(iV,iVi) {}
	ref_vector(valarray<T> &iV, slice is) : ref_matrix<m,n,T>(iV,gslice(is.start(),{is.size()},{is.stride()}))
	{
		//cout<<"ref_vector_source: "<<V<<endl;
		//cout<<"ref_matrix_source: "<<ref_matrix<m,n,T>::V<<endl;
		//cout<<"matrix_base_source: "<<matrix_base<m,n,T>::V<<endl;
	}
	//ref_vector(ref_matrix<m,n,T> &iM) : ref_matrix<m,n,T>(iM) {}
	
	void operator = (T it) { V[Vi] = it; }
	void operator = (array<T,m*n> iA) { for(int i=0; i<m*n; i++) V[Vi[i]] = iA[i]; }
	
	void operator = (ref_vector<m,n,T> iM) { V[Vi] = iM.v(); }
	void operator = (ref_vector<n,m,T> iM) { V[Vi] = iM.v(); }
	
	void operator = (matrix<m,n,T> iM) { V[Vi] = iM.v(); }
	
	T& operator [] (size_t ii)
	{
		return V[Vi[ii]];
	}
	//operator matrix<m,n,T> () { return ref_matrix<m,n,T>::operator matrix<m,n,T>(); }
};
template <class T> class ref_vector<1,1,T> : public ref_matrix<1,1,T>
{
protected:
	valarray<T> &V = ref_matrix<1,1,T>::rV;
	size_t x;
public:
	ref_vector(valarray<T> &iV, valarray<size_t> iVi) : ref_matrix<1,1,T>(iV,iVi) { x = iVi[0]; }
	ref_vector(valarray<T> &iV, slice is) : ref_matrix<1,1,T>(iV,gslice(is.start(),{is.size()},{is.stride()})), V(iV), x(0)
	{
		valarray<size_t> vix(V.size()), Vi;
		for(int i=0; i<V.size(); i++) vix[i] = i;
		x = valarray<size_t>(vix[is])[0];
	}
	ref_vector(valarray<T> &iV, size_t ix) : ref_matrix<1,1,T>(iV,x), V(iV), x(ix) {}
	
	void operator = (T it) { V[x] = it; }
	
	void operator = (ref_matrix<1,1,T> iM) { V[x] = iM; }
	
	T& operator [] (size_t ii) { return V[x]; }
	operator T& () { /*cout<<x<<": "<<V[x]<<endl;*/ return V[x]; }
};


template <size_t m,size_t n,class T> matrix<m,n,T> sum(matrix<m,n,valarray<T>> iM)
{
	matrix<m,n,T> oM;
	for(int i=0; i<m; i++) for(int j=0; j<n; j++) oM[i][j] = iM[i][j].sum();
	return oM;
}
template <size_t m,class T> matrix<m,1,T> sum(matrix<m,1,valarray<T>> iM)
{
	matrix<m,1,T> oM;
	for(int i=0; i<m; i++) oM[i] = iM[i].sum();
	return oM;
}
template <size_t n,class T> matrix<1,n,T> sum(matrix<1,n,valarray<T>> iM)
{
	matrix<1,n,T> oM;
	for(int i=0; i<n; i++)oM[i] = iM[i].sum();
	return oM;
}

template <size_t m,size_t n,class T = double> ostream& operator<< (ostream& os, matrix<m,n,T> iM)
{
	//os.precision(4);
	ios init(NULL);
	init.copyfmt(os);
	for(int i=0; i<m; i++)
	{
		os << '\t';
		for(int j=0; j<n; j++) os << setprecision(12) << setw(18) << iM[i][j] << " "; //<< '\t'
		os << endl;
	}
	os.copyfmt(init);
	return os;
}
template <size_t m,class T = double> ostream& operator<< (ostream& os, matrix<m,1,T> iM)
{
	//os.precision(4);
	ios init(NULL);
	init.copyfmt(os);
	for(int i=0; i<m; i++) { os << '\t' << setprecision(12) << setw(14) << iM[i] << " " << endl; }
	os.copyfmt(init);
	return os;
}
template <size_t n,class T = double> ostream& operator<< (ostream& os, matrix<1,n,T> iM)
{
	//os.precision(4);
	ios init(NULL);
	init.copyfmt(os);
	os << '\t';
	for(int j=0; j<n; j++) os << setprecision(12) << setw(14) << iM[j] << " "; //<< '\t'
	os << endl;
	os.copyfmt(init);
	return os;
}
template <size_t m,class T = double> ostream& operator<< (ostream& os, ref_select<m,T> iM)
{
	os << valarray<T>(iM); return os;
}
template <size_t m,size_t n,class T = double> ostream& operator<< (ostream& os, ref_matrix<m,n,T> iM)
{
	os << matrix<m,n,T>(iM); return os;
}
template <class T = double> ostream& operator<< (ostream& os, ref_matrix<1,1,T> iM)
{
	os << T(iM); return os;
}
