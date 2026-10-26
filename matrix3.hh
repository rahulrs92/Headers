#pragma once

#include <iostream>
#include <iomanip>
#include <array>
#include <initializer_list>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <valarray> // Preserved for std::slice/gslice and valarray element support

#include "pick.hh"

using namespace std;

// --- Helper: Convert slice/gslice to index arrays for std::array access ---
template<size_t N>
std::array<size_t, N> slice_to_array(slice s)
{
    std::array<size_t, N> idx;
    size_t curr = s.start();
    size_t str = s.stride();
    for(size_t i=0; i<N; ++i) { idx[i] = curr; curr += str; }
    return idx;
}

template<size_t N>
std::array<size_t, N> gslice_to_array(gslice gs)
{
    std::array<size_t, N> idx;
    valarray<size_t> sz = gs.size();
    valarray<size_t> st = gs.stride();
    size_t start = gs.start();
    for (size_t i = 0; i < N; ++i)
    {
        size_t temp = i, offset = 0;
        for (size_t d = sz.size(); d-- > 0; )
        {
            offset += (temp % sz[d]) * st[d];
            temp /= sz[d];
        }
        idx[i] = start + offset;
    }
    return idx;
}

// --- Forward Declarations ---
template <size_t m, size_t n, class T = double> class matrix_base;
template <size_t m, size_t n, class T> class ref_matrix;
template <size_t m, size_t n, class T> class ref_vector;
template <size_t m, class T> class ref_select;
template <size_t m, size_t n, class T = double> class matrix;

template <size_t m,class T = double> matrix<m,m,T> eye();
template <size_t m,class T = double> matrix<m,m,T> eye(T it);
template <size_t m,class T = double> matrix<m,m,T> eye(array<T,m> ita);

// --- Valarray IO ---
template <class T> ostream& operator<< (ostream& os, const valarray<T>& V)
{
	for(size_t i=0; i<V.size(); i++) os<<V[i]<<" "; return os;
}

// --- LU Decomposition for Inversion and Determinant ---
template <size_t N> class LUSolver
{
private:
    matrix<N, N> LU;
    size_t pivot[N];
    bool decomposed = false;

public:
    // 1. Factorize A (The O(N^3) step)
    LUSolver(matrix_base<N, N> A) : LU(A)
    {
        for (size_t i = 0; i < N; ++i) pivot[i] = i;

        for (size_t i = 0; i < N; ++i)
        {
            size_t max_row = i;
            double max_val = std::abs(LU(i, i));

            // Partial Pivoting for stability
            for (size_t k = i + 1; k < N; ++k)
            {
                if (std::abs(LU(k, i)) > max_val)
                {
                    max_val = std::abs(LU(k, i));
                    max_row = k;
                }
            }

            if (max_row != i)
            {
                for (size_t k = 0; k < N; ++k) std::swap(LU(i, k), LU(max_row, k));
                std::swap(pivot[i], pivot[max_row]);
            }

            // LU Decomposition (Doolittle's algorithm in-place)
            for (size_t j = i + 1; j < N; ++j)
            {
                LU(j, i) /= LU(i, i);
                for (size_t k = i + 1; k < N; ++k) LU(j, k) -= LU(j, i) * LU(i, k);
            }
        }
        decomposed = true;
    }

    // 2. Solve for a Matrix B (The O(K * N^2) step)
    template <size_t K>
    matrix<N, K> solve(const matrix_base<N, K>& B) const
    {
        matrix<N, K> X;

        for (size_t col = 0; col < K; ++col)
        {
            // Apply pivots and Forward Substitution (Ly = b)
            double y[N];
            for (size_t i = 0; i < N; ++i)
            {
                y[i] = B(pivot[i], col);
                for (size_t j = 0; j < i; ++j) y[i] -= LU(i, j) * y[j];
            }

            // Backward Substitution (Ux = y)
            for (int i = N - 1; i >= 0; --i)
            {
                double sum = 0;
                for (size_t j = i + 1; j < N; ++j) sum += LU(i, j) * X(j, col);
                X(i, col) = (y[i] - sum) / LU(i, i);
            }
        }
        return X;
    }
};

// --- Helper for isnan ---
template <class T> bool isnan(const valarray<T>& iv)
{
    for(auto a : iv) if(isnan(a)) { return 1; break; }
    return 0;
}

// --- Matrix Base Class ---
template <size_t m, size_t n, class T> class matrix_base
{
protected:
    std::array<T, m * n> V;

public:
    typedef T value_type;

    matrix_base()
    { 
        // Handles T=double (0.0) and T=valarray (empty valarray)
        if constexpr (!std::is_class<T>::value) V.fill(T(0)); 
    }
    
    // Supports scalar init: matrix<valarray> A(valarray(10));
    matrix_base(T iS) { V.fill(iS); }
    
    matrix_base(initializer_list<T> iV)
    { 
        size_t k = 0;
        for(auto& x : iV) { if(k<m*n) V[k++] = x; }
    }

    matrix_base(const matrix_base<m,n,T>& iM) { V = iM.V; }

    // Accessors
    std::array<T, m*n>& v() { return V; }
    const std::array<T, m*n>& v() const { return V; } // Const correctness added
    
    // Transpose
    matrix<n,m,T> tp() const
    {
        matrix<n,m,T> oM;
        for(size_t i=0; i<m; i++) for(size_t j=0; j<n; j++) oM(j, i) = (*this)(i, j);
        return oM;
    }
    
    void operator = (const matrix_base<m,n,T>& iM) { V = iM.V; }
    
    // Slicing helpers
    slice r(size_t iq) { return slice(iq*n,n,1); }
    slice c(size_t iq) { return slice(iq,m,n); }
    
    constexpr size_t M() { return m; }
    constexpr size_t N() { return n; }

    // Element Access
    T& operator () (size_t p, size_t q) { return V[n*p + q]; }
    const T& operator () (size_t p, size_t q) const { return V[n*p + q]; } // Const overload

    // Reference Views (Updated for std::array storage)
    ref_vector<1,n,T> operator [] (size_t q) { return ref_vector<1,n,T>(V.data(), r(q)); }
    ref_vector<m,1,T> operator () (size_t q) { return ref_vector<m,1,T>(V.data(), c(q)); }

    const matrix<1,n,T> operator [] (size_t q) const { return this->sub<1,n>(q,0); }
    const matrix<m,1,T> operator () (size_t q) const { return this->sub<m,1>(0,q); }
    
    // Sub-matrix Extraction
    template <int o, int p> 
    ref_matrix<(o<0 ? -o : o),(p<0 ? -p : p),T> sub(size_t ir=0, size_t ic=0)
    {
        constexpr size_t q = (o<0 ? -o : o), r = (p<0 ? -p : p);
        if(ir>=m || ic>=n) { cerr<<"sub error: origin out of range"<<endl; ir=0; ic=0; }
        
        size_t tr = (o<0 ? abs(int(o+ir)) : ir), tc = (p<0 ? abs(int(p+ic)) : ic);
        return ref_matrix<q,r,T>(V.data(), gslice(tr*n + tc,{q,r},{n,1}));
    }

    template <int o, int p> 
    matrix<(o<0 ? -o : o),(p<0 ? -p : p),T> sub(size_t ir=0, size_t ic=0) const
    {
        constexpr size_t q = (o<0 ? -o : o), r = (p<0 ? -p : p);
        if(ir>=m || ic>=n) { cerr<<"sub error: origin out of range"<<endl; ir=0; ic=0; }
        
        size_t tr = (o<0 ? abs(int(o+ir)) : ir), tc = (p<0 ? abs(int(p+ic)) : ic);
        auto Vi = gslice_to_array<q*r>(gslice(tr*n + tc,{q,r},{n,1}));
        matrix<q,r,T> oM;
        for(size_t i=0; i<q*r; i++) oM(i/r,i%r) = V[Vi[i]];
        return oM;
        // return ref_matrix<q,r,T>(V.data(), );
    }

    ref_matrix<m-1,n-1,T> minor(size_t ir=0, size_t ic=0)
    {
        if(ir>=m || ic>=n) { cerr<<"minor error: origin out of range"<<endl; ir=0; ic=0; }
        
        std::array<size_t, (m-1)*(n-1)> idx;
        size_t k=0;
        for(size_t i=0; i<m; i++) {
            if(i==ir) continue;
            for(size_t j=0; j<n; j++) {
                if(j==ic) continue;
                idx[k++] = i*n + j;
            }
        }
        return ref_matrix<m-1,n-1,T>(V.data(), idx);
    }
    
    // Row/Col Swaps
    matrix_base<m,n,T>& swap_rows(size_t r1,size_t r2)
	{
        for(size_t j=0; j<n; j++) swap(V[r1*n+j], V[r2*n+j]); return *this;
    }
    matrix_base<m,n,T>& swap_cols(size_t c1,size_t c2)
	{
        for(size_t i=0; i<m; i++) swap(V[i*n+c1], V[i*n+c2]); return *this;
    }

	// Clone
	matrix_base<m,n,T> clone()
	{
		matrix_base<m,n,T> oM;
		for(size_t i=0; i<m; i++) for(size_t j=0; j<n; j++) oM(i,j) = (*this)(i,j);
		return oM;
	}
    
    // --- Standard Arithmetic ---
    matrix<m,n,T> operator + () const { return *this; }
    matrix<m,n,T> operator + (T iS) const { matrix<m,n,T> R; for(size_t i=0; i<m*n; i++) R.v()[i] = V[i] + iS; return R; }
    matrix<m,n,T> operator + (const matrix_base<m,n,T>& iM) const
    { matrix<m,n,T> R; for(size_t i=0; i<m*n; i++) R.v()[i] = V[i] + iM.v()[i]; return R; }
    
    matrix<m,n,T> operator - () const { matrix<m,n,T> R; for(size_t i=0; i<m*n; i++) R.v()[i] = -V[i]; return R; }
    matrix<m,n,T> operator - (T iS) const { matrix<m,n,T> R; for(size_t i=0; i<m*n; i++) R.v()[i] = V[i] - iS; return R; }
    matrix<m,n,T> operator - (const matrix_base<m,n,T>& iM) const
    { matrix<m,n,T> R; for(size_t i=0; i<m*n; i++) R.v()[i] = V[i] - iM.v()[i]; return R; }
    
    void operator += (const matrix_base<m,n,T>& iM) { for(size_t i=0; i<m*n; i++) V[i] += iM.v()[i]; }
    void operator += (T iS) { for(size_t i=0; i<m*n; i++) V[i] += iS; }
    void operator -= (const matrix_base<m,n,T>& iM) { for(size_t i=0; i<m*n; i++) V[i] -= iM.v()[i]; }
    void operator -= (T iS) { for(size_t i=0; i<m*n; i++) V[i] -= iS; }
    void operator *= (T iS) { for(size_t i=0; i<m*n; i++) V[i] *= iS; }
    void operator /= (T iS) { for(size_t i=0; i<m*n; i++) V[i] /= iS; }
    
    bool operator == (const matrix_base<m,n,T>& iM) const { return V == iM.V; }
    bool operator != (const matrix_base<m,n,T>& iM) const { return V != iM.V; }
    
    matrix_base<m,n,bool> operator == (T it) { matrix_base<m,n,bool> R; for(size_t i=0; i<m*n; i++) R.v()[i]=(V[i]==it); return R; }
    matrix_base<m,n,bool> operator != (T it) { matrix_base<m,n,bool> R; for(size_t i=0; i<m*n; i++) R.v()[i]=(V[i]!=it); return R; }
    
    matrix<m,n,T> operator * (T iS) const { matrix<m,n,T> Mp; for(size_t i=0; i<m*n; i++) Mp.v()[i] = V[i] * iS; return Mp; }
    matrix<m,n,T> operator / (T iS) const { matrix<m,n,T> Mp; for(size_t i=0; i<m*n; i++) Mp.v()[i] = V[i] / iS; return Mp; }

    template <size_t o> matrix<m,o,T> operator * (const matrix_base<n,o,T>& iM) const
    {
        matrix<m,o,T> Mp = T(V[0] * 0); // Inits to 0
        for(size_t i=0; i<m; i++) for(size_t j=0; j<o; j++) for(size_t k=0; k<n; k++) 
        {
            T otv = (*this)(i,k) * iM(k,j);
            Mp(i,j) += otv;
            // if(isnan(otv)) cerr<<"matrix3: NaN error at "<<i<<", "<<j<<", "<<k<<": "<<(*this)(i,k)<<endl;
        }
        return Mp;
    }

    // The two cases of matrix division: "Right division" (A/B = A * B^-1) and "Left division" (A/B = B^-1 * A)
    // The cases are determined based on the dimensions of the divisor matrix (square only).
    // Square matrix dividends (m==n) are by default only right-divided.
    // Divisions are handled by LU decomposition for numerical stability
    matrix<m,n,T> operator / (const matrix_base<n,n,T>& iM) const
    {
        LUSolver<m> solver(iM.tp());
        return solver.solve((*this).tp()).tp();
    }
    template <size_t M=m, size_t N=n> auto operator / (const matrix_base<M,M,T>& iM) -> typename enable_if<(M!=N), matrix<M,N,T>>::type const
    {
        LUSolver<m> solver(iM);
        return solver.solve(*this);
    }
    
    T norm() const { T nrm = (*this)(0,0) * 0; for(size_t i=0; i<m*n; i++) nrm += V[i]*V[i]; return sqrt(nrm); }

    T trace() const { T tr = (*this)(0,0) * 0; for(size_t i=0; i<min(m,n); i++) tr += (*this)(i,i); return tr; }

    // --- RESTORED VALARRAY INTEROPERABILITY ---
    // Allows Matrix<double> + Matrix<valarray<double>> = Matrix<valarray<double>>
    
    matrix<m,n,valarray<T>> operator + (const matrix_base<m,n,valarray<T>>& ivM) const
    {
        matrix<m,n,valarray<T>> oM(ivM(0,0)); // Init with size of first valarray
        for(size_t i=0; i<m; i++) for(size_t j=0; j<n; j++) oM(i,j) = V[i*n + j] + ivM(i,j);
        return oM;
    }

    matrix<m,n,valarray<T>> operator - (const matrix_base<m,n,valarray<T>>& ivM) const
    {
        matrix<m,n,valarray<T>> oM(ivM(0,0)); 
        for(size_t i=0; i<m; i++) for(size_t j=0; j<n; j++) oM(i,j) = V[i*n + j] - ivM(i,j);
        return oM;
    }

    // Matrix<double> * Matrix<valarray<double>> -> Matrix<valarray<double>>
    template <size_t o> matrix<m,o,valarray<T>> operator * (const matrix_base<n,o,valarray<T>>& ivM) const
    {
        matrix<m,o,valarray<T>> oM(ivM(0,0)*T(0)); // Init zero matrix of correct valarray size
        for(size_t i=0; i<m; i++) for(size_t j=0; j<o; j++) for(size_t k=0; k<n; k++)
            oM(i,j) += V[i*n + k] * ivM(k,j);
        return oM;
    }

    /**
     * SVD Overload 1: Primary implementation for m >= n.
     * Uses One-Sided Jacobi Method.
     * A = U * diag(S) * V^T, where U is m x n, S is size n, V is n x n.
     */
    template <size_t M = m, size_t N = n>
    std::enable_if_t<(M >= N), void>
    SVD(matrix<m, n, T>& U, valarray<T>& S, matrix<n,n,T>& V, T tolerance = 1e-15)
    {
        static_assert(m >= n, "SVD implementation requires rows >= columns (m >= n)");

        U = (*this).clone();            // U starts as a copy of the original matrix
        S = valarray<T>(T(0), n);       // Initialize singular values to zero
        V = eye<n>();                   // Initialize V as Identity
        
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
                    // Compute the elements of the 2x2 rotation matrix
                    // We need dot products of columns i and j of U
                    T alpha = 0, beta = 0, gamma = 0;
                    for (size_t k = 0; k < m; k++)
                    {
                        alpha += U(k,i) * U(k,i);
                        beta  += U(k,j) * U(k,j);
                        gamma += U(k,i) * U(k,j);
                    }

                    if (abs(gamma) <= tolerance * sqrt(alpha * beta)) continue;

                    changed = true;
                    
                    // Jacobi rotation calculation
                    T zeta = (beta - alpha) / (2.0 * gamma);
                    T t = (zeta > 0 ? 1.0 : -1.0) / (abs(zeta) + sqrt(1.0 + zeta * zeta));
                    T cos_theta = 1.0 / sqrt(1.0 + t * t);
                    T sin_theta = t * cos_theta;

                    // Update columns of U and V
                    for (size_t k = 0; k < m; k++)
                    {
                        T u_ik = U(k,i);
                        T u_jk = U(k,j);
                        U(k,i) = cos_theta * u_ik - sin_theta * u_jk;
                        U(k,j) = sin_theta * u_ik + cos_theta * u_jk;
                    }
                    for (size_t k = 0; k < n; k++)
                    {
                        T v_ik = V(k,i);
                        T v_jk = V(k,j);
                        V(k,i) = cos_theta * v_ik - sin_theta * v_jk;
                        V(k,j) = sin_theta * v_ik + cos_theta * v_jk;
                    }
                }
            }
        }

        // After orthogonalization, singular values are the norms of columns of U
        for (size_t j = 0; j < n; j++) 
        {
            T norm = 0;
            for (size_t i = 0; i < m; i++) norm += U(i,j) * U(i,j);
            norm = sqrt(norm);
            S[j] = norm;

            // Normalize U columns to make them orthonormal
            if (norm > tolerance) for (size_t i = 0; i < m; i++) U(i,j) /= norm;
        }
    }

    /**
     * SVD Overload 2: Transpose delegation for n > m.
     * Computes SVD of A^T (n x m) using the m >= n method, then swaps U and V.
     * Output shapes: U is m x m, S is size m, V is n x m.
     */
    template <size_t M = m, size_t N = n>
    typename std::enable_if_t<(M < N), void>
    SVD(matrix<m, m, T>& U, valarray<T>& S, matrix<n, m, T>& V, T tolerance = 1e-15)
    {
        matrix<n, m, T> At = (*this).tp();
        
        matrix<n, m, T> U_temp;
        matrix<m, m, T> V_temp;

        // Calls SVD on A^T (n x m), which triggers Overload 1 because n >= m
        At.SVD(U_temp, S, V_temp, tolerance);

        // Swap U and V outputs (A = V_temp * diag(S) * U_temp^T)
        U = V_temp;  // size m x m
        V = U_temp;  // size n x m
    }

    // Full pseudo inversion using SVD for singular matrices
    matrix<n,m,T> pseudo_inverse(T tolerance = 1e-15) const
    {
        T norm_val = this->norm(); // Scale factor for numerical stability
        if (norm_val == T(0)) return matrix<n,m,T>(); // Guard against zero matrix
        
        T sf = norm_val / std::min(m, n);

        if constexpr (m >= n)
        {
            matrix<m,n,T> U;
            valarray<T> Sv;
            matrix<n,n,T> V, S_inv;

            (*this / sf).SVD(U, Sv, V);

            for (size_t i = 0; i < n; ++i) S_inv(i,i) = (std::abs(Sv[i]) > tolerance) ? (T(1) / Sv[i]) : T(0);

            matrix<n,m,T> pinv = V * S_inv * U.tp();
            return pinv / sf;
        } 
        else // m < n
        {
            matrix<m,m,T> U, S_inv;
            valarray<T> Sv;
            matrix<n,m,T> V;

            (*this / sf).SVD(U, Sv, V);

            for (size_t i = 0; i < m; ++i) S_inv(i,i) = (std::abs(Sv[i]) > tolerance) ? (T(1) / Sv[i]) : T(0);

            matrix<n,m,T> pinv = V * S_inv * U.tp();
            return pinv / sf;
        }
    }

};

// --- Free Functions for Matrix/Valarray mixing ---

template <size_t m,size_t n,class T> matrix<m,n,T> operator + (T it,matrix_base<m,n,T> iM) { return iM + it; }
template <size_t m,size_t n,class T> matrix<m,n,T> operator - (T it,matrix_base<m,n,T> iM) { return -(iM - it); }
template <size_t m,size_t n,class T> matrix<m,n,T> operator * (T it,matrix_base<m,n,T> iM) { return iM * it; }
// template <size_t m,size_t n,class T> matrix<m,n,T> operator / (T it,matrix_base<m,n,T> iM)
// {
//     matrix<m,n,T> R; for(size_t i=0; i<m*n; i++) R.v()[i] = it / iM.v()[i]; return R;
// }

// Valarray + Matrix<Valarray>
template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator + 
                                                            (valarray<T> iV,matrix_base<m,n,valarray<T>> iM) { return iM + iV; }
template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator - 
                                                            (valarray<T> iV,matrix_base<m,n,valarray<T>> iM) { return -(iM - iV); }
template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator * (valarray<T> iV,matrix_base<m,n,valarray<T>> iM) { return iM * iV; }

// Scalar + Matrix<Valarray>
template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator + (T it,matrix_base<m,n,valarray<T>> ivM)
{
    matrix<m,n,valarray<T>> oM(ivM);
    for(size_t i=0; i<m*n; i++) oM.v()[i] += it; 
    return oM;
}
template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator + (matrix_base<m,n,valarray<T>> ivM,T it) { return it + ivM; }
template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator - (T it,matrix_base<m,n,valarray<T>> ivM)
{
    matrix<m,n,valarray<T>> oM(ivM);
    for(size_t i=0; i<m*n; i++) oM.v()[i] = it - oM.v()[i]; 
    return oM;
}
template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator - (matrix_base<m,n,valarray<T>> ivM,T it) { return -(it - ivM); }
template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator * (T it,matrix_base<m,n,valarray<T>> ivM)
{
    matrix<m,n,valarray<T>> oM(ivM);
    for(size_t i=0; i<m*n; i++) oM.v()[i] *= it;
    return oM;
}
template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator * (matrix_base<m,n,valarray<T>> ivM,T it) { return it * ivM; }
// template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator / (T it,matrix_base<m,n,valarray<T>> ivM)
// {
//     matrix<m,n,valarray<T>> oM(ivM);
//     for(size_t i=0; i<m*n; i++) oM.v()[i] = it / oM.v()[i];
//     return oM;
// }

// Matrix<Valarray> + Matrix<Scalar> (Reverse of member function)
template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator + 
                                            (matrix_base<m,n,valarray<T>> ivM,matrix_base<m,n,T> iM) { return iM + ivM; }
template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> operator - 
                                            (matrix_base<m,n,valarray<T>> ivM,matrix_base<m,n,T> iM) { return -(iM - ivM); }

// Matrix<Valarray> * Matrix<Scalar>
template <size_t m,size_t n,size_t o,class T> matrix<m,o,valarray<T>> operator * (matrix_base<m,n,valarray<T>> ivM,matrix_base<n,o,T> iM)
{
    matrix<m,o,valarray<T>> oM(ivM(0,0)*T(0));
    for(size_t i=0; i<m; i++) for(size_t j=0; j<o; j++) for(size_t k=0; k<n; k++) 
        oM(i,j) += ivM(i,k) * iM(k,j);
    return oM;
}

// Pull functions (Slicing "Depth" of a valarray matrix)
template <size_t m,size_t n,class T> matrix<m,n,T> pull(matrix_base<m,n,valarray<T>> iM,size_t iq)
{
    matrix<m,n,T> oM;
    for(size_t i=0; i<m; i++) for(size_t j=0; j<n; j++) oM(i,j) = iM(i,j)[iq];
    return oM;
}

template <size_t m,size_t n,class T> matrix<m,n,valarray<T>> pull(matrix_base<m,n,valarray<T>> iM,valarray<size_t> iv)
{
    matrix<m,n,valarray<T>> oM(valarray<T>(iv.size()));
    for(size_t i=0; i<m; i++) for(size_t j=0; j<n; j++) oM(i,j) = iM(i,j)[iv];
    return oM;
}

template <size_t m,size_t n,size_t o,class T> matrix<m,n+o,T> cat(matrix_base<m,n,T> iM1,matrix_base<m,o,T> iM2)
{
    matrix<m,n+o,T> oM;
    for(size_t i=0; i<m; i++) for(size_t j=0; j<n; j++) oM(i,j) = iM1(i,j);
    for(size_t i=0; i<m; i++) for(size_t j=0; j<o; j++) oM(i,j+n) = iM2(i,j);
    return oM;
}

// Cross Product Specialization
template <class T> matrix<3,1,T> operator*(matrix<3,1,T> a, matrix<3,1,T> b)
{
    matrix<3,1,T> c;
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return c;
}

template <size_t m,class T> matrix<m,m,T> eye() { matrix<m,m,T> oM; for(size_t i=0; i<m; i++) oM(i,i) = T(1); return oM; }
template <size_t m,class T> matrix<m,m,T> eye(T it) { matrix<m,m,T> oM; for(size_t i=0; i<m; i++) oM(i,i) = it; return oM; }
template <size_t m,class T> matrix<m,m,T> eye(array<T,m> ita) { matrix<m,m,T> oM; for(size_t i=0; i<m; i++) oM(i,i) = ita[i]; return oM; }

// --- Concrete Matrix Class ---
template <size_t m,size_t n,class T> class matrix : public matrix_base<m,n,T>
{
public:
    using matrix_base<m,n,T>::matrix_base;
    matrix() : matrix_base<m,n,T>() {}
    matrix(const matrix_base<m,n,T>& iM) : matrix_base<m,n,T>(iM) {}
    matrix(T iS) : matrix_base<m,n,T>(iS) {}
    matrix(initializer_list<T> iV) : matrix_base<m,n,T>(iV) {}

};

// --- Square Matrix Specialization ---
template <size_t m,class T> class matrix<m,m,T> : public matrix_base<m,m,T>
{
protected:
    using matrix_base<m,m,T>::V;
public:
    using matrix_base<m,m,T>::matrix_base;
    
    matrix(const matrix_base<m,m,T>& iM) : matrix_base<m,m,T>(iM) {}

    T determinant()
    {
        matrix<m,m,T> L,D,R;
        LDR(L,D,R);
        T d = 1;
        for(size_t i=0; i<m; i++) d *= D.diag()[i];
        return d;
    }
    
    ref_select<m,T> diag() { return ref_select<m,T>(V.data(), slice(0,m,m+1)); }
    
    matrix<m,m,T> operator ^ (int ip)
    {
        matrix<m,m,T> oM((*this)(0,0) * 0), iM = (ip<0 ? inv() : *this);
        oM.diag() = pow(oM(0,0),0);
        for(int i=0; i<abs(ip); i++) oM = oM * iM;
        return oM;
    }

    matrix<m,m,T> divide_left(const matrix_base<m,m,T>& B)
    {
        LUSolver<m> solver(B);
        return solver.solve(*this);
    }
    
    void LDR(matrix<m,m,T> &oL,matrix<m,m,T> &oD,matrix<m,m,T> &oR)
    {
        matrix<m,m,T> L[m],D[m],R[m];
        L[0].diag() = 1; D[0] = *this; R[0].diag() = 1;
        for(size_t i=1; i<m; i++)
        {
            L[i].diag() = 1; R[i].diag() = 1;
            for(size_t j=i; j<m; j++)
            {
                L[i](j,i-1) = - (D[i-1](j,i-1) / D[i-1](i-1,i-1));
                R[i](i-1,j) = - (D[i-1](i-1,j) / D[i-1](i-1,i-1));
            }
            D[i] = L[i] * D[i-1] * R[i];
            L[i] = - L[i]; L[i].diag() = 1; L[i] = L[i-1] * L[i];
            R[i] = - R[i]; R[i].diag() = 1; R[i] = R[i] * R[i-1];
        }
        oL = L[m-1]; oD = D[m-1]; oR = R[m-1];
        matrix<m,m,T> Ld, Rd, Ldi, Rdi;
        Ld.diag() = oL.diag(); Rd.diag() = oR.diag();
        for(size_t i=0; i<m; i++) { Ldi(i,i) = 1/oL(i,i); Rdi(i,i) = 1/oR(i,i); }
        oL = oL * Ldi; oD = Ld * oD * Rd; oR = Rdi * oR;
    }

    matrix<m, m, T> inv()
    {
        matrix<m, m, T> A = (*this).clone(), I;
        for (size_t i = 0; i < m; i++) for (size_t j = 0; j < m; j++) I(i,j) = A(i,j) * 0;
        I.diag() = pow(I(0,0),0);
        for (size_t i = 0; i < m; i++)
		{
            T pivot = A(i, i);
            for (size_t j = 0; j < m; j++) { A(i, j) /= pivot; I(i, j) /= pivot; }
            for (size_t k = 0; k < m; k++)
			{
                if (k != i)
				{
                    T factor = A(k, i);
                    for (size_t j = 0; j < m; j++) { A(k, j) -= factor * A(i, j); I(k, j) -= factor * I(i, j); }
                }
            }
            // cout << "I(0,0): " << I(0,0) << endl;
        }
        return I;
    }

    void Eigens(valarray<T> &evl, matrix<m,m,T> &evt)
    {
        // Requires T to be scalar for logic, but implementation provided as requested
        // Note: Using std::array for temp storage would be better, but interface asks for valarray ref
        matrix<m,m,T> Q, Qt, E, nE, L, K, A;
        Q.diag() = 1; Qt.diag() = 1;
        for(size_t i=0; i<m; i++) { A(i) = (*this)(i); E(i) = (*this)(i); }
        
        while(1) {
            matrix<m,m,T> Qk, nQ;
            for(size_t i=0; i<m; i++) {
                Qk(i) = E(i);
                for(size_t j=0; j<i; j++) Qk(i) -= (Qk(j).tp()*E(i)) / (Qk(j).tp()*Qk(j)) * Qk(j);
                Qk(i) /= pow((Qk(i).tp()*Qk(i)), 0.5);
            }
            nE = Qk.tp() * E * Qk; nQ = Q * Qk;
            bool conv = true;
            for(size_t i=0; i<m; i++) if(std::abs(E.diag()[i]-nE.diag()[i]) > 1e-9) conv=false;
            if(conv) break;
            Q = nQ; E = nE;
        }
        evt = Q;
        evl.resize(m); for(size_t i=0; i<m; i++) evl[i] = E.diag()[i];
    }
};

template <size_t m,class T> T det(matrix<m,m,T> iM) { return iM.determinant(); }
template <size_t m,class T> T tr(matrix<m,m,T> iM) { T s=0; auto d=iM.diag(); for(size_t i=0; i<m; i++) s+=d[i]; return s; }

// --- Linear System Solver for division of vector specializations ---
template <size_t N>
matrix<N, 1> solve_linear_system(matrix<N, N> A, matrix<N, 1> b)
{
    // We work on copies of A and b (passed by value) to keep 
    // the original data intact while performing in-place elimination.
    
    // 1. Forward Elimination with Partial Pivoting
    for (size_t i = 0; i < N; ++i)
    {
        // Find the pivot row
        size_t pivot = i;
        for (size_t j = i + 1; j < N; ++j) if (std::abs(A(j, i)) > std::abs(A(pivot, i))) pivot = j;
       
        // Swap rows in A and b (only if necessary)
        if (pivot != i)
        {
            for (size_t k = i; k < N; ++k) std::swap(A(i, k), A(pivot, k));
            std::swap(b(i, 0), b(pivot, 0));
        }

        // Eliminate entries below the pivot
        for (size_t j = i + 1; j < N; ++j)
        {
            double factor = A(j, i) / A(i, i);
            for (size_t k = i + 1; k < N; ++k) A(j, k) -= factor * A(i, k);
            b(j, 0) -= factor * b(i, 0);
        }
    }

    // 2. Back Substitution
    matrix<N, 1> x;
    for (int i = N - 1; i >= 0; --i)
    {
        double sum = 0;
        for (size_t j = i + 1; j < N; ++j) sum += A(i, j) * x(j, 0);
        x(i, 0) = (b(i, 0) - sum) / A(i, i);
    }

    return x;
}

// --- Row Vector Specialization ---
template <size_t n,class T> class matrix<1,n,T> : public matrix_base<1,n,T>
{
protected:
    using matrix_base<1,n,T>::V;
public:
    using matrix_base<1,n,T>::matrix_base;
    
    matrix(const matrix_base<1,n,T>& iM) : matrix_base<1,n,T>(iM) {}
    matrix(ref_vector<1,n,T> irV); 
    
    void operator = (const matrix<1,n,T>& iM) { V = iM.V; }

    T& operator [] (size_t ii) { return V[ii]; }
	const T operator [] (size_t ii) const { return V[ii]; }

    matrix<1,n,T> operator / (const matrix<n,n,T>& iM)
    {
        matrix<n,1,T> col = (*this).tp();
        matrix<n,1,T> res_col = solve_linear_system(iM, col);
        return res_col.tp();
    }
    
    ref_vector<1,n-1,T> minor(size_t ic=0) { return matrix_base<1,n,T>::minor(0, ic); }
    
    T operator * (const matrix<n,1,T>& iM)
    {
        T ot = T(V[0] * 0); for(size_t k=0; k<n; k++) ot += V[k] * iM[k]; return ot;
    }
    valarray<T> operator * (const matrix<n,1,valarray<T>>& ivM)
    {
        valarray<T> oV(ivM[0].size()); for(size_t k=0; k<n; k++) oV += V[k] * ivM[k]; return oV;
    }
    template <size_t o> matrix<1,o,T> operator * (const matrix_base<n,o,T>& iM)
    {
        matrix<1,o,T> Mp = T(V[0] * 0); for(size_t i=0; i<o; i++) for(size_t k=0; k<n; k++) Mp[i] += V[k] * iM(k,i); return Mp;
    }
    template <size_t o> matrix<1,o,valarray<T>> operator * (const matrix_base<n,o,valarray<T>>& ivM)
    {
        matrix<1,o,valarray<T>> Mp(ivM(0,0)*T(0)); for(size_t i=0; i<o; i++) for(size_t k=0; k<n; k++) Mp[i] += V[k] * ivM(k,i); return Mp;
    }
};

template <size_t n,class T> valarray<T> operator * (const matrix<1,n,valarray<T>>& ivM,const matrix<n,1,T>& iM)
{
    valarray<T> oV(ivM[0].size()); for(size_t k=0; k<n; k++) oV += iM[k] * ivM[k]; return oV;
}

// --- Col Vector Specialization ---
template <size_t m,class T> class matrix<m,1,T> : public matrix_base<m,1,T>
{
protected:
    using matrix_base<m,1,T>::V;
public:
    using matrix_base<m,1,T>::matrix_base;
    
    matrix(const matrix_base<m,1,T>& iM) : matrix_base<m,1,T>(iM) {}
    matrix(ref_matrix<m,1,T> irV);
    
    void operator = (const matrix<m,1,T>& iM) { V = iM.V; }

    T& operator [] (size_t ii) { return V[ii]; }
	const T operator [] (size_t ii) const { return V[ii]; }

    matrix<m,1,T> operator / (const matrix<m,m,T>& iM) { return solve_linear_system(iM, *this); }
    matrix<m,1,T> operator / (T iS) { return (*this).matrix_base<m,1,T>::operator/(iS); }
    
    ref_vector<m-1,1,T> minor(size_t ir=0) { return matrix_base<m,1,T>::minor(ir, 0); }
};

// --- Sum functions for valarray matrices ---
template <size_t m,size_t n,class T> matrix<m,n,T> sum(const matrix<m,n,valarray<T>>& iM) {
    matrix<m,n,T> oM; for(size_t i=0; i<m; i++) for(size_t j=0; j<n; j++) oM(i,j) = iM(i,j).sum(); return oM;
}

// --- Reference Views (Impl) ---

template <size_t m,class T> class ref_select
{
protected:
    T* ptr;
    std::array<size_t, m> Vi;
public:
    ref_select(T* iPtr, slice is) : ptr(iPtr), Vi(slice_to_array<m>(is)) {}
    
    void operator = (T it) { for(size_t i=0; i<m; i++) ptr[Vi[i]] = it; }
    void operator = (const std::array<T,m>& iA) { for(size_t i=0; i<m; i++) ptr[Vi[i]] = iA[i]; }
    void operator = (const ref_select<m,T>& iR) { for(size_t i=0; i<m; i++) ptr[Vi[i]] = iR[i]; }
    // void operator = (const matrix<m,1,T>& iM) { for(size_t i=0; i<m; i++) ptr[Vi[i]] = iM[i]; }
    // void operator = (const matrix<1,m,T>& iM) { for(size_t i=0; i<m; i++) ptr[Vi[i]] = iM[i]; }
    
    T& operator [] (size_t ii) { return ptr[Vi[ii]]; }
    const T& operator [] (size_t ii) const { return ptr[Vi[ii]]; } // Const added

    void operator += (T it) { for(size_t i=0; i<m*n; i++) ptr[Vi[i]] += it; }
    void operator -= (T it) { for(size_t i=0; i<m*n; i++) ptr[Vi[i]] -= it; }
    void operator *= (T it) { for(size_t i=0; i<m*n; i++) ptr[Vi[i]] *= it; }
    void operator /= (T it) { for(size_t i=0; i<m*n; i++) ptr[Vi[i]] /= it; }
    
    operator matrix<m,1,T> () { matrix<m,1,T> oM; for(size_t i=0; i<m; i++) oM[i] = ptr[Vi[i]]; return oM; }
    operator valarray<T> () { valarray<T> V(m); for(size_t i=0; i<m; i++) V[i] = ptr[Vi[i]]; return V; }
};

template <size_t m,size_t n,class T> class ref_matrix : public matrix<m,n,T>
{
protected:
    T* ptr;
    std::array<size_t, m*n> Vi;
public:
    ref_matrix(T* iPtr, std::array<size_t, m*n> iVi) : matrix<m,n,T>(), ptr(iPtr), Vi(iVi)
    {
        for(size_t i=0; i<m*n; i++) this->V[i] = ptr[Vi[i]]; // Sync view
    }
    ref_matrix(T* iPtr, gslice is) : matrix<m,n,T>(), ptr(iPtr)
    {
        Vi = gslice_to_array<m*n>(is);
        for(size_t i=0; i<m*n; i++) this->V[i] = ptr[Vi[i]];
    }
    
    void operator = (T it) { for(size_t i=0; i<m*n; i++) ptr[Vi[i]] = it; }
    void operator = (const std::array<T,m*n>& iA) { for(size_t i=0; i<m*n; i++) ptr[Vi[i]] = iA[i]; }
    void operator = (const ref_matrix<m,n,T>& iM) { for(size_t i=0; i<m*n; i++) ptr[Vi[i]] = iM.v()[i]; }
    void operator = (const matrix<m,n,T>& iM) { for(size_t i=0; i<m*n; i++) ptr[Vi[i]] = iM.v()[i]; }

	// Add these to handle assignments from other views
    void operator = (const ref_select<m*n, T>& iS) { for (size_t i=0; i<m*n; i++) ptr[Vi[i]] = iS[i]; }
    
    ref_vector<1,n,T> operator [] (size_t ir)
    { 
        std::array<size_t, n> sub_idx;
        for(size_t k=0; k<n; k++) sub_idx[k] = Vi[ir*n + k];
        return ref_vector<1,n,T>(ptr, sub_idx);
    }
    ref_vector<m,1,T> operator () (size_t ic)
    { 
        std::array<size_t, m> sub_idx;
        for(size_t k=0; k<m; k++) sub_idx[k] = Vi[k*n + ic];
        return ref_vector<m,1,T>(ptr, sub_idx);
    }
    
    ref_select<(m>n ? n : m),T> diag()
    { 
        return ref_select<(m>n ? n : m),T>(ptr, slice(Vi[0], (m>n ? n : m), n+1)); 
    }
    
    void operator += (T it) { for(size_t i=0; i<m*n; i++) ptr[Vi[i]] += it; }
    void operator -= (T it) { for(size_t i=0; i<m*n; i++) ptr[Vi[i]] -= it; }
    void operator *= (T it) { for(size_t i=0; i<m*n; i++) ptr[Vi[i]] *= it; }
    void operator /= (T it) { for(size_t i=0; i<m*n; i++) ptr[Vi[i]] /= it; }
    void operator += (const matrix<m,n,T>& iM) { for(size_t i=0; i<m*n; i++) ptr[Vi[i]] += iM.v()[i]; }
    void operator -= (const matrix<m,n,T>& iM) { for(size_t i=0; i<m*n; i++) ptr[Vi[i]] -= iM.v()[i]; }
};

template <class T> class ref_matrix<1,1,T>
{
protected:
    T* ptr;
    size_t x;
public:
    ref_matrix(T* iPtr, std::array<size_t, 1> iVi) : ptr(iPtr), x(iVi[0]) {}
    ref_matrix(T* iPtr, size_t ix) : ptr(iPtr), x(ix) {}
    void operator = (T it) { ptr[x] = it; }
    operator T& () { return ptr[x]; }
};

template <size_t m,size_t n,class T> class ref_vector : public ref_matrix<m,n,T>
{
public:
    using ref_matrix<m,n,T>::ref_matrix;
    ref_vector(T* iPtr, slice is) : ref_matrix<m,n,T>(iPtr, slice_to_array<m*n>(is)) {}
    ref_vector(T* iPtr, std::array<size_t, m*n> idx) : ref_matrix<m,n,T>(iPtr, idx) {}

    void operator = (T it) { ref_matrix<m,n,T>::operator=(it); }
    void operator = (const matrix<m,n,T>& iM) { ref_matrix<m,n,T>::operator=(iM); }
    void operator = (const ref_select<n,T>& iM) { ref_matrix<m,n,T>::operator=(iM); }
    T& operator [] (size_t ii) { return this->ptr[this->Vi[ii]]; }
};

template <class T> class ref_vector<1,1,T> : public ref_matrix<1,1,T>
{
public:
    using ref_matrix<1,1,T>::ref_matrix;
    ref_vector(T* iPtr, slice is) : ref_matrix<1,1,T>(iPtr, slice_to_array<1>(is)[0]) {}
    void operator = (T it) { this->ptr[this->x] = it; }
    T& operator [] (size_t) { return this->ptr[this->x]; }
};

// --- Deferred Implementation of Constructors ---
template <size_t n,class T> matrix<1,n,T>::matrix(ref_vector<1,n,T> irV) : matrix_base<1,n,T>()
{
    for(size_t i=0; i<n; i++) V[i] = irV[i];
}
template <size_t m,class T> matrix<m,1,T>::matrix(ref_matrix<m,1,T> irV) : matrix_base<m,1,T>()
{
    for(size_t i=0; i<m; i++) V[i] = irV[i];
}

// --- Matrix IO ---
template <size_t m,size_t n,class T = double> ostream& operator<< (ostream& os, matrix<m,n,T> iM)
{
    ios init(NULL); init.copyfmt(os);
    for(size_t i=0; i<m; i++)
    {
        os << '\t';
        for(size_t j=0; j<n; j++) os << setprecision(12) << setw(18) << iM(i,j) << " ";
        os << endl;
    }
    os.copyfmt(init);
    return os;
}
template <size_t m,size_t n,class T = double> ostream& operator<< (ostream& os, ref_matrix<m,n,T> iM)
{
    os << matrix<m,n,T>(iM); return os;
}
template <size_t m,class T = double> ostream& operator<< (ostream& os, ref_select<m,T> iM)
{
    os << valarray<T>(iM); return os;
}


/*

BE#_H1
ABE#
BE#_H2

BE#

*/