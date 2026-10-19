#pragma once

#include<iostream>
#include<cmath>
#include<iostream>
#include<vector>
#include<valarray>
#include<array>
#include<algorithm>
using namespace std;

class Matrix
{
protected:
    int n=0;
    valarray<double> m;
    //Identity Matrix....
    Matrix Id() 
    {
        valarray<double> t(0.0,n*n);
        t[slice(0,n,n+1)]=1.0;
        return Matrix(t);
    }

public:
                    
    Matrix(vector<vector<double>> M)
    {
        int k=0;
        n=M.size();
        m = valarray<double>(n*n);
        for(int i=0;i<n;i++)
            for(int j=0;j<n;j++)
            {
                m[k]=M[i][j];
                k++;
            }
    }
    
    
    int getn() { return n; }

    Matrix(valarray<double> M)
    {   
        m=M;
        n=pow(m.size(),0.5);
        int k=0;        
    }
    Matrix() 
    {
        n=3;
        m={1,0,0,0,1,0,0,0,1};
    }
    
    //Accessing The elements.....
    double& operator () (int i, int j) 
    {
        return m[i*n+j];
    }
     
    
    //Matrix Addition.......
    Matrix operator + (Matrix M)
    {
        return Matrix(m+M.m);
    }
    
    //Matrix multiplication......
    Matrix operator * (double c) 
    {
       return Matrix(c*m); 
    }
    //Matrix multiplication......
    Matrix operator * (Matrix M)
    {
        int k=0;
        valarray<double> C(n*n);
        for(int i=0;i<n;i++)
            for(int j=0;j<n;j++)
            {
                C[k]=(valarray<double>(m[slice(i*n,n,1)])*valarray<double>(M.m[slice(j,n,n)])).sum();
                k++;                
            }
        return C;
    }
    
    //power of a Matrix....
    Matrix operator ^ (unsigned int p)
    {
        Matrix ind=Id();
        for(int i=0;i<p;i++)
            ind=ind**this;
        return ind;
    }
    
    //Trace of Matrix......
    double Tr()
    {
        return valarray<double>(m[slice(0,n,n+1)]).sum();
    }
   
    //Assignment operator.....
    Matrix operator = (vector<vector<double>> V)
    {
        return Matrix(V);
    }
    
    Matrix operator = (valarray<double> V) 
    {
        return Matrix(V);
    }
       
    //Faddeev-LeVerrier Algorithm.....
    valarray<double> Souriau() 
    {             
        valarray<double> a(n+1);
        Matrix I=Id();
        Matrix N=I;
        a[0]=1;
        for(int i=1;i<=n;i++)
        {
            auto tN=(*this)*N;
            a[i]=(-1.0/i)*(tN.Tr());
            N=tN+I*a[i];
        }
        return a;
    }       

    
    //Permutation operation........
    Matrix Permute(Matrix M,int i,int j)   
    { 
        valarray<double> t;
        int n=M.n;
        t=M.m[slice(i*n,n,1)];
        M.m[slice(i*n,n,1)]=M.m[slice(j*n,n,1)];
        M.m[slice(j*n,n,1)]=t;
        
        t=M.m[slice(i,n,n)];
        M.m[slice(i,n,n)]=M.m[slice(j,n,n)];
        M.m[slice(j,n,n)]=t; 
        
        return M;
    }  
    
    
    //C^-1AC........    
    Matrix C(Matrix A,int q)     
    {
       Matrix B=A; 
       int N=A.n-1,p=N-q;
       for(int i=0;i<=N;i++)
           for(int j=0;j<=N;j++)
               if(j==p)
                   B(i,j)=A(i,p)/A(p+1,p);
               else
                   B(i,j)=A(i,j)-A(i,p)*A(p+1,j)/A(p+1,p);
       double t2,t1; 
       for(int j=0;j<=N;j++)
       {
           t2=0;
           t1=0;
           if(j>=p && j<N)                                      
           {   
               for(int k=0;k<=p;k++)
               {
                   t1=t1+A(p+1,k)*B(k,j);
               }       
               B(p,j)=t1+A(p+1,j+1);
           }
           else
           {
               for(int k=0;k<=p;k++)
               {
                   t2=t2+A(p+1,k)*B(k,j);
               }
               B(p,j)=t2;
           }
       }
        for(int i=0;i<=N;i++)
        {
            for(int j=0;j<=N;j++) cout<<B(i,j)<<"\t";
            cout<<endl;
        }
        cout<<endl;
        

       return B;
    }
    
    
    //D1 D3.......
    array<valarray<double>,2> Diagonal(Matrix M,int i) 
    {
        return {M.m[gslice(0,{i,i},{M.n,1})],M.m[gslice(i*M.n+i,{M.n-i,M.n-i},{M.n,1})]};
    } 
    
    
    //check for pivot element......
    int Pivot(Matrix M,int r)    
    {
      int k;  
      for(k=0;k<r-1;k++)
          if(valarray<double>(M.m[slice(r*M.n,r-1,1)])[k]!=0)
          {
              return k;
              break;
          }      
      return k;  
    }
    
    //The Method of Danilevsky......
    vector<valarray<double>> Danilevsky(Matrix A)  
    {
      Matrix M=A;  
      vector<valarray<double>> D;  
      if(M.n==0)
      {
          cout<<"Not Applicable for Danilevsky"<<endl;
          return D;
      }
      if(M.n==1)
      {
          D.push_back({1.0,-m[0]});
          return D;
      }
      if(M.n==2)
      {
          D.push_back({1.0,-Matrix(m).Tr(),m[0]*m[3]-m[1]*m[2]});
          return D;
      }
      int t=0;  
      if(M.n>2)
      {
        for(int q=1;q<=M.n-1;q++)
        {
          t=Pivot(M,q);  
          if(M(M.n-q+1,M.n-q)!=0)
              M=C(M,q);              
          else if(t<q-1)
          {              
              M=Permute(M,q,t);
              M=C(M,q);
          }
          else if(t=q-1)
          {
              M=Diagonal(M,M.n-q+1)[0];
              D.push_back(Diagonal(M,M.n-q+1)[1][slice(0,M.n-q,1)]);
              vector<valarray<double>> D1=Danilevsky(M);
              for(int i=0;i<D1.size();i++)
                  D.push_back(D1[i]);
              return D;
          }
        }
        D.push_back(M.m[slice(0,M.n,1)]);
      }
      return D;
    }
    
    //Inverse.......
    valarray<double> Invert() 
    {
        valarray<double> a= Souriau();
        valarray<double> t(0.0,n*n);
        Matrix In= Matrix(t), M=Matrix(m); 
        for(int i=0;i<n;i++)
            In=In+(M^(n-i-1))*a[i];
        In=In*(-1.0/a[n]);
        return In.m;
    }
       
    //Determinant......
    double det() 
    {
      if(n==0)
      {
          cout<<"Not Applicable for Determinant"<<endl;
          return m[0];
      }
      if(n==1)
          return m[0];
      if(n==2)
          return m[0]*m[3]-m[1]*m[2];
      if(n>2)  
      {
          valarray<double> a=Souriau();
          return pow(-1.0,n)*a[n];          
      } 
    }
 
};




