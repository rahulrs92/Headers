#pragma once

#include<data_tools.hh>

#include <TGraph.h>
#include <TCanvas.h>
#include <TGraph2D.h>
#include <TGraph2DErrors.h>
#include <TF1.h>
#include <TF2.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TFitResult.h>
#include <TMultiGraph.h>

using namespace ROOT;
using namespace std;

template <typename T = double> class plot
{
protected:
	TGraph *G;
	vector<T> &x, &y;
	vector<T> vx, vy;
	size_t nr;
	
	void make_plot(string iT,string ilx,string ily)
	{
		G = new TGraph(nr,x.data(),y.data());//();//
		G->SetMarkerStyle(3);
		G->SetDrawOption("APL");
		G->SetTitle((iT + ";" + ilx + ";" + ily).data());
		G->GetYaxis()->SetTitleOffset(0.9);
	}
public:
	plot(Data<T>& iD,size_t ix,size_t iy) : x(iD[ix]), y(iD[iy]), nr(iD.order()[0])
	{
		make_plot(iD.title(),iD.label(ix),iD.label(iy));
	}
	
	plot(Data<T>& iD,size_t iy) : plot(iD,0,iy) {}
	plot(Data<T>& iD) : plot(iD,0,0) {}
	plot(valarray<T> ivx,valarray<T> ivy,string iT="Plot",string ilx="X",string ily="Y") : x(vx), y(vy), nr(ivx.size())
	{
		for(int i=0; i<ivx.size(); i++)
		{
			if(!isnan(ivx[i]) && !isnan(ivy[i])) { vx.push_back(ivx[i]); vy.push_back(ivy[i]); }
			else { nr--; }
		}
		make_plot(iT,ilx,ily);
	}
	plot(vector<T>& ivx,vector<T>& ivy,string iT="Plot",string ilx="X",string ily="Y") : x(ivx), y(ivy), nr(ivx.size())
	{
		make_plot(iT,ilx,ily);
	}
	plot(valarray<T> ivy,string iT="Plot",string ilx="X",string ily="Y") : x(vx), y(vy), nr(ivy.size())
	{
		for(int i=0; i<ivy.size(); i++)
		{
			if(!isnan(ivy[i])) { vx.push_back(i); vy.push_back(ivy[i]); }
			else { nr--; }
		}
		make_plot(iT,ilx,ily);
	}
	plot(vector<T>& ivy,string iT="Plot",string ilx="X",string ily="Y") : x(vx), y(ivy), nr(ivy.size())
	{
		for(int i=0; i<ivy.size(); i++) { vx.push_back(i); }
		make_plot(iT,ilx,ily);
	}
	
	TGraph* operator-> () { return G; }
	TGraph* get() { return G; }
	
	static pair<size_t,T> min(vector<T>& v)
	{
		T ov = v[0]; size_t oi = 0;
		for(int i=1; i<v.size(); i++) { ov = v[i]<ov ? v[i] : ov; oi = v[i]<ov ? i : oi; }
		return make_pair(oi,ov);
	}
	
	static pair<size_t,T> max(vector<T>& v)
	{
		T ov = v[0]; size_t oi = 0;
		for(int i=1; i<v.size(); i++) { ov = v[i]>ov ? v[i] : ov; oi = v[i]>ov ? i : oi; }
		return make_pair(oi,ov);
	}
	
	void axis_equal()
	{
		T ra_x = G->GetXaxis()->GetXmax() - G->GetXaxis()->GetXmin(); 
		T ra_y = G->GetYaxis()->GetXmax() - G->GetYaxis()->GetXmin(); 
		TAxis* axw = (ra_x > ra_y ? G->GetXaxis() : G->GetYaxis());
		TAxis* axn = (ra_x > ra_y ? G->GetYaxis() : G->GetXaxis());
		T rw = axw->GetXmax() - axw->GetXmin();
		T rn = axn->GetXmax() - axn->GetXmin();
		T mpw = (axw->GetXmax() + axw->GetXmin()) / 2;
		T mpn = (axn->GetXmax() + axn->GetXmin()) / 2;
		axn->SetLimits(mpn - rw/2, mpn + rw/2);
	}
	
	TGraph* axis(bool a)
	{
		TGraph* Ga = new TGraph();
		T mi = min(a ? y : x).second, ma = max(a ? y : x).second;
		Ga->AddPoint(a?0:mi, a?mi:0);
		Ga->AddPoint(0,0);
		Ga->AddPoint(a?0:ma, a?ma:0);
		Ga->SetDrawOption("APL");
		return Ga;
	}
};

template <typename T = double> class plot2 : public TObject
{
protected:
	TGraph2D *G;
	vector<T> &x, &y, &z;
	vector<T> vx, vy, vz;
	size_t nr;
	
	void make_plot(string iT,string ilx,string ily,string ilz)
	{
		G = new TGraph2D("plot","",nr,x.data(),y.data(),z.data());//();//
		//G->SetDrawOption("surf1");
		G->SetTitle((iT + ";" + ilx + ";" + ily + ";" + ilz).data());
		//g_no++;
	}
public:
	//static size_t g_no;
	
	plot2(Data<T>& iD,size_t ix,size_t iy,size_t iz) : x(iD[ix]), y(iD[iy]), z(iD[iz]), nr(iD.order()[0])
	{
		make_plot(iD.title(),iD.label(ix),iD.label(iy),iD.label(iz));
	}
	
	plot2(Data<T>& iD,size_t iz) : plot2(iD,0,0,iz) {}
	plot2(Data<T>& iD) : plot2(iD,0,0,0) {}
	plot2(valarray<T> ivx,valarray<T> ivy,valarray<T> ivz,
				string iT="Plot",string ilx="X",string ily="Y",string ilz="Z") : 
				x(vx), y(vy), z(vz), nr(ivz.size())
	{
		for(int i=0; i<ivx.size(); i++) { vx.push_back(ivx[i]); vy.push_back(ivy[i]); vz.push_back(ivz[i]); }
		make_plot(iT,ilx,ily,ilz);
	}
	
	plot2(vector<T>& ivx,vector<T>& ivy,vector<T>& ivz,
				string iT="Plot",string ilx="X",string ily="Y",string ilz="Z") : 
				x(ivx), y(ivy), z(ivz), nr(ivz.size()) { make_plot(iT,ilx,ily,ilz); }
	
	plot2(valarray<T> ivz,string iT="Plot",string ilx="X",string ily="Y",string ilz="Z") : x(vx), y(vy), z(vz), nr(ivz.size())
	{
		int n = pow(ivz.size(),0.5);
		for(int i=0; i<ivz.size(); i++) { vx.push_back(i/n); vy.push_back(i%n); vz.push_back(ivz[i]); }
		make_plot(iT,ilx,ily,ilz);
	}
	
	plot2(vector<T>& ivz,string iT="Plot",string ilx="X",string ily="Y",string ilz="Z") : x(vx), y(vy), z(ivz), nr(ivz.size())
	{
		int n = pow(ivz.size(),0.5);
		for(int i=0; i<ivz.size(); i++) { vx.push_back(i/n); vy.push_back(i%n); }
		make_plot(iT,ilx,ily,ilz);
	}
	
	TGraph2D* operator-> () { return G; }
	TGraph2D* get() { return G; }
	
	void axis_equal()
	{
		T ra_x = G->GetXaxis()->GetXmax() - G->GetXaxis()->GetXmin(); 
		T ra_y = G->GetYaxis()->GetXmax() - G->GetYaxis()->GetXmin(); 
		T ra_z = G->GetZaxis()->GetXmax() - G->GetZaxis()->GetXmin(); 
		TAxis* axw = ra_z > ra_x && ra_z > ra_y ? G->GetZaxis() : (ra_x > ra_y ? G->GetXaxis() : G->GetYaxis());
		TAxis* axi = ra_x < ra_z == ra_z < ra_y ? G->GetZaxis() : (ra_x > ra_y ? G->GetYaxis() : G->GetXaxis());
		TAxis* axn = ra_z < ra_x && ra_z < ra_y ? G->GetZaxis() : (ra_x > ra_y ? G->GetYaxis() : G->GetXaxis());
		T rw = axw->GetXmax() - axw->GetXmin();
		T mpw = (axw->GetXmax() + axw->GetXmin()) / 2;
		T mpi = (axi->GetXmax() + axi->GetXmin()) / 2;
		T mpn = (axn->GetXmax() + axn->GetXmin()) / 2;
		axi->SetLimits(mpi - rw/2, mpi + rw/2);
		axn->SetLimits(mpn - rw/2, mpn + rw/2);
	}
	
	TGraph2D* axis(size_t a)
	{
		TGraph2D* Ga = new TGraph2D();
		T mi = Data<T>::min(a==0 ? x : (a==1 ? y : z)).second;
		T ma = Data<T>::max(a==0 ? x : (a==1 ? y : z)).second;
		Ga->AddPoint(a==0?mi:0, a==1?mi:0, a>1?mi:0);
		Ga->AddPoint(0,0,0);
		Ga->AddPoint(a==0?ma:0, a==1?ma:0, a>1?ma:0);
		Ga->SetDrawOption("APL");
		return Ga;
	}
};


