#pragma once

#include<iostream>
#include<fstream>
#include<array>
#include<valarray>
#include<sstream>		// To count words in a line
#include<vector>
#include<map>
#include<limits>		// To define NaN
#include<cmath>			// To check for NaNs using isnan()
#include<iomanip>		// To set precision of the output streamer

#include<pick.hh>

using namespace std;

template <typename T> class selector
{
protected:
	int col;
	T val;
public:
	selector() {}
	selector(int c, T v) : col(c), val(v) {}
	int& column() { return col; }
	T& value() { return val; }
};

/* The class and static methods inside are designed to read a datafile and organize them.
 * The data in the text file must be pre-formatted as a table with fixed number of rows and
 * columns, with one row per line and the values delimited with one or more whitespace
 * characters.
 */

template <typename T = double> class Data
{
protected:
	vector<vector<T>> C;		//A vector of columns
													//C[i] will return column number i
	string Zs;							//An empty string (for referencing during exceptions)
	size_t nr,nc;						//Number of rows and columns
	string TT;							//Title of the dataset (default: Data)
	vector<string> L;				//Labels of each of the columns (default: Column #)

	vector<T> Ri;						//The row-index column
	vector<T> Ci;						//The column-index row
	row<T> aCi;							//Pointers to the contents of vector Ci
	//row<T> aRi;							//Pointers to the contents of vector Ri
	//row<T> arc;							//Pointers to a referred column
	row<T> arr;							//Pointers to a referred row

	vector<T> Zc;						//A zero column (for referencing during exceptions)
	vector<T> Zr;						//A zero column (for referencing during exceptions)
	T* aZc;									//Pointers to the contents of vector Zc
	T* aZr;									//Pointers to the contents of vector Zr

	map<string,size_t> l2c;

	array<bool,2> rng(size_t ii) { return {ii>0 && ii<=nr, ii>0 && ii<=nc}; }

	array<bool,2> rng(array<size_t,2> ia) { return {ia[0]>0 && ia[0]<=nr && ia[1]>0 && ia[1]<=nr && ia[0]<ia[1],
																									ia[0]>0 && ia[0]<=nc && ia[1]>0 && ia[1]<=nc && ia[0]<ia[1]}; }

	void c_vld(array<size_t,2>& ia)
	{
		ia = {(ia[0] > ia[1] ? ia[1] : ia[0]), (ia[0] > ia[1] ? ia[0] : ia[1])};
		ia = {(ia[0] < 1 ? 1 : ia[0]), (ia[1] > nc ? nc : ia[1])};
	}

	void r_vld(array<size_t,2>& ia)
	{
		ia = {(ia[0] > ia[1] ? ia[1] : ia[0]), (ia[0] > ia[1] ? ia[0] : ia[1])};
		ia = {(ia[0] < 1 ? 1 : ia[0]), (ia[1] > nr ? nr : ia[1])};
	}

	void update()
	{
		Ri = vector<T>(nr);
		Ci = vector<T>(nc);
		Zc = vector<T>(nr);
		Zr = vector<T>(nc);
		Zs = "";

		if constexpr (is_same<T,string>::value)
		{
			for(int i=0; i<nr; i++) { Ri[i] = to_string(i+1); Zc[i] = to_string(0); }
			for(int i=0; i<nc; i++) { Ci[i] = to_string(i+1); Zr[i] = to_string(0); }
		}
		else
		{
			for(int i=0; i<nr; i++) { Ri[i] = T(i+1); Zc[i] = T(0); }
			for(int i=0; i<nc; i++) { Ci[i] = T(i+1); Zr[i] = T(0); }
		}

		aCi.clear(); // IMPORTANT: Clear previous references first!
		aCi.set_row(0);
		for(int i=0; i<nc; i++) aCi.push_back(Ci[i]);
		//for(int i=0; i<=nr; i++) aRi.push_back(Ri[i]);
	}

public:
	static double constexpr NaN = std::numeric_limits<double>::quiet_NaN();

	Data()
	{
		nr = 0;
		nc = 0;
	}

	Data(const char* name)
	{
		read(name);
		// cout<<"YOOOO!"<<endl;
		TT = "Data";
		L = vector<string>(nc + 1,"Column ");
		L[0] = "Row no.";

		fstream f;
		string info = string(name) + "_info";
		f.open(info.data(), fstream::in);
		if(f.is_open())
		{
			cerr<<"info file found at: \n"<<info<<endl;
			string tts;
			for(int i=1; (i<=nc && !f.eof()); i++) { f>>L[i]; }
			f.close();
		}
		else for(int i=1; i<=nc; i++) L[i] = L[i] + to_string(i);
		update();
	}

	Data(string name) : Data(name.data()) {}

	Data(vector<vector<T>> iC)
	{
		C = iC;
		nc = C.size();
		nr = C[0].size();
		TT = "Data";
		L = vector<string>(nc + 1,"Column ");
		for(int i=0; i<=nc; i++) L[i] = L[i] + to_string(i);
		for(int i=0; i<=nc; i++) l2c[L[i]] = i;
		L[0] = "Row no.";
		update();
	}

	Data& operator = (const Data& iD)
	{
		if (this != &iD)
		{
			C = iD.C;
			nr = iD.nr;
			nc = iD.nc;
			TT = iD.TT;
			L = iD.L;
			update(); // Rebuilds Ci, Ri, Zc, Zr, AND resets aCi references to local memory
		}
	    return *this; // Fixes the trace trap
	}

	friend Data<T> operator + (const Data<T>& iD);
	
	vector<vector<T>> dataset() const { return C; }
	
	string& title() { return TT; }
	
	vector<string>& labels() { return L; }
	vector<string> labels() const { return L; }
	
	string& label(size_t ic) { return (ic>nc ? Zs : L[ic]); }
	
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
	
	/* This static method allows you to count the number of rows and columns
	 * in a datafile without creating any objects. To do this, you can simply
	 * call 'auto x = Data<type>::order("file.txt")'.
	 */
	static array<size_t,2> order(const char* name)
	{
		fstream f;
		try { f.open(name, fstream::in); }
		catch(exception X) { cerr<<"Error opening file: \n"<<X.what()<<endl; return {0,0}; }
		size_t N=0, M=0;
		string s;
		while(!f.eof())
		{
			getline(f,s);
			N++;
		}
		N = N - 1;
		f.close();
		f.open(name, fstream::in);
		getline(f,s);
		f.close();
		stringstream ss;
		ss.str(s);
		while(!ss.eof())
		{
			ss>>s;
			M++;
			// cout << M << " : " << s << endl;
		}
		return {N,M};		//{Number of rows, Number of Columns}
	}
	
	/* This static method allows you to read a datafile and store it in a two-dimensional
	 * vector, without creating any objects. Simply call 'auto x = Data<type>::read("file.txt",{nC,nR})'.
	 * In this configuration, x[i] will be a vector containing column number i, and x[i][j] refers
	 * to element j of that column.
	 */
	template <typename U = double> static vector<vector<U>> read(const char* name, array<size_t,2> order)
	{
		fstream f;
		//cout<<"Table order at template read: "<<order[0]<<","<<order[1]<<endl;
		vector<vector<U>> oC(order[1],vector<U>(order[0]));
		try { f.open(name, fstream::in); }
		catch(exception X) { cerr<<"Error opening file: \n"<<X.what()<<endl; return oC; }
		for(int i=0; i<order[0]; i++)
		{
			for(int j=0; j<order[1]; j++)
			{
				/* Here I use some stack exchange "Jugaad" to specify what to do
				 * in case of encountering different datatypes in U. Currently,
				 * I have dealt with int, float, double and string.
				 * When non-numeric values are encountered, float and double would
				 * store NaN, whereas int would store 0.
				 */
				string ts;
				if constexpr (is_same<U,double>::value || is_same<U,float>::value || is_same<U,int>::value)
				{
					double td;
					f>>td;
					if(f.fail())
					{
						f.clear();
						f>>ts;
						cerr<<"Found NaN at "<<i<<","<<j<<": "<<ts<<endl;
						td = NaN;
					}
					if constexpr (is_same<U,int>::value)
					{
						if(td != double(int(td)) && !isnan(td))
						cerr<<"Precise value floored at "<<i<<","<<j<<" from "<<td<<" to "<<int(td)<<"!"<<endl;
						oC[j][i] = (isnan(td) ? 0 : int(td));
					}
					else
					{
						oC[j][i] = td;
					}
				}
				else
				{
					f>>ts;
					oC[j][i] = ts;
				}
			}
		}
		f.close();
		return oC;
	}

	void read(const char* name)
	{
		auto o = order(name);
		nr = o[0]; nc = o[1];
		C = read<T>(name,o);

		/*
		fstream f;
		string info = string(name) + "_info";
		try { f.open(info.data(), fstream::in); }
		catch(exception X) { return; }

		cerr<<"info file found at: \n"<<info<<endl;
		string tts;
		for(int i=0; i<nc; i++) { f>>tts; cout<<tts<<endl; }
		f.close();
		*/
	}
	
	array<size_t,2> order() const { return {nr,nc}; }
	
	/* This function returns column number ic in U: a container of your choice.
	 * It is shown to work for vectors and valarrays, and will work if the size
	 * constructor in the first line (commented below), and referencing using []
	 * exists for the class.
	 */
	template <typename U = valarray<double>> U col(size_t ic)
	{
		U ou(nr);																		//most commonly sets the size of the container
		for(int i=0; i<nr; i++) ou[i] = (ic==0 ? i : (ic>nc ? 0 : C[ic-1][i]));
		return ou;
	}
	
	template <typename U = valarray<double>> U col(string il)
	{
		l2c.clear();
		for(int i=0; i<=nc; i++) l2c[L[i]] = i;
		return col<U>(l2c[il]);
	}
	
	column<T> operator [] (size_t ic) { return (ic==0 ? column<T>(Ri,ic) : column<T>(C[ic-1],ic)); }

	Data<T> operator [] (selector<T> s)
	{
		Data<T> oD;
		for (int i=0; i<nc; i++)
		{
			
		}
		return oD;
	}
	
	row<T> operator () (size_t ir)
	{
		if(ir==0) return aCi;
		arr.clear();
		arr.set_row(Ri[ir-1]);
		for(int i=0; i<nc; i++) arr.push_back(C[i][ir-1]);
		return (arr);
	}
	T& operator () (size_t ic, size_t ir)
	{
		return (ic==0 ? Ri[ir-1] : (ic>nc ? Zr[ir-1] : C[ic-1][ir-1]));
	}
	valarray<T> operator () (size_t ic, array<size_t,2> iar)
	{
		r_vld(iar);
		valarray<T> ov(T(0),iar[1]-iar[0]+1);
		for(int i=iar[0]; i<=iar[1]; i++) ov[i-iar[0]] = (ic==0 ? Ri : C[ic-1])[i];
		return ov;
	}
	
	valarray<T> operator () (array<size_t,2> iac, size_t ir)
	{
		c_vld(iac);
		valarray<T> ov(T(0),iac[1]-iac[0]+1);
		for(int i=iac[0]; i<=iac[1]; i++) ov[i-iac[0]] = (ir==0 ? Ci[i] : C[i][ir-1]);
		return ov;
	}
	
	Data<T> operator () (array<size_t,2> iac, array<size_t,2> iar)
	{
		c_vld(iac);r_vld(iar);
		cout<<"inputs: "<<iac[0]<<", "<<iac[1]<<", "<<iar[0]<<", "<<iar[1]<<endl;
		cout<<"vector size: "<<(iac[1]+1)-iac[0]<<", "<<(iar[1]+1)-iar[0]<<endl;
		vector<vector<T>> oC(iac[1]-iac[0]+1,vector<T>(iar[1]-iar[0]+1,0));
		for(int i=iac[0]; i<=iac[1]; i++) 
			for(int j=iar[0]; j<=iar[1]; j++)
			{
				cout<<"loop parameters: "<<i-iac[0]<<", "<<j-iar[0]<<", "<<i<<", "<<j<<", "<<(*this)[i][j]<<endl;
				oC[i-iac[0]][j-iar[0]] = (*this)[i][j];
			}
		/*
		cout<<"Size of oC: "<<oC[1].size()<<endl;
		*/
		return Data<T>(oC);
	}
	
	/* This function allows the user to add two datasets, essentially
	 * placing the columns of 'iD' next to the columns of this object,
	 * thus making a larger combined dataset.
	 * It also handles the exception if the operands have different number
	 * of rows. The returned set acquires the lower among the number of
	 * rows of the operands, creating the "largest possible rectangle"
	 * from the combined dataset. However, this exception should be
	 * avoided as the operation would result in the loss of some data.
	 */ 
	Data<T> operator + (const Data<T>& iD)
	{
		if(iD.order()[0] == 0 || iD.order()[1] == 0) return *this;
		
		Data<T> oD = *this;
		size_t inr = iD.order()[0];
		size_t inc = iD.order()[1];
		
		size_t onc = nc + inc;
		size_t onr = (inr>=nr ? nr : inr);
		
		auto iC = iD.dataset(), oC = C;
		auto iL = iD.labels(), oL = L;
		
		/* Here I "trim" the columns of the dataset with longer columns
		*/
		for(int i=nr; i<inr; i++) for(int j=0; j<inc; j++) iC[j].pop_back();
		for(int i=inr; i<nr; i++) for(int j=0; j<nc; j++) oC[j].pop_back();
	
		/* Here I append the columns of the output dataset with the ones
		* from the input dataset.
		*/
		for(int i=0; i<onc-nc; i++)
		{
			oC.push_back(iC[i]);
			oL.push_back(iL[i+1]);
		}
		
		/* Here I use my friend privileges to set the properties
		* of the output Data
		*/
		oD.C = oC; oD.L = oL; oD.nr = onr; oD.nc = onc;
		oD.TT = TT + " + " + iD.TT;
		oD.update();

		// cout<<"Adding datasets: "<<TT<<" + "<<iD.TT<<endl;
		// cout<<"Resulting dataset: "<<oD.order()[0]<<" x "<<oD.order()[1]<<endl;
		// cout<<"Resulting dataset labels: "<<endl;
		// for(int i=1; i<=oD.order()[1]; i++) cout<<oD.labels()[i]<<" ";
		// cout<<endl;
		
		return oD;
	}
	
	/* This function allows the user to merge two datasets, essentially
	 * placing the rows of 'iD' below the rows of this object, thus making
	 * a larger combined dataset.
	 * It also handles the exception if the operands have different number
	 * of columns. The returned set acquires the lower among the number of
	 * columns of the operands, creating the "largest possible rectangle"
	 * from the combined dataset. However, this exception should be avoided
	 * as the operation would result in the loss of some data.
	 */ 
	Data<T> operator / (const Data& iD)
	{
		if(iD.order()[0] == 0 || iD.order()[1] == 0) return *this;
		
		Data<T> oD = *this;
		size_t inr = iD.order()[0];
		size_t inc = iD.order()[1];
		
		size_t onr = nr + inr;
		size_t onc = (inc>=nc ? nc : inc);
		
		auto iC = iD.dataset(), oC = C;
		auto iL = iD.labels(), oL = L;
		
		/* Here I delete the columns of the dataset with more columns
		 */
		for(int i=nc; i<inc; i++) iC.pop_back();
		for(int i=inc; i<nc; i++) oC.pop_back();
		
		/* Here I append the columns of the output dataset with the ones
		 * from the input dataset.
		 */
		for(int i=0; i<onr-nr; i++) for(int j=0; j<onc; j++) oC[j].push_back(iC[j][i]);
		
		/* Here I use my friend privileges to set the properties
		 * of the output Data
		 */
		oD.C = oC; oD.L = oL; oD.nr = onr; oD.nc = onc;
		oD.TT = TT + " / " + iD.TT;
		oD.update();
		
		return oD;
	}

	selector<T> operator == (pair<size_t,T> iD)
	{
		return selector<T>(iD.first, iD.second);
	}

	selector<T> operator == (pair<string,T> iD)
	{
		return selector<T>(labels()[iD.first], iD.second);
	}

	void clean_nans()
	{
		bool f_nan = 0;
		size_t i_nan = 0;
		for(int i=0; i<nr; i++)
		{
			for(int j=0; j<nc; j++)
			{
				//cout<<"i: "<<i<<", j: "<<j<<", C: "<<C[j][i];
				if(isnan(C[j][i])) { f_nan = 1; i_nan = i; break; }
			}
			//cout<<endl;
			if(f_nan) break;
		}
		if(f_nan)
		{
			//cout<<"i_nan: "<<i_nan<<endl;
			for(int j=0; j<nc; j++) C[j].erase(C[j].begin() + i_nan);
			nr--;
			update();
			clean_nans();
		}
	}
	
	
	vector<T> convert_to_mesh(size_t ix, size_t iy, size_t iz)
	{
		// Extract the x, y, and z columns
		vector<T>& x_col = C[ix - 1];	// Adjusting for 0-based indexing
		vector<T>& y_col = C[iy - 1];
		vector<T>& z_col = C[iz - 1];
		
		// Get unique x and y values and sort them
		vector<T> unique_x(x_col.begin(), x_col.end());
		vector<T> unique_y(y_col.begin(), y_col.end());
		
		sort(unique_x.begin(), unique_x.end());
		sort(unique_y.begin(), unique_y.end());
		
		unique_x.erase(unique(unique_x.begin(), unique_x.end()), unique_x.end());
		unique_y.erase(unique(unique_y.begin(), unique_y.end()), unique_y.end());
		
		size_t nx = unique_x.size();
		size_t ny = unique_y.size();
		
		// Helper lambda to convert double to string without trailing zeros
		auto double_to_string = [](T val) -> string
		{
			stringstream ss;
			ss << fixed << setprecision(15) << val;
			string str = ss.str();
			// Remove trailing zeros
			str.erase(str.find_last_not_of('0') + 1, string::npos);
			// Remove trailing decimal point if no fractional part
			if(str.back() == '.') str.pop_back();
			return str;
		};
		
		// Create a map for quick lookup: (x, y) -> z_index
		map<pair<T, T>, size_t> coord_to_index;
		for(size_t i = 0; i < x_col.size(); i++)
		{
			coord_to_index[{x_col[i], y_col[i]}] = i;
		}
		
		// Create the new matrix: ny rows (for each y), nx columns (for each x)
		vector<vector<T>> mesh(nx, vector<T>(ny, T(0)));
		
		for(size_t j = 0; j < ny; j++)
		{
			for(size_t i = 0; i < nx; i++)
			{
				auto it = coord_to_index.find({unique_x[i], unique_y[j]});
				if(it != coord_to_index.end())
				{
					mesh[i][j] = z_col[it->second];
				}
			}
		}
		
		// Create new labels: column labels are x coordinates
		vector<string> new_labels(nx + 1);
		new_labels[0] = "Y \\ X";
		for(size_t i = 0; i < nx; i++)
		{
			new_labels[i + 1] = double_to_string(unique_x[i]);
		}
		
		// Update the Data object
		C = mesh;
		L = new_labels;
		nc = nx;
		nr = ny;
		TT = "Mesh";
		update();

		// Return the unique y values for reference
		return unique_y;		
	}
	
};


template <typename T = double> ostream& operator<< (ostream& os, Data<T> iD)
{
	//os.precision(4);
	ios init(NULL);
	init.copyfmt(os);
	//os << iD.labels()[0];
	//for(int j=1; j<=iD.order()[1]; j++) os << setprecision(13) << setw(15) << iD.labels()[j];
	//os << endl;
	os << setprecision(5) << setw(15) << "Row\\Col";
	for(int i=0; i<=iD.order()[0]; i++)
	{
		if(i!=0) os << setprecision(5) << setw(15) << iD[0][i];//(i)[0];//
		os << '\t';
		for(int j=1; j<=iD.order()[1]; j++) os << setprecision(5) << setw(15) << iD[j][i]; //(i)[j];//<< '\t'  //
		os << endl;
	}
	os.copyfmt(init);
	//cout<<"YES"<<endl;
	return os;
}

