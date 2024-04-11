
// Import libraries and header files
#include <string.h>
#include <sys/resource.h>
#include <map>
#include <regex>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <sstream>

// Bring the FBGLPK namespace into the current scope
using namespace FBGLPK;

// Declare static variables
static double* Vars; // Pointer to an array of doubles
static double FBAtime = -1;
static double Flag = -1;

double rate = 0;

static map <string, string> FBAreact; // Map that stores key-value pairs of strings
static map <string, string> FBAplace; // Map that stores key-value pairs of strings

static vector<vector<double>*> M; // 2D vector of pointers to doubles

static bool init = false; // Boolean variable initialized to false

static int timeindex;
static int Msize;

// method 1:
void read_map_string_int(string fname, unordered_map<string,int>& m) {
  ifstream f (fname);
  string line;
  if(f.is_open()) {
    cout << "#### " << fname << "####" << endl;
    int j = 1;
    while (getline(f,line))
    {
      line.erase(remove( line.begin(), line.end(), '\"' ),line.end());
      m.insert(pair<string,int>(line,j));
      //cout << line << ";" << j << endl;
      ++j;
    }
    f.close();
  }
  else
  {
    std::cerr<<"\nUnable to open " << fname << ": file do not exists\n";
    exit(EXIT_FAILURE);
  }
}

// method 2:
void read_map_string_double(string fname, unordered_map<string,double>& m) {
  ifstream f (fname);
  string line;
  if(f.is_open()) {
    size_t pos = 0, length = 0;
    cout << "#### " << fname << "####" << endl;
    int j = 1;
    while (getline(f,line)) {
      line.erase(remove( line.begin(), line.end(), '\"' ),line.end());
      pos = 0;
      // read rates
      length = line.length();
      
      pos = line.find(',');
      if( pos == string::npos)
        pos = length;
      m.insert(pair<string,double>(line.substr(0,pos) , stod(line.substr(pos+1,length))) );
      cout <<line.substr(0,pos) << ": " << stod(line.substr(pos+1,length)) << " ";
      cout << endl;
      ++j;
    }
    f.close();
  }
  else {
    std::cerr<<"\nUnable to open " << fname << ": file do not exists\n";
    exit(EXIT_FAILURE);
  }
}

// method 3:
void read_constant(string fname, double& rate) {
  ifstream f (fname);
  string line;
  if(f.is_open()) {
    int i = 1;
    while (getline(f,line)) {
      switch(i) {
      case 1:
        rate = stod(line);
        //cout << "p" << i << ": " << line << "\t" << p1 << endl;
        break;
      }
      ++i;
    }
    f.close();
  }
  else {
    std::cerr<<"\nUnable to open " << fname << ": file do not exists\": file do not exists\n";
    exit(EXIT_FAILURE);
  }
}

// method 4:
void read_vector_vector(string fname, vector<vector<double>*>& v) {
  ifstream f (fname);
  string line;
  if(f.is_open())
  {
    //cout << "#### " << fname << " ####" << endl;
    while (getline(f,line))
    {
      size_t pos = 0, c_pos = 0, length = 0;
      int i = 1;
      vector<double>* vec = new vector<double>;
      // read rates
      length = line.length();
      do
      {
        pos = line.find(',',c_pos);
        if( pos == string::npos)
          pos = length;
        //cout << "***" << line.substr(c_pos,pos - c_pos) << endl;
        vec -> push_back(stod(line.substr(c_pos,pos - c_pos)));
        cout << to_string(i) << ": " << to_string(vec -> at(vec -> size() -1)) << " | ";
        c_pos = pos+1;
        ++i;
      }
      while(pos != length);
      cout << endl;
      v.push_back(vec);
    }
  }
}

// method 5:
void read_vector(string fname, vector<double> V) {
  string line;
  ifstream inFile(fname);
  //cout << "#### " << fname << " ####" << endl;
  while (getline(inFile, line)) {
    V.push_back(stod(line));
    //cout << "***" << stod(line) << endl;
  }
}

// method 6:
// the function "search" takes two parameters: 
// a double value named time and a 2D vector of pointers to doubles named matrix. 

// The purpose "search" is to search for the right time epoch in the matrix vector
//that corresponds to the given time value.
void search(double time, vector<vector<double>*> matrix) {
  // Initialize variables
  int idx = 0; // Current row index
  bool fnd_epoch = false; // Flag indicating if the right time epoch has been found
  unsigned int iter = 0; // Current iteration count
  
  // Loop over rows of the matrix until the right time epoch is found or the end of the matrix is reached
  while(iter < (matrix.size()-1) && !fnd_epoch) {
    // Check if the time value is within the time epoch defined by the current row and the next row
    if(matrix[iter+1] -> at(0) > time && time >= matrix[iter] -> at(0)) {
      // If it is, set the timeindex variable to the current row index and set the fnd_epoch flag to true
      timeindex = idx;
      fnd_epoch = true;
    }
    // If the end of the matrix is reached and the time value is greater than the last time epoch
    // then set the timeindex variable to the last row index and set the fnd_epoch flag to true
    else if (iter == (matrix.size()-2)) {
      cout << "time is " << time << endl;
      cout << "M["<< iter+1 <<"] =" << matrix[iter+1] -> at(0) << endl;
      if(matrix[iter+1] -> at(0) < time) {
        timeindex = idx;
        fnd_epoch = true;
      }
    }
    // If the time value is not within the current time epoch
    // then increment the row index and iteration count and continue to the next row
    else {
      ++idx;
      ++iter;
    }
  }
  
  // If the right time epoch is not found, throw an exception
  if(!fnd_epoch) {
    throw std::invalid_argument( "rE: unable to find a contact matrix for time > " + 
                                 to_string(time) + "! Please, kill me and check the input parameters" );
  }
}

// method 7:
double trunc(double value, double decimal) {
  const double multiplier = std::pow(10.0, decimal);
  return std::floor(value * multiplier) / multiplier;
  //return ((unsigned long int)(value * multiplier)) / multiplier;
}

// data strucure inizialization
void init_data_structures(const struct InfTr* Trans, map <string,int>& NumTrans)  {
  
  read_vector_vector("./M", M);
  Msize = M.size();
  
  FBAreact["EX_lcts_e_in"] = "EX_lcts_e";
  FBAreact["EX_glc_e_in"] = "EX_glc_D_e";
  
  FBAplace["EX_lcts_e_in"] = "lcts_e";
  FBAplace["EX_glc_e_in"] = "glc_e";
  
  init = true;
  Flag = 1;
  
}

double FBA(double *Value,
           vector<class FBGLPK::LPprob>& vec_fluxb,
           map <string,int>& NumTrans,
           map <string,int>& NumPlaces,
           const vector<string> & NameTrans,
           const struct InfTr* Trans,
           const int T,
           const double& time) {
  
  bool FBAmarking = 0;
  // truncation floating
  double decimalTrunc = 12;
  
  if(Flag == -1) init_data_structures(Trans, NumTrans);
  
  cout << "FBAtime: " << FBAtime << endl;
  cout << "glc_e(t = " << time << ") = " << Value[NumPlaces["glc_e"]] << " (mmol)" << endl;
  cout << "lcts_e(t = " << time << ") = " << Value[NumPlaces["lcts_e"]] << " (mmol)" << endl;
  
  if(FBAtime != time) {
    
    for (map<string, string>::iterator p = FBAreact.begin(); p != FBAreact.end(); ++p) {
      
      int index = vec_fluxb[0].fromNametoid(p->second);
      string TypeBound = "GLP_DB";
      
      double Ub = vec_fluxb[0].getUpBounds(index);
      double Lb = vec_fluxb[0].getLwBounds(index);
      
      // Check if the first element of the map is equal to reaction for Glucose uptake
      double Met = trunc(Value[NumPlaces[FBAplace[p->first]]], decimalTrunc);
      Lb = -Met;
      if(abs(Lb) > Met){
        if(-Met >= vec_fluxb[0].getLwBounds(index)){
          Lb = -Met;
        } else {
          Lb = vec_fluxb[0].getLwBounds(index);
        }
      }
      
      cout << "Trans: " <<  p->first << ", Met (mmol): " << Met << ";" << endl;
      cout << "Trans: " <<  p->first << ", Lb (mmol/gDW*h): " << Lb << ";" << endl;
      cout << "Trans: " <<  p->first << ", Ub (mmol/gDW*h): " << Ub << ";" << endl;
      
      double Lbtr = trunc(Lb, decimalTrunc);
      double Ubtr = trunc(Ub, decimalTrunc);
      
      vec_fluxb[0].update_bound(index, TypeBound, Lbtr, Ubtr);
      
    }
    
    vec_fluxb[0].solve();
    Vars = vec_fluxb[0].getVariables();
    
    FBAtime = time;
    
  }
  
  bool Out = 0;
  bool In = 0;
  
  string str = NameTrans[T];
  
  if(str.find("_out") != string::npos){
    str = std::regex_replace(str, std::regex("_out"), "_in");
    Out = 1;
    } else {
      In = 1;
      }
    
    int index = vec_fluxb[0].fromNametoid(FBAreact.find(str) -> second);
    
    cout <<"\nSolution:\n\n";
    
    double glc = trunc(Value[NumPlaces.find("glc_e") -> second], decimalTrunc);
    double lcts = trunc(Value[NumPlaces.find("lcts_e") -> second], decimalTrunc);
    
    rate = Vars[index];
    
    if((Out) && (rate > 0))
      rate = rate;
    else if((Out) && (rate < 0))
      rate = 0;
    else if((In) && (rate > 0))
      rate = 0;
    else if((In) && (rate < 0))
      rate = -rate;
    
    cout << "Firing transition: " << NameTrans[T] << endl;
    cout << "transition associated to: " << FBAreact.find(str) -> second << endl;
    cout << "flux estimated from fba = " << Vars[index] << " (mmol/gDW*h)" << endl;
    
    cout << "index: " <<  index << endl;
    
    return(rate);
  
}

// CarbonFlux transition
double CarbonFlux(double *Value,
                  vector<class FBGLPK::LPprob>& vec_fluxb,
                  map <string,int>& NumTrans,
                  map <string,int>& NumPlaces,
                  const vector<string> & NameTrans,
                  const struct InfTr* Trans,
                  const int T,
                  const double& time) {
  
  double rate = 0;
  
  if(init == false) {
    void init_data_structures(const struct InfTr* Trans, map <string,int>& NumTrans);
    }
  
  if(timeindex > 0 && timeindex < Msize) {
    // timeindex == 0: the Cpp code is in the first iteration
    // NB: the second IF loop scans upward and downward indexes
    
    // cout << "IF (search) loop condition: " << time << "<=" 
    // << M[timeindex-1]-> at(0) << " || " << time  
    // << ">" <<  M[timeindex]-> at(0) << endl;
    
    if(time < M[timeindex]-> at(0) || time >= M[timeindex+1]-> at(0)) {
      // if simulation time (time) gets through the time interval already crossed
      // during previous iterations, then the IF loop called search() and rates updated
      search(time, M);
    }
  }
  
  else if(timeindex == 0){
    if(time >= M[timeindex + 1]-> at(0)) {
      search(time, M);
      }
    }
  else{
    if(time < M[timeindex]-> at(0)) {
      search(time, M);
      }
    }
  
  if(NameTrans[T] == "Tglc") {
    rate = M[timeindex]-> at(1);
    }
  else if(NameTrans[T] == "Tlcts"){
    rate = M[timeindex]-> at(2);
    }
  
  double intensity = 1.0;
  for (unsigned int k=0; k<Trans[T].InPlaces.size(); k++) {
    intensity *= pow(Value[Trans[T].InPlaces[k].Id],Trans[T].InPlaces[k].Card);
    }
  
  rate = rate*intensity;
  
  cout << "rate_" << NameTrans[T] << " = " << rate << endl;
  
  return(rate);
  
}
