//#include <stdafx.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <ctime>
#include <algorithm>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <iomanip>
#include <cstdlib>
#include <dirent.h>


/*
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rstat.h>
*/

#include </usr/local/include/gsl/gsl_linalg.h>
#include </usr/local/include/gsl/gsl_vector.h>
#include </usr/local/include/gsl/gsl_blas.h>
#include </usr/local/include/gsl/gsl_randist.h>
#include </usr/local/include/gsl/gsl_rng.h>
#include </usr/local/include/gsl/gsl_statistics.h>
#include </usr/local/include/gsl/gsl_sf_exp.h>
#include </usr/local/include/gsl/gsl_sf_pow_int.h>
#include </usr/local/include/gsl/gsl_sf_log.h>
#include </usr/local/include/gsl/gsl_statistics_double.h>
#include </usr/local/include/gsl/gsl_multimin.h>
//#include </usr/local/include/gsl/gsl_rstat.h>
//#include </usr/local/Cellar/gsl/1.16/include/gsl/gsl_multimin.h>




using namespace std;

string MAChome="/Users/johann/Documents/SYMBA/";
string MACwork="/Users/admin/Documents/GNT/SYMBA/";
string MACwork2="/Users/johanngnt/Documents/GNT/SYMBA/";
string PCwork="C:\\Users\\lussange\\Desktop\\";
string Machine=MAChome;
string Postdoc2022="/Users/johann/Desktop/Work/Postdoc2022/";


int DayTicks = 1;
int Week = 5*DayTicks;
int Month = 21*DayTicks;
int Year = 252*DayTicks;
int Year2 = 286*DayTicks;



// Initialize a brand new RNG state, with unique seed
gsl_rng * make_rng()
{
    gsl_rng_default_seed = static_cast<unsigned long>(time(NULL)) + rand();
    const gsl_rng_type * T = gsl_rng_ranlux389;
    gsl_rng * r = gsl_rng_alloc (T);
    return r;
}







// This converts an int to a string
string IntToString (int A) {
    std::ostringstream ostr;
    ostr << A;
    std::string NewA = ostr.str();
    return NewA;
};



// Trunks any double number with a number "Digits" of Significant digits
double DigitTrunk (double x, int Digits, string FloorOrCeil) {
    double res=0;
    //ofstream outputDigitTrunk(Machine + "DigitTrunk.txt", ofstream::app);
    /*
     double x=0; double y=0;
     gsl_rng * r = make_rng();
     gsl_rng_set(r, static_cast<unsigned long int>(time(0)/10)); // gsl_rng_set(const gsl_rng * r, unsigned long int s)
     x=100*gsl_rng_uniform(r);
     y=100*gsl_rng_uniform(r);
     outputV << "x=" << x << ", y=" << y << endl;
     outputV.close();
     */
    // Trunk is the interger version of x, and Residue is its decimal rest
    int Trunk=0; double Residue=0;
    if (FloorOrCeil=="Floor") {
        Trunk=int(floor(x));
        Residue=x-Trunk;
        Residue*=pow(10, Digits);
        Residue=floor(Residue);
        Residue/=pow(10, Digits);
        res=Trunk+Residue;
    }
    else if (FloorOrCeil=="Ceil") {
        Trunk=int(ceil(x));
        Residue=x-(Trunk-1);
        Residue*=pow(10, Digits);
        Residue=ceil(Residue);
        Residue/=pow(10, Digits);
        res=Trunk+Residue;
    }
    else if (FloorOrCeil=="Int") {
        Trunk=int(x);
        Residue=x-(Trunk-1);
        Residue*=pow(10, Digits);
        Residue=double(int(Residue));
        Residue/=pow(10, Digits);
        res=Trunk+Residue;
    };
    //outputDigitTrunk.precision(16);
    //outputDigitTrunk << "x=" << x << ", Digits=" << Digits << ", FloorOrCeil=" << FloorOrCeil <<" : result=" << res << endl;
    //outputDigitTrunk.close();
    if (abs(res-ceil(res))>0) {cout << "DIGITTRUNK() ISSUE : res=" << res << endl;};
    return res;
};


// Generating random variables for specific time delays
vector<int> STLRandomInt (int Length, string S, string Plot) {
    //gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937); // or gsl_rng_default instead of gsl_rng_mt19937
    gsl_rng * r = make_rng();
    gsl_rng_set(r, static_cast<unsigned long int>(time(0))); //gsl_rng_set(const gsl_rng * r, unsigned long int s)
    vector<int> V;
    srand (unsigned (time(0)));
    
    if (S=="Uniform") {for (int t=0; t<Length; t++) {V.push_back(int(1000*gsl_rng_uniform (r)));};};
    if (S=="Rand") {for (int t=0; t<Length; t++) {V.push_back(int(rand()));};};
    
    if (Plot=="On") {
        ofstream outputV(Machine + "V.txt", ofstream::app);
        for (int t=0; t<Length; t++) {outputV << V[t] << ", ";};
        outputV << endl <<endl;
        outputV.close();
    };
    
    gsl_rng_free (r);
    return V;
};




// Generating random variables for specific time delays
vector<double> STLRandom (int Length, string S, string Plot) {
    //gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937); // or gsl_rng_default instead of gsl_rng_mt19937
    gsl_rng * r = make_rng();
    gsl_rng_set(r, static_cast<unsigned long int>(time(0))); //gsl_rng_set(const gsl_rng * r, unsigned long int s)
    //gsl_rng * r2 = gsl_rng_alloc (gsl_rng_mt19937); // or gsl_rng_default instead of gsl_rng_mt19937
    gsl_rng * r2 = make_rng();
    gsl_rng_set(r2, static_cast<unsigned long int>(time(0))); //gsl_rng_set(const gsl_rng * r, unsigned long int s)
    vector<double> V;
    if (S=="Poisson") {for (int t=0; t<Length; t++) {V.push_back((gsl_ran_poisson (r, 1))*(gsl_ran_poisson (r, 1)));};};
    if (S=="PoissonOriginal") {for (int t=0; t<Length; t++) {V.push_back(gsl_ran_poisson (r, 1));};};
    if (S=="PoissonNew") {for (int t=0; t<Length; t++) {V.push_back(abs(10*(1+gsl_ran_gaussian (r, 1)+gsl_ran_poisson (r, 1))));};};
    if (S=="PoissonJump") {for (int t=0; t<Length; t++) {V.push_back(gsl_ran_poisson (r, 10));};};
    if (S=="StandardNormal") {for (int t=0; t<Length; t++) {V.push_back(gsl_ran_gaussian (r, 1));};};
    if (S=="StandardNormal2") {for (int t=0; t<Length; t++) {V.push_back(gsl_ran_gaussian (r2, 1));};};
    if (S=="NormalNew") {for (int t=0; t<Length; t++) {double x = 0.5+0.2*gsl_ran_gaussian (r, 1); if (x<0.05) {x=gsl_rng_uniform (r);}; if (x>0.95) {x=gsl_rng_uniform (r);}; V.push_back(x);};}; // This is our NEBLearningRate
    if (S=="Uniform") {for (int t=0; t<Length; t++) {V.push_back(gsl_rng_uniform (r));};};
    if (S=="UniformPercent") {for (int t=0; t<Length; t++) {V.push_back(100*(gsl_rng_uniform (r)));};};
    if (S=="UniformPercentFloor") {for (int t=0; t<Length; t++) {V.push_back(floor(100*(gsl_rng_uniform (r))));};};
    if (S=="C") {for (int t=0; t<Length; t++) {V.push_back(rand());};};
    
    if (Plot=="On") {
        ofstream outputV(Machine + "V.txt", ofstream::app);
        for (int t=0; t<Length; t++) {outputV << V[t] << ", ";};
        outputV << endl <<endl;
        outputV.close();
    };
    
    gsl_rng_free (r);
    return V;
};




// Function plotting an STL vector
void PlotSTL(vector<double> V, string Name, string XL) {
    // Randomization of output file name
    string A2 = Machine;
    string B2=Name;
    if (Name=="000") {int B = int(1000*((STLRandom(5, "Uniform", "NoPlot"))[2])); B2=IntToString(B);};
    A2+=B2;
    const char * Title2 = A2.c_str();
    ofstream outputV(Title2, ofstream::app);
    for (int t=0; t<int(V.size()); t++) {
        double x=V[t];
        if (XL=="XL") {x=int(x);};
        if (t==(int(V.size()))-1) {outputV << x; break;};
        outputV << x << ",";
    };
    outputV << endl << endl;
    outputV.close();
};




// Function plotting an STL vector
void PlotSTLInt(vector<int> V, string Name) {
    // Randomization of output file name
    string A2 = Machine;
    string B2=Name;
    if (Name=="000") {int B = int(1000*((STLRandom(5, "Uniform", "NoPlot"))[2])); B2=IntToString(B);};
    A2+=B2;
    const char * Title2 = A2.c_str();
    ofstream outputV(Title2, ofstream::app);
    for (int t=0; t<(int(V.size())); t++) {
        if (t==(int(V.size()))-1) {outputV << V[t]; break;};
        outputV << V[t] << endl;
    };
    outputV << endl << endl;
    outputV.close();
};




// Function plotting an STL vector
void PlotSTLMarket(vector<double> V, string Name) {
    // Randomization of output file name
    string A2 = Machine;
    string B2=Name;
    if (Name=="000") {int B = int(1000*((STLRandom(5, "Uniform", "NoPlot"))[2])); B2=IntToString(B);};
    A2+=B2;
    const char * Title2 = A2.c_str();
    ofstream outputV(Title2, ofstream::app);
    for (int t=0; t<(int(V.size())); t++) {
        int x = int(V[t]);
        if (t==(int(V.size()))-1) {outputV << x; break;};
        outputV << x << endl;
    };
    outputV << endl << endl;
    outputV.close();
};




// Function plotting a GSL matrix
void PlotGSLMatrix(gsl_matrix * M, string Name, int F) {
    // Randomization of output file name
    string A2 = Machine;
    string B2=Name;
    if (Name=="000") {int B = int(1000*((STLRandom(5, "Uniform", "NoPlot"))[2])); B2=IntToString(B);};
    A2+=B2;
    const char * Title2 = A2.c_str();
    ofstream outputM(Title2, ofstream::app);
    for (int t=0; t<int(M->size2); t++) {
        for (int j=0; j<int(M->size1); j++) {
            //int x = int(floor(F*gsl_matrix_get(M,j,t)));
            double x = F*gsl_matrix_get(M,j,t);
            if (j==(int(M->size1))-1) {outputM << fixed << setprecision(10) << x << endl; continue;};
            outputM << fixed << setprecision(10) << x << "\t";
        };
    };
    outputM << endl << endl;
    outputM.close();
};




// This compute the mean and variance from an STL vector
vector<double> StaticMoments (vector<double> V, int Start, int End, int Delta) {
    vector<double> Result;
    double Mean=0;
    double Variance=0;
    double Skewness=0;
    double Kurtosis=0;
    double RealizedVariance=0;
    double RealizedVolatility=0;
    double SSXlag=0; double SSXX=0; double SSXX2=0;
    //int SeriesSize=int(V.size()); // Original
    //if (End>=SeriesSize) {End=SeriesSize-1;}; // New
    if (Start>=End) {Start=End;}; // Original
    //if (Start<0) {Start=0;}; // New
    
    for (int i=Start; i<=End; i++) {
        Mean+=V[i];
        if (i-Delta>=0) {
            RealizedVariance+=(log(V[i]/V[i-Delta]))*(log(V[i]/V[i-Delta]));
            RealizedVolatility+=(log(V[i]/V[i-Delta]))*(log(V[i]/V[i-Delta]));
        };
    };
    Mean=Mean/(End-Start+1); // Mean computation
    RealizedVolatility=sqrt(252*RealizedVolatility/(End-Start+1));
    for (int i=Start; i<=End; i++) {
        double x=V[i] - Mean;
        Variance+=x*x;
        Skewness+=x*x*x;
        Kurtosis+=x*x*x*x;
    };
    Variance/=End-Start+1; // Variance computation
    Skewness=(Skewness/(End-Start+1))*(1/(Variance*sqrt(Variance)));
    Kurtosis=(Kurtosis/(End-Start+1))*(1/(Variance*Variance));
    
    for (int i=Start; i<=End-Delta; i++) {
        double x=V[i] - Mean;
        double y=V[i+Delta] - Mean;
        SSXlag+=x*y;
        SSXX+=x*x;
        SSXX2+=y*y;
    };
    double AutoCorr=SSXlag/((sqrt(SSXX))*(sqrt(SSXX2))); // Function returning the p-lag autocorrelation of a gsl_matrix, between time steps begin and end
    
    Result.push_back (Mean);
    Result.push_back (Variance);
    Result.push_back (Skewness);
    Result.push_back (Kurtosis); // Distributions with kurtosis<3 (>3) are platykurtic (leptokurtic)
    Result.push_back (sqrt(Variance));
    Result.push_back (RealizedVolatility);
    Result.push_back (AutoCorr);
    return Result;
};





// This computes the mean and variance for a given lag at each time step of a given time series
vector<vector<double> > DynamicMoments (vector<double> V, int Start, int End, int Lag, int Delta) {
    vector<vector<double> > Result;
    vector<double> ResultMean, ResultVariance, ResultSkewness, ResultKurtosis, ResultStdDev, TimeSeriesPlus3Sigma, TimeSeriesMinus3Sigma, ResultRealizedVolatility, ResultAutoCorr, ResultSquarredLogReturns1, ResultSquarredLogReturns2, ResultSquarredLogReturns3;
    if (End >= int(V.size())) {End=int(V.size()) -1;}; // Original
    Lag-=1;
    //if (Start>=End) {Start=End;}; // New
    //if (Start<0) {Start=0;}; // New
    
    double SquarredLogReturns1=0;
    double SquarredLogReturns2=0;
    double SquarredLogReturns3=0;
    for (int i=End; i>=Start; i--) {
        if (i<Start+Lag) {Lag-=1;};
        vector<double> SM = StaticMoments(V, i-Lag, i, Delta);
        double m = SM[0];
        double v = SM[1];
        double s = SM[2];
        double k = SM[3];
        double sd = SM[4];
        double rv = SM[5];
        double ac = SM[6];
        if (i-1>=0) {SquarredLogReturns1=(log(V[i]/V[i-1]))*(log(V[i]/V[i-1]));};
        if (i-Week>=0) {SquarredLogReturns2=(log(V[i]/V[i-Week]))*(log(V[i]/V[i-Week]));};
        if (i-2*Week>=0) {SquarredLogReturns3=(log(V[i]/V[i-2*Week]))*(log(V[i]/V[i-2*Week]));};
        ResultMean.push_back(m);
        ResultVariance.push_back(v);
        ResultSkewness.push_back(s);
        ResultKurtosis.push_back(k);
        ResultStdDev.push_back(sd);
        ResultRealizedVolatility.push_back(rv);
        ResultAutoCorr.push_back(ac);
        ResultSquarredLogReturns1.push_back(SquarredLogReturns1);
        ResultSquarredLogReturns2.push_back(SquarredLogReturns2);
        ResultSquarredLogReturns3.push_back(SquarredLogReturns3);
    };
    reverse(ResultMean.begin(), ResultMean.end());
    reverse(ResultVariance.begin(), ResultVariance.end());
    reverse(ResultSkewness.begin(), ResultSkewness.end());
    reverse(ResultKurtosis.begin(), ResultKurtosis.end());
    reverse(ResultStdDev.begin(), ResultStdDev.end());
    reverse(ResultRealizedVolatility.begin(), ResultRealizedVolatility.end());
    reverse(ResultAutoCorr.begin(), ResultAutoCorr.end());
    reverse(ResultSquarredLogReturns1.begin(), ResultSquarredLogReturns1.end());
    reverse(ResultSquarredLogReturns2.begin(), ResultSquarredLogReturns2.end());
    reverse(ResultSquarredLogReturns3.begin(), ResultSquarredLogReturns3.end());
    
    for (int i=Start; i<=End; i++) {
        TimeSeriesPlus3Sigma.push_back(V[i]+3*ResultStdDev[i-Start]);
        TimeSeriesMinus3Sigma.push_back(V[i]-3*ResultStdDev[i-Start]);
    };
    
    Result.push_back(ResultMean);
    Result.push_back(ResultVariance);
    Result.push_back(ResultSkewness);
    Result.push_back(ResultKurtosis);
    Result.push_back(ResultStdDev);
    Result.push_back(TimeSeriesPlus3Sigma);
    Result.push_back(TimeSeriesMinus3Sigma);
    Result.push_back(ResultRealizedVolatility);
    Result.push_back(ResultAutoCorr);
    Result.push_back(ResultSquarredLogReturns1);
    Result.push_back(ResultSquarredLogReturns2);
    Result.push_back(ResultSquarredLogReturns3);
    
    return Result;
};




//  AC-logreturns of [t-∆, t] & [t-2∆, t-∆] for ∆={3, w, 2w, 3w, m}
double MALAutoCorrelation (gsl_matrix * X, int t, int p, int Stock) {
    double Mean1=0; double Mean2=0; double SSXY=0; double SSXX=0; double SSYY=0;
    int Start1=max(t-p, 0); int Denominator1=max(t-Start1+1, 1);
    int Start2=max(t-2*p, 0); int Denominator2=max(Start1-Start2+1, 1);
    for (int i=Start1; i<=t; i++) {Mean1+=gsl_matrix_get (X, Stock, i)/Denominator1;};
    for (int i=Start2; i<=Start1; i++) {Mean2+=gsl_matrix_get (X, Stock, i)/Denominator2;};
    for (int i=Start1; i<=t; i++) {
        double x=gsl_matrix_get (X, Stock, i) - Mean1;
        double y=gsl_matrix_get (X, Stock, max(i-p, 0)) - Mean2;
        SSXX+=x*x;
        SSYY+=y*y;
        SSXY+=x*y;
    };
    double Result=SSXY/((sqrt(SSXX))*(sqrt(SSYY))); // Autocorrelation between intervals [t-p, t] and [t-2p, t-p]
    if (sqrt(SSXX)*sqrt(SSYY)==0) {Result=0;};
    if ((t-p<0) || (!Result)) {Result=0;}; // Failsafe
    return Result;
};


//  AC-logreturns of [t-∆, t] & [t-∆-∂,t-∂] for ∂={1, 3, w, 2w, 3w, m}
double MALAutoCorrelationBlend (gsl_matrix * X, int t, int p, int Shift, int Stock) {
    double Mean1=0; double Mean2=0; double SSXY=0; double SSXX=0; double SSYY=0;
    int Start1=max(t-p, 0); int Start2=max(t-p-Shift, 0);
    int Denominator=max(t-Start1+1, 1);
    for (int i=Start1; i<=t; i++) {Mean1+=gsl_matrix_get (X, Stock, i)/Denominator;};
    for (int i=Start2; i<=max(t-Shift, 0); i++) {Mean2+=gsl_matrix_get (X, Stock, i)/Denominator;};
    for (int i=Start1; i<=t; i++) {
        double x=gsl_matrix_get (X, Stock, i) - Mean1;
        double y=gsl_matrix_get (X, Stock, max(i-Shift, 0)) - Mean2;
        SSXX+=x*x;
        SSYY+=y*y;
        SSXY+=x*y;
    };
    double Result=SSXY/((sqrt(SSXX))*(sqrt(SSYY))); // Autocorrelation between intervals [t-p, t] and [t-2p, t-p]
    if (sqrt(SSXX)*sqrt(SSYY)==0) {Result=0;};
    if ((t-p<0) || (!Result)) {Result=0;}; // Failsafe
    return Result;
};


double MALVolatility (gsl_matrix * X, int t, int p, int Stock) {
    double Mean=0; double Variance=0;
    int Start=max(t-p, 0); int Denominator=max(t-Start+1, 1);
    for (int i=Start; i<=t; i++) {Mean+=gsl_matrix_get (X, Stock, i)/Denominator;};
    for (int i=Start; i<=t; i++) {double x=gsl_matrix_get (X, Stock, i) - Mean; Variance+=x*x;};
    double Result=sqrt(Variance/Denominator)/Mean; // Standard deviation on the interval [t-p, t]
    if ((t-p<0) || (!Result)) {Result=0;}; // Failsafe
    //if (p==2*Week) {cout << "Denominator=" << Denominator << endl;};
    return Result;
};


double MALVolatility2 (gsl_matrix * X, int t1, int t2, int Stock) {
    double Mean=0; double Variance=0;
    int Start=max(t1, 0); int Denominator=max(t2-Start+1, 1);
    for (int i=Start; i<=t2; i++) {Mean+=gsl_matrix_get (X, Stock, i)/Denominator;};
    for (int i=Start; i<=t2; i++) {double x=gsl_matrix_get (X, Stock, i) - Mean; Variance+=x*x;};
    double Result=sqrt(Variance/Denominator)/Mean; // Standard deviation on the interval [t-p, t]
    if ((t1<0) || (!Result)) {Result=0;}; // Failsafe
    //if (p==2*Week) {cout << "Denominator=" << Denominator << endl;};
    return Result;
};




double StaticAutoCorrelation (gsl_matrix * X, int Start, int End, int p, int Stock, int j) { // BBB
    double Mean=0; double SSXlag=0; double SSXX=0; double SSXX2=0;
    for (int i=Start; i<=End; i++) {Mean+=gsl_matrix_get (X, Stock+8*j, i);}; Mean=Mean/(End-Start+1); // Mean of MarketBidAskTrue!!!
    for (int i=Start; i<=End-p; i++) {
        double x=gsl_matrix_get (X, Stock+8*j, i) - Mean;
        double y=gsl_matrix_get (X, Stock+8*j, i+p) - Mean;
        SSXlag+=x*y;
        SSXX+=x*x;
        SSXX2+=y*y;
    };
    double Result=SSXlag/((sqrt(SSXX))*(sqrt(SSXX2))); // p-lag autocorrelation of a gsl_matrix, between Start and End
    return Result;
};






gsl_matrix * DynamicAutoCorrelation (gsl_matrix * X, int Lag, int Stock, int j) {
    gsl_matrix * Result = gsl_matrix_alloc (5, int(X->size2));
    for (int i=0; i<int(X->size2)-Lag; i++) {
        gsl_matrix_set (Result, 0+6*j, i, StaticAutoCorrelation (X, i, i+Lag, DayTicks, Stock, j));
        gsl_matrix_set (Result, 1+6*j, i, StaticAutoCorrelation (X, i, i+Lag, 2*DayTicks, Stock, j));
        gsl_matrix_set (Result, 2+6*j, i, StaticAutoCorrelation (X, i, i+Lag, Week, Stock, j));
        gsl_matrix_set (Result, 3+6*j, i, StaticAutoCorrelation (X, i, i+Lag, 2*Week, Stock, j));
        gsl_matrix_set (Result, 4+6*j, i, StaticAutoCorrelation (X, i, i+Lag, Month, Stock, j));
    }
    return Result;
};






gsl_matrix * ProgressiveAutoCorrelation (gsl_matrix * X, int Stock, int j) {
    gsl_matrix * Result = gsl_matrix_calloc (5, int(X->size2));
    for (int i=1; i<int(X->size2); i++) {
        if (i>DayTicks) {gsl_matrix_set (Result, 0+6*j, i, StaticAutoCorrelation (X, 0, i, DayTicks, Stock, j));};
        if (i>2*DayTicks) {gsl_matrix_set (Result, 1+6*j, i, StaticAutoCorrelation (X, 0, i, 2*DayTicks, Stock, j));};
        if (i>Week) {gsl_matrix_set (Result, 2+6*j, i, StaticAutoCorrelation (X, 0, i, Week, Stock, j));};
        if (i>2*Week) {gsl_matrix_set (Result, 3+6*j, i, StaticAutoCorrelation (X, 0, i, 2*Week, Stock, j));};
        if (i>Month) {gsl_matrix_set (Result, 4+6*j, i, StaticAutoCorrelation (X, 0, i, Month, Stock, j));};
    }
    return Result;
};






gsl_matrix * ExtensiveAutoCorrelation (gsl_matrix * X, int Stock, int j) {
    gsl_matrix * Result = gsl_matrix_calloc (5, int(X->size2));
    for (int i=1; i<int(X->size2); i++) {
        gsl_matrix_set (Result, 0+6*j, i, StaticAutoCorrelation (X, 0, int(X->size2)-1, i, Stock, j));
    }
    return Result;
};






double StaticCorrelation (gsl_matrix * X, gsl_matrix * Y, int XStart, int XEnd, int p, int Stock, int j) {
    double XMean=0; double YMean=0; double SSXY=0; double SSXX=0; double SSYY=0;
    for (int i=XStart; i<=XEnd; i++) {XMean+=gsl_matrix_get (X, Stock+6*j, i)/(XEnd-XStart+1);}; // Mean of MarketBidAskTrue!!!
    for (int i=XStart+p; i<=XEnd+p; i++) {YMean+=gsl_matrix_get (Y, Stock+6*j, i)/(XEnd-XStart+1);}; // Mean of MarketBidAskTrue!!!
    for (int i=XStart; i<=XEnd; i++) {
        double x=gsl_matrix_get (X, Stock+6*j, i) - XMean;
        double y=gsl_matrix_get (Y, Stock+6*j, i+p) - YMean;
        SSXY+=x*y;
        SSXX+=x*x;
        SSYY+=y*y;
    };
    return SSXY/((sqrt(SSXX))*(sqrt(SSYY)));
};





gsl_matrix * CFMCorrelation (gsl_matrix * X, gsl_matrix * Y, int XStart, int XEnd, int Lag, int Stock, int j) {
    gsl_matrix * Result = gsl_matrix_calloc (1, Lag);
    for (int i=0; i<Lag; i++) {gsl_matrix_set (Result, 0, i, StaticCorrelation (X, Y, XStart, XEnd, i, Stock, j));} // Be careful that the Lag+XEnd <= int(X->size2) !!!
    return Result;
};






// Function returning the p-lag autocorrelation of a gsl_matrix, between time steps begin and end
gsl_matrix * AutocorrelationLagP (gsl_matrix * X) {
    gsl_matrix * Res=gsl_matrix_alloc (int(X->size1), int(X->size2));
    for (int j=0; j<int(X->size1); j++) {
        for (int p=0; p<int(X->size2); p++) {
            double SSXlag=0; double SSXX1=0; double SSXX2=0; double Mean=0;
            
            for (int t=0; t<int(X->size2); t++) {Mean+=gsl_matrix_get (X, j, t);};
            Mean/=int(X->size2); // Mean computation
            
            for (int t=0; t<int(X->size2)-p; t++) {
                double x=gsl_matrix_get (X, j, t) - Mean;
                double y=gsl_matrix_get (X, j, t+p) - Mean;
                SSXlag+=x*y;
                SSXX1+=x*x;
                SSXX2+=y*y;
            }; // closes t loop
            gsl_matrix_set (Res, j, p, SSXlag/((sqrt(SSXX1))*(sqrt(SSXX2)))); // Function returning the p-lag autocorrelation
        }; // closes p loop
    }; // closes j loop
    
    return Res;
}








// Function returning the Pearson correlation between two vectors
double StaticCorrelation (vector<double> X, vector<double> Y, int Start, int End) {
    double Res; double SSXY=0; double SSXX=0; double SSYY=0; double MeanX=0; double MeanY=0;
    int Size=min(int(X.size()), int(Y.size()));
    if (Start<0) {Start=0;};
    if (End>Size) {End=Size;};
    if (End==0) {End=Size;};
    for (int i=Start; i<=End; i++) {MeanX+=X[i]; MeanY+=Y[i];}; MeanX/=End-Start+1; MeanY/=End-Start+1; // Mean
    
    for (int i=Start; i<=End; i++) {SSXY+=(X[i] - MeanX)*(Y[i] - MeanY);}
    for (int i=Start; i<=End; i++) {SSXX+=(X[i] - MeanX)*(X[i] - MeanX);}
    for (int i=Start; i<=End; i++) {SSYY+=(Y[i] - MeanY)*(Y[i] - MeanY);}
    Res=SSXY/(sqrt(SSXX*SSYY));
    
    return Res;
};








vector<double> DynamicCorrelation (vector<double> X, vector<double> Y, int Lag) {
    vector<double> Res;
    int Size=min(int(X.size()), int(Y.size()));
    for (int k=0; k<Lag; k++) {Res.push_back(0);};
    for (int k=0; k<Size-Lag; k++) {
        Res.push_back(StaticCorrelation(X,Y,k,k+Lag));
    };
    
    return Res;
};







vector<vector<double> > GSLMatrixToSTLMatrix (gsl_matrix * M) {
    vector<vector<double> > STLM;
    vector<double> TempLine;
    for (int i=0; i<int(M->size1); i++) {
        for (int j=0; j<int(M->size2); j++) {
            //vector<double> Temp; V.push_back(Temp); double x=gsl_matrix_get(M,i,j); if (x!=0) {(V[i]).push_back(x);};
            double x=gsl_matrix_get(M,i,j);
            TempLine.push_back(x);
        };
        STLM.push_back(TempLine);
        TempLine.clear();
    };
    return STLM;
};



// Function to shuffle a given STL vector V={0,1,2,3,...,N} (including 0, excluding N)
int myrandom (int i) {srand (unsigned (time(0))); return rand()%i;};
vector<int> Shuffle(int N)
{
    vector<int> V;
    srand ( unsigned ( time(0) ) );
    for (int i=0; i<N; ++i) V.push_back(i);
    random_shuffle ( V.begin(), V.end() );
    random_shuffle ( V.begin(), V.end(), myrandom);
    return V;
};



vector<int> Shuffle2(int N) // Including N
{
    vector<int> V;
    srand ( unsigned ( time(0) ) );
    for (int i=0; i<=N; ++i) V.push_back(i);
    random_shuffle ( V.begin(), V.end() );
    random_shuffle ( V.begin(), V.end(), myrandom);
    return V;
};




// Fail for Tool=2 and Lag=0 because Start=
double BinaryProjection (gsl_matrix * ReflexiveValues, int t, int Tool, int Lag, int Future) {
    double Result=0; int Start=t-(Lag+1)*Future; if (Start<0) {Start=0;}; int Past=t-Start;
    double Mean0=0; for (int i=Start; i<t-Past/2; i++) {Mean0+=gsl_matrix_get(ReflexiveValues, 0, i)/(Past/2);};
    double Mean1=0; for (int i=t-Past/2; i<t; i++) {Mean1+=gsl_matrix_get(ReflexiveValues, 0, i)/(Past/2);};
    if (Tool==0) {Result=gsl_matrix_get(ReflexiveValues, 0, t)-(Mean1-Mean0);} // counter-trend following
    else if (Tool==1) {Result=0.5*(Mean1+Mean0);} // mean reversion
    else if (Tool==2) {Result=gsl_matrix_get(ReflexiveValues, 0, t)+(Mean1-Mean0);} // trend following
    if (Result<0.001) {Result=0.001;}; // Failsafe
    //if (Result<1.0) {cout << "Pt=" << gsl_matrix_get(ReflexiveValues, 0, t) << ", Mean0=" << Mean0 << ", Mean1=" << Mean1 << ", Tool=" << Tool << ", Future=" << Future << ", Lag=" << Lag << ", Result=" << Result << endl;}; // Failsafe
    return Result;
};





// KEVIN
// Fonction A(i) dans la log-likelihood
void function_A(vector<double>* t ,double beta, vector<double>& A ){
    size_t n = (*t).size();
    A[0] = 0;
    for (size_t i = 1; i < n; i++) {A[i] = A[i-1]*exp(-beta*( (*t)[i] - (*t)[i-1] )) + exp(-beta*( (*t)[i] - (*t)[i-1] ));}
};

// Fonction log-likelihood
double function_f(const gsl_vector* v, void* params){
    double alpha = gsl_vector_get(v,0);
    double beta = gsl_vector_get(v,1);
    double mu = gsl_vector_get(v,2);
    vector<double>* t = (vector<double>*) params;
    size_t n = (*t).size();
    vector<double> A (n);
    function_A(t,beta,A);
    const double b = alpha/beta;
    double r = mu*(*t)[n-1]-b*(A[n-1] + 1);
    for (size_t i = 0; i < n; i++) {r = r -log(mu + alpha*A[i]) + b;}
    return r;
};

// Fonction B(i) dans la log-likelihood
void function_B( vector<double>* t,double beta, vector<double>& B){
    size_t n = (*t).size();
    B[0] = 0;
    for (size_t i = 1; i < n; i++) {
        double tmp = 0;
        for (size_t j = 0; j < i; j++) {tmp = tmp + ( (*t)[i] - (*t)[j])*exp(-beta*((*t)[i]-(*t)[j]));}
        B[i] = tmp;
    };
};

// Gradient de la log-likelihood
void function_df(const gsl_vector* v, void* params, gsl_vector* df){
    double alpha = gsl_vector_get(v,0);
    double beta = gsl_vector_get(v,1);
    double mu = gsl_vector_get(v,2);
    vector<double> *t = (vector<double> *) params;
    size_t n = (*t).size();
    vector<double> A (n);
    vector<double> B (n);
    function_A(t, beta, A);
    function_B(t, beta, B);
    double dmu = (*t)[n-1];
    double dalpha = 0;
    double dbeta = 0;
    for (size_t i = 0; i < n; i++) {dmu = dmu - 1/( mu + alpha*A[i]);}
    for (size_t i = 0; i < n; i++) {dalpha = dalpha - (A[i]/(mu + alpha*A[i])) + (1/beta);}
    dalpha = dalpha - (1/beta)*(A[n-1] + 1 );
    const double C1 = alpha/beta;
    const double C2 = alpha/pow(beta,2);
    for (size_t i = 0; i < n; i++) {dbeta = dbeta + (alpha*B[i]/(mu + alpha*A[i])) - C2 ;}
    dbeta = dbeta + C1*B[n-1] + C2*A[n-1]+C2;
    gsl_vector_set(df,0,dalpha);
    gsl_vector_set(df,1,dbeta);
    gsl_vector_set(df, 2,dmu);
}

// Fonction fdf, calcule la fonction f et le gradient de f
void function_fdf(const gsl_vector* v, void* params, double* f, gsl_vector* df) {
    *f = function_f(v,params);
    function_df(v, params, df );
}

// Minimizer de la log-likelihood
int EstimationOzaki(std::vector<double>* t, double* Init, const double stepsize, const size_t NbIteration, double* R){
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = 0; // NULL
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;
    size_t iter = 0;
    int status;
    double size;
    // Starting point
    x = gsl_vector_alloc (3);
    gsl_vector_set (x, 0, Init[0]);
    gsl_vector_set (x, 1, Init[1]);
    gsl_vector_set (x, 2, Init[2]);
    // Set initial step sizes to 1
    ss = gsl_vector_alloc (3);
    gsl_vector_set_all (ss, stepsize);
    // Initialize method and iterate
    minex_func.n = 3;
    minex_func.f = & function_f;
    minex_func.params =  t;
    s = gsl_multimin_fminimizer_alloc (T, 3);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 0.0003);
    }
    while ((status == GSL_CONTINUE) && (iter < NbIteration));
    R[0] = gsl_vector_get(s->x,0);
    R[1] = gsl_vector_get(s->x,1);
    R[2] = gsl_vector_get(s->x,2);
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
    return status;
}

// Calcul de lambda(t)
double ExpIntensityState(const double mu, const double alpha, const double beta,std::vector<double>* t, const double rescale ,const size_t n){
    double res = mu;
    const double T = (double) n/rescale;
    for (vector<double>::iterator it = (*t).begin(); (it != (*t).end()) && (*it <= T); it++) {res = res + alpha*exp(-beta*(T-*it));};
    return res;
}

// Fonction d'integration
void IntegrateExpIntensity(const double mu, const double alpha, const double beta, std::vector<double>* t, vector<double>& I){
    size_t n = (*t).size();
    vector<double> A(n);
    function_A(t,beta,A);
    vector<double> tmp;
    I.push_back(0);
    for (size_t i = 1; i < n; i++) {
        double res = mu*(t[0])[i];
        res = res - (alpha/beta)*A[i] + i;
        tmp.push_back(res);
    };
    for (size_t i = 1; i < n; i++) {
        I.push_back(tmp[i]-tmp[i-1]);
    };
};

/*
// Fonction retournant l'intensité du processus de hawkes; t = les closing prices dans l'ordre chronologique ; pct = seuil (0.1); NbIteration = le nombre d'iteration du minimiseur; R = les parametres du processus de Hawkes; I = goodness of fit
vector<double>* Intensity(std::vector<double>* t, const double rescale, const double pct , const size_t NbIteration, double* R, vector<double>& I){
    size_t n = (*t).size();
    std::vector<double>* res = new vector<double>[2]; // JJJ10
    std::vector<double> jump;
    double logr[n-1];
    for (size_t i = 0; i < n-1; i++) {
        logr[i] = log((*t)[i+1]/(*t)[i]);
    };
    gsl_rstat_quantile_workspace* quantile1 = gsl_rstat_quantile_alloc(pct); // POSTDOC2022
    gsl_rstat_quantile_workspace* quantile2 = gsl_rstat_quantile_alloc(1-pct); // POSTDOC2022
    for (size_t i = 0; i < n-1; i++) {
        gsl_rstat_quantile_add(logr[i], quantile1);
        gsl_rstat_quantile_add(logr[i], quantile2);
    };
    double bsup = gsl_rstat_quantile_get(quantile2);
    double binf = gsl_rstat_quantile_get(quantile1);
    for (size_t i = 1; i < n-1; i++) {
        double tmp = logr[i];
        if ( tmp > bsup ) {jump.push_back((i-1)/rescale);};
        if ( tmp < binf) {jump.push_back((i-1)/rescale);};
    };
    double Init[3];
    Init[0] = 0.6;
    Init[1] = 0.4;
    Init[2] = 0.5;
    double stepsize = 0.001;
    EstimationOzaki(&jump, Init, stepsize, NbIteration, R);
    for (size_t i = 0; i < n-1; i++) {
        double tmp = ExpIntensityState( R[2], R[0], R[1], &jump, rescale, i );
        (res[0]).push_back(tmp);
        (res[1]).push_back(i/rescale);
    };
    IntegrateExpIntensity(R[2], R[0], R[1], &jump, I);
    return res;
}

// Hawkes process historical intensity
gsl_matrix * HistoricalIntensity (gsl_matrix * TimeSeries, int k) {
    int T = int(TimeSeries->size2);
    gsl_matrix * Res = gsl_matrix_calloc (1, T);
    vector<double> V, GoF; for (int t=0; t<T; t++) {V.push_back(gsl_matrix_get (TimeSeries, k, t)); GoF.push_back(1);};
    vector<double> * pV= & V;
    const double Rescale=10; // Rescale factor
    const double Threshold=0.1; // Threshold (0.1)
    const size_t NumIteration=100000; // Number of iterations of the minimizer (stops before once a minimum is found)
    double R; double * pR = & R; // Hawkes process parameters (does not need to be initialized)
    vector<double> & pGoF = GoF; // Goodness of fit
    vector<double> * pRes = Intensity (pV, Rescale, Threshold, NumIteration, pR, pGoF);
    for (int t=0; t<T; t++) {gsl_matrix_set (Res, 0, t, (*pRes)[t]);}; // Hawkes process intensity
    gsl_matrix_set (Res, 0, T-1, gsl_matrix_get (Res, 0, T-2));
    delete[] pRes; // JJJ10 void gsl_matrix_free(gsl_matrix * m)
    return Res;
}; // closes HistoricalIntensity()

// Hawkes process spot intensity
gsl_matrix * SpotIntensity (gsl_matrix * TimeSeries, int k) {
    gsl_matrix * Res = gsl_matrix_calloc (1, int(TimeSeries->size2));
    for (int T=5; T<int(TimeSeries->size2); T++) {
        vector<double> V, GoF; for (int t=0; t<T; t++) {V.push_back(gsl_matrix_get (TimeSeries, k, t)); GoF.push_back(1);};
        vector<double> * pV= & V;
        const double Rescale=10; // Rescale factor
        const double Threshold=0.1; // Threshold (0.1)
        const size_t NumIteration=100000; // Number of iterations of the minimizer (stops before once a minimum is found)
        double R; double * pR = & R; // Hawkes process parameters (does not need to be initialized)
        vector<double> & pGoF = GoF; // Goodness of fit
        vector<double> * pRes = Intensity (pV, Rescale, Threshold, NumIteration, pR, pGoF);
        gsl_matrix_set (Res, 0, T, (*pRes)[T-2]); // Hawkes process intensity
    }; // closes T loop
    return Res;
}; // closes SpotIntensity()

// Hawkes process spot intensity with first 100 time steps at mean
gsl_matrix * SpotIntensity2 (gsl_matrix * TimeSeries, int k) {
    gsl_matrix * Res = gsl_matrix_calloc (1, int(TimeSeries->size2));
    for (int T=5; T<int(TimeSeries->size2); T++) {
        vector<double> V, GoF; for (int t=0; t<T; t++) {V.push_back(gsl_matrix_get (TimeSeries, k, t)); GoF.push_back(1);};
        vector<double> * pV= & V;
        const double Rescale=10; // Rescale factor
        const double Threshold=0.1; // Threshold (0.1)
        const size_t NumIteration=1000000; // Number of iterations of the minimizer (stops before once a minimum is found)
        double R; double * pR = & R; // Hawkes process parameters (does not need to be initialized)
        vector<double> & pGoF = GoF; // Goodness of fit
        vector<double> * pRes = Intensity (pV, Rescale, Threshold, NumIteration, pR, pGoF);
        gsl_matrix_set (Res, 0, T, (*pRes)[T-2]); // Hawkes process intensity
    }; // closes T loop
    double Mean=0;
    for (int T=100; T<int(TimeSeries->size2); T++) {Mean+=gsl_matrix_get (Res, 0, T)/(int(TimeSeries->size2)-100);};
    for (int T=0; T<100; T++) {gsl_matrix_set (Res, 0, T, Mean);};
    return Res;
}; // closes SpotIntensity2()
*/


// ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE ***
// ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE ***
// ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE ***
// ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE ***


// Takes as input the file constructed by ImportMacro (). Then the file must be updated by hand as such:
// LIBOR rates taken from the FED in StLouis https://fred.stlouisfed.org/categories/33003/downloaddata
// Exchange rates taken from the Bank of England at http://www.bankofengland.co.uk/boeapps/iadb/Rates.asp?TD=12&TM=Dec&TY=2017&
gsl_matrix * Macroeconomics () {
    double FXFee=1-0.3*0.01; // FX fees from IG
    string Path=Machine + "Symba/CSV/MACROECONOMICS.csv";
    ifstream FileInput(Path); string Line; int LineNb=0;
    vector<double> V1, V2, V3, V4, V5;
    while (getline (FileInput, Line)) {LineNb++;
        istringstream LineStream(Line); string Item; int ItemNb=0;
        while (getline (LineStream, Item, ',')) {ItemNb++;
            if ((ItemNb==1) && (LineNb>1)) {V1.push_back(atof(Item.c_str()));};
            if ((ItemNb==2) && (LineNb>1)) {V2.push_back(atof(Item.c_str()));};
            if ((ItemNb==3) && (LineNb>1)) {V3.push_back(atof(Item.c_str()));};
            if ((ItemNb==4) && (LineNb>1)) {V4.push_back(atof(Item.c_str()));};
            if ((ItemNb==5) && (LineNb>1)) {V5.push_back(atof(Item.c_str()));};
        };
    };
    int Size=int(V1.size()); FileInput.close();
    gsl_matrix * Result = gsl_matrix_calloc (6, Size);
    for (int i=0; i<Size; i++) {
        gsl_matrix_set (Result, 0, i, V1[i]); // Dates
        gsl_matrix_set (Result, 1, i, V2[i]); // LIBOR in percent and based on USD
        gsl_matrix_set (Result, 2, i, V3[i]); // LIBOR in percent and based on GBP
        gsl_matrix_set (Result, 3, i, V4[i]*FXFee); // 1 EUR in USD
        gsl_matrix_set (Result, 4, i, V5[i]*FXFee); // 1 EUR in GBP
        gsl_matrix_set (Result, 5, i, 1); // No currency conversion
    };
    /*
     for (int i=1; i<Size; i++) {
     if ((abs(gsl_matrix_get(Result, 1, i))>0.1) || (abs(gsl_matrix_get(Result, 1, i))==0)) {gsl_matrix_set(Result, 1, i, gsl_matrix_get(Result, 1, i-1));};
     if ((abs(gsl_matrix_get(Result, 2, i))>0.1) || (abs(gsl_matrix_get(Result, 2, i))==0)) {gsl_matrix_set(Result, 2, i, gsl_matrix_get(Result, 2, i-1));};
     if ((abs(gsl_matrix_get(Result, 3, i))>2) || (abs(gsl_matrix_get(Result, 3, i))==0)) {gsl_matrix_set(Result, 3, i, gsl_matrix_get(Result, 3, i-1));};
     if ((abs(gsl_matrix_get(Result, 4, i))>2) || (abs(gsl_matrix_get(Result, 4, i))==0)) {gsl_matrix_set(Result, 4, i, gsl_matrix_get(Result, 4, i-1));};
     };
     */
    //PlotGSLMatrix (Result, "Macroeconomics.csv", 1);
    return Result;
};





// Barclay Equity L/S Index between 1.1.2007 and 1.1.2018 (132 months, one month being 21.75=22 days)
void StockIndices () {
    int Size=2872;
    string Path=Machine + "Symba/CSV/BarclayHFI.csv"; vector<double> V;
    ifstream FileInput(Path); string Line; int LineNb=0;
    while (getline (FileInput, Line)) {LineNb++;
        istringstream LineStream(Line); string Item; int ItemNb=0;
        while (getline (LineStream, Item, '\r')) {ItemNb++; V.push_back(atof(Item.c_str()));};
    };
    // Output on 132 months
    gsl_matrix * Result = gsl_matrix_calloc (1, 132); gsl_matrix_set (Result, 0, 0, 1);
    for (int i=1; i<132; i++) {gsl_matrix_set (Result, 0, i, gsl_matrix_get (Result, 0, i-1)*(V[i-1]+100)/100.0);}
    // Output on 2520 days
    gsl_matrix * Result2 = gsl_matrix_calloc (2, Size);
    for (int i=0; i<Size; i++) {gsl_matrix_set (Result2, 0, i, gsl_matrix_get (Result, 0, i/22));}
    
    Path=Machine + "Symba/CSV/NasdaqComposite.csv"; vector<double> W;
    ifstream FileInput2(Path); string Line2; LineNb=0;
    while (getline (FileInput2, Line2)) {LineNb++;
        istringstream LineStream(Line2); string Item; int ItemNb=0;
        while (getline (LineStream, Item, '\r')) {ItemNb++; W.push_back(atof(Item.c_str()));};
    };
    vector<int> J = Shuffle(int(W.size())-10);
    int Diff=Size-int(W.size());
    vector<double> N;
    for (int i=0; i<int(W.size()); i++) {
        N.push_back(W[i]);
        for (int k=0; k<Diff; k++) {
            if (J[k]==i) {N.push_back(W[i]);};
        };
    }
    for (int i=0; i<Size; i++) {
        gsl_matrix_set (Result2, 1, i, N[i]/N[0]);
    }
    PlotGSLMatrix (Result2, "StockIndices.csv", 1); // Column 0 is Barclay, column 1 is Nasdaq
};






class Share {
public:
    string Symbol, Title, Exchange, File, Country, Currency, Sector, InPF;
    string temp_price, temp_volume, temp_date;
    vector<double> close_prices;
    vector<double> volumes;
    vector<string> dates;
    gsl_matrix * Data; // Prices & cashflows
    
    
    void Gen (string Name, gsl_matrix * Macroeconomics, int FirstDate) {
        //int Numeraire=5;
        //if (Currency=="USD") {Numeraire=3;} else if (Currency=="GBP") {Numeraire=4;}; // This is switched off so as to study local market microstructure apart from FX perturbations // FXX
        vector<double> LocalDates, LocalOpen, LocalClose, LocalHigh, LocalLow, LocalVolumes;
        ifstream FileInput(Name); string Line; int LineNb=0;
        while (getline (FileInput, Line)) {LineNb++;
            istringstream LineStream(Line); string Item; int ItemNb=0;
            while (getline (LineStream, Item, ',')) {ItemNb++;
                if ((LineNb>1) && (ItemNb==3)) {LocalDates.push_back(atof(Item.c_str()));};
                if ((LineNb>1) && (ItemNb==5)) {LocalOpen.push_back(atof(Item.c_str()));};
                if ((LineNb>1) && (ItemNb==8)) {LocalClose.push_back(atof(Item.c_str()));};
                if ((LineNb>1) && (ItemNb==9)) {LocalVolumes.push_back(atof(Item.c_str()));};
                if ((LineNb>1) && (ItemNb==6)) {LocalHigh.push_back(atof(Item.c_str()));};
                if ((LineNb>1) && (ItemNb==7)) {LocalLow.push_back(atof(Item.c_str()));};
            };
        };
        FileInput.close();
        // Computing the data size from FirstDate to today
        int Size=0; bool Cond=0;
        for (int i=0; i<int(Macroeconomics->size2); i++) {
            if (gsl_matrix_get (Macroeconomics, 0, i)==FirstDate) {Cond=1;};
            if (Cond==1) {Size+=1;};
        };
        // Populating gsl_matrix * Data
        int PastLag=120+Year; // Steps before FirstDate necessary for BinaryProjector() to start at FirstDate
        int FirstStep=0; // STL vector index where data matches FirstDate
        int VectorSize=int(LocalDates.size());
        Data = gsl_matrix_calloc (6, Size+PastLag);
        Cond=0;
        for (int i=0; i<Size; i++) {gsl_matrix_set(Data, 0, i+PastLag, gsl_matrix_get (Macroeconomics, 0, i));}; //Dates
        for (int i=0; i<Size; i++) {
            for (int k=0; k<VectorSize; k++) {
                if ((LocalDates[k]==gsl_matrix_get (Macroeconomics, 0, i)) && Cond==0) {FirstStep=k; Cond=1;};
                if (LocalDates[k]==gsl_matrix_get (Macroeconomics, 0, i)) {
                    gsl_matrix_set(Data, 1, i+PastLag, LocalOpen[k]); //Open in €
                    gsl_matrix_set(Data, 2, i+PastLag, LocalClose[k]); //Close in €
                    gsl_matrix_set(Data, 3, i+PastLag, LocalVolumes[k]); // Daily volume
                    gsl_matrix_set(Data, 4, i+PastLag, LocalHigh[k]); //High in €
                    gsl_matrix_set(Data, 5, i+PastLag, LocalLow[k]); //Low in €
                    break;
                };
            }; // closes k loop
        }; // closes i loop
        // Stiching so as to match Macroeconomics and other exchanges data
        for (int i=1; i<Size; i++) {
            if (gsl_matrix_get(Data, 1, i+PastLag)<0.0000001) {gsl_matrix_set(Data, 1, i+PastLag, gsl_matrix_get(Data, 1, i+PastLag-1));};
            if (gsl_matrix_get(Data, 2, i+PastLag)<0.0000001) {gsl_matrix_set(Data, 2, i+PastLag, gsl_matrix_get(Data, 2, i+PastLag-1));};
            if (gsl_matrix_get(Data, 3, i+PastLag)<0.0000001) {gsl_matrix_set(Data, 3, i+PastLag, gsl_matrix_get(Data, 3, i+PastLag-1));};
            if (gsl_matrix_get(Data, 4, i+PastLag)<0.0000001) {gsl_matrix_set(Data, 4, i+PastLag, gsl_matrix_get(Data, 4, i+PastLag-1));};
            if (gsl_matrix_get(Data, 5, i+PastLag)<0.0000001) {gsl_matrix_set(Data, 5, i+PastLag, gsl_matrix_get(Data, 5, i+PastLag-1));};
        };
        // Populating the region 0<t<PastLag for BinaryProjector
        if (FirstStep>PastLag) {
            for (int k=FirstStep-PastLag; k<FirstStep; k++) {
                gsl_matrix_set(Data, 0, k-FirstStep+PastLag, LocalDates[k]);
                gsl_matrix_set(Data, 1, k-FirstStep+PastLag, LocalOpen[k]);
                gsl_matrix_set(Data, 2, k-FirstStep+PastLag, LocalClose[k]);
                gsl_matrix_set(Data, 3, k-FirstStep+PastLag, LocalVolumes[k]);
                gsl_matrix_set(Data, 4, k-FirstStep+PastLag, LocalHigh[k]);
                gsl_matrix_set(Data, 5, k-FirstStep+PastLag, LocalLow[k]);
            };
        };
        
        // Stitching stock splits
        for (int t=1; t<Size+PastLag; t++) {
            double Split=1;
            if (gsl_matrix_get(Data, 2, t-1)>0.00001) {Split=gsl_matrix_get(Data, 2, t)/gsl_matrix_get(Data, 2, t-1);};
            if (((Split>=4) || (Split<=0.25)) && (Split>0)) {
                for (int i=t; i<Size+PastLag; i++) {gsl_matrix_set(Data, 2, i, gsl_matrix_get(Data, 2, i)/Split);};
            }; // if
        }; // t loop
        
        
    }; // closes Gen()
    
    void Gen_trunk_old (string Name, gsl_matrix * Macroeconomics, int FirstDate, int LastDate) {
        //int Numeraire=5;
        //if (Currency=="USD") {Numeraire=3;} else if (Currency=="GBP") {Numeraire=4;}; // This is switched off so as to study local market microstructure apart from FX perturbations // FXX
        vector<double> LocalDates, LocalOpen, LocalClose, LocalHigh, LocalLow, LocalVolumes;
        ifstream FileInput(Name); string Line; int LineNb=0;
        while (getline (FileInput, Line)) {LineNb++;
            istringstream LineStream(Line); string Item; int ItemNb=0;
            while (getline (LineStream, Item, ',')) {ItemNb++;
                if ((LineNb>1) && (ItemNb==3)) {LocalDates.push_back(atof(Item.c_str()));};
                if ((LineNb>1) && (ItemNb==5)) {LocalOpen.push_back(atof(Item.c_str()));};
                if ((LineNb>1) && (ItemNb==8)) {LocalClose.push_back(atof(Item.c_str()));};
                if ((LineNb>1) && (ItemNb==9)) {LocalVolumes.push_back(atof(Item.c_str()));};
                if ((LineNb>1) && (ItemNb==6)) {LocalHigh.push_back(atof(Item.c_str()));};
                if ((LineNb>1) && (ItemNb==7)) {LocalLow.push_back(atof(Item.c_str()));};
            };
        };
        FileInput.close();
        // Computing the data size from FirstDate to today
        int Size=0; bool Cond=0;
        /*
        for (int i=0; i<int(Macroeconomics->size2); i++) {
            if (gsl_matrix_get (Macroeconomics, 0, i)==FirstDate) {Cond=1;};
            if (Cond==1) {Size+=1;};
        };
        */
        /// TRUNK_ZZZ
        /*
        for (int i = 0; i < int(Macroeconomics->size2); i++) {
                int CurrentDate = gsl_matrix_get(Macroeconomics, 0, i);
                if (CurrentDate == FirstDate) Cond = 1;
                if (Cond && CurrentDate <= LastDate) Size += 1;
                if (CurrentDate > LastDate) break;
        };
        */
        
        
        // Populating gsl_matrix * Data
        int PastLag=120+Year; // Steps before FirstDate necessary for BinaryProjector() to start at FirstDate
        int FirstStep=0; // STL vector index where data matches FirstDate
        int VectorSize=int(LocalDates.size());
        Data = gsl_matrix_calloc (6, Size+PastLag);
        Cond=0;
        for (int i=0; i<Size; i++) {gsl_matrix_set(Data, 0, i+PastLag, gsl_matrix_get (Macroeconomics, 0, i));}; //Dates
        for (int i=0; i<Size; i++) {
            for (int k=0; k<VectorSize; k++) {
                if ((LocalDates[k]==gsl_matrix_get (Macroeconomics, 0, i)) && Cond==0) {FirstStep=k; Cond=1;};
                if (LocalDates[k]==gsl_matrix_get (Macroeconomics, 0, i)) {
                    gsl_matrix_set(Data, 1, i+PastLag, LocalOpen[k]); //Open in €
                    gsl_matrix_set(Data, 2, i+PastLag, LocalClose[k]); //Close in €
                    gsl_matrix_set(Data, 3, i+PastLag, LocalVolumes[k]); // Daily volume
                    gsl_matrix_set(Data, 4, i+PastLag, LocalHigh[k]); //High in €
                    gsl_matrix_set(Data, 5, i+PastLag, LocalLow[k]); //Low in €
                    break;
                };
            }; // closes k loop
        }; // closes i loop
        // Stiching so as to match Macroeconomics and other exchanges data
        for (int i=1; i<Size; i++) {
            if (gsl_matrix_get(Data, 1, i+PastLag)<0.0000001) {gsl_matrix_set(Data, 1, i+PastLag, gsl_matrix_get(Data, 1, i+PastLag-1));};
            if (gsl_matrix_get(Data, 2, i+PastLag)<0.0000001) {gsl_matrix_set(Data, 2, i+PastLag, gsl_matrix_get(Data, 2, i+PastLag-1));};
            if (gsl_matrix_get(Data, 3, i+PastLag)<0.0000001) {gsl_matrix_set(Data, 3, i+PastLag, gsl_matrix_get(Data, 3, i+PastLag-1));};
            if (gsl_matrix_get(Data, 4, i+PastLag)<0.0000001) {gsl_matrix_set(Data, 4, i+PastLag, gsl_matrix_get(Data, 4, i+PastLag-1));};
            if (gsl_matrix_get(Data, 5, i+PastLag)<0.0000001) {gsl_matrix_set(Data, 5, i+PastLag, gsl_matrix_get(Data, 5, i+PastLag-1));};
        };
        // Populating the region 0<t<PastLag for BinaryProjector
        if (FirstStep>PastLag) {
            for (int k=FirstStep-PastLag; k<FirstStep; k++) {
                gsl_matrix_set(Data, 0, k-FirstStep+PastLag, LocalDates[k]);
                gsl_matrix_set(Data, 1, k-FirstStep+PastLag, LocalOpen[k]);
                gsl_matrix_set(Data, 2, k-FirstStep+PastLag, LocalClose[k]);
                gsl_matrix_set(Data, 3, k-FirstStep+PastLag, LocalVolumes[k]);
                gsl_matrix_set(Data, 4, k-FirstStep+PastLag, LocalHigh[k]);
                gsl_matrix_set(Data, 5, k-FirstStep+PastLag, LocalLow[k]);
            };
        };
        
        // Stitching stock splits
        for (int t=1; t<Size+PastLag; t++) {
            double Split=1;
            if (gsl_matrix_get(Data, 2, t-1)>0.00001) {Split=gsl_matrix_get(Data, 2, t)/gsl_matrix_get(Data, 2, t-1);};
            if (((Split>=4) || (Split<=0.25)) && (Split>0)) {
                for (int i=t; i<Size+PastLag; i++) {gsl_matrix_set(Data, 2, i, gsl_matrix_get(Data, 2, i)/Split);};
            }; // if
        }; // t loop
        
        
    }; // closes Gen_trunk_old()
    
    
    void Gen_trunk (string Name, gsl_matrix * Macroeconomics, int FirstDate, int LastDate) {
        vector<double> LocalDates, LocalOpen, LocalClose, LocalHigh, LocalLow, LocalVolumes;
        ifstream FileInput(Name); string Line; int LineNb=0;
        while (getline (FileInput, Line)) {LineNb++;
            istringstream LineStream(Line); string Item; int ItemNb=0;
            while (getline (LineStream, Item, ',')) {ItemNb++;
                if ((LineNb>1) && (ItemNb==3)) {LocalDates.push_back(atof(Item.c_str()));};
                if ((LineNb>1) && (ItemNb==5)) {LocalOpen.push_back(atof(Item.c_str()));};
                if ((LineNb>1) && (ItemNb==8)) {LocalClose.push_back(atof(Item.c_str()));};
                if ((LineNb>1) && (ItemNb==9)) {LocalVolumes.push_back(atof(Item.c_str()));};
                if ((LineNb>1) && (ItemNb==6)) {LocalHigh.push_back(atof(Item.c_str()));};
                if ((LineNb>1) && (ItemNb==7)) {LocalLow.push_back(atof(Item.c_str()));};
            };
        };
        FileInput.close();
        /*
        for (int i = 0; i < 5 && i < LocalDates.size(); i++) { // Adjust '5' to print more or fewer lines
            std::cout << "Initial data in Gen_trunk for " << Name << ": Date: " << LocalDates[i] << ", Open: " << LocalOpen[i] << std::endl;
        }
        */
        // Computing the data size from FirstDate to LastDate
        int Size=0; bool Cond=0; int LastDate_index = 0;
        int Macroeconomics_FirstDate_index = 0; int LocalDates_FirstDate_index = 0;
        for (int i = 0; i < int(Macroeconomics->size2); i++) {
            int CurrentDate = gsl_matrix_get(Macroeconomics, 0, i);
            if (CurrentDate == FirstDate) {
                Cond = 1;
                Macroeconomics_FirstDate_index = i;
                //std::cout << "FirstDate found in Macroeconomics: " << CurrentDate << std::endl;
            }
            if (Cond && CurrentDate <= LastDate) {
                Size++;
                //std::cout << "CurrentDate: " << CurrentDate << std::endl;
            }
            if (CurrentDate > LastDate) {
                LastDate_index = CurrentDate-1;
                //std::cout << "Reached date beyond LastDate: " << CurrentDate << std::endl;
                break;
            }
        };
        for (int k = 0; k < int(LocalDates.size()); k++) {
            if (LocalDates[k] == FirstDate) {LocalDates_FirstDate_index = k; break;}
        };
        //std::cout << "Size calculated for date range: " << Size << std::endl;

        // Populating gsl_matrix * Data
        //int PastLag=120+Year; // Adjust 'Year' as needed for your implementation
        Data = gsl_matrix_calloc (6, Size);
        //std::cout << "Data matrix allocated with size: " << Data->size1 << "x" << Data->size2 << std::endl;
        //int FirstStep=0; // STL vector index where data matches FirstDate
        Cond=0;
        for (int i=0; i<Size; i++) {
            gsl_matrix_set(Data, 0, i, gsl_matrix_get(Macroeconomics, 0, i + Macroeconomics_FirstDate_index)); //Dates
        }

        for (int i=0; i<Size; i++) {
            gsl_matrix_set(Data, 1, i, LocalOpen[i + LocalDates_FirstDate_index]); //Open in €
            gsl_matrix_set(Data, 2, i, LocalClose[i + LocalDates_FirstDate_index]); //Close in €
            gsl_matrix_set(Data, 3, i, LocalVolumes[i + LocalDates_FirstDate_index]); // Daily volume
            gsl_matrix_set(Data, 4, i, LocalHigh[i + LocalDates_FirstDate_index]); //High in €
            gsl_matrix_set(Data, 5, i, LocalLow[i + LocalDates_FirstDate_index]); //Low in €
            //std::cout << "Setting data for date: " << LocalDates[k] << " Open: " << LocalOpen[k] << " Close: " << LocalClose[k] << std::endl;
        }

        // Stiching so as to match Macroeconomics and other exchanges data
        for (int i=1; i<Size; i++) {
            if (gsl_matrix_get(Data, 2, i) < 0.0000001) {
                gsl_matrix_set(Data, 2, i, gsl_matrix_get(Data, 2, i-1));
                //std::cout << "No data found for date: " << gsl_matrix_get(Data, 0, i+PastLag) << ", using previous day's data." << std::endl;
            }
            // Repeat this if-check for the other data points (Close, Volumes, High, Low)
        }

        // Stitching stock splits
        for (int t = 1; t < Size; t++) {
            // Make sure we do not go beyond the matrix size when accessing 't-1'
            if (t-1 >= Data->size2) {
                //std::cerr << "Index " << t-1 << " is out of bounds for the Data matrix with size " << Data->size2 << std::endl;
                continue; // Skip this iteration if out of bounds
            }

            double Split = 1;
            if (gsl_matrix_get(Data, 2, t-1) > 0.00001) {
                Split = gsl_matrix_get(Data, 2, t) / gsl_matrix_get(Data, 2, t-1);
            }
            if (((Split >= 4) || (Split <= 0.25)) && (Split > 0)) {
                // When adjusting for a split, do not go beyond the matrix size
                for (int i = t; i < Size; i++) {
                    if (i >= Data->size2) {
                        //std::cerr << "Index " << i << " is out of bounds for the Data matrix with size " << Data->size2 << std::endl;
                        break; // Exit the loop if we're about to go out of bounds
                    }
                    gsl_matrix_set(Data, 2, i, gsl_matrix_get(Data, 2, i) / Split);
                }
                //std::cout << "Adjusting for stock split on date: " << gsl_matrix_get(Data, 0, t) << " Split ratio: " << Split << std::endl;
            }
        }

        std::cout << "Final dataset size: " << Data->size2 << " for " << Name << std::endl;
        
        std::cout << "Final data in Gen_trunk for " << Name << std::endl;
        for (int i = 0; i < 5; i++) { // Print some of the final loaded data
            std::cout << std::setprecision(10) << "Date: " << gsl_matrix_get(Data, 0, i) << ", Data: " << gsl_matrix_get(Data, 2, i) << std::endl;
        }
        
    }

    
    
    
    
    // Broker fees of IG are downloaded at https://www.ig.com/fr/actions/conditions-actions
    // This outputs files LSEnew.txt, NYSEnew.txt, NASDAQnew.txt which are the stocks offered by IG. These files must be formatted to Legacy MAC OS (CR) text format, and can in turn be used as direct data feed
    void ProcessBrokerData () {
        InPF="Out";
        vector<vector<string> > W;
        string Path=Machine + "Symba/CSV/Tiered Margin_cfd.csv";
        ifstream FileInput(Path); string Line; int LineNb=0;
        vector<string> V1, V2, V3, V4;
        while (getline (FileInput, Line)) {LineNb++;
            istringstream LineStream(Line); string Item; int ItemNb=0;
            while (getline (LineStream, Item, '\t')) {ItemNb++;
                if ((ItemNb==1) && (LineNb>5)) {V1.push_back(Item.c_str());}; // Symbol
                if ((ItemNb==2) && (LineNb>5)) {V2.push_back(Item.c_str());}; // Title
                if ((ItemNb==3) && (LineNb>5)) {V3.push_back(Item.c_str());}; // Country
                if ((ItemNb==5) && (LineNb>5)) {V4.push_back(Item.c_str());}; // Can go short?
            };
        };
        // Reuters RIC: .L (London Stock Exchange), .O (NASDAQ), .N (NYSE), .P (Nyse ARCA), .PK (OTC Market Group), .OB (?), .K (New York Consolidated), .A (American Stock Exchange)
        string Name=Machine + "Symba/CSV/" + Exchange + "/" + File;
        string Suffix=".L";
        if (Exchange=="NYSE") {Suffix=".N";} // Some stocks will taken as .K on the New York Consolidated!
        else if (Exchange=="NASDAQ") {Suffix=".O";};
        W.push_back(V1); W.push_back(V2); W.push_back(V3); W.push_back(V4);
        ifstream FileInput2(Name); string Line2; LineNb=0;
        while (getline (FileInput2, Line2)) {LineNb++;
            istringstream LineStream(Line2); string Item; int ItemNb=0;
            while (getline (LineStream, Item, ',')) {ItemNb++;
                if ((LineNb==2) && (ItemNb==1)) {
                    for (int i=0; i<int(W[0].size()); i++) {
                        if (((Item.c_str()==W[0][i]) || (Item.c_str()+Suffix==W[0][i])) && (W[2][i]==Country) && (W[3][i]=="Yes")) {
                            Symbol=W[0][i]; Title=W[1][i]; InPF="In";
                        };
                    }; // closes i loop
                } // closes if
            }; // closes while
        }; // closes while
        if (InPF=="In") {
            ofstream outputLSE(Machine + "Symba/CSV/LSEnew.txt", ofstream::app);
            ofstream outputNYSE(Machine + "Symba/CSV/NYSEnew.txt", ofstream::app);
            ofstream outputNASDAQ(Machine + "Symba/CSV/NASDAQnew.txt", ofstream::app);
            if (Exchange=="LSE") {outputLSE << File << ',' << Symbol << ',' << Title << ',' << Exchange << ',' << Country << ',' << Currency << endl;}
            if (Exchange=="NYSE") {outputNYSE << File << ',' << Symbol << ',' << Title << ',' << Exchange << ',' << Country << ',' << Currency << endl;}
            if (Exchange=="NASDAQ") {outputNASDAQ << File << ',' << Symbol << ',' << Title << ',' << Exchange << ',' << Country << ',' << Currency << endl;}
        };
        
    }; // closes ProcessBrokerData()
    
    
}; // closes Share class


 
vector<Share> PortfolioGenerator (string Name, int FirstDate) {
    string ProcessedData="Yes";
    vector<Share> PF, PFFinal; int TickStart=int(time(0));
    string DPath=Machine + "Symba/CSV/" + Name + ".txt";
    ifstream DFileInput(DPath);string DLine; vector<string> VExchanges, VCountries, VCurrencies;
    while (getline (DFileInput, DLine, '\r')) {
        istringstream LineStream(DLine); string DItem; int DItemNb=0;
        while (getline (LineStream, DItem, ',')) {DItemNb++;
            if (DItemNb==1) {VExchanges.push_back(DItem);}
            else if (DItemNb==2) {VCountries.push_back(DItem);}
            else if (DItemNb==3) {VCurrencies.push_back(DItem);};
        };
    };
    int ExchangeSize=int(VExchanges.size()); DFileInput.close();
    gsl_matrix * Macro = Macroeconomics();
    for (int m=0; m<ExchangeSize; m++) {
        string FPath=Machine + "Symba/CSV/" + VExchanges[m] + "new.txt";
        vector<string> VFiles, VSymbol, VTitles;
        if (ProcessedData=="No") {FPath=Machine + "Symba/CSV/" + VExchanges[m] + ".txt";};
        ifstream FFileInput(FPath); string FLine;
        if (ProcessedData=="No") {
            while (getline (FFileInput, FLine, '\r')) {VFiles.push_back(FLine);};
        }
        else if (ProcessedData=="Yes") {
            while (getline (FFileInput, FLine, '\r')) {
                istringstream LineStream(FLine); string DItem; int DItemNb=0;
                while (getline (LineStream, DItem, ',')) {DItemNb++;
                    if (DItemNb==1) {VFiles.push_back(DItem);}
                    else if (DItemNb==2) {VSymbol.push_back(DItem);}
                    else if (DItemNb==3) {VTitles.push_back(DItem);};
                };
            };
        };
        FFileInput.close();
        
        int StockSize=int(VFiles.size());
        for (int s=0; s<StockSize; s++) {
            Share X; X.Exchange=VExchanges[m]; X.File=VFiles[s]; X.Country=VCountries[m]; X.Currency=VCurrencies[m];
            if (ProcessedData=="Yes") {X.Symbol=VSymbol[s]; X.Title=VTitles[s];}
            else if (ProcessedData=="No") {X.ProcessBrokerData();};
            //cout << X.InPF << endl;
            if ((X.InPF=="In") && (ProcessedData=="No")) {
                X.Gen(Machine + "Symba/CSV/" + VExchanges[m] + "/" + VFiles[s], Macro, FirstDate);
                PF.push_back(X);
                cout << VExchanges[m] << " (" << floor(s*100.0/StockSize) << "%): " << X.Symbol << ", " << X.Title << " (" << X.Country << ", " << X.Currency << ")" << endl;
            };
            if (ProcessedData=="Yes") {
                X.Gen(Machine + "Symba/CSV/" + VExchanges[m] + "/" + VFiles[s], Macro, FirstDate);
                PF.push_back(X);
                cout << VExchanges[m] << " (" << floor(s*100.0/StockSize) << "%): " << X.Symbol << ", " << X.Title << " (" << X.Country << ", " << X.Currency << ")" << endl;
            };
            
        }; // closes s loop
    }; // closes m loop
    
    // Considering older shares since Jan 2007 only
    for (int s=0; s<int(PF.size()); s++) {
        int Count=0;
        for (int t=0; t<int(PF[s].Data->size2); t++) {
            if (gsl_matrix_get(PF[s].Data, 2, t)<0.0000001) {Count+=1;};
        };
        if (Count<Month) {PFFinal.push_back(PF[s]);};
    };
    int PFSize=int(PFFinal.size());
    //int PFSize=int(PF.size());
    cout << "****************" << endl << "PF.size()=" << PFSize << endl << "Loading time: " << int(time(0)) - TickStart << "s" << endl << "****************" << endl << endl;
    
    return PFFinal;
}; // closes PortfolioGenerator()



vector<Share> PortfolioGenerator_trunk (string Name, int FirstDate, int LastDate) {
    string ProcessedData="Yes";
    vector<Share> PF, PFFinal; int TickStart=int(time(0));
    string DPath=Machine + "Symba/CSV/" + Name + ".txt";
    ifstream DFileInput(DPath);string DLine; vector<string> VExchanges, VCountries, VCurrencies;
    while (getline (DFileInput, DLine, '\r')) {
        istringstream LineStream(DLine); string DItem; int DItemNb=0;
        while (getline (LineStream, DItem, ',')) {DItemNb++;
            if (DItemNb==1) {VExchanges.push_back(DItem);}
            else if (DItemNb==2) {VCountries.push_back(DItem);}
            else if (DItemNb==3) {VCurrencies.push_back(DItem);};
        };
    };
    int ExchangeSize=int(VExchanges.size()); DFileInput.close();
    gsl_matrix * Macro = Macroeconomics();
    for (int m=0; m<ExchangeSize; m++) {
        string FPath=Machine + "Symba/CSV/" + VExchanges[m] + "new.txt";
        vector<string> VFiles, VSymbol, VTitles;
        if (ProcessedData=="No") {FPath=Machine + "Symba/CSV/" + VExchanges[m] + ".txt";};
        ifstream FFileInput(FPath); string FLine;
        if (ProcessedData=="No") {
            while (getline (FFileInput, FLine, '\r')) {VFiles.push_back(FLine);};
        }
        else if (ProcessedData=="Yes") {
            while (getline (FFileInput, FLine, '\r')) {
                istringstream LineStream(FLine); string DItem; int DItemNb=0;
                while (getline (LineStream, DItem, ',')) {DItemNb++;
                    if (DItemNb==1) {VFiles.push_back(DItem);}
                    else if (DItemNb==2) {VSymbol.push_back(DItem);}
                    else if (DItemNb==3) {VTitles.push_back(DItem);};
                };
            };
        };
        FFileInput.close();
        
        int StockSize=int(VFiles.size());
        for (int s=0; s<StockSize; s++) {
            Share X; X.Exchange=VExchanges[m]; X.File=VFiles[s]; X.Country=VCountries[m]; X.Currency=VCurrencies[m];
            if (ProcessedData=="Yes") {X.Symbol=VSymbol[s]; X.Title=VTitles[s];}
            else if (ProcessedData=="No") {X.ProcessBrokerData();};
            //cout << X.InPF << endl;
            if ((X.InPF=="In") && (ProcessedData=="No")) {
                X.Gen_trunk(Machine + "Symba/CSV/" + VExchanges[m] + "/" + VFiles[s], Macro, FirstDate, LastDate);
                PF.push_back(X);
                cout << VExchanges[m] << " (" << floor(s*100.0/StockSize) << "%): " << X.Symbol << ", " << X.Title << " (" << X.Country << ", " << X.Currency << ")" << endl;
            };
            if (ProcessedData=="Yes") {
                X.Gen_trunk(Machine + "Symba/CSV/" + VExchanges[m] + "/" + VFiles[s], Macro, FirstDate, LastDate);
                PF.push_back(X);
                cout << VExchanges[m] << " (" << floor(s*100.0/StockSize) << "%): " << X.Symbol << ", " << X.Title << " (" << X.Country << ", " << X.Currency << ")" << endl;
            };
            
        }; // closes s loop
    }; // closes m loop
    
    // Considering older shares since FirstDate only
    for (int s=0; s<int(PF.size()); s++) {
        int Count=0;
        for (int t=0; t<int(PF[s].Data->size2); t++) {
            if (gsl_matrix_get(PF[s].Data, 2, t)<0.0000001) {Count+=1;};
        };
        if (Count<Month) {PFFinal.push_back(PF[s]);};
    };
    int PFSize=int(PFFinal.size());
    //int PFSize=int(PF.size());
    cout << "****************" << endl << "PF.size()=" << PFSize << endl << "Loading time: " << int(time(0)) - TickStart << "s" << endl << "****************" << endl << endl;
    
    return PFFinal;
}; // closes PortfolioGenerator()



vector<Share> PortfolioGenerator_trunk_old (string Name, int FirstDate, int LastDate) {
    string ProcessedData="Yes";
    vector<Share> PF, PFFinal; int TickStart=int(time(0));
    string DPath=Machine + "Symba/CSV/" + Name + ".txt";
    ifstream DFileInput(DPath);string DLine; vector<string> VExchanges, VCountries, VCurrencies;
    while (getline (DFileInput, DLine, '\r')) {
        istringstream LineStream(DLine); string DItem; int DItemNb=0;
        while (getline (LineStream, DItem, ',')) {DItemNb++;
            if (DItemNb==1) {VExchanges.push_back(DItem);}
            else if (DItemNb==2) {VCountries.push_back(DItem);}
            else if (DItemNb==3) {VCurrencies.push_back(DItem);};
        };
    };
    int ExchangeSize=int(VExchanges.size()); DFileInput.close();
    gsl_matrix * Macro = Macroeconomics();
    for (int m=0; m<ExchangeSize; m++) {
        string FPath=Machine + "Symba/CSV/" + VExchanges[m] + "new.txt";
        vector<string> VFiles, VSymbol, VTitles;
        if (ProcessedData=="No") {FPath=Machine + "Symba/CSV/" + VExchanges[m] + ".txt";};
        ifstream FFileInput(FPath); string FLine;
        if (ProcessedData=="No") {
            while (getline (FFileInput, FLine, '\r')) {VFiles.push_back(FLine);};
        }
        else if (ProcessedData=="Yes") {
            while (getline (FFileInput, FLine, '\r')) {
                istringstream LineStream(FLine); string DItem; int DItemNb=0;
                while (getline (LineStream, DItem, ',')) {DItemNb++;
                    if (DItemNb==1) {VFiles.push_back(DItem);}
                    else if (DItemNb==2) {VSymbol.push_back(DItem);}
                    else if (DItemNb==3) {VTitles.push_back(DItem);};
                };
            };
        };
        FFileInput.close();
        
        int StockSize=int(VFiles.size());
        
        // TRUNK_ZZZ
        for (int s=0; s<StockSize; s++) {
            Share X; X.Exchange=VExchanges[m]; X.File=VFiles[s]; X.Country=VCountries[m]; X.Currency=VCurrencies[m];
            if (ProcessedData=="Yes") {X.Symbol=VSymbol[s]; X.Title=VTitles[s];}
            else if (ProcessedData=="No") {X.ProcessBrokerData();};
            //cout << X.InPF << endl;
            if ((X.InPF=="In") && (ProcessedData=="No")) {
                X.Gen_trunk(Machine + "Symba/CSV/" + VExchanges[m] + "/" + VFiles[s], Macro, FirstDate, LastDate);
                PF.push_back(X);
                cout << VExchanges[m] << " (" << floor(s*100.0/StockSize) << "%): " << X.Symbol << ", " << X.Title << " (" << X.Country << ", " << X.Currency << ")" << endl;
            };
            if (ProcessedData=="Yes") {
                X.Gen_trunk(Machine + "Symba/CSV/" + VExchanges[m] + "/" + VFiles[s], Macro, FirstDate, LastDate);
                PF.push_back(X);
                cout << VExchanges[m] << " (" << floor(s*100.0/StockSize) << "%): " << X.Symbol << ", " << X.Title << " (" << X.Country << ", " << X.Currency << ")" << endl;
            };
        }; // closes s loop
        
        /*
        // Considering older shares since Jan 2007 only
        for (int s=0; s<int(PF.size()); s++) {
            int Count=0;
            for (int t=0; t<int(PF[s].Data->size2); t++) {
                if (gsl_matrix_get(PF[s].Data, 2, t)<0.0000001) {Count+=1;};
            };
            if (Count<Month) {PFFinal.push_back(PF[s]);};
        };
        */
    }; // closes m loop
    
    std::cout << "PF.size() after filtering: " << PF.size() << std::endl;


    
    int PFSize = int(PF.size());
    cout << "****************" << endl << "PF.size()=" << PFSize << endl << "Loading time: " << int(time(0)) - TickStart << "s" << endl << "****************" << endl << endl;

    return PFFinal;

    //return PF;
}; // closes PortfolioGenerator_trunk()




gsl_matrix * PortfolioOutput (vector<Share> PF, int Size) {
    if (Size==0) {Size=int(PF.size());};
    if (PF.empty()) {
        std::cerr << "Error: PF vector is empty. No data to process." << std::endl;
        // Handle the error appropriately, perhaps return an empty gsl_matrix or nullptr
        return nullptr;
    }
    string A = Machine+"PF.csv"; ofstream outputSpecs(A.c_str(), ofstream::app);
    vector<int> J=Shuffle(int(PF.size()));
    gsl_matrix * Result = gsl_matrix_calloc (Size+1, PF[0].Data->size2-381);
    outputSpecs << fixed << "Date (yyyymmdd),"; for (int i=1; i<Size+1; i++) {outputSpecs << PF[J[i-1]].Symbol << " (" << PF[J[i-1]].Exchange << "),";}; outputSpecs << endl;
    //outputSpecs << fixed << "Date,"; for (int i=1; i<Size+1; i++) {outputSpecs << PF[J[i-1]].Exchange <<",";}; outputSpecs << endl;
    outputSpecs.close();
    for (int t=0; t<PF[0].Data->size2-381; t++) {gsl_matrix_set(Result, 0, t, gsl_matrix_get(PF[0].Data, 0, t+381));}; // Dates
    for (int i=1; i<Size+1; i++) {
        for (int t=0; t<PF[0].Data->size2-381; t++) {
            gsl_matrix_set(Result, i, t, gsl_matrix_get(PF[J[i-1]].Data, 2, t+381)); // Close prices
        };
    };
    PlotGSLMatrix(Result, "PF.csv", 1);
    return Result;
};




//POSTDOC2022a
vector<string> getFileNames(const string& directory) {
    vector<string> files;
    DIR* dir;
    struct dirent* ent;
    if ((dir = opendir (directory.c_str())) != NULL) {
        while ((ent = readdir (dir)) != NULL) {
            files.push_back(ent->d_name);
        }
        closedir (dir);
    } else {
        perror ("");
    }
    return files;
}

vector<Share> readStockData(const string& directory) {
    vector<Share> data;
    vector<string> files = getFileNames(directory);
    
    for (const auto& file : files) {
        ifstream fin(directory + "/" + file);
        string line, temp;
        if (getline(fin, line)) { // Skip header line
            Share share;
            while (getline(fin, line)) {
                stringstream s(line);
                getline(s, share.Symbol, ';');
                getline(s, temp, ';'); // Skip interval
                getline(s, share.Title, ';');
                getline(s, share.temp_date, ';');
                getline(s, temp, ';'); // Skip open
                getline(s, temp, ';'); // Skip high
                getline(s, temp, ';'); // Skip low
                getline(s, share.temp_price, ';');
                getline(s, share.temp_volume, ';');
                replace(share.temp_price.begin(), share.temp_price.end(), ',', '.');
                replace(share.volumes.begin(), share.volumes.end(), ',', '.');
                //replace(share.dates.begin(), share.dates.end(), '/', ' ');
                //replace(share.dates.begin(), share.dates.end(), ':', ' ');
                share.close_prices.push_back(atof(share.temp_price.c_str()));
                share.volumes.push_back(atof(share.temp_volume.c_str()));
                share.dates.push_back(share.temp_date);
                /*
                double close;
                s >> close;
                share.close_prices.push_back(close);
                //if (share.Symbol=="TEM.L") {cout << "share.Symbol " << share.Symbol << ", price " << close << endl;};
                double volume;
                if (s >> volume) {
                    share.volumes.push_back(volume);
                }
                else {
                    // If 'volume' is not available for a given line,
                    // append a default value to the 'volumes' vector.
                    share.volumes.push_back(0);
                }
                */
            }
            share.File = file;
            int size = int(share.close_prices.size());
            if (size > 0) {
                share.Data = gsl_matrix_alloc(6, size);
                for(int t = 0; t < size; t++) {
                    gsl_matrix_set(share.Data, 2, t, share.close_prices[t]);
                    gsl_matrix_set(share.Data, 3, t, share.volumes[t]);
                    //if (share.Symbol=="TEM.L") {cout << share.Symbol << ", price " << gsl_matrix_get(share.Data, i, 0) << ", volume " << gsl_matrix_get(share.Data, i, 1) << endl;};
                }
                share.Exchange="LSE";
                share.Country="UK";
                share.Currency="GBP";
                share.Sector="All";
                share.InPF="In";
                
                /*
                 vector<double> close_prices;
                 vector<double> volumes;
                 vector<string> dates;
                 gsl_matrix * Data; // Prices & cashflows
                */
                
                if (share.close_prices.size()>=1010) { // If enough points to fill 27/09/2018 to 27/09/2022
                    share.close_prices.erase(share.close_prices.begin(), share.close_prices.begin() + share.close_prices.size()-1010);
                    share.volumes.erase(share.volumes.begin(), share.volumes.begin() + share.volumes.size()-1010);
                    share.dates.erase(share.dates.begin(), share.dates.begin() + share.dates.size()-1010);
                    gsl_matrix_free(share.Data);
                    share.Data = gsl_matrix_alloc(6, 1010);
                    for(int t = 0; t < 1010; t++) {
                        gsl_matrix_set(share.Data, 2, t, share.close_prices[t]);
                        gsl_matrix_set(share.Data, 3, t, share.volumes[t]);
                    };
                    cout << share.dates[0] << " : ";
                    if (share.dates[0]=="27/09/2018 00:00:00") {data.push_back(share); cout << share.Symbol << " InPF!" << endl;};
                };
                
                //if (share.close_prices.size()==2527) {data.push_back(share);}; // Keep only stocks that are continuously trading from 27/09/2012 to 27/09/2022
                //data.push_back(share);
            }
        }
        fin.close();
        
    }
    //cout << "XXX = " << data[0].Symbol << ", " << data[0].Title << ", " << data[0].close_prices[5] << ", " << data[0].volumes[5] << ", "  << gsl_matrix_get(data[0].Data, 5, 0) << ", " << gsl_matrix_get(data[0].Data, 5, 1) << endl;
    return data;
}


gsl_matrix* PortfolioOutput2(vector<Share> PF, int Size) {
    int sizePF = int(PF.size());
    
    /*
    gsl_matrix* Mpricesizes = gsl_matrix_calloc(1, sizePF);
    for (int i = 0; i < sizePF ; i++) {
        gsl_matrix_set(Mpricesizes, 0, i, int(PF[i].close_prices.size()));
    };
    PlotGSLMatrix(Mpricesizes, "Mpricesizes.csv", 1);
    exit(0); // 431 shares have a size of 2527 trading days
    */
    
    if (Size == 0) {Size = sizePF;}
    string A = Machine+"PF2.csv";
    ofstream outputSpecs(A.c_str(), ofstream::app);
    vector<int> J = Shuffle2(int(PF.size()));
    gsl_matrix* Result = gsl_matrix_calloc(2*(Size + 1)+1, 4000);
    
    outputSpecs << "Date (yyyymmdd), ";
    for (int i = 0; i < Size; i++) {
        outputSpecs << PF[J[i]].Symbol << " (price), " << PF[J[i]].Symbol << " (volume), ";
        //cout << "i=" << i << ", J[i - 1]=" << J[i - 1] << ", PF[J[i - 1]].Symbol=" << PF[J[i - 1]].Symbol << endl;
    }
    outputSpecs << endl;
    
    
    //cout << "PF[10].Data->size1=" << PF[10].Data->size1 << "PF[10].Data->size2=" << PF[10].Data->size2 << endl;
    //int sizedates=2527;
    int sizedates=int(PF[J[0]].dates.size());
    for (int t = 0; t < sizedates ; t++) {
        //if (t<sizedates) {outputSpecs << PF[J[0]].dates[t] << ", ";};
        outputSpecs << PF[J[0]].dates[t] << ", ";
        for (int i = 0; i < Size; i++) {
            outputSpecs << PF[J[i]].close_prices[t] << ", " << PF[J[i]].volumes[t] << ", ";
            /*
            if ((t<int((PF[J[i-1]].close_prices).size())) && (t<int((PF[J[i-1]].volumes).size()))) {
                outputSpecs << PF[J[i-1]].close_prices[t] << ", " << PF[J[i-1]].volumes[t] << ", ";
            }
            else {outputSpecs << ", " << ", ";};
            */
            //gsl_matrix_set(Result, i, t, PF[J[i-1]].close_prices[t]); // Close prices
            //gsl_matrix_set(Result, i+1, t, PF[J[i-1]].volumes[t]); // Volumes
            //cout << "i=" << i << ", t=" << t << ", gsl_matrix_get(PF[J[i-1]].Data, 1, t)=" << gsl_matrix_get(PF[J[i-1]].Data, 1, t) << endl;
        }
        outputSpecs << endl;
    }
    
    //PlotGSLMatrix(Result, "PF2.csv", 1); // This line is commented out because it depends on your implementation of PlotGSLMatrix.
     
    
    outputSpecs.close();
    return Result;
}



class Coin {
public:
    string Symbol, Title, Exchange, File, Country, Currency, Sector, InPF;
    string temp_price, temp_volume, temp_date;
    vector<double> close_prices_BTC;
    vector<double> volumes_BTC;
    vector<double> close_prices_BUSD;
    vector<double> volumes_BUSD;
    vector<double> close_prices_ETH;
    vector<double> volumes_ETH;
    vector<double> close_prices_EUR;
    vector<double> volumes_EUR;
    vector<double> close_prices_GBP;
    vector<double> volumes_GBP;
    vector<double> close_prices_USDT;
    vector<double> volumes_USDT;
    vector<string> dates;
    gsl_matrix * Data; // Prices & cashflows
};



vector<Coin> readStockData_coin(const string& directory) {
    vector<Coin> data;
    vector<string> files = getFileNames(directory);
    
    for (const auto& file : files) {
        ifstream fin(directory + "/" + file);
        string line, temp;
        if (getline(fin, line)) { // Skip header line
            Coin coin;
            coin.Symbol = file.substr(0, file.find('.'));
            coin.Title = coin.Symbol; // For now...
            coin.File = file;
            coin.Country = "UK"; // For now...
            coin.Currency = "GBP"; // For now...
            coin.Exchange = "BINANCE";
            coin.Sector="All";
            coin.InPF="In";
            while (getline(fin, line)) {
                stringstream s(line);
                
                getline(s, coin.temp_date, ',');
                coin.dates.push_back(coin.temp_date);
                
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                
                getline(s, coin.temp_price, ',');
                //replace(coin.temp_price.begin(), coin.temp_price.end(), ',', '.');
                coin.close_prices_BTC.push_back(atof(coin.temp_price.c_str()));
                getline(s, coin.temp_volume, ',');
                //replace(coin.temp_volume.begin(), coin.temp_volume.end(), ',', '.');
                coin.volumes_BTC.push_back(atof(coin.temp_volume.c_str()));

                getline(s, coin.temp_price, ',');
                //replace(coin.temp_price.begin(), coin.temp_price.end(), ',', '.');
                coin.close_prices_BUSD.push_back(atof(coin.temp_price.c_str()));
                getline(s, coin.temp_volume, ',');
                //replace(coin.temp_volume.begin(), coin.temp_volume.end(), ',', '.');
                coin.volumes_BUSD.push_back(atof(coin.temp_volume.c_str()));

                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                
                getline(s, coin.temp_price, ',');
                //replace(coin.temp_price.begin(), coin.temp_price.end(), ',', '.');
                coin.close_prices_ETH.push_back(atof(coin.temp_price.c_str()));
                getline(s, coin.temp_volume, ',');
                //replace(coin.temp_volume.begin(), coin.temp_volume.end(), ',', '.');
                coin.volumes_ETH.push_back(atof(coin.temp_volume.c_str()));

                getline(s, coin.temp_price, ',');
                //replace(coin.temp_price.begin(), coin.temp_price.end(), ',', '.');
                coin.close_prices_EUR.push_back(atof(coin.temp_price.c_str()));
                getline(s, coin.temp_volume, ',');
                //replace(coin.temp_volume.begin(), coin.temp_volume.end(), ',', '.');
                coin.volumes_EUR.push_back(atof(coin.temp_volume.c_str()));

                getline(s, coin.temp_price, ',');
                //replace(coin.temp_price.begin(), coin.temp_price.end(), ',', '.');
                coin.close_prices_GBP.push_back(atof(coin.temp_price.c_str()));
                getline(s, coin.temp_volume, ',');
                //replace(coin.temp_volume.begin(), coin.temp_volume.end(), ',', '.');
                coin.volumes_GBP.push_back(atof(coin.temp_volume.c_str()));

                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval

                getline(s, coin.temp_price, ',');
                //replace(coin.temp_price.begin(), coin.temp_price.end(), ',', '.');
                coin.close_prices_USDT.push_back(atof(coin.temp_price.c_str()));
                getline(s, coin.temp_volume, ',');
                //replace(coin.temp_volume.begin(), coin.temp_volume.end(), ',', '.');
                coin.volumes_USDT.push_back(atof(coin.temp_volume.c_str()));

                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                getline(s, temp, ','); // Skip interval
                
            }

            int size = 1462; // From 27/09/2018 to 27/09/2022
            //int size = 1097; // From 27/09/2019 to 27/09/2022
            if (int(coin.close_prices_BTC.size()) > 0) {
                
                int trunk_start = 440; // BINANCE: 440 time steps to 27/09/2018 (excluded), 805 to 27/09/2019 (excluded)
                coin.dates.erase(coin.dates.begin(), coin.dates.begin() + trunk_start);
                coin.close_prices_BTC.erase(coin.close_prices_BTC.begin(), coin.close_prices_BTC.begin() + trunk_start);
                coin.volumes_BTC.erase(coin.volumes_BTC.begin(), coin.volumes_BTC.begin() + trunk_start);
                coin.close_prices_BUSD.erase(coin.close_prices_BUSD.begin(), coin.close_prices_BUSD.begin() + trunk_start);
                coin.volumes_BUSD.erase(coin.volumes_BUSD.begin(), coin.volumes_BUSD.begin() + trunk_start);
                coin.close_prices_ETH.erase(coin.close_prices_ETH.begin(), coin.close_prices_ETH.begin() + trunk_start);
                coin.volumes_ETH.erase(coin.volumes_ETH.begin(), coin.volumes_ETH.begin() + trunk_start);
                coin.close_prices_EUR.erase(coin.close_prices_EUR.begin(), coin.close_prices_EUR.begin() + trunk_start);
                coin.volumes_EUR.erase(coin.volumes_EUR.begin(), coin.volumes_EUR.begin() + trunk_start);
                coin.close_prices_GBP.erase(coin.close_prices_GBP.begin(), coin.close_prices_GBP.begin() + trunk_start);
                coin.volumes_GBP.erase(coin.volumes_GBP.begin(), coin.volumes_GBP.begin() + trunk_start);
                coin.close_prices_USDT.erase(coin.close_prices_USDT.begin(), coin.close_prices_USDT.begin() + trunk_start);
                coin.volumes_USDT.erase(coin.volumes_USDT.begin(), coin.volumes_USDT.begin() + trunk_start);

                int trunk_end = 235; // BINANCE: 235 time steps from 27/09/2022 (excluded)
                coin.dates.erase(coin.dates.end() - trunk_end, coin.dates.end());
                coin.close_prices_BTC.erase(coin.close_prices_BTC.end() - trunk_end, coin.close_prices_BTC.end());
                coin.volumes_BTC.erase(coin.volumes_BTC.end() - trunk_end, coin.volumes_BTC.end());
                coin.close_prices_BUSD.erase(coin.close_prices_BUSD.end() - trunk_end, coin.close_prices_BUSD.end());
                coin.volumes_BUSD.erase(coin.volumes_BUSD.end() - trunk_end, coin.volumes_BUSD.end());
                coin.close_prices_ETH.erase(coin.close_prices_ETH.end() - trunk_end, coin.close_prices_ETH.end());
                coin.volumes_ETH.erase(coin.volumes_ETH.end() - trunk_end, coin.volumes_ETH.end());
                coin.close_prices_EUR.erase(coin.close_prices_EUR.end() - trunk_end, coin.close_prices_EUR.end());
                coin.volumes_EUR.erase(coin.volumes_EUR.end() - trunk_end, coin.volumes_EUR.end());
                coin.close_prices_GBP.erase(coin.close_prices_GBP.end() - trunk_end, coin.close_prices_GBP.end());
                coin.volumes_GBP.erase(coin.volumes_GBP.end() - trunk_end, coin.volumes_GBP.end());
                coin.close_prices_USDT.erase(coin.close_prices_USDT.end() - trunk_end, coin.close_prices_USDT.end());
                coin.volumes_USDT.erase(coin.volumes_USDT.end() - trunk_end, coin.volumes_USDT.end());

                coin.Data = gsl_matrix_alloc(6, size);
                for(int t = 0; t < size; t++) {
                    if (coin.close_prices_BTC[t]<0) {coin.InPF="Out";};
                    //if (coin.close_prices_USDT[t]<0) {coin.InPF="Out";}; //FFFKKK
                    gsl_matrix_set(coin.Data, 0, t, coin.close_prices_USDT[t]);
                    gsl_matrix_set(coin.Data, 1, t, coin.volumes_USDT[t]);
                    
                    gsl_matrix_set(coin.Data, 2, t, coin.close_prices_BTC[t]);
                    gsl_matrix_set(coin.Data, 3, t, coin.volumes_BTC[t]);
                    
                    gsl_matrix_set(coin.Data, 4, t, coin.close_prices_EUR[t]);
                    gsl_matrix_set(coin.Data, 5, t, coin.volumes_EUR[t]);
                }

                if (coin.InPF=="In") {data.push_back(coin);}; // Keep only stocks that are continuously traded
                //data.push_back(coin);
            }
        }
        fin.close();
        
    }
    //cout << "XXX = " << data[0].Symbol << ", " << data[0].Title << ", " << data[0].close_prices[5] << ", " << data[0].volumes[5] << ", "  << gsl_matrix_get(data[0].Data, 5, 0) << ", " << gsl_matrix_get(data[0].Data, 5, 1) << endl;
    return data;
}


gsl_matrix* PortfolioOutput_coin(vector<Coin> PF, int Size) {
    int sizePF = int(PF.size());
    
    /*
    gsl_matrix* Mpricesizes = gsl_matrix_calloc(1, sizePF);
    for (int i = 0; i < sizePF ; i++) {
        gsl_matrix_set(Mpricesizes, 0, i, int(PF[i].close_prices_BTC.size()));
    };
    PlotGSLMatrix(Mpricesizes, "Mpricesizes_coin.csv", 1);
    //exit(0);
    */
    
    if (Size == 0) {Size = sizePF;}
    string A = Machine+"PF_coin.csv";
    ofstream outputSpecs(A.c_str(), ofstream::app);
    vector<int> J = Shuffle2(int(PF.size()));
    gsl_matrix* Result = gsl_matrix_calloc(2*(Size + 1)+1, 4000);
    
    outputSpecs << "Date (yyyymmdd), ";
    for (int i = 0; i < Size; i++) {
        cout << "i=" << i << ", J[i]=" << J[i] << ", PF[J[i]].Symbol=" << PF[J[i]].Symbol << endl;
        outputSpecs << PF[J[i]].Symbol << " (price), " << PF[J[i]].Symbol << " (volume), ";
    }
    outputSpecs << endl;
    
    //int sizedates=2527;
    int sizedates=int(PF[J[0]].dates.size());
    for (int t = 0; t < sizedates ; t++) {
        outputSpecs << PF[J[0]].dates[t] << ", ";
        for (int i = 0; i < Size; i++) {
            outputSpecs << PF[J[i]].close_prices_BTC[t] << ", " << PF[J[i]].volumes_BTC[t] << ", ";
            //outputSpecs << PF[J[i]].close_prices_USDT[t] << ", " << PF[J[i]].volumes_USDT[t] << ", "; //FFFKKK
        }
        outputSpecs << endl;
    }
    
    outputSpecs.close();
    return Result;
}







vector<vector<Share> > PFTrainingTesting (vector<Share> PF, int TrainSize) {
    int Size=int(PF.size());
    vector<int> J=Shuffle(Size);
    vector<vector<Share> > PFVec;
    vector<Share> PFTraining, PFTesting;
    for (int k=0; k<TrainSize; k++) {PFTraining.push_back(PF[J[k]]);}; PFVec.push_back(PFTraining);
    for (int k=TrainSize; k<Size; k++) {PFTesting.push_back(PF[J[k]]);}; PFVec.push_back(PFTesting);
    return PFVec;
};



vector<vector<Coin> > PFTrainingTesting_coin (vector<Coin> PF, int TrainSize) {
    int Size=int(PF.size());
    vector<int> J=Shuffle(Size);
    vector<vector<Coin> > PFVec;
    vector<Coin> PFTraining, PFTesting;
    for (int k=0; k<TrainSize; k++) {PFTraining.push_back(PF[J[k]]);}; PFVec.push_back(PFTraining);
    for (int k=TrainSize; k<Size; k++) {PFTesting.push_back(PF[J[k]]);}; PFVec.push_back(PFTesting);
    return PFVec;
};


vector<vector<gsl_matrix *> > PortfolioMultiOutput (vector<Share> PF, int Size) {
    if (Size==0) {Size=int(PF.size());};
    string A = Machine+"PF_PortfolioMultiOutput.csv"; ofstream outputSpecs(A.c_str(), ofstream::app);
    vector<int> J=Shuffle(int(PF.size()));
    gsl_matrix * Result = gsl_matrix_calloc (Size+1, PF[0].Data->size2-381);
    gsl_matrix * Result2 = gsl_matrix_calloc (Size+1, PF[0].Data->size2-381);
    outputSpecs << fixed << "Date (yyyymmdd),"; for (int i=1; i<Size+1; i++) {outputSpecs << PF[J[i-1]].Symbol << " (" << PF[J[i-1]].Exchange << "),";}; outputSpecs << endl;
    //outputSpecs << fixed << "Date,"; for (int i=1; i<Size+1; i++) {outputSpecs << PF[J[i-1]].Exchange <<",";}; outputSpecs << endl;
    outputSpecs.close();
    for (int t=0; t<PF[0].Data->size2-381; t++) {// Dates
        gsl_matrix_set(Result, 0, t, gsl_matrix_get(PF[0].Data, 0, t+381));
        gsl_matrix_set(Result2, 0, t, gsl_matrix_get(PF[0].Data, 0, t+381));
    };
    for (int i=1; i<Size+1; i++) {
        for (int t=0; t<PF[0].Data->size2-381; t++) {
            gsl_matrix_set(Result, i, t, gsl_matrix_get(PF[J[i-1]].Data, 2, t+381)); // Close prices
            gsl_matrix_set(Result2, i, t, gsl_matrix_get(PF[J[i-1]].Data, 3, t+381)); // Volumes
        };
        //cout << PF[J[i-1]].Exchange << ": " << PF[J[i-1]].Symbol << ", " << PF[J[i-1]].Title << " selected" << endl;
    };
    PlotGSLMatrix(Result, "PF_PortfolioMultiOutput.csv", 1);
    //outputMoments << "Log-return" << '\t' << "AC-1w Log-return" << '\t' << "AC-2w Log-return" << '\t' << "AC-m Log-return" << '\t' << "AC-3m Log-return" << '\t' << "AC-6m Log-return" << '\t' << "AC-y Log-return" << '\t' << "Abs-log-return" << '\t' << "AC-1w Abs-log-return" << '\t' << "AC-2w Abs-log-return" << '\t' << "AC-m Abs-log-return" << '\t' << "AC-3m Abs-log-return" << '\t' << "AC-6m Abs-log-return" << '\t' << "AC-y Abs-log-return" << '\t' << "w-Volatility" << '\t' << "AC-w w-Volatility" << '\t' << "2w-Volatility" << '\t' << "AC-2w 2w-Volatility" << '\t' << "m-Volatility" << '\t' << "AC-m m-Volatility" << '\t' << "3m-Volatility" << '\t' << "AC-3m 3m-Volatility" << '\t' << "6m-Volatility" << '\t' << "AC-6m 6m-Volatility" << '\t' << "y-Volatility" << '\t' << "AC-y y-Volatility" << '\t' << "Volumes(bsp)" << '\t' << "AC-1w Volume" << '\t' << "AC-2w Volume" << '\t' << "AC-m Volume" << '\t' << "AC-3m Volume" << '\t' << "AC-6m Volume" << '\t' << "AC-y Volume" << endl;
    vector<vector<gsl_matrix *> > FinalResult;
    int T=int(PF[0].Data->size2-381);
    for (int m=1; m<Size+1; m++) {
        vector<gsl_matrix *> Temp;
        gsl_matrix * Moments = gsl_matrix_calloc (64, T);
        gsl_matrix * Prices = gsl_matrix_calloc (1, T); for (int t=0; t<T; t++) {gsl_matrix_set (Prices, 0, t, gsl_matrix_get (Result, m, t));};
        //PlotGSLMatrix(Prices, "Prices.csv", 1);
        for (int t=1; t<T; t++) {
            if (gsl_matrix_get(Result, m, t-1)!=0) {gsl_matrix_set(Moments, 0, t, log(gsl_matrix_get(Result, m, t)/gsl_matrix_get(Result, m, t-1)));} // Log-return
            else {gsl_matrix_set(Moments, 0, t, 0);}; // Log-return failsafe
            gsl_matrix_set(Moments, 7, t, abs(gsl_matrix_get(Moments, 0, t))); // Abs log-return
            //if (gsl_matrix_get(Result2, m, t-1)!=0) {gsl_matrix_set(Moments, 26, t, gsl_matrix_get(Result2, m, t)/gsl_matrix_get(Result2, m, t-1));} // Volumes-return
            //else {gsl_matrix_set(Moments, 26, t, 1);}; // Volumes-return failsafe
            gsl_matrix_set(Moments, 26, t, gsl_matrix_get(Result2, m, t)); // Volumes
        };
        double Crash=0; double TrendCounter=0;
        for (int t=1; t<T; t++) {
            gsl_matrix_set(Moments, 1, t, MALAutoCorrelation (Moments, t, Week, 0)); // AC of log-returns at lag of Week
            gsl_matrix_set(Moments, 2, t, MALAutoCorrelation (Moments, t, 2*Week, 0)); // AC of log-returns at lag of 2*Week
            gsl_matrix_set(Moments, 3, t, MALAutoCorrelation (Moments, t, Month, 0)); // AC of log-returns at lag of Month
            gsl_matrix_set(Moments, 4, t, MALAutoCorrelation (Moments, t, 3*Month, 0)); // AC of log-returns at lag of 3*Month
            gsl_matrix_set(Moments, 5, t, MALAutoCorrelation (Moments, t, 6*Month, 0)); // AC of log-returns at lag of 6*Month
            gsl_matrix_set(Moments, 6, t, MALAutoCorrelation (Moments, t, Year, 0)); // AC of log-returns at lag of Year
            gsl_matrix_set(Moments, 8, t, abs(MALAutoCorrelation (Moments, t, Week, 7))); // AC of abs log-returns at lag of Week
            gsl_matrix_set(Moments, 9, t, abs(MALAutoCorrelation (Moments, t, 2*Week, 7))); // AC of abs log-returns at lag of 2*Week
            gsl_matrix_set(Moments, 10, t, abs(MALAutoCorrelation (Moments, t, Month, 7))); // AC of abs log-returns at lag of Month
            gsl_matrix_set(Moments, 11, t, abs(MALAutoCorrelation (Moments, t, 3*Month, 7))); // AC of abs log-returns at lag of 3*Month
            gsl_matrix_set(Moments, 12, t, abs(MALAutoCorrelation (Moments, t, 6*Month, 7))); // AC of abs log-returns at lag of 6*Month
            gsl_matrix_set(Moments, 13, t, abs(MALAutoCorrelation (Moments, t, Year, 7))); // AC of abs log-returns at lag of Year
            gsl_matrix_set(Moments, 14, t, MALVolatility (Result, t, Week, m)); // Volatility at lag of Week
            gsl_matrix_set(Moments, 16, t, MALVolatility (Result, t, 2*Week, m)); // Volatility at lag of 2*Week
            gsl_matrix_set(Moments, 18, t, MALVolatility (Result, t, Month, m)); // Volatility at lag of Month
            gsl_matrix_set(Moments, 20, t, MALVolatility (Result, t, 3*Month, m)); // Volatility at lag of 3*Month
            gsl_matrix_set(Moments, 22, t, MALVolatility (Result, t, 6*Month, m)); // Volatility at lag of 6*Month
            gsl_matrix_set(Moments, 24, t, MALVolatility (Result, t, Year, m)); // Volatility at lag of Year
            gsl_matrix_set(Moments, 27, t, MALAutoCorrelation (Moments, t, Week, 26)); // AC of volumes at lag of Week
            gsl_matrix_set(Moments, 28, t, MALAutoCorrelation (Moments, t, 2*Week, 26)); // AC of volumes at lag of 2*Week
            gsl_matrix_set(Moments, 29, t, MALAutoCorrelation (Moments, t, Month, 26)); // AC of volumes at lag of Month
            gsl_matrix_set(Moments, 30, t, MALAutoCorrelation (Moments, t, 3*Month, 26)); // AC of volumes at lag of 3*Month
            gsl_matrix_set(Moments, 31, t, MALAutoCorrelation (Moments, t, 6*Month, 26)); // AC of volumes at lag of 6*Month
            gsl_matrix_set(Moments, 32, t, MALAutoCorrelation (Moments, t, Year, 26)); // AC of volumes at lag of Year
            
            // AC-logreturns of [t-∆, t] & [t-2∆, t-∆] for ∆={w, 2w, 3w, m}
            gsl_matrix_set(Moments, 33, t, MALAutoCorrelation (Moments, t, 3*Week, 0)); // AC of log-returns at lag of Week
            // AC-logreturns of [t-∆, t] & [t-∆-∂,t-∂] for ∂={1, 2, 3, 4, w} and ∆={w, 2w, 3w, m}
            gsl_matrix_set(Moments, 34, t, MALAutoCorrelationBlend (Moments, t, Week, 1, 0)); // AC of log-returns at lag of Week shifted by 1
            gsl_matrix_set(Moments, 35, t, MALAutoCorrelationBlend (Moments, t, Week, 2, 0)); // AC of log-returns at lag of Week shifted by 2
            gsl_matrix_set(Moments, 36, t, MALAutoCorrelationBlend (Moments, t, Week, 3, 0)); // AC of log-returns at lag of Week shifted by 3
            gsl_matrix_set(Moments, 37, t, MALAutoCorrelationBlend (Moments, t, Week, 4, 0)); // AC of log-returns at lag of Week shifted by 4
            gsl_matrix_set(Moments, 38, t, MALAutoCorrelationBlend (Moments, t, Week, Week, 0)); // AC of log-returns at lag of Week shifted by Week
            
            gsl_matrix_set(Moments, 39, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*1, 0)); // AC of log-returns at lag of 2Week shifted by 2
            gsl_matrix_set(Moments, 40, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*2, 0)); // AC of log-returns at lag of 2Week shifted by 4
            gsl_matrix_set(Moments, 41, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*3, 0)); // AC of log-returns at lag of 2Week shifted by 6
            gsl_matrix_set(Moments, 42, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*4, 0)); // AC of log-returns at lag of 2Week shifted by 8
            gsl_matrix_set(Moments, 43, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*Week, 0)); // AC of log-returns at lag of 2Week shifted by 2Week
            
            gsl_matrix_set(Moments, 44, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*1, 0)); // AC of log-returns at lag of 3Week shifted by 3
            gsl_matrix_set(Moments, 45, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*2, 0)); // AC of log-returns at lag of 3Week shifted by 6
            gsl_matrix_set(Moments, 46, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*3, 0)); // AC of log-returns at lag of 3Week shifted by 9
            gsl_matrix_set(Moments, 47, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*4, 0)); // AC of log-returns at lag of 3Week shifted by 12
            gsl_matrix_set(Moments, 48, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*Week, 0)); // AC of log-returns at lag of 3Week shifted by 3Week
            
            gsl_matrix_set(Moments, 49, t, MALAutoCorrelationBlend (Moments, t, Month, 4*1, 0)); // AC of log-returns at lag of Month shifted by 4
            gsl_matrix_set(Moments, 50, t, MALAutoCorrelationBlend (Moments, t, Month, 4*2, 0)); // AC of log-returns at lag of Month shifted by 8
            gsl_matrix_set(Moments, 51, t, MALAutoCorrelationBlend (Moments, t, Month, 4*3, 0)); // AC of log-returns at lag of Month shifted by 12
            gsl_matrix_set(Moments, 52, t, MALAutoCorrelationBlend (Moments, t, Month, 4*4, 0)); // AC of log-returns at lag of Month shifted by 16
            gsl_matrix_set(Moments, 53, t, MALAutoCorrelationBlend (Moments, t, Month, Month, 0)); // AC of log-returns at lag of Month shifted by Month
            
            //Systemics
            gsl_matrix_set(Moments, 54, t, gsl_matrix_get(Result, m, t)); // Market prices
            gsl_matrix_set(Moments, 55, t, -100+100*gsl_matrix_get(Result, m, t)/gsl_matrix_get(Result, m, t-1)); // Percentage of returns -100+100*P(t)/P(t-1)
            gsl_matrix_set(Moments, 56, t, gsl_matrix_get(Result2, m, t)); // Volumes
            gsl_matrix_set(Moments, 57, t, gsl_matrix_get(Moments, 14, t)); // w-volatility
            gsl_matrix_set(Moments, 58, t, gsl_matrix_get(Moments, 18, t)); // m-volatility
            gsl_matrix_set(Moments, 59, t, gsl_matrix_get(Moments, 22, t)); // 6m-volatility
            if (t>2*Week) {
                double Mean1=0; for (int k=t-2*Week; k<t-Week; k++) {Mean1+=gsl_matrix_get (Result, m, k)/Week;};
                double Mean2=0; for (int k=t-Week; k<t; k++) {Mean2+=gsl_matrix_get (Result, m, k)/Week;};
                gsl_matrix_set (Moments, 60, t, -100+100*Mean2/Mean1); // -100+100*<P(t-w,t)/P(t-2w,t-w)> Allows to see crashes (<-20% in the distribution)
            };
            if (gsl_matrix_get (Moments, 60, t)<=-20) {Crash=1;} else {Crash=0;}; gsl_matrix_set (Moments, 61, t, Crash); // Crash
            if ((gsl_matrix_get (Moments, 60, t-1)>0) && (gsl_matrix_get (Moments, 60, t)<0)) {TrendCounter=0;}; // Trend reversion
            if ((gsl_matrix_get (Moments, 60, t-1)<0) && (gsl_matrix_get (Moments, 60, t)>0)) {TrendCounter=0;}; // Trend reversion
            if ((gsl_matrix_get (Moments, 60, t-1)>0) && (gsl_matrix_get (Moments, 60, t)>0)) {TrendCounter+=1;}; // Trend reversion
            if ((gsl_matrix_get (Moments, 60, t-1)<0) && (gsl_matrix_get (Moments, 60, t)<0)) {TrendCounter-=1;}; // Trend reversion
            gsl_matrix_set (Moments, 62, t, TrendCounter);
            
        };
        for (int t=1; t<T; t++) {
            gsl_matrix_set(Moments, 15, t, MALAutoCorrelation (Moments, t, Week, 14)); // AC of volatility at lag of Week for lag of Week
            gsl_matrix_set(Moments, 17, t, MALAutoCorrelation (Moments, t, 2*Week, 16)); // AC of volatility at lag of 2*Week for lag of 2*Week
            gsl_matrix_set(Moments, 19, t, MALAutoCorrelation (Moments, t, Month, 18)); // AC of volatility at lag of Month for lag of Month
            gsl_matrix_set(Moments, 21, t, MALAutoCorrelation (Moments, t, 3*Month, 20)); // AC of volatility at lag of 3*Month for lag of 3*Month
            gsl_matrix_set(Moments, 23, t, MALAutoCorrelation (Moments, t, 6*Month, 22)); // AC of volatility at lag of 6*Month for lag of 6*Month
            gsl_matrix_set(Moments, 25, t, MALAutoCorrelation (Moments, t, Year, 24)); // AC of volatility at lag of Year for lag of Year
        };
        
        vector<double> V, GoF; for (int t=0; t<T; t++) {V.push_back(gsl_matrix_get (Prices, 0, t)); GoF.push_back(1);}; // FFF2
        //vector<double> * pV= & V; // POSTDOC2022
        // const double Rescale=10; // Rescale factor // POSTDOC2022
        //const double Threshold=0.1; // Threshold (0.1) // POSTDOC2022
        //const size_t NumIteration=10000000; // Number of iterations of the minimizer (stops before once a minimum is found)
        // double R; double * pR = & R; // Hawkes process parameters (does not need to be initialized) // POSTDOC2022
        //vector<double> & pGoF = GoF; // Goodness of fit // POSTDOC2022
        //vector<double> * pRes = Intensity (pV, Rescale, Threshold, NumIteration, pR, pGoF); // POSTDOC2022
        /*
        vector<double> * pRes = 0;
        for (int t=0; t<T; t++) {gsl_matrix_set (Moments, 63, t, (*pRes)[t]);}; // Hawkes process intensity in column 9
        for (int t=0; t<T; t++) { // Hawkes process intensity in column 9
            if ((*pRes)[t]>10) { // Maximum value of Hawkes at 10
                gsl_matrix_set (Moments, 63, t, 10);
            }
            else if ((*pRes)[t]<0) { // Minimum value of Hawkes at 0
                gsl_matrix_set (Moments, 63, t, 0);
            };
        };
        */
        //gsl_matrix * Hawkes = HistoricalIntensity (Prices, 0);
        //for (int t=1; t<T; t++) {gsl_matrix_set(Moments, 42, t, gsl_matrix_get (Hawkes, 0, t));};
        
        //delete pV; delete pR; delete pRes;
        //pV=0; pR=0; pRes=0;
        
        PlotGSLMatrix(Moments, "MomentsReal_PortfolioMultiOutput.xls", 1);
        ofstream outputLog(Machine+"MomentsReal_PortfolioMultiOutput.xls", ofstream::app); outputLog << endl;
        Temp.push_back(Moments);
        FinalResult.push_back(Temp);
        //Temp.erase(Temp.begin(), Temp.end());
    };
    
    
    return FinalResult;
};



vector<vector<gsl_matrix *> > PortfolioMultiOutput_trunk (vector<Share> PF, int Size) {
    if (Size==0) {Size=int(PF.size());};
    string A = Machine+"PF_PortfolioMultiOutput_trunk.csv"; ofstream outputSpecs(A.c_str(), ofstream::app);
    vector<int> J=Shuffle(int(PF.size()));
    gsl_matrix * Result = gsl_matrix_calloc (Size+1, PF[0].Data->size2);
    gsl_matrix * Result2 = gsl_matrix_calloc (Size+1, PF[0].Data->size2);
    outputSpecs << fixed << "Date (yyyymmdd),"; for (int i=1; i<Size+1; i++) {outputSpecs << PF[J[i-1]].Symbol << " (" << PF[J[i-1]].Exchange << "),";}; outputSpecs << endl;
    //outputSpecs << fixed << "Date,"; for (int i=1; i<Size+1; i++) {outputSpecs << PF[J[i-1]].Exchange <<",";}; outputSpecs << endl;
    outputSpecs.close();
    for (int t=0; t<PF[0].Data->size2; t++) {// Dates
        gsl_matrix_set(Result, 0, t, gsl_matrix_get(PF[0].Data, 0, t));
        gsl_matrix_set(Result2, 0, t, gsl_matrix_get(PF[0].Data, 0, t));
    };
    for (int i=1; i<Size+1; i++) {
        for (int t=0; t<PF[0].Data->size2; t++) {
            gsl_matrix_set(Result, i, t, gsl_matrix_get(PF[J[i-1]].Data, 2, t)); // Close prices
            gsl_matrix_set(Result2, i, t, gsl_matrix_get(PF[J[i-1]].Data, 3, t)); // Volumes
        };
        //cout << PF[J[i-1]].Exchange << ": " << PF[J[i-1]].Symbol << ", " << PF[J[i-1]].Title << " selected" << endl;
    };
    PlotGSLMatrix(Result, "PF_PortfolioMultiOutput_trunk.csv", 1);
    //outputMoments << "Log-return" << '\t' << "AC-1w Log-return" << '\t' << "AC-2w Log-return" << '\t' << "AC-m Log-return" << '\t' << "AC-3m Log-return" << '\t' << "AC-6m Log-return" << '\t' << "AC-y Log-return" << '\t' << "Abs-log-return" << '\t' << "AC-1w Abs-log-return" << '\t' << "AC-2w Abs-log-return" << '\t' << "AC-m Abs-log-return" << '\t' << "AC-3m Abs-log-return" << '\t' << "AC-6m Abs-log-return" << '\t' << "AC-y Abs-log-return" << '\t' << "w-Volatility" << '\t' << "AC-w w-Volatility" << '\t' << "2w-Volatility" << '\t' << "AC-2w 2w-Volatility" << '\t' << "m-Volatility" << '\t' << "AC-m m-Volatility" << '\t' << "3m-Volatility" << '\t' << "AC-3m 3m-Volatility" << '\t' << "6m-Volatility" << '\t' << "AC-6m 6m-Volatility" << '\t' << "y-Volatility" << '\t' << "AC-y y-Volatility" << '\t' << "Volumes(bsp)" << '\t' << "AC-1w Volume" << '\t' << "AC-2w Volume" << '\t' << "AC-m Volume" << '\t' << "AC-3m Volume" << '\t' << "AC-6m Volume" << '\t' << "AC-y Volume" << endl;
    vector<vector<gsl_matrix *> > FinalResult;
    int T=int(PF[0].Data->size2);
    for (int m=1; m<Size+1; m++) {
        vector<gsl_matrix *> Temp;
        gsl_matrix * Moments = gsl_matrix_calloc (64, T);
        gsl_matrix * Prices = gsl_matrix_calloc (1, T); for (int t=0; t<T; t++) {gsl_matrix_set (Prices, 0, t, gsl_matrix_get (Result, m, t));};
        //PlotGSLMatrix(Prices, "Prices.csv", 1);
        for (int t=1; t<T; t++) {
            if (gsl_matrix_get(Result, m, t-1)!=0) {gsl_matrix_set(Moments, 0, t, log(gsl_matrix_get(Result, m, t)/gsl_matrix_get(Result, m, t-1)));} // Log-return
            else {gsl_matrix_set(Moments, 0, t, 0);}; // Log-return failsafe
            gsl_matrix_set(Moments, 7, t, abs(gsl_matrix_get(Moments, 0, t))); // Abs log-return
            //if (gsl_matrix_get(Result2, m, t-1)!=0) {gsl_matrix_set(Moments, 26, t, gsl_matrix_get(Result2, m, t)/gsl_matrix_get(Result2, m, t-1));} // Volumes-return
            //else {gsl_matrix_set(Moments, 26, t, 1);}; // Volumes-return failsafe
            gsl_matrix_set(Moments, 26, t, gsl_matrix_get(Result2, m, t)); // Volumes
        };
        double Crash=0; double TrendCounter=0;
        for (int t=1; t<T; t++) {
            gsl_matrix_set(Moments, 1, t, MALAutoCorrelation (Moments, t, Week, 0)); // AC of log-returns at lag of Week
            gsl_matrix_set(Moments, 2, t, MALAutoCorrelation (Moments, t, 2*Week, 0)); // AC of log-returns at lag of 2*Week
            gsl_matrix_set(Moments, 3, t, MALAutoCorrelation (Moments, t, Month, 0)); // AC of log-returns at lag of Month
            gsl_matrix_set(Moments, 4, t, MALAutoCorrelation (Moments, t, 3*Month, 0)); // AC of log-returns at lag of 3*Month
            gsl_matrix_set(Moments, 5, t, MALAutoCorrelation (Moments, t, 6*Month, 0)); // AC of log-returns at lag of 6*Month
            gsl_matrix_set(Moments, 6, t, MALAutoCorrelation (Moments, t, Year, 0)); // AC of log-returns at lag of Year
            gsl_matrix_set(Moments, 8, t, abs(MALAutoCorrelation (Moments, t, Week, 7))); // AC of abs log-returns at lag of Week
            gsl_matrix_set(Moments, 9, t, abs(MALAutoCorrelation (Moments, t, 2*Week, 7))); // AC of abs log-returns at lag of 2*Week
            gsl_matrix_set(Moments, 10, t, abs(MALAutoCorrelation (Moments, t, Month, 7))); // AC of abs log-returns at lag of Month
            gsl_matrix_set(Moments, 11, t, abs(MALAutoCorrelation (Moments, t, 3*Month, 7))); // AC of abs log-returns at lag of 3*Month
            gsl_matrix_set(Moments, 12, t, abs(MALAutoCorrelation (Moments, t, 6*Month, 7))); // AC of abs log-returns at lag of 6*Month
            gsl_matrix_set(Moments, 13, t, abs(MALAutoCorrelation (Moments, t, Year, 7))); // AC of abs log-returns at lag of Year
            gsl_matrix_set(Moments, 14, t, MALVolatility (Result, t, Week, m)); // Volatility at lag of Week
            gsl_matrix_set(Moments, 16, t, MALVolatility (Result, t, 2*Week, m)); // Volatility at lag of 2*Week
            gsl_matrix_set(Moments, 18, t, MALVolatility (Result, t, Month, m)); // Volatility at lag of Month
            gsl_matrix_set(Moments, 20, t, MALVolatility (Result, t, 3*Month, m)); // Volatility at lag of 3*Month
            gsl_matrix_set(Moments, 22, t, MALVolatility (Result, t, 6*Month, m)); // Volatility at lag of 6*Month
            gsl_matrix_set(Moments, 24, t, MALVolatility (Result, t, Year, m)); // Volatility at lag of Year
            gsl_matrix_set(Moments, 27, t, MALAutoCorrelation (Moments, t, Week, 26)); // AC of volumes at lag of Week
            gsl_matrix_set(Moments, 28, t, MALAutoCorrelation (Moments, t, 2*Week, 26)); // AC of volumes at lag of 2*Week
            gsl_matrix_set(Moments, 29, t, MALAutoCorrelation (Moments, t, Month, 26)); // AC of volumes at lag of Month
            gsl_matrix_set(Moments, 30, t, MALAutoCorrelation (Moments, t, 3*Month, 26)); // AC of volumes at lag of 3*Month
            gsl_matrix_set(Moments, 31, t, MALAutoCorrelation (Moments, t, 6*Month, 26)); // AC of volumes at lag of 6*Month
            gsl_matrix_set(Moments, 32, t, MALAutoCorrelation (Moments, t, Year, 26)); // AC of volumes at lag of Year
            
            // AC-logreturns of [t-∆, t] & [t-2∆, t-∆] for ∆={w, 2w, 3w, m}
            gsl_matrix_set(Moments, 33, t, MALAutoCorrelation (Moments, t, 3*Week, 0)); // AC of log-returns at lag of Week
            // AC-logreturns of [t-∆, t] & [t-∆-∂,t-∂] for ∂={1, 2, 3, 4, w} and ∆={w, 2w, 3w, m}
            gsl_matrix_set(Moments, 34, t, MALAutoCorrelationBlend (Moments, t, Week, 1, 0)); // AC of log-returns at lag of Week shifted by 1
            gsl_matrix_set(Moments, 35, t, MALAutoCorrelationBlend (Moments, t, Week, 2, 0)); // AC of log-returns at lag of Week shifted by 2
            gsl_matrix_set(Moments, 36, t, MALAutoCorrelationBlend (Moments, t, Week, 3, 0)); // AC of log-returns at lag of Week shifted by 3
            gsl_matrix_set(Moments, 37, t, MALAutoCorrelationBlend (Moments, t, Week, 4, 0)); // AC of log-returns at lag of Week shifted by 4
            gsl_matrix_set(Moments, 38, t, MALAutoCorrelationBlend (Moments, t, Week, Week, 0)); // AC of log-returns at lag of Week shifted by Week
            
            gsl_matrix_set(Moments, 39, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*1, 0)); // AC of log-returns at lag of 2Week shifted by 2
            gsl_matrix_set(Moments, 40, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*2, 0)); // AC of log-returns at lag of 2Week shifted by 4
            gsl_matrix_set(Moments, 41, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*3, 0)); // AC of log-returns at lag of 2Week shifted by 6
            gsl_matrix_set(Moments, 42, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*4, 0)); // AC of log-returns at lag of 2Week shifted by 8
            gsl_matrix_set(Moments, 43, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*Week, 0)); // AC of log-returns at lag of 2Week shifted by 2Week
            
            gsl_matrix_set(Moments, 44, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*1, 0)); // AC of log-returns at lag of 3Week shifted by 3
            gsl_matrix_set(Moments, 45, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*2, 0)); // AC of log-returns at lag of 3Week shifted by 6
            gsl_matrix_set(Moments, 46, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*3, 0)); // AC of log-returns at lag of 3Week shifted by 9
            gsl_matrix_set(Moments, 47, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*4, 0)); // AC of log-returns at lag of 3Week shifted by 12
            gsl_matrix_set(Moments, 48, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*Week, 0)); // AC of log-returns at lag of 3Week shifted by 3Week
            
            gsl_matrix_set(Moments, 49, t, MALAutoCorrelationBlend (Moments, t, Month, 4*1, 0)); // AC of log-returns at lag of Month shifted by 4
            gsl_matrix_set(Moments, 50, t, MALAutoCorrelationBlend (Moments, t, Month, 4*2, 0)); // AC of log-returns at lag of Month shifted by 8
            gsl_matrix_set(Moments, 51, t, MALAutoCorrelationBlend (Moments, t, Month, 4*3, 0)); // AC of log-returns at lag of Month shifted by 12
            gsl_matrix_set(Moments, 52, t, MALAutoCorrelationBlend (Moments, t, Month, 4*4, 0)); // AC of log-returns at lag of Month shifted by 16
            gsl_matrix_set(Moments, 53, t, MALAutoCorrelationBlend (Moments, t, Month, Month, 0)); // AC of log-returns at lag of Month shifted by Month
            
            //Systemics
            gsl_matrix_set(Moments, 54, t, gsl_matrix_get(Result, m, t)); // Market prices
            gsl_matrix_set(Moments, 55, t, -100+100*gsl_matrix_get(Result, m, t)/gsl_matrix_get(Result, m, t-1)); // Percentage of returns -100+100*P(t)/P(t-1)
            gsl_matrix_set(Moments, 56, t, gsl_matrix_get(Result2, m, t)); // Volumes
            gsl_matrix_set(Moments, 57, t, gsl_matrix_get(Moments, 14, t)); // w-volatility
            gsl_matrix_set(Moments, 58, t, gsl_matrix_get(Moments, 18, t)); // m-volatility
            gsl_matrix_set(Moments, 59, t, gsl_matrix_get(Moments, 22, t)); // 6m-volatility
            if (t>2*Week) {
                double Mean1=0; for (int k=t-2*Week; k<t-Week; k++) {Mean1+=gsl_matrix_get (Result, m, k)/Week;};
                double Mean2=0; for (int k=t-Week; k<t; k++) {Mean2+=gsl_matrix_get (Result, m, k)/Week;};
                gsl_matrix_set (Moments, 60, t, -100+100*Mean2/Mean1); // -100+100*<P(t-w,t)/P(t-2w,t-w)> Allows to see crashes (<-20% in the distribution)
            };
            if (gsl_matrix_get (Moments, 60, t)<=-20) {Crash=1;} else {Crash=0;}; gsl_matrix_set (Moments, 61, t, Crash); // Crash
            if ((gsl_matrix_get (Moments, 60, t-1)>0) && (gsl_matrix_get (Moments, 60, t)<0)) {TrendCounter=0;}; // Trend reversion
            if ((gsl_matrix_get (Moments, 60, t-1)<0) && (gsl_matrix_get (Moments, 60, t)>0)) {TrendCounter=0;}; // Trend reversion
            if ((gsl_matrix_get (Moments, 60, t-1)>0) && (gsl_matrix_get (Moments, 60, t)>0)) {TrendCounter+=1;}; // Trend reversion
            if ((gsl_matrix_get (Moments, 60, t-1)<0) && (gsl_matrix_get (Moments, 60, t)<0)) {TrendCounter-=1;}; // Trend reversion
            gsl_matrix_set (Moments, 62, t, TrendCounter);
            
        };
        for (int t=1; t<T; t++) {
            gsl_matrix_set(Moments, 15, t, MALAutoCorrelation (Moments, t, Week, 14)); // AC of volatility at lag of Week for lag of Week
            gsl_matrix_set(Moments, 17, t, MALAutoCorrelation (Moments, t, 2*Week, 16)); // AC of volatility at lag of 2*Week for lag of 2*Week
            gsl_matrix_set(Moments, 19, t, MALAutoCorrelation (Moments, t, Month, 18)); // AC of volatility at lag of Month for lag of Month
            gsl_matrix_set(Moments, 21, t, MALAutoCorrelation (Moments, t, 3*Month, 20)); // AC of volatility at lag of 3*Month for lag of 3*Month
            gsl_matrix_set(Moments, 23, t, MALAutoCorrelation (Moments, t, 6*Month, 22)); // AC of volatility at lag of 6*Month for lag of 6*Month
            gsl_matrix_set(Moments, 25, t, MALAutoCorrelation (Moments, t, Year, 24)); // AC of volatility at lag of Year for lag of Year
        };
        
        PlotGSLMatrix(Moments, "MomentsReal_PortfolioMultiOutput.xls", 1);
        ofstream outputLog(Machine+"MomentsReal_PortfolioMultiOutput.xls", ofstream::app); outputLog << endl;
        Temp.push_back(Moments);
        FinalResult.push_back(Temp);
        //Temp.erase(Temp.begin(), Temp.end());
    };
    
    return FinalResult;
};



vector<vector<gsl_matrix *> > PortfolioMultiOutput2 (vector<Share> PF, int Size) {
    if (Size==0) {Size=int(PF.size());};
    string A = Machine+"PF_PortfolioMultiOutput2.csv"; ofstream outputSpecs(A.c_str(), ofstream::app);
    vector<int> J=Shuffle(int(PF.size()));
    gsl_matrix * Result = gsl_matrix_calloc (Size+1, PF[0].Data->size2);
    gsl_matrix * Result2 = gsl_matrix_calloc (Size+1, PF[0].Data->size2);
    outputSpecs << fixed << "Date (yyyymmdd),"; for (int i=1; i<Size+1; i++) {outputSpecs << PF[J[i-1]].Symbol << " (" << PF[J[i-1]].Exchange << "),";}; outputSpecs << endl;
    //outputSpecs << fixed << "Date,"; for (int i=1; i<Size+1; i++) {outputSpecs << PF[J[i-1]].Exchange <<",";}; outputSpecs << endl;
    outputSpecs.close();
    for (int t=0; t<PF[0].Data->size2; t++) {// Dates
        gsl_matrix_set(Result, 0, t, gsl_matrix_get(PF[0].Data, 0, t));
        gsl_matrix_set(Result2, 0, t, gsl_matrix_get(PF[0].Data, 0, t));
    };
    for (int i=1; i<Size+1; i++) {
        for (int t=0; t<PF[0].Data->size2; t++) {
            gsl_matrix_set(Result, i, t, gsl_matrix_get(PF[J[i-1]].Data, 2, t)); // Close prices
            gsl_matrix_set(Result2, i, t, gsl_matrix_get(PF[J[i-1]].Data, 3, t)); // Volumes
        };
        //cout << PF[J[i-1]].Exchange << ": " << PF[J[i-1]].Symbol << ", " << PF[J[i-1]].Title << " selected" << endl;
    };
    PlotGSLMatrix(Result, "PF_PortfolioMultiOutput2.csv", 1);
    //outputMoments << "Log-return" << '\t' << "AC-1w Log-return" << '\t' << "AC-2w Log-return" << '\t' << "AC-m Log-return" << '\t' << "AC-3m Log-return" << '\t' << "AC-6m Log-return" << '\t' << "AC-y Log-return" << '\t' << "Abs-log-return" << '\t' << "AC-1w Abs-log-return" << '\t' << "AC-2w Abs-log-return" << '\t' << "AC-m Abs-log-return" << '\t' << "AC-3m Abs-log-return" << '\t' << "AC-6m Abs-log-return" << '\t' << "AC-y Abs-log-return" << '\t' << "w-Volatility" << '\t' << "AC-w w-Volatility" << '\t' << "2w-Volatility" << '\t' << "AC-2w 2w-Volatility" << '\t' << "m-Volatility" << '\t' << "AC-m m-Volatility" << '\t' << "3m-Volatility" << '\t' << "AC-3m 3m-Volatility" << '\t' << "6m-Volatility" << '\t' << "AC-6m 6m-Volatility" << '\t' << "y-Volatility" << '\t' << "AC-y y-Volatility" << '\t' << "Volumes(bsp)" << '\t' << "AC-1w Volume" << '\t' << "AC-2w Volume" << '\t' << "AC-m Volume" << '\t' << "AC-3m Volume" << '\t' << "AC-6m Volume" << '\t' << "AC-y Volume" << endl;
    vector<vector<gsl_matrix *> > FinalResult;
    int T=int(PF[0].Data->size2);
    for (int m=1; m<Size+1; m++) {
        vector<gsl_matrix *> Temp;
        gsl_matrix * Moments = gsl_matrix_calloc (64, T);
        gsl_matrix * Prices = gsl_matrix_calloc (1, T); for (int t=0; t<T; t++) {gsl_matrix_set (Prices, 0, t, gsl_matrix_get (Result, m, t));};
        //PlotGSLMatrix(Prices, "Prices.csv", 1);
        for (int t=1; t<T; t++) {
            if (gsl_matrix_get(Result, m, t-1)!=0) {gsl_matrix_set(Moments, 0, t, log(gsl_matrix_get(Result, m, t)/gsl_matrix_get(Result, m, t-1)));} // Log-return
            else {gsl_matrix_set(Moments, 0, t, 0);}; // Log-return failsafe
            gsl_matrix_set(Moments, 7, t, abs(gsl_matrix_get(Moments, 0, t))); // Abs log-return
            //if (gsl_matrix_get(Result2, m, t-1)!=0) {gsl_matrix_set(Moments, 26, t, gsl_matrix_get(Result2, m, t)/gsl_matrix_get(Result2, m, t-1));} // Volumes-return
            //else {gsl_matrix_set(Moments, 26, t, 1);}; // Volumes-return failsafe
            gsl_matrix_set(Moments, 26, t, gsl_matrix_get(Result2, m, t)); // Volumes
        };
        double Crash=0; double TrendCounter=0;
        for (int t=1; t<T; t++) {
            gsl_matrix_set(Moments, 1, t, MALAutoCorrelation (Moments, t, Week, 0)); // AC of log-returns at lag of Week
            gsl_matrix_set(Moments, 2, t, MALAutoCorrelation (Moments, t, 2*Week, 0)); // AC of log-returns at lag of 2*Week
            gsl_matrix_set(Moments, 3, t, MALAutoCorrelation (Moments, t, Month, 0)); // AC of log-returns at lag of Month
            gsl_matrix_set(Moments, 4, t, MALAutoCorrelation (Moments, t, 3*Month, 0)); // AC of log-returns at lag of 3*Month
            gsl_matrix_set(Moments, 5, t, MALAutoCorrelation (Moments, t, 6*Month, 0)); // AC of log-returns at lag of 6*Month
            gsl_matrix_set(Moments, 6, t, MALAutoCorrelation (Moments, t, Year, 0)); // AC of log-returns at lag of Year
            gsl_matrix_set(Moments, 8, t, abs(MALAutoCorrelation (Moments, t, Week, 7))); // AC of abs log-returns at lag of Week
            gsl_matrix_set(Moments, 9, t, abs(MALAutoCorrelation (Moments, t, 2*Week, 7))); // AC of abs log-returns at lag of 2*Week
            gsl_matrix_set(Moments, 10, t, abs(MALAutoCorrelation (Moments, t, Month, 7))); // AC of abs log-returns at lag of Month
            gsl_matrix_set(Moments, 11, t, abs(MALAutoCorrelation (Moments, t, 3*Month, 7))); // AC of abs log-returns at lag of 3*Month
            gsl_matrix_set(Moments, 12, t, abs(MALAutoCorrelation (Moments, t, 6*Month, 7))); // AC of abs log-returns at lag of 6*Month
            gsl_matrix_set(Moments, 13, t, abs(MALAutoCorrelation (Moments, t, Year, 7))); // AC of abs log-returns at lag of Year
            gsl_matrix_set(Moments, 14, t, MALVolatility (Result, t, Week, m)); // Volatility at lag of Week
            gsl_matrix_set(Moments, 16, t, MALVolatility (Result, t, 2*Week, m)); // Volatility at lag of 2*Week
            gsl_matrix_set(Moments, 18, t, MALVolatility (Result, t, Month, m)); // Volatility at lag of Month
            gsl_matrix_set(Moments, 20, t, MALVolatility (Result, t, 3*Month, m)); // Volatility at lag of 3*Month
            gsl_matrix_set(Moments, 22, t, MALVolatility (Result, t, 6*Month, m)); // Volatility at lag of 6*Month
            gsl_matrix_set(Moments, 24, t, MALVolatility (Result, t, Year, m)); // Volatility at lag of Year
            gsl_matrix_set(Moments, 27, t, MALAutoCorrelation (Moments, t, Week, 26)); // AC of volumes at lag of Week
            gsl_matrix_set(Moments, 28, t, MALAutoCorrelation (Moments, t, 2*Week, 26)); // AC of volumes at lag of 2*Week
            gsl_matrix_set(Moments, 29, t, MALAutoCorrelation (Moments, t, Month, 26)); // AC of volumes at lag of Month
            gsl_matrix_set(Moments, 30, t, MALAutoCorrelation (Moments, t, 3*Month, 26)); // AC of volumes at lag of 3*Month
            gsl_matrix_set(Moments, 31, t, MALAutoCorrelation (Moments, t, 6*Month, 26)); // AC of volumes at lag of 6*Month
            gsl_matrix_set(Moments, 32, t, MALAutoCorrelation (Moments, t, Year, 26)); // AC of volumes at lag of Year
            
            // AC-logreturns of [t-∆, t] & [t-2∆, t-∆] for ∆={w, 2w, 3w, m}
            gsl_matrix_set(Moments, 33, t, MALAutoCorrelation (Moments, t, 3*Week, 0)); // AC of log-returns at lag of Week
            // AC-logreturns of [t-∆, t] & [t-∆-∂,t-∂] for ∂={1, 2, 3, 4, w} and ∆={w, 2w, 3w, m}
            gsl_matrix_set(Moments, 34, t, MALAutoCorrelationBlend (Moments, t, Week, 1, 0)); // AC of log-returns at lag of Week shifted by 1
            gsl_matrix_set(Moments, 35, t, MALAutoCorrelationBlend (Moments, t, Week, 2, 0)); // AC of log-returns at lag of Week shifted by 2
            gsl_matrix_set(Moments, 36, t, MALAutoCorrelationBlend (Moments, t, Week, 3, 0)); // AC of log-returns at lag of Week shifted by 3
            gsl_matrix_set(Moments, 37, t, MALAutoCorrelationBlend (Moments, t, Week, 4, 0)); // AC of log-returns at lag of Week shifted by 4
            gsl_matrix_set(Moments, 38, t, MALAutoCorrelationBlend (Moments, t, Week, Week, 0)); // AC of log-returns at lag of Week shifted by Week
            
            gsl_matrix_set(Moments, 39, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*1, 0)); // AC of log-returns at lag of 2Week shifted by 2
            gsl_matrix_set(Moments, 40, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*2, 0)); // AC of log-returns at lag of 2Week shifted by 4
            gsl_matrix_set(Moments, 41, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*3, 0)); // AC of log-returns at lag of 2Week shifted by 6
            gsl_matrix_set(Moments, 42, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*4, 0)); // AC of log-returns at lag of 2Week shifted by 8
            gsl_matrix_set(Moments, 43, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*Week, 0)); // AC of log-returns at lag of 2Week shifted by 2Week
            
            gsl_matrix_set(Moments, 44, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*1, 0)); // AC of log-returns at lag of 3Week shifted by 3
            gsl_matrix_set(Moments, 45, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*2, 0)); // AC of log-returns at lag of 3Week shifted by 6
            gsl_matrix_set(Moments, 46, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*3, 0)); // AC of log-returns at lag of 3Week shifted by 9
            gsl_matrix_set(Moments, 47, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*4, 0)); // AC of log-returns at lag of 3Week shifted by 12
            gsl_matrix_set(Moments, 48, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*Week, 0)); // AC of log-returns at lag of 3Week shifted by 3Week
            
            gsl_matrix_set(Moments, 49, t, MALAutoCorrelationBlend (Moments, t, Month, 4*1, 0)); // AC of log-returns at lag of Month shifted by 4
            gsl_matrix_set(Moments, 50, t, MALAutoCorrelationBlend (Moments, t, Month, 4*2, 0)); // AC of log-returns at lag of Month shifted by 8
            gsl_matrix_set(Moments, 51, t, MALAutoCorrelationBlend (Moments, t, Month, 4*3, 0)); // AC of log-returns at lag of Month shifted by 12
            gsl_matrix_set(Moments, 52, t, MALAutoCorrelationBlend (Moments, t, Month, 4*4, 0)); // AC of log-returns at lag of Month shifted by 16
            gsl_matrix_set(Moments, 53, t, MALAutoCorrelationBlend (Moments, t, Month, Month, 0)); // AC of log-returns at lag of Month shifted by Month
            
            //Systemics
            gsl_matrix_set(Moments, 54, t, gsl_matrix_get(Result, m, t)); // Market prices
            gsl_matrix_set(Moments, 55, t, -100+100*gsl_matrix_get(Result, m, t)/gsl_matrix_get(Result, m, t-1)); // Percentage of returns -100+100*P(t)/P(t-1)
            gsl_matrix_set(Moments, 56, t, gsl_matrix_get(Result2, m, t)); // Volumes
            gsl_matrix_set(Moments, 57, t, gsl_matrix_get(Moments, 14, t)); // w-volatility
            gsl_matrix_set(Moments, 58, t, gsl_matrix_get(Moments, 18, t)); // m-volatility
            gsl_matrix_set(Moments, 59, t, gsl_matrix_get(Moments, 22, t)); // 6m-volatility
            if (t>2*Week) {
                double Mean1=0; for (int k=t-2*Week; k<t-Week; k++) {Mean1+=gsl_matrix_get (Result, m, k)/Week;};
                double Mean2=0; for (int k=t-Week; k<t; k++) {Mean2+=gsl_matrix_get (Result, m, k)/Week;};
                gsl_matrix_set (Moments, 60, t, -100+100*Mean2/Mean1); // -100+100*<P(t-w,t)/P(t-2w,t-w)> Allows to see crashes (<-20% in the distribution)
            };
            if (gsl_matrix_get (Moments, 60, t)<=-20) {Crash=1;} else {Crash=0;}; gsl_matrix_set (Moments, 61, t, Crash); // Crash
            if ((gsl_matrix_get (Moments, 60, t-1)>0) && (gsl_matrix_get (Moments, 60, t)<0)) {TrendCounter=0;}; // Trend reversion
            if ((gsl_matrix_get (Moments, 60, t-1)<0) && (gsl_matrix_get (Moments, 60, t)>0)) {TrendCounter=0;}; // Trend reversion
            if ((gsl_matrix_get (Moments, 60, t-1)>0) && (gsl_matrix_get (Moments, 60, t)>0)) {TrendCounter+=1;}; // Trend reversion
            if ((gsl_matrix_get (Moments, 60, t-1)<0) && (gsl_matrix_get (Moments, 60, t)<0)) {TrendCounter-=1;}; // Trend reversion
            gsl_matrix_set (Moments, 62, t, TrendCounter);
            
        };
        for (int t=1; t<T; t++) {
            gsl_matrix_set(Moments, 15, t, MALAutoCorrelation (Moments, t, Week, 14)); // AC of volatility at lag of Week for lag of Week
            gsl_matrix_set(Moments, 17, t, MALAutoCorrelation (Moments, t, 2*Week, 16)); // AC of volatility at lag of 2*Week for lag of 2*Week
            gsl_matrix_set(Moments, 19, t, MALAutoCorrelation (Moments, t, Month, 18)); // AC of volatility at lag of Month for lag of Month
            gsl_matrix_set(Moments, 21, t, MALAutoCorrelation (Moments, t, 3*Month, 20)); // AC of volatility at lag of 3*Month for lag of 3*Month
            gsl_matrix_set(Moments, 23, t, MALAutoCorrelation (Moments, t, 6*Month, 22)); // AC of volatility at lag of 6*Month for lag of 6*Month
            gsl_matrix_set(Moments, 25, t, MALAutoCorrelation (Moments, t, Year, 24)); // AC of volatility at lag of Year for lag of Year
        };
        
        vector<double> V, GoF; for (int t=0; t<T; t++) {V.push_back(gsl_matrix_get (Prices, 0, t)); GoF.push_back(1);};
        PlotGSLMatrix(Moments, "MomentsReal_PortfolioMultiOutput2.xls", 1);
        ofstream outputLog(Machine+"MomentsReal_PortfolioMultiOutput2.xls", ofstream::app); outputLog << endl;
        Temp.push_back(Moments);
        FinalResult.push_back(Temp);
        //Temp.erase(Temp.begin(), Temp.end());
    };
    
    
    return FinalResult;
};






vector<vector<gsl_matrix *> > PortfolioMultiOutput_coin (vector<Coin> PF, int Size) {
    if (Size==0) {Size=int(PF.size());};
    string A = Machine+"PF_PortfolioMultiOutput_coin.csv"; ofstream outputSpecs(A.c_str(), ofstream::app);
    vector<int> J=Shuffle(int(PF.size()));
    gsl_matrix * Result = gsl_matrix_calloc (Size+1, PF[0].Data->size2);
    gsl_matrix * Result2 = gsl_matrix_calloc (Size+1, PF[0].Data->size2);
    outputSpecs << fixed << "Date (yyyymmdd),"; for (int i=0; i<Size; i++) {outputSpecs << PF[J[i]].Symbol << " (" << PF[J[i]].Exchange << "),";}; outputSpecs << endl;
    //outputSpecs << fixed << "Date,"; for (int i=1; i<Size+1; i++) {outputSpecs << PF[J[i-1]].Exchange <<",";}; outputSpecs << endl;
    outputSpecs.close();
    for (int t=0; t<PF[0].Data->size2; t++) {// Dates
        gsl_matrix_set(Result, 0, t, t+1);
        gsl_matrix_set(Result2, 0, t, t+1);
    };
    for (int i=1; i<Size+1; i++) {
        for (int t=0; t<PF[0].Data->size2; t++) {
            gsl_matrix_set(Result, i, t, gsl_matrix_get(PF[J[i-1]].Data, 2, t)); // BTC close prices
            gsl_matrix_set(Result2, i, t, gsl_matrix_get(PF[J[i-1]].Data, 3, t)); // BTC volumes
            //gsl_matrix_set(Result, i, t, gsl_matrix_get(PF[J[i-1]].Data, 0, t)); // USDT close prices //FFFKKK
            //gsl_matrix_set(Result2, i, t, gsl_matrix_get(PF[J[i-1]].Data, 1, t)); // USDT volumes //FFFKKK
        };
        //cout << PF[J[i-1]].Exchange << ": " << PF[J[i-1]].Symbol << ", " << PF[J[i-1]].Title << " selected" << endl;
    };
    PlotGSLMatrix(Result, "PF_PortfolioMultiOutput_coin.csv", 1);
    //outputMoments << "Log-return" << '\t' << "AC-1w Log-return" << '\t' << "AC-2w Log-return" << '\t' << "AC-m Log-return" << '\t' << "AC-3m Log-return" << '\t' << "AC-6m Log-return" << '\t' << "AC-y Log-return" << '\t' << "Abs-log-return" << '\t' << "AC-1w Abs-log-return" << '\t' << "AC-2w Abs-log-return" << '\t' << "AC-m Abs-log-return" << '\t' << "AC-3m Abs-log-return" << '\t' << "AC-6m Abs-log-return" << '\t' << "AC-y Abs-log-return" << '\t' << "w-Volatility" << '\t' << "AC-w w-Volatility" << '\t' << "2w-Volatility" << '\t' << "AC-2w 2w-Volatility" << '\t' << "m-Volatility" << '\t' << "AC-m m-Volatility" << '\t' << "3m-Volatility" << '\t' << "AC-3m 3m-Volatility" << '\t' << "6m-Volatility" << '\t' << "AC-6m 6m-Volatility" << '\t' << "y-Volatility" << '\t' << "AC-y y-Volatility" << '\t' << "Volumes(bsp)" << '\t' << "AC-1w Volume" << '\t' << "AC-2w Volume" << '\t' << "AC-m Volume" << '\t' << "AC-3m Volume" << '\t' << "AC-6m Volume" << '\t' << "AC-y Volume" << endl;
    vector<vector<gsl_matrix *> > FinalResult;
    int T=int(PF[0].Data->size2);
    for (int m=1; m<Size+1; m++) {
        vector<gsl_matrix *> Temp;
        gsl_matrix * Moments = gsl_matrix_calloc (64, T);
        gsl_matrix * Prices = gsl_matrix_calloc (1, T); for (int t=0; t<T; t++) {gsl_matrix_set (Prices, 0, t, gsl_matrix_get (Result, m, t));};
        //PlotGSLMatrix(Prices, "Prices.csv", 1);
        for (int t=1; t<T; t++) {
            if (gsl_matrix_get(Result, m, t-1)!=0) {gsl_matrix_set(Moments, 0, t, log(gsl_matrix_get(Result, m, t)/gsl_matrix_get(Result, m, t-1)));} // Log-return
            else {gsl_matrix_set(Moments, 0, t, 0);}; // Log-return failsafe
            gsl_matrix_set(Moments, 7, t, abs(gsl_matrix_get(Moments, 0, t))); // Abs log-return
            //if (gsl_matrix_get(Result2, m, t-1)!=0) {gsl_matrix_set(Moments, 26, t, gsl_matrix_get(Result2, m, t)/gsl_matrix_get(Result2, m, t-1));} // Volumes-return
            //else {gsl_matrix_set(Moments, 26, t, 1);}; // Volumes-return failsafe
            gsl_matrix_set(Moments, 26, t, gsl_matrix_get(Result2, m, t)); // Volumes
        };
        double Crash=0; double TrendCounter=0;
        for (int t=1; t<T; t++) {
            gsl_matrix_set(Moments, 1, t, MALAutoCorrelation (Moments, t, Week, 0)); // AC of log-returns at lag of Week
            gsl_matrix_set(Moments, 2, t, MALAutoCorrelation (Moments, t, 2*Week, 0)); // AC of log-returns at lag of 2*Week
            gsl_matrix_set(Moments, 3, t, MALAutoCorrelation (Moments, t, Month, 0)); // AC of log-returns at lag of Month
            gsl_matrix_set(Moments, 4, t, MALAutoCorrelation (Moments, t, 3*Month, 0)); // AC of log-returns at lag of 3*Month
            gsl_matrix_set(Moments, 5, t, MALAutoCorrelation (Moments, t, 6*Month, 0)); // AC of log-returns at lag of 6*Month
            gsl_matrix_set(Moments, 6, t, MALAutoCorrelation (Moments, t, Year, 0)); // AC of log-returns at lag of Year
            gsl_matrix_set(Moments, 8, t, abs(MALAutoCorrelation (Moments, t, Week, 7))); // AC of abs log-returns at lag of Week
            gsl_matrix_set(Moments, 9, t, abs(MALAutoCorrelation (Moments, t, 2*Week, 7))); // AC of abs log-returns at lag of 2*Week
            gsl_matrix_set(Moments, 10, t, abs(MALAutoCorrelation (Moments, t, Month, 7))); // AC of abs log-returns at lag of Month
            gsl_matrix_set(Moments, 11, t, abs(MALAutoCorrelation (Moments, t, 3*Month, 7))); // AC of abs log-returns at lag of 3*Month
            gsl_matrix_set(Moments, 12, t, abs(MALAutoCorrelation (Moments, t, 6*Month, 7))); // AC of abs log-returns at lag of 6*Month
            gsl_matrix_set(Moments, 13, t, abs(MALAutoCorrelation (Moments, t, Year, 7))); // AC of abs log-returns at lag of Year
            gsl_matrix_set(Moments, 14, t, MALVolatility (Result, t, Week, m)); // Volatility at lag of Week
            gsl_matrix_set(Moments, 16, t, MALVolatility (Result, t, 2*Week, m)); // Volatility at lag of 2*Week
            gsl_matrix_set(Moments, 18, t, MALVolatility (Result, t, Month, m)); // Volatility at lag of Month
            gsl_matrix_set(Moments, 20, t, MALVolatility (Result, t, 3*Month, m)); // Volatility at lag of 3*Month
            gsl_matrix_set(Moments, 22, t, MALVolatility (Result, t, 6*Month, m)); // Volatility at lag of 6*Month
            gsl_matrix_set(Moments, 24, t, MALVolatility (Result, t, Year, m)); // Volatility at lag of Year
            gsl_matrix_set(Moments, 27, t, MALAutoCorrelation (Moments, t, Week, 26)); // AC of volumes at lag of Week
            gsl_matrix_set(Moments, 28, t, MALAutoCorrelation (Moments, t, 2*Week, 26)); // AC of volumes at lag of 2*Week
            gsl_matrix_set(Moments, 29, t, MALAutoCorrelation (Moments, t, Month, 26)); // AC of volumes at lag of Month
            gsl_matrix_set(Moments, 30, t, MALAutoCorrelation (Moments, t, 3*Month, 26)); // AC of volumes at lag of 3*Month
            gsl_matrix_set(Moments, 31, t, MALAutoCorrelation (Moments, t, 6*Month, 26)); // AC of volumes at lag of 6*Month
            gsl_matrix_set(Moments, 32, t, MALAutoCorrelation (Moments, t, Year, 26)); // AC of volumes at lag of Year
            
            // AC-logreturns of [t-∆, t] & [t-2∆, t-∆] for ∆={w, 2w, 3w, m}
            gsl_matrix_set(Moments, 33, t, MALAutoCorrelation (Moments, t, 3*Week, 0)); // AC of log-returns at lag of Week
            // AC-logreturns of [t-∆, t] & [t-∆-∂,t-∂] for ∂={1, 2, 3, 4, w} and ∆={w, 2w, 3w, m}
            gsl_matrix_set(Moments, 34, t, MALAutoCorrelationBlend (Moments, t, Week, 1, 0)); // AC of log-returns at lag of Week shifted by 1
            gsl_matrix_set(Moments, 35, t, MALAutoCorrelationBlend (Moments, t, Week, 2, 0)); // AC of log-returns at lag of Week shifted by 2
            gsl_matrix_set(Moments, 36, t, MALAutoCorrelationBlend (Moments, t, Week, 3, 0)); // AC of log-returns at lag of Week shifted by 3
            gsl_matrix_set(Moments, 37, t, MALAutoCorrelationBlend (Moments, t, Week, 4, 0)); // AC of log-returns at lag of Week shifted by 4
            gsl_matrix_set(Moments, 38, t, MALAutoCorrelationBlend (Moments, t, Week, Week, 0)); // AC of log-returns at lag of Week shifted by Week
            
            gsl_matrix_set(Moments, 39, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*1, 0)); // AC of log-returns at lag of 2Week shifted by 2
            gsl_matrix_set(Moments, 40, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*2, 0)); // AC of log-returns at lag of 2Week shifted by 4
            gsl_matrix_set(Moments, 41, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*3, 0)); // AC of log-returns at lag of 2Week shifted by 6
            gsl_matrix_set(Moments, 42, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*4, 0)); // AC of log-returns at lag of 2Week shifted by 8
            gsl_matrix_set(Moments, 43, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*Week, 0)); // AC of log-returns at lag of 2Week shifted by 2Week
            
            gsl_matrix_set(Moments, 44, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*1, 0)); // AC of log-returns at lag of 3Week shifted by 3
            gsl_matrix_set(Moments, 45, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*2, 0)); // AC of log-returns at lag of 3Week shifted by 6
            gsl_matrix_set(Moments, 46, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*3, 0)); // AC of log-returns at lag of 3Week shifted by 9
            gsl_matrix_set(Moments, 47, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*4, 0)); // AC of log-returns at lag of 3Week shifted by 12
            gsl_matrix_set(Moments, 48, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*Week, 0)); // AC of log-returns at lag of 3Week shifted by 3Week
            
            gsl_matrix_set(Moments, 49, t, MALAutoCorrelationBlend (Moments, t, Month, 4*1, 0)); // AC of log-returns at lag of Month shifted by 4
            gsl_matrix_set(Moments, 50, t, MALAutoCorrelationBlend (Moments, t, Month, 4*2, 0)); // AC of log-returns at lag of Month shifted by 8
            gsl_matrix_set(Moments, 51, t, MALAutoCorrelationBlend (Moments, t, Month, 4*3, 0)); // AC of log-returns at lag of Month shifted by 12
            gsl_matrix_set(Moments, 52, t, MALAutoCorrelationBlend (Moments, t, Month, 4*4, 0)); // AC of log-returns at lag of Month shifted by 16
            gsl_matrix_set(Moments, 53, t, MALAutoCorrelationBlend (Moments, t, Month, Month, 0)); // AC of log-returns at lag of Month shifted by Month
            
            //Systemics
            gsl_matrix_set(Moments, 54, t, gsl_matrix_get(Result, m, t)); // Market prices
            gsl_matrix_set(Moments, 55, t, -100+100*gsl_matrix_get(Result, m, t)/gsl_matrix_get(Result, m, t-1)); // Percentage of returns -100+100*P(t)/P(t-1)
            gsl_matrix_set(Moments, 56, t, gsl_matrix_get(Result2, m, t)); // Volumes
            gsl_matrix_set(Moments, 57, t, gsl_matrix_get(Moments, 14, t)); // w-volatility
            gsl_matrix_set(Moments, 58, t, gsl_matrix_get(Moments, 18, t)); // m-volatility
            gsl_matrix_set(Moments, 59, t, gsl_matrix_get(Moments, 22, t)); // 6m-volatility
            if (t>2*Week) {
                double Mean1=0; for (int k=t-2*Week; k<t-Week; k++) {Mean1+=gsl_matrix_get (Result, m, k)/Week;};
                double Mean2=0; for (int k=t-Week; k<t; k++) {Mean2+=gsl_matrix_get (Result, m, k)/Week;};
                gsl_matrix_set (Moments, 60, t, -100+100*Mean2/Mean1); // -100+100*<P(t-w,t)/P(t-2w,t-w)> Allows to see crashes (<-20% in the distribution)
            };
            if (gsl_matrix_get (Moments, 60, t)<=-20) {Crash=1;} else {Crash=0;}; gsl_matrix_set (Moments, 61, t, Crash); // Crash
            if ((gsl_matrix_get (Moments, 60, t-1)>0) && (gsl_matrix_get (Moments, 60, t)<0)) {TrendCounter=0;}; // Trend reversion
            if ((gsl_matrix_get (Moments, 60, t-1)<0) && (gsl_matrix_get (Moments, 60, t)>0)) {TrendCounter=0;}; // Trend reversion
            if ((gsl_matrix_get (Moments, 60, t-1)>0) && (gsl_matrix_get (Moments, 60, t)>0)) {TrendCounter+=1;}; // Trend reversion
            if ((gsl_matrix_get (Moments, 60, t-1)<0) && (gsl_matrix_get (Moments, 60, t)<0)) {TrendCounter-=1;}; // Trend reversion
            gsl_matrix_set (Moments, 62, t, TrendCounter);
            
        };
        for (int t=1; t<T; t++) {
            gsl_matrix_set(Moments, 15, t, MALAutoCorrelation (Moments, t, Week, 14)); // AC of volatility at lag of Week for lag of Week
            gsl_matrix_set(Moments, 17, t, MALAutoCorrelation (Moments, t, 2*Week, 16)); // AC of volatility at lag of 2*Week for lag of 2*Week
            gsl_matrix_set(Moments, 19, t, MALAutoCorrelation (Moments, t, Month, 18)); // AC of volatility at lag of Month for lag of Month
            gsl_matrix_set(Moments, 21, t, MALAutoCorrelation (Moments, t, 3*Month, 20)); // AC of volatility at lag of 3*Month for lag of 3*Month
            gsl_matrix_set(Moments, 23, t, MALAutoCorrelation (Moments, t, 6*Month, 22)); // AC of volatility at lag of 6*Month for lag of 6*Month
            gsl_matrix_set(Moments, 25, t, MALAutoCorrelation (Moments, t, Year, 24)); // AC of volatility at lag of Year for lag of Year
        };
        
        vector<double> V, GoF; for (int t=0; t<T; t++) {V.push_back(gsl_matrix_get (Prices, 0, t)); GoF.push_back(1);};
        PlotGSLMatrix(Moments, "MomentsReal_PortfolioMultiOutput_coin.xls", 1);
        ofstream outputLog(Machine+"MomentsReal_PortfolioMultiOutput_coin.xls", ofstream::app); outputLog << endl;
        Temp.push_back(Moments);
        FinalResult.push_back(Temp);
        //Temp.erase(Temp.begin(), Temp.end());
    };
    
    
    return FinalResult;
};



// END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE ***
// END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE ***
// END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE ***
// END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE ***
// END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE ***







class Stock {
public:
    int StockName, StockQuantity, StockQuantityFirst;
    double StockBidValue, StockAskValue;
};


class Agent {
public:
    // GENERAL VARIABLES
    int AgentName; // Name of agent
    double RFA, RBA, RFAFirst, RBAFirst, NAVOnJanFirst; // Risk-free asset (Zero Coupon Bond, ZCB), Risk-based assets (total stock holdings)
    int Bankruptcy, LiquidationFloor, TradingWindow; // Bankruptcy parameters and trading window
    vector<Stock> Stocks; // Collection of stocks
    gsl_matrix * BiasedValues; // J*T matrix built via CointegratedWalk()
    vector<int> TradingWindowClock; // Clock counter reinitialized at each trading publication
    vector<int> ExitHorizons; // Times at which the agent must exit its position
    vector<vector<double> > Qs; // Qs[j][i] is the arg max_a Q(s,a) of state s of 5-index i for stock j for model-based RL (updated every TradingWindow)
    
    // COINTEGRATION VARIABLES
    vector<double> Leash, LeashVol, Accuracy; // One for each j stock. Accuracy is the prob. of bubble burst, on average 0.2,0.75,1,2 per year for Accuracy=13,11,10,9 resp.
    
    // BIAS VARIABLES
    double Gesture; // Commercial gesture in percent added (or subtracted!) to the Bid, subtracted (or added!) to the Ask
    double Versatility; // How often (number of time steps) on average does the agent change his weight
    double Reflexivity; // Reflexivity coefficient
    
    // RL VARIABLE
    //double DiscountFactor; // For the ARMA forecast
    int Future, History; // History is the past lag over which the agent performs comparison methods for FSIndex, FAIndexReal, FAIndexVirtual, ForecastReal, ForecastVirtual, FResultDisReal, FResultDisVirtual, AND TSIndex, TAIndexReal, QuantitiesReal, TransactionPriceReal, TAIndexVirtual, QuantitiesVirtual, TransactionPriceVirtual, TResultDisReal, TResultDisVirtual and also FResultReal, FResultVirtual, TResultReal, TResultVirtual (formerly HistoryRes).
    int Rough, Smooth, Reflexive; // Forecast states
    int Tool, Lag, Weight;// Forecast actions
    int Mu, Sig, RFALevel, RBALevel, Liquid; // Trading states
    int Quant, Pinch; // Trading actions
    int QTradeReal, QTradeVirtual; // Trading quantity to file in the OB
    int Exploration; // Epsilon-Greedy (for Exploration=1), Softmax (for Exploration=2), Pursuit (for Exploration=3), Off-Policy Watkin's Q-learning (for Exploration=4), SAM (none)
    int OffPolicy, OffPolicy2; // 0: no off-policy, 1: off-policy
    double Epsilon; // Epsilon-greedy method parameter (probability value below which exploration is performed)
    //int Update; // Action-value functions Q(s,a) and policies Pi(s,a) updating methods
    
    // Neuroeconomics biases
    int NEB; // Given the value {0,1,2,3,4,5,6} if belonging to class {no NEB, NEB124, NEB356, NEB89, NEB123456, NEB12489, NEB35689}.
    double Human, NEBLossAversion, NEBPositivity, NEBNegativity, NEBFear, NEBFear2, NEBGreed, NEBLearningRate;
    
    // RL INDICES // DDD
    int FS, FA, TS, TA; // Number of states S and actions A
    int FSDim[3], FADim[3], FQDim[6], FPiDim[6], TSDim[5], TADim[2], TQDim[7], TPiDim[7]; // Dimensions specifications
    //int FSDim[2], FADim[2], FQDim[4], FPiDim[4], TSDim[5], TADim[2], TQDim[7], TPiDim[7]; // Dimensions specifications // DDD2
    int FSInd[3], FAIndReal[3], FAIndVirtual[3], FAInd5[3], FQIndReal[6], FQIndVirtual[6], FQInd5[6], FPiInd[6];
    int TSInd[5], TAIndReal[2], TAIndVirtual[2], TAInd5[2], TQIndReal[7], TQIndVirtual[7], TQInd5[7], TPiInd[7];
    int FTenPiId[6], TTenPiId[7]; // To select an exploratory policy
    // int FTenPiId[4], TTenPiId[7]; // To select an exploratory policy // DDD2
    
    // MMM
    vector<vector<double> > dRoughVec, dSmoothVec, dMuPosVec, dMuNegVec, dSigVec, LiquidPercentVec, dFResultVec, dTResultVec; // Past values to determine median & percentiles and discrete score via History // MMM
    vector<vector<int> > FSIndex, TSIndex, FAIndexReal, FAIndexVirtual, TAIndexReal, TAIndexVirtual; // Real and virtual actions taken (identified by its vector index) at time t // MMM
    vector<vector<double> > FPi, TPi, FQ, TQ; // Action-value function Q(s)=E[R|s,a] as a SxA matrix // MMM
    vector<vector<double> > ForecastReal, ForecastVirtual, Forecast5; // MMM
    vector<vector<int> > FNumberA, TNumberA; // The number of times a was performed in s for SAM // MMM
    vector<vector<int> > QuantitiesReal, QuantitiesVirtual; // Quantities of stock at each time step that are willed to be transacted // MMM
    vector<vector<double> > TransactionPriceReal, TransactionPriceVirtual; // Transaction prices // MMM
    
    // LOG // MMM
    vector<vector<double> > LvolLog, SvolLog; // MMM
    vector<vector<int> > RoughLog, SmoothLog, ReflexiveLog, ToolRealLog, LagRealLog, WeightRealLog, ToolVirtualLog, LagVirtualLog, WeightVirtualLog, MuLog, SigLog, RFALog, RBALog, LiquidityLog, PinchRealLog, QuantRealLog, PinchVirtualLog, QuantVirtualLog, FResultDisReal, FResultDisVirtual, TResultDisReal, TResultDisVirtual; // MMM
    
    
    
    // Returns the index of a vector<double> Tensor (decomposed from a tensor of dimensions vector<int> Dimensions) corresponding to the vector<int> Indices. This can be used to update the corresponding element in Tensor, or to simply access it
    int get_tensor_index (const int IndSize, int Indices[], int Dimensions[]) {
        int Res=0;
        int ResTemp;
        for (int i=0; i<IndSize; i++) {
            ResTemp=Indices[i];
            for (int j=i+1; j<IndSize; j++) {
                ResTemp*=Dimensions[j];
            };
            Res+=ResTemp;
        };
        return Res;
    }; // closes method get_tensor_index()
    
    // Returns the D-dim tensor coordinates corresponding to the Index of the STL vector representation of that tensor
    vector<int> get_tensor_coord (int Index, int DimSize, int Dimensions[]) {
        //ofstream outputLog("/Users/admin/Documents/GNT/SYMBA/SimLog.txt", ofstream::app);
        vector<int> Res;
        int Multi1, Multi2, Sum;
        for (int k=0; k<DimSize; k++) {
            Res.push_back(0);
            Multi1=1; Sum=0;
            for (int n=k+1; n<DimSize; n++) {Multi1*=Dimensions[n];};
            for (int i=0; i<k; i++) {
                if (k==0) {break;};
                Multi2=1;
                for (int j=i+1; j<DimSize; j++) {Multi2*=Dimensions[j];};
                Sum+=Res[i]*Multi2;
            };
            //outputLog << "Index=" << Index << ", Sum=" << Sum << ", Multi1=" << Multi1 << endl;
            Res[k]=(Index - Sum)/Multi1;
        };
        return Res;
    }; // closes method get_coord()
    
    double Capital() { // Computation of total capital amount of portfolio
        // The NAV or Net Asset Value is given by the sum of the RFA and the collection of stocks valued at market price
        double Temp=0;
        for (int i=0; i<int(Stocks.size()); i++) {Temp += ((Stocks[i]).StockAskValue)*((Stocks[i]).StockQuantity);};
        return RFA+Temp;
    };
    
    double StockHoldings() { // Computation of stock holdings
        double Temp=0;
        for (int i=0; i<int(Stocks.size()); i++) {Temp += ((Stocks[i]).StockAskValue)*((Stocks[i]).StockQuantity);};
        return Temp;
    };
    
    void Liquidation(int i, vector<Agent> &Market) {
        Market[i].RFA=0; Market[i].RBA=0; // Clearing all his bond holdings
        for (int j=0; j<int(Market[i].Stocks.size()); j++) {Market[i].Stocks[j].StockQuantity = 0;}; // Clearing all the agent stocks holdings
        Market[i].Bankruptcy=0; // The agent gets entirely liquidated
    };
    
    
    
    
    
    
    
    
    
    
    
    // RL algorithm for forecasting and placing order in the OB (SS6)
    void RL (int j, int t, double Rate, gsl_matrix * ReflexiveValues, double VSpread, double LiquidPercent, int Time, int NumberOfStocks, string TradingFrequencyCond, string Plot, string VersatilityCondition, double MarketPerformance, int TimeSinceJan1st, int LearningPhase, int LeaderAgent, int LeaderQuant, int LeaderPinch, int ClusterLimit, int Trunk) {
        
        
        /*
         F(s): Rough={0,1,2,3}, Smooth={0,1,2,3}, Reflexive={0,1,2}
         F(a): Tool={0,1,2}, Lag={0,1,2}, Weight={0,1,2}
         T(s): Mu={0,1,2,3,4,5,6}, Sig={0,1,2}, RFA={0,1,2}, RBA={0,1,2}, Liquid={0,1,2,3}
         T(a): Quant={0,1,2,3,4,5,6}, Pinch={0,1,2}
         */
        ofstream outputLog(Machine+"SimLog.txt", ofstream::app);
        // ofstream outputDebug("/Users/admin/Documents/GNT/SYMBA/Debug.txt", ofstream::app);
        double Seed=1000*t*(RFA+RBA+AgentName);
        //double Seed=1000*t*gsl_matrix_get(ReflexiveValues, j, t);
        //ERROR RNG
        //#define GSL_DLL
        //#define DGSL_DLL
        //#define DWIN32
        //gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937); // or gsl_rng_default instead of gsl_rng_mt19937
        gsl_rng * r = make_rng();
        gsl_rng_set(r, static_cast<unsigned long int>(time(0)+Seed)); //gsl_rng_set(const gsl_rng * r, unsigned long int s)
        vector<double> VRan; for (int i=0; i<36; i++) {VRan.push_back(gsl_rng_uniform (r));}; // Generating a vector of random numbers without seeding issues
        //#undef GSL_DLL
        //#undef DGSL_DLL
        //#undef DWIN32
        //delete r; r=0; // Delete and de-allocate pointer
        if (Plot=="On") {
            outputLog << "***************RL() FOR AGENT " << AgentName << " STOCK " << j << " AT TIME t=" << t << "**************" << endl;
            outputLog << "ReflexiveValues[]={"; for (int h=0; h<=t; h++) {outputLog << gsl_matrix_get(ReflexiveValues, j, h); if (h<t) {outputLog << ", ";};};
            outputLog << "}" << endl;
            outputLog << "Seed=" << Seed << endl;
            outputLog << "VRan[]={" ; for (int u=0; u<int(VRan.size()); u++) {if (u<int(VRan.size())-1) {outputLog << VRan[u] << ", ";} else {outputLog << VRan[u];};}; outputLog << "}" << endl;
            outputLog << "At step ABM I.6, the bid and ask values are first defined (as resp. the time-discounted minimum and maximum values between present price and its forecast). Then at step ABM II.3, the gesture on bid-ask spread (''Pinch'') is incorporated in these values. Then at step ABM II.3, the actual quantity to short (QTradeReal<0) or long (QTradeReal>0) is also updated in class agent member (''Quant''). Notice that these quantities are already formated to the real numbers to trade, and as such are ready for transaction (except the sign)." << endl;
            
            
            /*************************************ABM I.1******************************************/
            // Selecting the trading window (Future) according to present time
            outputLog << endl << "    ABM I.1: Selecting the trading window (Future) according to present time" << endl;
            outputLog << "Future=" << Future << ", History=" << History << endl;
            
            
            /*************************************ABM I.2******************************************/
            outputLog << endl << "    ABM I.2: Defining present state s" << endl;
        };// closes Plot condition
        // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 1.1 computed..." << endl;
        // Get the state s
        double Mean=0; double Variance=0; int Start=t-(3*Future); if (Start<0) {Start=0;}; int End=t; // BOLLINGER ISSUE
        for (int p=Start; p<=End; p++) {Mean+=gsl_matrix_get(ReflexiveValues, j, p);}; Mean=Mean/(End-Start+1); // Mean
        for (int p=Start; p<=End; p++) {Variance+=abs(gsl_matrix_get(ReflexiveValues, j, p) - Mean);};
        double Lvol=Variance/(End-Start+1); // Variance Lvol
        if (Lvol<0.001) {Lvol=0.001;}; // Failsafe
        if (Plot=="On") {outputLog << "For t=[" << Start << "," << End << "], Mean=" << Mean << ", Variance=" << Variance << "  " << Lvol << endl;
            LvolLog[j].push_back(Lvol);}; // LLL
        Mean=0; Variance=0; Start=t-(Future)-1; if (Start<0) {Start=0;}; End=t; // BOLLINGER ISSUE
        for (int p=Start; p<=End; p++) {Mean+=gsl_matrix_get(ReflexiveValues, j, p);}; Mean=Mean/(End-Start+1); // Mean
        for (int p=Start; p<=End; p++) {Variance+=abs(gsl_matrix_get(ReflexiveValues, j, p) - Mean);};
        double Svol=Variance/(End-Start+1); // Variance Svol
        if (Svol<0.001) {Svol=0.001;}; // Failsafe
        if (Plot=="On") {outputLog << "For t=[" << Start << "," << End << "], Mean=" << Mean << ", Variance=" << Variance << " => Svol=" << Svol << endl;
            SvolLog[j].push_back(Svol);}; // LLL
        
        double dRough=Svol/Lvol; double dRoughPercentile=0; int RoughVecSize=int(dRoughVec[j].size());
        if (Lvol==0) {dRough=0;};
        if (t==1) {dRoughVec[j].push_back(dRough); dRoughPercentile=0.5;}
        else { // sorting the vector dRoughVec by ascending values
            for (int i=0; i<RoughVecSize; i++) {
                if (dRough<dRoughVec[j][i]) {dRoughVec[j].insert(dRoughVec[j].begin()+i, dRough*(1 + 0.0001*VRan[29])); dRoughPercentile=i*1.0/RoughVecSize; break;};
                if (i==RoughVecSize-1) {dRoughVec[j].push_back(dRough*(1 + 0.0001*VRan[29])); dRoughPercentile=1; break;};
            }; // closes for loop
        }; // closes else
        
        /*
         if (dRoughPercentile<0.25) {Rough=0;}
         else if ((dRoughPercentile>=0.25) && (dRoughPercentile<0.5)) {Rough=1;}
         else if ((dRoughPercentile>=0.5) && (dRoughPercentile<0.75)) {Rough=2;}
         else if (dRoughPercentile>=0.75) {Rough=3;};
         */
        
        if (dRoughPercentile<0.25) {Rough=0;}
        else if ((dRoughPercentile>=0.25) && (dRoughPercentile<0.75)) {Rough=1;}
        else if (dRoughPercentile>=0.75) {Rough=2;};
        
        if ((Lvol==0) || (t<=3)) {Rough=1;};
        if (int(dRoughVec[j].size()) - Future - History >0) {dRoughVec[j].erase (dRoughVec[j].begin());}; // OOO
        if (Plot=="On") {RoughLog[j].push_back(Rough);}; // LLL
        
        double dSmooth=Lvol/(3*Future); double dSmoothPercentile=0; int SmoothVecSize=int(dSmoothVec[j].size());
        if (3*Future==0) {dSmooth=0;};
        if (t==1) {dSmoothVec[j].push_back(dSmooth); dSmoothPercentile=0.5;}
        else { // sorting the vector SmoothVec by ascending values
            for (int i=0; i<SmoothVecSize; i++) {
                if (dSmooth<dSmoothVec[j][i]) {dSmoothVec[j].insert(dSmoothVec[j].begin()+i, dSmooth*(1 + 0.0001*VRan[29])); dSmoothPercentile=i*1.0/SmoothVecSize; break;};
                if (i==SmoothVecSize-1) {dSmoothVec[j].push_back(dSmooth*(1 + 0.0001*VRan[29])); dSmoothPercentile=1; break;};
            }; // closes for loop
        }; // closes else
        /*
         if (dSmoothPercentile<0.25) {Smooth=0;}
         else if ((dSmoothPercentile>=0.25) && (dSmoothPercentile<0.5)) {Smooth=1;}
         else if ((dSmoothPercentile>=0.5) && (dSmoothPercentile<0.75)) {Smooth=2;}
         else if (dSmoothPercentile>=0.75) {Smooth=3;};
         */
        
        if (dSmoothPercentile<0.25) {Smooth=0;}
        else if ((dSmoothPercentile>=0.25) && (dSmoothPercentile<0.75)) {Smooth=1;}
        else if (dSmoothPercentile>=0.75) {Smooth=2;};
        
        if ((t<=3) || (3*Future==0)) {Smooth=1;};
        if (int(dSmoothVec[j].size()) - Future - History >0) {dSmoothVec[j].erase (dSmoothVec[j].begin());}; // OOO
        if (Plot=="On") {SmoothLog[j].push_back(Smooth);}; // LLL
        
        // DDD
        double dReflexive=0; Start=t-3*Future; if (Start<0) {Start=0;};
        for (int p=Start; p<=t; p++) {dReflexive+=100*abs(gsl_matrix_get(ReflexiveValues, j, p)-gsl_matrix_get(BiasedValues, j, p))/(gsl_matrix_get(ReflexiveValues, j, p)*(t-Start+1));};
        Reflexive=0; // Reflexive={0,1,2}={[0%,10%[, [10%,30%[, [30%,+\infty[} for dReflexive=<100*abs(Pp-Biased)/Pp> over p[t-3*Future,t]
        if ((dReflexive>=0) && (dReflexive<10)) {Reflexive=0;}
        else if ((dReflexive>=10) && (dReflexive<30)) {Reflexive=1;}
        else if (dReflexive>=30) {Reflexive=2;};
        if (Plot=="On") {ReflexiveLog[j].push_back(Reflexive);}; // LLL
        
        if (Plot=="On") {
            outputLog << "dRough=Svol/Lvol=" << dRough << " and Rough=" << Rough << endl;
            outputLog << "dSmooth=Lvol/(3*Future)=" << dSmooth << " and Smooth=" << Smooth << endl;
            outputLog << "dReflexive=mean <100*(P(t)-B(t))/P(t)> over t=[t-3*Future, t]=" << dReflexive << " and Reflexive=" << Reflexive << endl; // DDD
            outputLog << "P(t)=" << gsl_matrix_get(ReflexiveValues, j, t) << ", B(t)=" << gsl_matrix_get(BiasedValues, j, t) << endl;
        }; // closes Plot condition
        FSInd[0]=Rough; FSInd[1]=Smooth; FSInd[2]=Reflexive; // DDD
        int STensorIndex = get_tensor_index (3, FSInd, FSDim); // Vector index in the tensor of all possible RL states // DDD
        FSIndex[j].push_back(STensorIndex); // Recording state s_t (ABM#1 for Forecast)
        if (int(FSIndex[j].size()) - Future - 1 >0) {FSIndex[j].erase (FSIndex[j].begin());}; // Updates by checking FSIndex[int(XXX.size()) - Future - 1] //OOO
        if (Plot=="On") {
            outputLog << "STensorIndex=" << STensorIndex << endl;
            outputLog << "FSIndex[]={" ; for (int u=0; u<int(FSIndex[j].size()); u++) {if (u<int(FSIndex[j].size())-1) {outputLog << FSIndex[j][u] << ", ";} else {outputLog << FSIndex[j][u];};}; outputLog << "}" << endl;
            outputLog << "FSDim[]={" ; for (int u=0; u<3; u++) {if (u<3-1) {outputLog << FSDim[u] << ", ";} else {outputLog << FSDim[u];};}; outputLog << "}" << endl;
            outputLog << "FSInd[]={" ; for (int u=0; u<3; u++) {if (u<3-1) {outputLog << FSInd[u] << ", ";} else {outputLog << FSInd[u];};}; outputLog << "}" << endl;
            
            
            /*************************************ABM I.3******************************************/
            outputLog << endl << "    ABM I.3: Initializing the start policy pi_0" << endl;
        }; // closes Plot condition
        // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 1.2 computed..." << endl;
        // Initializing the start policy $\pi_0$ out of which an an action will be selected (we initialize such that in all states s we take Tool=0, Lag=0)
        if (t==1) {
            FPiInd[0]=0; FPiInd[1]=0; FPiInd[2]=0; FPiInd[3]=0; FPiInd[4]=0; FPiInd[5]=0; // DDD
            for (int i=0; i<FS*FA; i++) {FPi[j].push_back(1.0/FA);};
            if (Plot=="On") {outputLog << "Policy initialized... all " << FA << " possible actions are given equiprobability 1/" << FA << "=" << 1.0/FA << endl;};
        }; // HHH All actions are equiprobable
        if (t>1) {if (Plot=="On") {outputLog << "Policy already initialized..." << endl;};};
        
        
        /*************************************ABM I.4******************************************/
        if (Plot=="On") {outputLog << endl << "    ABM I.4: Using current policy pi to select real action a=(Tool,Lag,Reflexive) in present state s" << endl;};
        // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 1.3 computed..." << endl;
        // We now use the current policy to perform an action $a$ in our state $s$. Note: There is no need to use \verb?PiStack?, all we need is when we have a vector of action each with a probability in $[0,1]$, we generate a random number from a uniform distribution: if that number is larger than the first probability it is subtracted from it, and if that new number is still larger than the second probability, we do likewise until we arrive to the probability corresponding to action $a$ to be picked up. This is because the sum of all these probabilities is $1$.
        FTenPiId[0]=Rough; FTenPiId[1]=Smooth; FTenPiId[2]=Reflexive; FTenPiId[3]=0; FTenPiId[4]=0; FTenPiId[5]=0; // DDD
        //FTenPiId[0]=Rough; FTenPiId[1]=Smooth; FTenPiId[2]=0; FTenPiId[3]=0; // DDD2
        int Fid = 0;
        Fid = get_tensor_index (6, FTenPiId, FPiDim); // Vector index of pi-tensor representation corresponding to present s and a=(0,0) // DDD
        //Fid = get_tensor_index (4, FTenPiId, FPiDim); // Vector index of pi-tensor representation corresponding to present s and a=(0,0) // DDD2
        if (Plot=="On") {
            outputLog << "FTenPiId[]={"; for (int u=0; u<6; u++) {if (u<6-1) {outputLog << FTenPiId[u] << ", ";} else {outputLog << FTenPiId[u];};}; outputLog << "}" << endl; // DDD
            //outputLog << "FTenPiId[]={"; for (int u=0; u<4; u++) {if (u<4-1) {outputLog << FTenPiId[u] << ", ";} else {outputLog << FTenPiId[u];};}; outputLog << "}" << endl; // DDD2
            outputLog << "Fid=" << Fid << ", FPi[j].size()=" << FPi[j].size() << endl;
        }; // closes Plot condition
        
        double Fix=VRan[0];
        double Fix2=VRan[1];
        double Act=0; int Pos=Fid;
        if (Plot=="On") {outputLog << "Fix=" << Fix << ", Fix2=" << Fix2 << endl;};
        
        if (Exploration==1) { // Epsilon-greedy: Choosing always the action with highest probability, but once in a while (with a small probability epsilon) choosing a (uniform) random one (which is a problem: the worst potential action can be taken as likely as the best, with in some cases devastating returns).
            if (Plot=="On") {outputLog << "Epsilon-greedy method selected" << endl;}; // DDD1
            if (Fix<Epsilon) { // Exploration
                Pos=Fid + int(floor(FA*Fix2));
                Act=FPi[j][Pos];
                if (Plot=="On") {outputLog << "Exploration: Chose random action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
            }
            else { // Exploitation
                for (int n=Fid; n<Fid+FA; n++) { // Running over the actions a associated with that state s
                    if (Plot=="On") {outputLog << "For n=" << n << ", FPi[n]=" << floor(100*FPi[j][n]) << "%" << endl;};
                    if (FPi[j][n]>=Act) {Pos=n; Act=FPi[j][Pos];}; // We try and find the action with greatest probability (greedy move)
                }; // closes for
                vector<int> PosVec; // Indices of equiprobable actions
                for (int n=Fid; n<Fid+FA; n++) {if (FPi[j][n]==Act) {PosVec.push_back(n);};}; // We record all the indices of equiprobable actions
                Pos=PosVec[int(floor(VRan[23]*PosVec.size()))]; // We choose a random indice among these!
                PosVec.erase(PosVec.begin(),PosVec.end());
                if (Plot=="On") {outputLog << "Exploitation epsilon-greedy: Chose greediest action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
            }; // closes else
        } // closes if
        
        else if (Exploration==2) { // Softmax: Same as above except the random choice in now not (uniform) equiprobable for all actions, but graded according to their estimated value. One method (using a Gibbs/Boltzmann distribution $e^{Q_t(a) / \tau} / \sum_{b=1}^n e^{Q_t(b) / \tau}$, with $\tau$ the temperature) allows to even shift this grading all the way back to epsilon greedy methods (when $\tau \rightarrow 0$).
            if (Plot=="On") {outputLog << "Softmax method selected" << endl;};
            if (Fix<Epsilon) { // Exploration
                Pos=Fid + int(floor(FA*Fix2));
                Act=FPi[j][Pos];
                if (Plot=="On") {outputLog << "Exploration: Chose random action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
            }
            else { // Exploitation
                vector<int> N = Shuffle(FA);
                if (Plot=="On") {outputLog << "From n=" << Fid << " to " << Fid+FA-1 << ": Fix2=" << Fix2 << endl;}; // DDD1
                for (int n=Fid; n<Fid+FA; n++) { // Running over the actions a associated with that state s
                    //if (Plot=="On") {outputLog << "For n=" << n << ", FPi[n]=" << floor(100*FPi[n]) << "%" << endl;};
                    int Newn = N[n-Fid]+Fid;
                    Fix2-=FPi[j][Newn];
                    if (Plot=="On") {outputLog << "For n=" << Newn << ", FPi[j][Newn]=" << floor(100*FPi[j][Newn]) << "%: Fix2-=FPi[j][Newn]=" << Fix2 << endl;};
                    if (Fix2<=0) {Pos=Newn; Act=FPi[j][Pos]; break;}; // We try and find a new action according to its graded probability
                }; // closes for
                if (Plot=="On") {outputLog << "Exploitation softmax: Chose greediest action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
            }; // closes else
        } // closes else if
        
        else if (Exploration==3) { // Pursuit: The probability of selecting the greedy action $a_{t+1}=a_{t+1}^{\ast}$ at next time step is incremented a fraction $\beta$ of the way toward $1$, while the probabilities of selecting the other actions are decremented toward $0$, respectively shown by the following two expressions: $\pi_{t+1}(a_{t+1}^{\ast}) = \pi_{t}(a_{t+1}^{\ast}) + \beta \left[ 1- \pi_{t}(a_{t+1}^{\ast}) \right]$ and $\pi_{t+1}(a) = \pi_{t}(a) + \beta \left[ 0- \pi_{t}(a) \right] , \hspace{5mm} \forall a\neq a_{t+1}^{\ast}$.
            if (Plot=="On") {outputLog << "Pursuit method selected" << endl;}; // DDD1
            if (Fix<Epsilon*(1-(double(t)/Time))) { // Exploration
                Pos=Fid + int(floor(FA*Fix2));
                Act=FPi[j][Pos];
                if (Plot=="On") {outputLog << "Exploration: Chose random action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
            }
            else { // Exploitation
                for (int n=Fid; n<Fid+FA; n++) { // Running over the actions a associated with that state s
                    if (Plot=="On") {outputLog << "For n=" << n << ", FPi[j][n]=" << floor(100*FPi[j][n]) << "%" << endl;};
                    if (FPi[j][n]>=Act) {Pos=n; Act=FPi[j][Pos];}; // We try and find the action with greatest probability (greedy move)
                }; // closes for
                vector<int> PosVec; // Indices of equiprobable actions
                for (int n=Fid; n<Fid+FA; n++) {if (FPi[j][n]==Act) {PosVec.push_back(n);};}; // We record all the indices of equiprobable actions
                Pos=PosVec[int(floor(VRan[23]*PosVec.size()))]; // We choose a random indice among these!
                PosVec.erase(PosVec.begin(),PosVec.end());
                if (Plot=="On") {outputLog << "Exploitation pursuit: Chose greediest action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
            }; // closes else
        } // closes else if
        
        else if (Exploration==4) { // This is for off-policy Watkin's Q-learning
            if (Plot=="On") {outputLog << "Off-policy Watkin's Q-learning method selected" << endl;}; // DDD1
            for (int n=Fid; n<Fid+FA; n++) { // Running over the actions a associated with that state s
                if (Plot=="On") {outputLog << "For n=" << n << ", FPi[j][n]=" << floor(100*FPi[j][n]) << "%" << endl;};
                if (FPi[j][n]>=Act) {Pos=n; Act=FPi[j][Pos];}; // We try and find the action with greatest probability (greedy move)
            }; // closes for
            vector<int> PosVec; // Indices of equiprobable actions
            for (int n=Fid; n<Fid+FA; n++) {if (FPi[j][n]==Act) {PosVec.push_back(n);};}; // We record all the indices of equiprobable actions
            Pos=PosVec[int(floor(VRan[23]*PosVec.size()))]; // We choose a random indice among these!
            PosVec.erase(PosVec.begin(),PosVec.end());
            if (Plot=="On") {outputLog << "Exploitation Off-policy: Chose greediest action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
        }; // closes else if
        
        int FormerWeight=Weight;
        vector<int> TenPiCoo = get_tensor_coord (Pos, 6, FPiDim); // Getting the tensor coordinates from that index
        //vector<int> TenPiCoo = get_tensor_coord (Pos, 4, FPiDim); // Getting the tensor coordinates from that index // DDD2
        Tool=TenPiCoo[3]; // Selection of first action from policy
        Lag=TenPiCoo[4]; // Selection of second action from policy
        Weight=TenPiCoo[5]; // Selection of third action from policy
        
        //if (Tool==0) {Tool+=1;}; // we change Poly to Linear UUU
        
        // Tuning the frequency of changing Weight according to Versatility
        if ((VRan[35] > 1.0/Versatility) && (t>1)) { // Case where we do not want to switch but maintain the previous choices for Weight
            Weight=FormerWeight;
            FTenPiId[0]=Rough; FTenPiId[1]=Smooth; FTenPiId[2]=Reflexive; FTenPiId[3]=Tool; FTenPiId[4]=Lag; FTenPiId[5]=Weight;
            Pos = get_tensor_index (6, FTenPiId, FPiDim); // Index representation of present s and a with former reflexivity weight
            TenPiCoo = get_tensor_coord (Pos, 6, FPiDim); // Getting the tensor coordinates from that index
            if (Plot=="On") {outputLog << "Weight overridden to FormerWeight=" << FormerWeight << " since VRan[35]=" << VRan[35] << ">" << "1.0/Versatility=" << 1.0/Versatility << endl;};
        }; // closes if
        
        // Turning Versatility condition on or off
        if (VersatilityCondition=="Off") {
            Weight=1; // Always fixed at Reflexivity
            FTenPiId[0]=Rough; FTenPiId[1]=Smooth; FTenPiId[2]=Reflexive; FTenPiId[3]=Tool; FTenPiId[4]=Lag; FTenPiId[5]=Weight;
            Pos = get_tensor_index (6, FTenPiId, FPiDim); // Index representation of present s and a with reflexivity weight (Weight=1)
            TenPiCoo = get_tensor_coord (Pos, 6, FPiDim); // Getting the tensor coordinates from that index
            if (Plot=="On") {outputLog << "Weight overridden to 1 (i.e to the natural Reflexivity of the agent) since VersatilityCondition is turned off" << endl;};
        }; // closes if
        
        if (Plot=="On") {outputLog << "TenPiCoo[]={" ; for (int u=0; u<int(TenPiCoo.size()); u++) {if (u<int(TenPiCoo.size())-1) {outputLog << TenPiCoo[u] << ", ";} else {outputLog << TenPiCoo[u];};}; outputLog << "}, ";
            outputLog << "Tool=TenPiCoo[3]=" << Tool << ", ";
            outputLog << "Lag=TenPiCoo[4]=" << Lag << ", ";
            outputLog << "Weight=TenPiCoo[5]=" << Weight << endl;
        }; // closes Plot condition
        // Defining the real action as tensor index to take according to above policy which selected Tool and Lag
        //Lag+=1; // To make sure the value 0 means 1, 1 means 2, and 2 means 3, resp.
        FAIndReal[0]=Tool; FAIndReal[1]=Lag; FAIndReal[2]=Weight;
        int ARealTensorIndex = get_tensor_index (3, FAIndReal, FADim); // Vector index in the tensor of all possible RL actions
        FAIndexReal[j].push_back(ARealTensorIndex); // Recording real action a_t (ABM#2 for Forecast)
        if (int(FAIndexReal[j].size()) - Future - 1 >0) {FAIndexReal[j].erase (FAIndexReal[j].begin());}; // OOO
        if (Plot=="On") {
            outputLog << "FAIndReal[0]=Tool=" << Tool << ", FAIndReal[1]=Lag=" << Lag << ", FAIndReal[2]=Weight=" << Weight << ", ";
            outputLog << "ARealTensorIndex=" << ARealTensorIndex << ", ";
            outputLog << "FAIndexReal[j][]={" ; for (int u=0; u<int(FAIndexReal[j].size()); u++) {if (u<int(FAIndexReal[j].size())-1) {outputLog << FAIndexReal[j][u] << ", ";} else {outputLog << FAIndexReal[j][u];};}; outputLog << "}" << endl;
            if (Plot=="On") {ToolRealLog[j].push_back(Tool); LagRealLog[j].push_back(Lag); WeightRealLog[j].push_back(Weight);}; // LLL
        }; // closes Plot condition
        // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 1.4 computed..." << endl;
        
        /*************************************ABM I.5******************************************/
        int VirtualTool=0; int VirtualLag=0; int VirtualWeight=0;
        if (Exploration==4) {
            if (Plot=="On") {    outputLog << endl << "    ABM I.5: Randomly selecting a virtual action a=(VirtualTool,VirtualLag)" << endl;};
            // Defining the virtual action as tensor index to take (not based on any exploration policy)
            VirtualTool=int(1000*VRan[2])%3; // 0, 1, or 2
            VirtualLag=int(1000*VRan[3])%3; // 1, 2, or 3
            VirtualWeight=int(1000*VRan[33])%3; // 1, 2, or 3
            for (int i=0; i<100; i++) {int VT=int(1000*VRan[4])%3; if (VT!=Tool) {VirtualTool=VT; break;};}; // Making sure virtual action not like real one
            for (int i=0; i<100; i++) {int VL=int(1000*VRan[5])%3; if (VL!=Lag) {VirtualLag=VL; break;};}; // Making sure virtual action not like real one
            for (int i=0; i<100; i++) {int VL=int(1000*VRan[34])%3; if (VL!=Weight) {VirtualWeight=VL; break;};}; // Making sure virtual action not like real one
            FAIndVirtual[0]=VirtualTool; FAIndVirtual[1]=VirtualLag; FAIndVirtual[2]=VirtualWeight; // Turlututu
            int AVirtualTensorIndex = get_tensor_index (3, FAIndVirtual, FADim); // Vector index in the tensor of all possible RL actions
            FAIndexVirtual[j].push_back(AVirtualTensorIndex); // Recording real action a_t (ABM#2 for Forecast)
            if (int(FAIndexVirtual[j].size()) - Future - 1 >0) {FAIndexVirtual[j].erase (FAIndexVirtual[j].begin());}; // OOO
            if (Plot=="On") {
                outputLog << "FAIndVirtual[0]=VirtualTool=" << VirtualTool << ", FAIndVirtual[1]=VirtualLag=" << VirtualLag << ", FAIndVirtual[2]=VirtualWeight=" << VirtualWeight << ", ";
                outputLog << "AVirtualTensorIndex=" << AVirtualTensorIndex << endl;
                outputLog << "FAIndexVirtual[j][]={" ; for (int u=0; u<int(FAIndexVirtual[j].size()); u++) {if (u<int(FAIndexVirtual[j].size())-1) {outputLog << FAIndexVirtual[j][u] << ", ";} else {outputLog << FAIndexVirtual[j][u];};}; outputLog << "}" << endl;
                if (Plot=="On") {ToolVirtualLog[j].push_back(VirtualTool); LagVirtualLog[j].push_back(VirtualLag); WeightVirtualLog[j].push_back(VirtualWeight);}; // LLL
            }; // closes Plot condition
        }; // closes Exploration==4 condition
        // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 1.5 computed..." << endl;
        
        
        /*************************************ABM I.6******************************************/
        if (Plot=="On") {outputLog << endl << "    ABM I.6: Forecasting the time series Future time steps ahead" << endl;};
        // F() must work on Reflexive not Reference, and the result of this is backed up and used for evaluation of the result of the whole F(). Then once this is done, we weight-average it with the value at time t+Future (not a forecast thereof!) of BiasedValues the cointegrated time series, and produce Mini and Maxi that will be input to T().
        double dForecastReal=gsl_matrix_get(ReflexiveValues, j, t); double Arma=0; double Linear=0; double Polynomial=0;
        double MeanGen=0; double MeanStep=0; double y0=0; double y1=0; double Poly=1; double DiscountSum=0; int Count=0;
        vector<double> Vec;
        if (Plot=="On") {outputLog << "We are in the case Tool=" << Tool;};
        if (t<(Lag+1)*Future) { // If not enough steps backwards everyone (Arma, Linear, Polynomial) is given the Mean of past reference values
            for (int n=0; n<=t; n++) {MeanGen+=gsl_matrix_get(ReflexiveValues, j, n);};
            dForecastReal=MeanGen/(t+1);
            if (Plot=="On") {outputLog << " and t=" << t << "<" << "(Lag+1)*Future-1=" << (Lag+1)*Future-1 << endl;};
            if (Plot=="On") {outputLog << "MeanGen=" << MeanGen << ", dForecastReal=" << dForecastReal << endl;};
        } // closes if
        else { // If enough steps backwards computing the mean over every (Lag+1) past array
            if (Plot=="On") {outputLog << " and t=" << t << ">=" << "(Lag+1)*Future-1=" << (Lag+1)*Future-1 << endl;};
            dForecastReal=BinaryProjection (ReflexiveValues, t, Tool, Lag, Future);
        }; // closes else
        //if (dForecastReal<0) {dForecastReal=0.01;}; // Just in case because of Polynomial easily having exaggerated expectations
        if (dForecastReal<=0) {dForecastReal=gsl_matrix_get(ReflexiveValues, j, t);}; // Failsafe for Polynomial exaggerated expectations
        if (Plot=="On") {outputLog << "dForecastReal=" << dForecastReal << endl;};
        //dForecastReal*=1/(1-(Future)*Rate/Year); // Time-discounting
        if (Plot=="On") {outputLog << "Time-discounted dForecastReal=" << dForecastReal << endl;};
        
        double dForecastVirtual=gsl_matrix_get(ReflexiveValues, j, t);
        if (Exploration==4) {
            Arma=0; Linear=0; Polynomial=0;
            MeanGen=0; MeanStep=0; y0=0; y1=0; Poly=1; DiscountSum=0; Count=0;
            vector<double> VecVirtual;
            if (Plot=="On") {outputLog << "We are in the case VirtualTool=" << VirtualTool;};
            if (t<(VirtualLag+1)*Future) { // If not enough steps backwards everyone (Arma, Linear, Polynomial) is given the Mean of past reference values
                for (int n=0; n<=t; n++) {MeanGen+=gsl_matrix_get(ReflexiveValues, j, n)/(t+1);};
                dForecastVirtual=MeanGen;
                if (Plot=="On") {outputLog << " and t=" << t << "<" << "(VirtualLag+1)*Future-1=" << (VirtualLag+1)*Future-1 << endl;};
            } // closes if
            else { // If enough steps backwards computing the mean over every (VirtualLag+1) past array
                dForecastVirtual=BinaryProjection (ReflexiveValues, t, VirtualTool, VirtualLag, Future);
            }; // closes else
            //if (dForecastVirtual<0) {dForecastVirtual=0.01;}; // Just in case because of Polynomial easily having exaggerated expectations
            if (dForecastVirtual<=0) {dForecastVirtual=gsl_matrix_get(ReflexiveValues, j, t);}; // Failsafe for Polynomial exaggerated expectations
            if (Plot=="On") {outputLog << "dForecastVirtual=" << dForecastVirtual << endl;};
            //dForecastVirtual*=1/(1-(Future)*Rate/Year); // Time-discounting
            if (Plot=="On") {outputLog << " Time-discounted dForecastVirtual=" << dForecastVirtual << endl;};
        }; // closes Exploration==4 condition
        
        //XXXExploration5
        if ((t-Future>=0) && (t%(2+Future/Month)==0)) {
            OffPolicy2=OffPolicy;
            if (Plot=="On") {outputLog << "We are in the case OffPolicy2=OffPolicy=" << OffPolicy << endl;};
        }
        else {
            OffPolicy2=0;
            if (Plot=="On") {outputLog << "We are in the case OffPolicy2=0" << endl;};
        };
        double dForecast5=gsl_matrix_get(ReflexiveValues, j, t); double MinimumCostFunction=9999999;
        int OptimalTool5=0; int OptimalLag5=0; int OptimalRWeight5=0;
        if (OffPolicy2==1) {
            for (int k1=0; k1<3; k1++) { // Tool
                for (int k2=0; k2<3; k2++) { // Lag
                    for (int k3=0; k3<3; k3++) { // RWeight
                        dForecast5=BinaryProjection (ReflexiveValues, t-Future, k1, k2, Future);
                        double RWeight5=0;
                        if (Reflexivity<=0.5) {
                            if (k3==0) {RWeight5=0;}
                            else if (k3==1) {RWeight5=100*Reflexivity;}
                            else if (k3==2) {RWeight5=200*Reflexivity;};
                        }
                        else {
                            if (k3==0) {RWeight5=100*(2*Reflexivity-1);}
                            else if (k3==1) {RWeight5=100*Reflexivity;}
                            else if (k3==2) {RWeight5=100;};
                        };
                        double XBid = gsl_matrix_get (BiasedValues, j , t-Future);
                        double YBid = dForecast5;
                        double b = RWeight5;
                        double a = abs(100-b);
                        double ForecastReference5=(XBid*a + YBid*b)/(a+b);
                        double CostFunction=100*abs(gsl_matrix_get(ReflexiveValues, j, t)-ForecastReference5)/gsl_matrix_get(ReflexiveValues, j, t); // Comparison between what F() would have forecasted with k1, k2, k3 and present realization
                        if (MinimumCostFunction>=CostFunction) {MinimumCostFunction=CostFunction; OptimalTool5=k1; OptimalLag5=k2; OptimalRWeight5=k3;}
                        if (Plot=="On") {outputLog << "For Tool=" << k1 << ", Lag="<< k2 << ", RWeight="<< k3 << ": dForecast5=$" << dForecast5 << ", ForecastReference5=$" << ForecastReference5 << ", P(t)=$" << gsl_matrix_get(ReflexiveValues, j, t) << ", CostFunction=" << CostFunction << "%" << endl;};
                    }; // closes k3
                }; // closes k2
            }; // closes k1
            if (Plot=="On") {outputLog << "MinimumCostFunction=" << MinimumCostFunction << " for OptimalTool5=" << OptimalTool5 << " (0: anti-trend, 1: mean-reversion, 2: trend), OptimalLag5=" << OptimalLag5 << ", OptimalRWeight5=" << OptimalRWeight5 << endl;};
        }; // closes OffPolicy2==1 condition
        
        // Computation of Mini and Maxi : minimum = (a<b) ? a : b;
        double RWeightReal=100*Reflexivity; double RWeightVirtual=100*Reflexivity;
        // Transforming Weight to real reflexivity coefficients RWeightReal
        if (Reflexivity<=0.5) {
            if (Weight==0) {RWeightReal=0;}
            else if (Weight==1) {RWeightReal=100*Reflexivity;}
            else if (Weight==2) {RWeightReal=200*Reflexivity;};
        }
        else {
            if (Weight==0) {RWeightReal=100*(2*Reflexivity-1);}
            else if (Weight==1) {RWeightReal=100*Reflexivity;}
            else if (Weight==2) {RWeightReal=100;};
        };
        //double XBid = gsl_matrix_get (BiasedValues, j , Time-1); if (t<Time-Future) {XBid = gsl_matrix_get (BiasedValues, j , t+Future);}; // Failsafe QQQ
        double XBid = gsl_matrix_get (BiasedValues, j , t); // Changed to present time (more rigorous approach)
        double YBid = dForecastReal;
        double b = RWeightReal;
        double a = abs(100-b);
        double ForecastReference=(XBid*a + YBid*b)/(a+b); // Now taking into account the dynamic reflexivity via Weight
        ForecastReal[j].push_back(ForecastReference); // Pushing the value into the vector
        if (int(ForecastReal[j].size()) - Future - 1 >0) {ForecastReal[j].erase (ForecastReal[j].begin());};
        if (Plot=="On") {
            outputLog << "Forecasted with ";
            if (Tool==0) {outputLog << "Anti-trend following";}
            else if (Tool==1) {outputLog << "Mean reversion";}
            else if (Tool==2) {outputLog << "Trend following";};
            outputLog << " ((Lag+1)*Future=" << (Lag+1)*Future << ", Future=" << Future << ")" << endl;
            outputLog << "ForecastReal[j][]={" ; for (int u=0; u<int(ForecastReal[j].size()); u++) {if (u<int(ForecastReal[j].size())-1) {outputLog << ForecastReal[j][u] << ", ";} else {outputLog << ForecastReal[j][u];};}; outputLog << "}" << endl;
        }; // closes Plot condition
        
        double ForecastReferenceVirtual=ForecastReference;
        if (Exploration==4) {
            if (Reflexivity<=0.5) {
                if (VirtualWeight==0) {RWeightVirtual=0;}
                else if (VirtualWeight==1) {RWeightVirtual=100*Reflexivity;}
                else if (VirtualWeight==2) {RWeightVirtual=200*Reflexivity;};
            }
            else {
                if (VirtualWeight==0) {RWeightVirtual=100*(2*Reflexivity-1);}
                else if (VirtualWeight==1) {RWeightVirtual=100*Reflexivity;}
                else if (VirtualWeight==2) {RWeightVirtual=100;};
            };
            //if (VirtualWeight==0) {RWeightVirtual=0*Reflexivity;} else if (VirtualWeight==1) {RWeightVirtual=50*Reflexivity;} else if (VirtualWeight==2) {RWeightVirtual=100*Reflexivity;}; // Transforming VirtualWeight to real reflexivity coefficients RWeightVirtual
            YBid = dForecastVirtual;
            b = RWeightVirtual;
            a = abs(100-b);
            ForecastReferenceVirtual=(XBid*a + YBid*b)/(a+b); // Now taking into account the dynamic reflexivity via Weight
            ForecastVirtual[j].push_back(ForecastReferenceVirtual); // Pushing the value into the vector
            if (int(ForecastVirtual[j].size()) - Future - 1 >0) {ForecastVirtual[j].erase (ForecastVirtual[j].begin());};
            if (Plot=="On") {outputLog << " to " << dForecastVirtual;};
            if (Plot=="On") {
                outputLog << ", forecasted with ";
                if (VirtualTool==0) {outputLog << "Anti-trend following";}
                else if (VirtualTool==1) {outputLog << "Mean reversion";}
                else if (VirtualTool==2) {outputLog << "Trend following";};
                outputLog << " ((VirtualLag+1)*Future=" << (VirtualLag+1)*Future << ", Future=" << Future << ")" << endl;
                outputLog << "ForecastVirtual[j][]={" ; for (int u=0; u<int(ForecastVirtual[j].size()); u++) {if (u<int(ForecastVirtual[j].size())-1) {outputLog << ForecastVirtual[j][u] << ", ";} else {outputLog << ForecastVirtual[j][u];};}; outputLog << "}" << endl;
            }; // closes Plot condition
        }; // closes Exploration==4 condition
        
        double ReflexiveVal_t = gsl_matrix_get (ReflexiveValues, j, t);
        //double ReflexiveVal2_t=ReflexiveVal_t*(1+Future*Rate/Year); // Time-discounting
        double ReflexiveVal2_t=ReflexiveVal_t;
        double Mini=ReflexiveVal2_t; double Maxi=ReflexiveVal2_t; int MiniTime=t; int MaxiTime=t+Future;
        if (ForecastReference<= ReflexiveVal2_t) {Mini=ForecastReference; Maxi=ReflexiveVal2_t; MiniTime=t+Future; MaxiTime=t;}
        else if (ForecastReference > ReflexiveVal2_t) {Maxi=ForecastReference; Mini=ReflexiveVal2_t; MaxiTime=t+Future; MiniTime=t;};
        Stocks[j].StockBidValue = Mini;
        Stocks[j].StockAskValue = Maxi;
        if (Trunk>-1) {
            Stocks[j].StockBidValue=DigitTrunk (Mini, Trunk, "Ceil");
            Stocks[j].StockAskValue=DigitTrunk (Maxi, Trunk, "Ceil");
        };
        
        if (Plot=="On") {
            outputLog << endl << "Reflexivity=" << Reflexivity << ", Weight=" << Weight << "=> RWeightReal=" << RWeightReal << "% as reflexivity coefficient, and Reflexive=" << Reflexive << endl;
            if (Exploration==4) {outputLog << "VirtualWeight=" << VirtualWeight << "=> RWeightVirtual=" << RWeightVirtual << "% as reflexivity coefficient..." << endl;};
            outputLog << endl << "Updating Bid and Ask:" << endl;
            outputLog <<"dForecastReal=" << dForecastReal << ", time-discounted P(t)=" << ReflexiveVal2_t << ", t=" << t << ", B(t)=" << gsl_matrix_get (BiasedValues, j , t) << ", ForecastReference=" << ForecastReference << endl;
            if (Exploration==4) {outputLog <<"dForecastVirtual=" << dForecastVirtual << ", time-discounted P(t)=" << ReflexiveVal2_t << ", t=" << t << ", ForecastReferenceVirtual=" << ForecastReferenceVirtual << endl;};
            outputLog << "Mini=Stocks[j].StockBidValue=" << Mini << ", MiniTime=" << MiniTime << endl;
            outputLog << "Maxi=Stocks[j].StockAskValue=" << Maxi << ", MaxiTime=" << MaxiTime << endl;
            
            
            
            /*************************************ABM I.7******************************************/
            outputLog << endl << "    ABM I.7: Calculating the return if t>Future as the percentage of difference between forecast and realization" << endl;
        }; // closes Plot condition
        // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 1.6 computed..." << endl;
        
        if (t==1) {
            dFResultVec[j].push_back(0);
            for (int i=0; i<FS*FA; i++) {FQ[j].push_back(0); FNumberA[j].push_back(0);};
            if (Plot=="On") {outputLog << "Initialization of FQ[], FNumberA[], by Agent " << AgentName << ", t=" << t << endl;};
        }; // closes if
        // We now calculate the return at time $t-T$ as the difference between the forecast and its realization
        if (t-Future<=0) {if (Plot=="On") {outputLog << "We are not in the case t>Future (t=" << t << ", Future=" << Future << ")" << endl;};}
        else if (t-Future>0) {
            if (Plot=="On") {
                outputLog << "We are in the case t=" << t << ">" << Future << "=Future (ForecastReal.size()=" << ForecastReal.size() << ")" << endl;
                outputLog << "At time t-Future=" << t-Future << ", real forecast was $" << ForecastReal[j][max(0,int(ForecastReal[j].size()) - Future - 1)] << " and at present time t=" << t << ", market price is $" << ReflexiveVal_t << endl;
                if (Exploration==4) {outputLog << " and virtual forecast was $" << ForecastVirtual[j][max(0,int(ForecastVirtual[j].size()) - Future - 1)] << endl;};
            }; // closes Plot condition
            double dFResultReal=0; double dFResultVirtual=0;
            //dFResultReal=-100*abs(ReflexiveVal_t - ForecastReal[int(ForecastReal.size()) - Future - 1])/ReflexiveVal_t; // minus the percentage of the absolute value of miss compared to market price.
            //dFResultVirtual=-100*abs(ReflexiveVal_t - ForecastVirtual[int(ForecastVirtual.size()) - Future - 1])/ReflexiveVal_t; // minus the percentage of the absolute value of miss compared to market price.
            dFResultReal=100*abs(gsl_matrix_get(ReflexiveValues, j, t) - ForecastReal[j][max(0,int(ForecastReal[j].size()) - Future - 1)])/gsl_matrix_get(ReflexiveValues, j, t); // minus the percentage of the absolute value of miss compared to market price.
            if (Exploration==4) {dFResultVirtual=100*abs(gsl_matrix_get(ReflexiveValues, j, t) - ForecastVirtual[j][max(0,int(ForecastVirtual[j].size()) - Future - 1)])/gsl_matrix_get(ReflexiveValues, j, t);}; // minus the percentage of the absolute value of miss compared to market price.
            //FResultReal.push_back(dFResultReal); // Recording result of forecast (ABM#6 for Forecast)
            //if (int(FResultReal.size()) - Future - 1 > 0) {FResultReal.erase (FResultReal.begin());}; // OOO
            //FResultVirtual.push_back(dFResultVirtual); // Recording result of forecast (ABM#6 for Forecast)
            //if (int(FResultVirtual.size()) - Future - 1 > 0) {FResultVirtual.erase (FResultVirtual.begin());}; // OOO
            
            // Finding the optimal action at time t-Future
            //if (OffPolicy2==1) // XXXExploration5
            
            
            
            if (Plot=="On") {
                outputLog << "So dFResultReal=" << dFResultReal << "%" << endl;
                if (Exploration==4) {outputLog << "dFResultVirtual=" << dFResultVirtual << "%" << endl;};
                //outputLog << "FResultReal[]={" ; for (int u=0; u<int(FResultReal.size()); u++) {if (u<int(FResultReal.size())-1) {outputLog << FResultReal[u] << ", ";} else {outputLog << FResultReal[u];};}; outputLog << "}" << endl;
                //outputLog << "FResultVirtual[]={" ; for (int u=0; u<int(FResultVirtual.size()); u++) {if (u<int(FResultVirtual.size())-1) {outputLog << FResultVirtual[u] << ", ";} else {outputLog << FResultVirtual[u];};}; outputLog << "}" << endl;
                
                
                /*************************************ABM I.8******************************************/
                outputLog << endl << "    ABM I.8: Returns above now discretized according to past returns (cf. reinforcement comparison & actor-critic methods)" << endl;
            }; // closes Plot condition
            // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 1.7 computed..." << endl;
            
            // We now discretise this result by giving it the following score -2, -1, +1, or +2 depending on where it stands in all past rewards
            
            double dFResultPercentile=0; int dFResultDisReal=0; int dFResultVecSize=int(dFResultVec[j].size());
            if (dFResultReal<0.001) {dFResultReal=0.001;}; // Failsafe
            for (int i=0; i<dFResultVecSize; i++) { // sorting the vector SmoothVec by ascending values
                if (dFResultReal<dFResultVec[j][i]) {dFResultVec[j].insert(dFResultVec[j].begin()+i, dFResultReal*(1 + 0.0001*VRan[25])); dFResultPercentile=i*1.0/dFResultVecSize; break;};
                if (i==dFResultVecSize-1) {dFResultVec[j].push_back(dFResultReal*(1 + 0.0001*VRan[25])); dFResultPercentile=1; break;};
            }; // closes for loop
            if (dFResultPercentile<0.05) {dFResultDisReal=+4;}
            else if ((dFResultPercentile>=0.05) && (dFResultPercentile<0.25)) {dFResultDisReal=2;}
            else if ((dFResultPercentile>=0.25) && (dFResultPercentile<0.5)) {dFResultDisReal=1;}
            else if ((dFResultPercentile>=0.5) && (dFResultPercentile<0.75)) {dFResultDisReal=-1;}
            else if ((dFResultPercentile>=0.75) && (dFResultPercentile<0.95)) {dFResultDisReal=-2;}
            else if (dFResultPercentile>=0.95) {dFResultDisReal=-4;};
            if (dFResultVecSize==1) {dFResultDisReal=1;}; // In the beginning it means nothing
            if ((dFResultReal<5) && (dFResultDisReal<0)) {dFResultDisReal=1;}; // Making sure a forecast smaller than 5% mismatch is not accounted as a negative reward XYZ333
            
            if (dFResultVecSize - Future - History >0) {dFResultVec[j].erase (dFResultVec[j].begin());}; // OOO
            FResultDisReal[j].push_back(dFResultDisReal); // LLL
            //if (int(FResultDisReal.size()) - Future - 1 >0) {FResultDisReal.erase (FResultDisReal.begin());}; // OOO
            
            
            if (Plot=="On") {
                //outputLog << "FRMean=" << FRMean << ", FRVariance=" << FRVariance << ", FRStdDev=" << FRStdDev << endl;
                //outputLog << "BBMinus=" << BBMinus << ", BBPlus=" << BBPlus << endl;
                outputLog << "dFResultDisReal=" << dFResultDisReal << " (dFResultReal=" << dFResultReal << ")" << endl;
                outputLog << "FResultDisReal[j][]={" ; for (int u=0; u<int(FResultDisReal[j].size()); u++) {if (u<int(FResultDisReal[j].size())-1) {outputLog << FResultDisReal[j][u] << ", ";} else {outputLog << FResultDisReal[j][u];};}; outputLog << "}" << endl;
            }; // closes Plot condition
            
            int dFResultDisVirtual=0;
            if (Exploration==4) {
                dFResultPercentile=0; dFResultVecSize=int(dFResultVec[j].size());
                if (dFResultVirtual<0.001) {dFResultVirtual=0.001;}; // Failsafe
                for (int i=0; i<dFResultVecSize; i++) { // sorting the vector by ascending values
                    if (dFResultVirtual<dFResultVec[j][i]) {dFResultVec[j].insert(dFResultVec[j].begin()+i, dFResultVirtual*(1 + 0.0001*VRan[26])); dFResultPercentile=i*1.0/dFResultVecSize; break;};
                    if (i==dFResultVecSize-1) {dFResultVec[j].push_back(dFResultVirtual*(1 + 0.0001*VRan[26])); dFResultPercentile=1; break;};
                }; // closes for loop
                if (dFResultPercentile<0.05) {dFResultDisVirtual=4;}
                else if ((dFResultPercentile>=0.05) && (dFResultPercentile<0.25)) {dFResultDisVirtual=2;}
                else if ((dFResultPercentile>=0.25) && (dFResultPercentile<0.5)) {dFResultDisVirtual=1;}
                else if ((dFResultPercentile>=0.5) && (dFResultPercentile<0.75)) {dFResultDisVirtual=-1;}
                else if ((dFResultPercentile>=0.75) && (dFResultPercentile<0.95)) {dFResultDisVirtual=-2;}
                else if (dFResultPercentile>=0.95) {dFResultDisVirtual=-4;};
                if (dFResultVecSize==1) {dFResultDisVirtual=1;}; // In the beginning it means nothing
                if (dFResultVecSize - Future - History >0) {dFResultVec[j].erase (dFResultVec[j].begin());}; // OOO
                //if (int(FResultDisVirtual.size()) - Future - 1 >0) {FResultDisVirtual.erase (FResultDisVirtual.begin());}; // OOO
                FResultDisVirtual[j].push_back(dFResultDisVirtual); // LLL
                if (Plot=="On") {
                    //outputLog << "FRMeanVirtual=" << FRMeanVirtual << ", FRVarianceVirtual=" << FRVarianceVirtual << ", FRStdDevVirtual=" << FRStdDevVirtual << endl;
                    //outputLog << "BBMinusVirtual=" << BBMinusVirtual << ", BBPlusVirtual=" << BBPlusVirtual << endl;
                    outputLog << "dFResultDisVirtual=" << dFResultDisVirtual << " (dFResultVirtual=" << dFResultVirtual << ")" << endl;
                    outputLog << "FResultDisVirtual[j][]={" ; for (int u=0; u<int(FResultDisVirtual[j].size()); u++) {if (u<int(FResultDisVirtual[j].size())-1) {outputLog << FResultDisVirtual[j][u] << ", ";} else {outputLog << FResultDisVirtual[j][u];};}; outputLog << "}" << endl;
                }; // closes Plot condition
            }; // closes if
            
            if ((NEBPositivity>VRan[17]) && (dFResultDisReal>0)) {
                double FirstdFResultDisReal=dFResultDisReal;
                dFResultDisReal+=1; // Positivity bias: better recalling pleasant memories than unpleasant ones
                if (Plot=="On") {outputLog << "NEBPositivity=" << NEBPositivity << " : positivity bias switched on => dFResultDisReal=" << FirstdFResultDisReal << "->" << dFResultDisReal << endl;};
            };
            
            if ((NEBNegativity>VRan[17]) && (dFResultDisReal<0)) {
                double FirstdFResultDisReal=dFResultDisReal;
                dFResultDisReal+=1; // Positivity bias: better recalling pleasant memories than unpleasant ones
                if (Plot=="On") {outputLog << "NEBNegativity=" << NEBNegativity << " : negativity bias switched on => dFResultDisReal=" << FirstdFResultDisReal << "->" << dFResultDisReal << endl;};
            };
            
            int dFResultDis5=4;
            if (OffPolicy2==1) {
                if (Plot=="On") {outputLog << "dFResultDis5=" << 4 << " (always) corresponding to a mismatch MinimumCostFunction=100*abs(P(t)-ˆP(t-T))/P(t)=" << MinimumCostFunction << "%" << endl;};
            };
            
            /*************************************ABM I.9******************************************/
            if (Plot=="On") {outputLog << endl << "    ABM I.9: Updating Q(s,a) via SAM according to these discretized returns" << endl;};
            // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 1.8 computed..." << endl;
            
            // We now update the action-value function $Q(s,a)$ based on this return via the SAM or Sample Average Method (ABM#7)
            vector<int> OldStates = get_tensor_coord (FSIndex[j][max(0,int(FSIndex[j].size()) - Future - 1)], 3, FSDim); // Coord in s-tensor of past state $s_{t-T}$
            vector<int> OldRealActions = get_tensor_coord (FAIndexReal[j][max(0,int(FAIndexReal[j].size()) - Future - 1)], 3, FADim); //Coord in a-tensor of past action $a_{t-T}$
            vector<int> OldVirtualActions, OldActions5;
            int xxx=max(0, int(FResultDisReal[j].size()) - Future - 1);
            dFResultDisReal=FResultDisReal[j][xxx]; // Past value t-Future time steps ago that will be used to update both Q() and pi() // TTT1
            if (Exploration==4) {
                OldVirtualActions = get_tensor_coord (FAIndexVirtual[j][max(0,int(FAIndexVirtual[j].size()) - Future - 1)], 3, FADim); //Coord in a-tensor of past action $a_{t-T}$
                dFResultDisVirtual=FResultDisVirtual[j][max(0,int(FResultDisVirtual[j].size()) - Future - 1)]; // Past value t-Future time steps ago that will be used to update both Q() and pi() // TTT1
            };
            FAInd5[0]=OptimalTool5; FAInd5[1]=OptimalLag5; FAInd5[2]=OptimalRWeight5; int ATensorIndex5 = 0;
            if (OffPolicy2==1) {
                ATensorIndex5 = get_tensor_index (3, FAInd5, FADim);
                OldActions5 = get_tensor_coord (ATensorIndex5, 3, FADim);
            };
            
            if (Plot=="On") {
                outputLog << "OldStates[]={" ; for (int u=0; u<int(OldStates.size()); u++) {if (u<int(OldStates.size())-1) {outputLog << OldStates[u] << ", ";} else {outputLog << OldStates[u];};}; outputLog << "}" << endl;
                outputLog << "OldRealActions[]={" ; for (int u=0; u<int(OldRealActions.size()); u++) {if (u<int(OldRealActions.size())-1) {outputLog << OldRealActions[u] << ", ";} else {outputLog << OldRealActions[u];};}; outputLog << "}" << endl;
                if (Exploration==4) {outputLog << "OldVirtualActions[]={" ; for (int u=0; u<int(OldVirtualActions.size()); u++) {if (u<int(OldVirtualActions.size())-1) {outputLog << OldVirtualActions[u] << ", ";} else {outputLog << OldVirtualActions[u];};}; outputLog << "}" << endl;};
                if (OffPolicy2==1) {outputLog << "ATensorIndex5=" << ATensorIndex5 << endl;};
                outputLog << "FQ[], FNumberA[], initiated with size " << FNumberA.size() << " (=FS*FA)" << endl;
                outputLog << "We are at t=Future=" << t << endl;
            }; // closes Plot condition
            FQIndReal[0]=OldStates[0]; FQIndReal[1]=OldStates[1]; FQIndReal[2]=OldStates[2];
            FQIndReal[3]=OldRealActions[0]; FQIndReal[4]=OldRealActions[1]; FQIndReal[5]=OldRealActions[2];
            int QRealTensorIndex = 0;
            QRealTensorIndex = get_tensor_index (6, FQIndReal, FQDim); // Vector index of Q-tensor corresponding to $s_{t-T}$ x $a_{t-T}$
            int QVirtualTensorIndex=0;
            if (Exploration==4) {
                FQIndVirtual[0]=OldStates[0]; FQIndVirtual[1]=OldStates[1]; FQIndVirtual[2]=OldStates[2]; FQIndVirtual[3]=OldVirtualActions[0]; FQIndVirtual[4]=OldVirtualActions[1]; FQIndVirtual[5]=OldVirtualActions[2];
                QVirtualTensorIndex = get_tensor_index (6, FQIndVirtual, FQDim); // Vector index of Q-tensor corresponding to $s_{t-T}$ x $a_{t-T}$
            };
            int QTensorIndex5=0;
            if (OffPolicy2==1) {
                FQInd5[0]=OldStates[0]; FQInd5[1]=OldStates[1]; FQInd5[2]=OldStates[2]; FQInd5[3]=OldActions5[0]; FQInd5[4]=OldActions5[1]; FQInd5[5]=OldActions5[2];
                QTensorIndex5 = get_tensor_index (6, FQInd5, FQDim); // Vector index of Q-tensor corresponding to $s_{t-T}$ x $a_{t-T}$
            };
            if (Plot=="On") {
                outputLog << "FQIndReal[]={" ; for (int u=0; u<6; u++) {if (u<6-1) {outputLog << FQIndReal[u] << ", ";} else {outputLog << FQIndReal[u];};}; outputLog << "}" << endl;
                if (Exploration==4) {outputLog << "FQIndVirtual[]={" ; for (int u=0; u<6; u++) {if (u<6-1) {outputLog << FQIndVirtual[u] << ", ";} else {outputLog << FQIndVirtual[u];};}; outputLog << "}" << endl;};
                outputLog << "FQDim[]={" ; for (int u=0; u<6; u++) {if (u<6-1) {outputLog << FQDim[u] << ", ";} else {outputLog << FQDim[u];};}; outputLog << "}" << endl;
                outputLog << "QRealTensorIndex=" << QRealTensorIndex << endl;
                if (Exploration==4) {outputLog << "QVirtualTensorIndex=" << QVirtualTensorIndex << endl;};
                if (OffPolicy2==1) {outputLog << "QTensorIndex5=" << QTensorIndex5 << endl;};
            }; // closes Plot condition
            
            FNumberA[j][QRealTensorIndex]+=1; // Adding this action as taken in the ActionNumber tensor so as to count for SAM
            if (FNumberA[j][QRealTensorIndex]==1) {FQ[j][QRealTensorIndex]=dFResultReal;} // No SAM if not yet populated (i.e =0)
            else {FQ[j][QRealTensorIndex]=(FQ[j][QRealTensorIndex]*(FNumberA[j][QRealTensorIndex]-1) + dFResultReal)/FNumberA[j][QRealTensorIndex];}; // SAM
            if (Plot=="On") {
                outputLog << "Number of times this real action was previously taken in that former state: FNumberA[j][QRealTensorIndex]-1=" << FNumberA[j][QRealTensorIndex]-1 << ", so FQ[j][QRealTensorIndex]=" << FQ[j][QRealTensorIndex] << endl;
            }; // closes plot condition
            
            if (Exploration==4) {
                FNumberA[j][QVirtualTensorIndex]+=1; // Adding this action as taken in the ActionNumber tensor so as to count for SAM
                if (FNumberA[j][QVirtualTensorIndex]==1) {FQ[j][QVirtualTensorIndex]=dFResultVirtual;} // No SAM if not yet populated
                else {FQ[j][QVirtualTensorIndex]=(FQ[j][QVirtualTensorIndex]*(FNumberA[j][QVirtualTensorIndex]-1) + dFResultVirtual)/FNumberA[j][QVirtualTensorIndex];}; // SAM
                if (Plot=="On") {
                    outputLog << "Number of times this virtual action was previously taken in that former state: FNumberA[j][QVirtualTensorIndex]-1=" << FNumberA[j][QVirtualTensorIndex]-1 << ", so FQ[j][QVirtualTensorIndex]=" << FQ[j][QVirtualTensorIndex] << endl;
                }; // closes plot condition
            }; // closes if
            
            if (OffPolicy2==1) {
                FNumberA[j][QTensorIndex5]+=1; // Adding this action as taken in the ActionNumber tensor so as to count for SAM
                if (FNumberA[j][QTensorIndex5]==1) {FQ[j][QTensorIndex5]=MinimumCostFunction;} // No SAM if not yet populated
                else {FQ[j][QTensorIndex5]=(FQ[j][QTensorIndex5]*(FNumberA[j][QTensorIndex5]-1) + MinimumCostFunction)/FNumberA[j][QTensorIndex5];}; // SAM
                if (Plot=="On") {
                    outputLog << "Number of times this Exploration5 action was previously taken in that former state: FNumberA[j][QTensorIndex5]-1=" << FNumberA[j][QTensorIndex5]-1 << ", so FQ[j][QTensorIndex5]=" << FQ[j][QTensorIndex5] << endl;
                }; // closes plot condition
            }; // closes if
            
            
            
            /*************************************ABM I.10******************************************/
            if (Plot=="On") {outputLog << endl << "    ABM I.10: Updating pi(s,a) according to these discretized returns" << endl;};
            // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 1.9 computed..." << endl;
            
            // We now update existing policy based on this discretized results FResultDisReal and FResultDisVirtual (ABM#8)
            FSInd[0]=OldStates[0]; FSInd[1]=OldStates[1]; FSInd[2]=OldStates[2];
            int OldSTensorIndex = get_tensor_index (3, FSInd, FSDim); // Vector index of the s-tensor associated with $s_{t-T}$
            int OldSCount=0;
            for (int i=0; i<int(FSIndex[j].size())-Future; i++) {if (FSIndex[j][i]==OldSTensorIndex) {OldSCount+=1;};}; // Nb of times $s_{t-T}$ seen
            if (Plot=="On") {
                outputLog << "State at time t-Future: FSInd[0]=OldStates[0]=" << FSInd[0] << ", FSInd[1]=OldStates[1]=" << FSInd[1] << ", FSInd[2]=OldStates[2]=" << FSInd[2] << endl;
                outputLog << "OldSTensorIndex=" << OldSTensorIndex << endl;
                outputLog << "Number of times state at FSIndex.size()-Future encountered: OldSCount=" << OldSCount << endl;
            }; // closes Plot condition
            FQIndReal[0]=OldStates[0]; FQIndReal[1]=OldStates[1]; FQIndReal[2]=OldStates[2]; // Using results from above on Q
            FQIndReal[3]=0; FQIndReal[4]=0; FQIndReal[5]=0; // As start of the Pi vector index to span
            int PiSTensorIndex = get_tensor_index (6, FQIndReal, FQDim); // Vector index of the Q-tensor associated with $s_{t-T}$
            if (Plot=="On") {
                outputLog << "FQIndReal[0]=OldStates[0]=" << FQIndReal[0] << ", FQIndReal[1]=OldStates[1]=" << FQIndReal[1] << ", FQIndReal[2]=OldStates[2]=" << FQIndReal[2] << endl;
                outputLog << "FQIndReal[3]=" << FQIndReal[3] << ", FQIndReal[4]=" << FQIndReal[4] << ", FQIndReal[5]=" << FQIndReal[5] << ", which is PiSTensorIndex=" << PiSTensorIndex << " from which we span all actions to update their probability" << endl;
                //outputLog << "State s was encountered OldSCount=" << OldSCount << " times" << endl;
                outputLog << "Results was dFResultDisReal=" << dFResultDisReal << endl;
                if (Exploration==4) {outputLog << "And dFResultDisVirtual=" << dFResultDisVirtual << endl;};
                if (OffPolicy2==1) {outputLog << "And dFResultDis5=" << dFResultDis5 << endl;};
            }; // closes Plot condition
            double Sum=0;
            for (int i=PiSTensorIndex; i<PiSTensorIndex+FA; i++) {
                //if (Plot=="On") {outputLog << "Old FPi[j][" << i << "]=" << floor(100*FPi[j][i]) << "%" << endl;};
                Sum+=FPi[j][i];
            };
            if (Plot=="On") {
                outputLog << "Sum=" << 100*Sum << endl;
                if ((int(100*Sum)<99) || (int(100*Sum)>101)) {outputLog << "ErrorOriginalF[]" << endl;};
            };
            Sum=0;
            
            //if (1-NEB8>VRan[19]) {
            // First for real
            for (int i=PiSTensorIndex; i<PiSTensorIndex+FA; i++) {
                if (Plot=="On") {outputLog << "Old real FPi[j][" << i << "]=" << 100*FPi[j][i] << "% => ";};
                if ((i==QRealTensorIndex) && (dFResultDisReal>=0)) {
                    for (int f=0; f<dFResultDisReal; f++) {FPi[j][i]=FPi[j][i]*(1-NEBLearningRate/FNumberA[j][QRealTensorIndex]) + NEBLearningRate/FNumberA[j][QRealTensorIndex];}; // DDD3
                } // closes if
                else if ((i!=QRealTensorIndex) && (dFResultDisReal>=0)) {
                    for (int f=0; f<dFResultDisReal; f++) {FPi[j][i]=FPi[j][i]*(1-NEBLearningRate/FNumberA[j][QRealTensorIndex]);};
                } // closes else if
                else if ((i==QRealTensorIndex) && (dFResultDisReal<0)) {
                    for (int f=0; f<abs(dFResultDisReal); f++) {FPi[j][i]=FPi[j][i]*(1-NEBLearningRate/FNumberA[j][QRealTensorIndex]);};
                } // closes else if
                else if ((i!=QRealTensorIndex) && (dFResultDisReal<0)) {
                    for (int f=0; f<abs(dFResultDisReal); f++) {FPi[j][i]=FPi[j][i]*(1-NEBLearningRate/FNumberA[j][QRealTensorIndex]) + NEBLearningRate/(FA-1)/FNumberA[j][QRealTensorIndex];};
                }; // closes else if
                if (Plot=="On") {outputLog << "New real FPi[j][" << i << "]=" << 100*FPi[j][i] << "%" << endl;};
                Sum+=FPi[j][i];
            }; // closes i loop
            if (Plot=="On") {
                outputLog << "(QRealTensorIndex=" << QRealTensorIndex << " was the real action update), FNumberA[j][QRealTensorIndex]=" << FNumberA[j][QRealTensorIndex] << endl;
                outputLog << "Sum=" << 100*Sum << endl;
                if ((int(100*Sum)<99) || (int(100*Sum)>101)) {outputLog << "ErrorRealF[]: For OldSCount=" << OldSCount << ", dFResultDisReal=" << dFResultDisReal << endl;};
            }; // closes Plot condition
            Sum=0;
            //}; // closes if
            // Second for virtual
            //if ((Exploration==4) && (1-NEB8>VRan[19])) { // This is for off-policy Watkin's Q-learning
            if (Exploration==4) { // This is for off-policy Watkin's Q-learning
                for (int i=PiSTensorIndex; i<PiSTensorIndex+FA; i++) {
                    if (Plot=="On") {outputLog << "Old virtual FPi[j][" << i << "]=" << 100*FPi[j][i] << "% => ";};
                    if ((i==QVirtualTensorIndex) && (dFResultDisVirtual>=0)) {
                        for (int f=0; f<dFResultDisVirtual; f++) {FPi[j][i]=FPi[j][i]*(1-NEBLearningRate/FNumberA[j][QVirtualTensorIndex]) + NEBLearningRate/FNumberA[j][QVirtualTensorIndex];};
                    } // closes if
                    else if ((i!=QVirtualTensorIndex) && (dFResultDisVirtual>=0)) {
                        for (int f=0; f<dFResultDisVirtual; f++) {FPi[j][i]=FPi[j][i]*(1-NEBLearningRate/FNumberA[j][QVirtualTensorIndex]);};
                    } // closes else if
                    else if ((i==QVirtualTensorIndex) && (dFResultDisVirtual<0)) {
                        for (int f=0; f<abs(dFResultDisVirtual); f++) {FPi[j][i]=FPi[j][i]*(1-NEBLearningRate/FNumberA[j][QVirtualTensorIndex]);};
                    } // closes else if
                    else if ((i!=QVirtualTensorIndex) && (dFResultDisVirtual<0)) {
                        for (int f=0; f<abs(dFResultDisVirtual); f++) {FPi[j][i]=FPi[j][i]*(1-NEBLearningRate/FNumberA[j][QVirtualTensorIndex]) + NEBLearningRate/(FA-1)/FNumberA[j][QVirtualTensorIndex];};
                    }; // closes else if
                    if (Plot=="On") {outputLog << "New virtual FPi[j][" << i << "]=" << 100*FPi[j][i] << "%" << endl;};
                    Sum+=FPi[j][i];
                }; // closes i loop
                if (Plot=="On") {
                    outputLog << "(QVirtualTensorIndex=" << QVirtualTensorIndex << " was the real action update), FNumberA[j][QVirtualTensorIndex]=" << FNumberA[j][QVirtualTensorIndex] << endl;
                    outputLog << "Sum=" << 100*Sum << endl; if ((int(100*Sum)<99) || (int(100*Sum)>101)) {outputLog << "ErrorVirtualF[]: For OldSCount=" << OldSCount << ", dFResultDisVirtual=" << dFResultDisVirtual << endl;};
                }; // closes Plot condition
                Sum=0;
            }; // closes if Exploration
            // Second for Exploration5
            if (OffPolicy2==1) { // This is for off-policy Watkin's Q-learning
                for (int i=PiSTensorIndex; i<PiSTensorIndex+FA; i++) {
                    if (Plot=="On") {outputLog << "Old Exploration5 FPi[j][" << i << "]=" << 100*FPi[j][i] << "% => ";};
                    if (i==QTensorIndex5) {
                        for (int f=0; f<dFResultDis5; f++) {FPi[j][i]=FPi[j][i]*(1-NEBLearningRate/FNumberA[j][QTensorIndex5]) + NEBLearningRate/FNumberA[j][QTensorIndex5];};
                    } // closes if
                    else if (i!=QTensorIndex5) {
                        for (int f=0; f<dFResultDis5; f++) {FPi[j][i]=FPi[j][i]*(1-NEBLearningRate/FNumberA[j][QTensorIndex5]);};
                    } // closes else if
                    if (Plot=="On") {outputLog << "New Exploration5 FPi[j][" << i << "]=" << 100*FPi[j][i] << "%" << endl;};
                    Sum+=FPi[j][i];
                }; // closes i loop
                if (Plot=="On") {
                    outputLog << "(QTensorIndex5=" << QTensorIndex5 << " was the real action update), FNumberA[j][QTensorIndex5]=" << FNumberA[j][QTensorIndex5] << endl;
                    outputLog << "Sum=" << 100*Sum << endl; if ((int(100*Sum)<99) || (int(100*Sum)>101)) {outputLog << "Error5F[]: For OldSCount=" << OldSCount << ", dFResultDis5=" << dFResultDis5 << endl;};
                }; // closes Plot condition
                Sum=0;
            }; // closes if Exploration
        }; // closes if (t-Future>0)
        // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 1.10 computed..." << endl;
        
        
        
        
        
        
        /*******************************************************************************/
        /*******************************************************************************/
        /*******************************************************************************/
        /***********************************RL TRADING**********************************/
        /*******************************************************************************/
        /*******************************************************************************/
        /*******************************************************************************/
        
        /*************************************ABM II.1**************************************************/
        if (Plot=="On") {outputLog << endl << endl << "    ABM 2.1: Defining present state s" << endl;};
        // Griding growth forecast results for RL states s
        double dMu=100*(ForecastReference-gsl_matrix_get(ReflexiveValues, j, t))/(gsl_matrix_get(ReflexiveValues, j, t));
        int dMuPosVecSize=int(dMuPosVec[j].size()); int dMuNegVecSize=int(dMuNegVec[j].size());
        double dMuPosPercentile=0; double dMuNegPercentile=0;
        if (t==1) {dMuPosVec[j].push_back(0); dMuPosPercentile=0.5; dMuNegVec[j].push_back(0); dMuNegPercentile=0.5;}
        else if ((t>1) && (dMu>=0)) {
            for (int i=0; i<dMuPosVecSize; i++) {
                if (dMu<dMuPosVec[j][i]) {dMuPosVec[j].insert(dMuPosVec[j].begin()+i, dMu*(1 + 0.0001*VRan[30])); dMuPosPercentile=i*1.0/dMuPosVecSize; break;};
                if (i==dMuPosVecSize-1) {dMuPosVec[j].push_back(dMu*(1 + 0.0001*VRan[30])); dMuPosPercentile=1; break;};
            }; // closes for loop
            
            /*
             if (dMuPosPercentile<0.05) {Mu=3;}
             else if ((dMuPosPercentile>=0.05) && (dMuPosPercentile<0.33)) {Mu=4;}
             else if ((dMuPosPercentile>=0.33) && (dMuPosPercentile<0.67)) {Mu=5;}
             else if (dMuPosPercentile>=0.67) {Mu=6;};
             */
            
            if (dMuPosPercentile<0.05) {Mu=1;}
            else if (dMuPosPercentile>=0.05) {Mu=2;};
            // NEB4 implementation
            //if ((NEB4p>VRan[16]) && (NEB4==0)) {if (dMuPosPercentile>=0.025) {Mu=2;};}; // Regressive conservatism: overestimating high values & likelihoods while underestimating low values & likelihoods
            //if ((NEB4p>VRan[16]) && (NEB4==1)) {if (dMuPosPercentile<=0.075) {Mu=1;};}; // Regressive conservatism: underestimating high values & likelihoods while overestimating low values & likelihoods
            
        } // closes else if
        else if ((t>1) && (dMu<0)) {
            for (int i=0; i<dMuNegVecSize; i++) {
                if (dMu<dMuNegVec[j][i]) {dMuNegVec[j].insert(dMuNegVec[j].begin()+i, dMu*(1 + 0.0001*VRan[30])); dMuNegPercentile=i*1.0/dMuNegVecSize; break;};
                if (i==dMuNegVecSize-1) {dMuNegVec[j].push_back(dMu*(1 + 0.0001*VRan[30])); dMuNegPercentile=1; break;};
            }; // closes for loop
            
            /*
             if (dMuNegPercentile<0.33) {Mu=0;}
             else if ((dMuNegPercentile>=0.33) && (dMuNegPercentile<0.67)) {Mu=1;}
             else if ((dMuNegPercentile>=0.67) && (dMuNegPercentile<0.95)) {Mu=2;}
             else if (dMuNegPercentile>=0.95) {Mu=3;};
             */
            
            if (dMuNegPercentile<=0.95) {Mu=0;}
            else if (dMuNegPercentile>0.95) {Mu=1;};
            // NEB4 implementation
            //if ((NEB4p>VRan[16]) && (NEB4==0)) {if (dMuNegPercentile<=0.975) {Mu=0;}}; // Regressive conservatism: overestimating high values & likelihoods while underestimating low values & likelihoods
            //if ((NEB4p>VRan[16]) && (NEB4==1)) {if (dMuNegPercentile>=0.925) {Mu=1;};}; // Regressive conservatism: underestimating high values & likelihoods while overestimating low values & likelihoods
            
        }; // closes else if
        if (int(dMuPosVec[j].size()) - Future - History >0) {dMuPosVec[j].erase (dMuPosVec[j].begin());}; // OOO
        if (int(dMuNegVec[j].size()) - Future - History >0) {dMuNegVec[j].erase (dMuNegVec[j].begin());}; // OOO
        
        // NEB10 implementation
        //if ((NEB10>VRan[18]) && (Mu<2)) {Mu+=1;
        //    if (Plot=="On") {outputLog << "NEB10=" << NEB10 << " => DelayDiscounting : Mu=" << Mu-1 << "->" << Mu << endl;};
        //};
        
        
        // Now griding volatility forecast results for RL states s
        //double dSig=100*Lvol/(gsl_matrix_get(ReferenceValues, j, t));
        double dSig=Lvol;
        int dSigVecSize=int(dSigVec[j].size()); double dSigPercentile=0;
        if (t==1) {dSigVec[j].push_back(dSig); dSigPercentile=0.5;}
        else {
            for (int i=0; i<dSigVecSize; i++) {
                if (dSig<dSigVec[j][i]) {dSigVec[j].insert(dSigVec[j].begin()+i, dSig*(1 + 0.0001*VRan[31])); dSigPercentile=i*1.0/dSigVecSize; break;};
                if (i==dSigVecSize-1) {dSigVec[j].push_back(dSig*(1 + 0.0001*VRan[31])); dSigPercentile=1; break;};
            }; // closes for loop
            if (dSigPercentile<0.33) {Sig=0;}
            else if ((dSigPercentile>=0.33) && (dSigPercentile<0.67)) {Sig=1;}
            else if (dSigPercentile>=0.67) {Sig=2;};
        } // closes else
        if (int(dSigVec[j].size()) - Future - History >0) {dSigVec[j].erase (dSigVec[j].begin());}; // OOO
        
        
        
        
        if (t==1) {RFAFirst=RFA;}; // To make sure that it is the market price at t=1 that acts as standard
        if ((RFAFirst==0) && (t>=2*Week)) {RFAFirst=RFA;}; // Failsafe if agent has RFAFirst=0
        /*
         if ((RFA>=0) && (RFA<0.25*RFAFirst)) {RFALevel=0;}
         else if ((RFA>=0.25*RFAFirst) && (RFA<0.6*RFAFirst)) {RFALevel=1;}
         else if (RFA>=0.6*RFAFirst) {RFALevel=2;};
         */
        if ((RFA>=0) && (RFA<0.6*RFAFirst)) {RFALevel=0;}
        else if (RFA>=0.6*RFAFirst) {RFALevel=1;};
        
        // NEB10 implementation
        //if (NEB10>VRan[18]) {RFALevel=1;
        //if (Plot=="On") {outputLog << "NEB10=" << NEB10 << " => DelayDiscounting : RFALevel set at 1" << endl;};
        //};
        
        
        if (Plot=="On") {RFALog[j].push_back(RFALevel);}; // LLL
        // Now griding RBA since it is a parameter for shorting stocks
        RBA=StockHoldings();
        if (t==1) {RBAFirst=RBA;}; // To make sure that it is the market price at t=1 that acts as standard
        if ((RBAFirst==0) && (t>=2*Week)) {RBAFirst=RBA;}; // Failsafe if agent has Q=0 stock holdings and hence RBAFirst=0
        /*
         if ((RBA>=0) && (RBA<0.25*RBAFirst)) {RBALevel=0;}
         else if ((RBA>=0.25*RBAFirst) && (RBA<0.6*RBAFirst)) {RBALevel=1;}
         else if (RBA>=0.6*RBAFirst) {RBALevel=2;};
         */
        if ((RBA>=0) && (RBA<0.6*RBAFirst)) {RBALevel=0;}
        else if (RBA>=0.6*RBAFirst) {RBALevel=1;};
        
        // NEB10 implementation
        //if (NEB10>VRan[18]) {RBALevel=1;
        //if (Plot=="On") {outputLog << "NEB10=" << NEB10 << " => DelayDiscounting : RBALevel set at 1" << endl;};
        //};
        
        if (Plot=="On") {RBALog[j].push_back(RBALevel);}; // LLL
        // Now griding Liquid as it is a parameter for the Law of Supply and Offer and Gesture parameter
        int LiquidPercentVecSize=int(LiquidPercentVec[j].size()); double LiquidPercentPercentile=0;
        if (t==1) {LiquidPercentVec[j].push_back(LiquidPercent); LiquidPercentPercentile=0.5;}
        else {
            for (int i=0; i<LiquidPercentVecSize; i++) { // Important: LiquidPercentVec built out of non-zero values of LiquidPercent!
                if ((LiquidPercent<LiquidPercentVec[j][i]) && (LiquidPercent!=0)) {LiquidPercentVec[j].insert(LiquidPercentVec[j].begin()+i, LiquidPercent*(1 + 0.0001*VRan[32])); LiquidPercentPercentile=i*1.0/LiquidPercentVecSize; break;};
                if ((i==LiquidPercentVecSize-1) && (LiquidPercent!=0)) {LiquidPercentVec[j].push_back(LiquidPercent*(1 + 0.0001*VRan[32])); LiquidPercentPercentile=1; break;};
            }; // closes for loop
            /*
             if (LiquidPercentPercentile<0.33) {Liquid=1;}
             else if ((LiquidPercentPercentile>=0.33) && (LiquidPercentPercentile<0.67)) {Liquid=2;}
             else if (LiquidPercentPercentile>=0.67) {Liquid=3;};
             if (LiquidPercent==0) {Liquid=0;};
             */
            if (LiquidPercentPercentile<0.33) {Liquid=1;}
            else if (LiquidPercentPercentile>=0.33) {Liquid=2;};
            if (LiquidPercent==0) {Liquid=0;};
        } // closes else
        if (int(LiquidPercentVec[j].size()) - Future - History >0) {LiquidPercentVec[j].erase (LiquidPercentVec[j].begin());}; // OOO
        
        
        if (t==1) {Mu=1; Sig=1; Liquid=0;};
        TSInd[0]=Mu; TSInd[1]=Sig; TSInd[2]=RFALevel; TSInd[3]=RBALevel; TSInd[4]=Liquid;
        int TSid=get_tensor_index (5, TSInd, TSDim);// This gives the state index tensor for our second RL algorithm
        TSIndex[j].push_back(TSid); // Back up
        if (int(TSIndex[j].size()) - Future - 1 >0) {TSIndex[j].erase (TSIndex[j].begin());}; // OOO
        if (Plot=="On") {
            outputLog << "dMu=" << dMu << ", dSig=" << dSig << ", RFA=" << RFA << " (RFAFirst=" << RFAFirst << "), RBA=" << RBA << " (RBAFirst=" << RBAFirst << "), LiquidPercent=" << LiquidPercent << endl;
            outputLog << "Mu=" << Mu << ", Sig=" << Sig << ", RFALevel=" << RFALevel << ", RBALevel=" << RBALevel << ", Liquid=" << Liquid << endl;
            outputLog << "TSDim[]={" ; for (int u=0; u<5; u++) {if (u<5-1) {outputLog << TSDim[u] << ", ";} else {outputLog << TSDim[u];};}; outputLog << "}" << endl;
            outputLog << "TSInd[]={" ; for (int u=0; u<5; u++) {if (u<5-1) {outputLog << TSInd[u] << ", ";} else {outputLog << TSInd[u];};}; outputLog << "}" << endl;
            outputLog << "TSid=" << TSid << endl;
            outputLog << "TSIndex[j][]={" ; for (int u=0; u<int(TSIndex[j].size()); u++) {if (u<int(TSIndex[j].size())-1) {outputLog << TSIndex[j][u] << ", ";} else {outputLog << TSIndex[j][u];};}; outputLog << "}" << endl;
            MuLog[j].push_back(Mu); // LLL
            SigLog[j].push_back(Sig); // LLL
            LiquidityLog[j].push_back(Liquid); // LLL
            
            
            /*************************************ABM II.2******************************************/
            outputLog << endl << "    ABM 2.2: Initializing the start policy pi_0" << endl;
        }; // closes Plot condition
        // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 2.1 computed..." << endl;
        
        if (t==1) {
            TPiInd[0]=1; TPiInd[1]=0; TPiInd[2]=1; TPiInd[3]=1; TPiInd[4]=1; // State s
            TPiInd[5]=1; TPiInd[6]=1; // Actions a
            for (int i=0; i<TS*TA; i++) {TPi[j].push_back(1.0/TA);}; // HHH : All actions have equiprobability
            if (Plot=="On") {outputLog << "Policy initialized... all " << TA << " possible actions are given equiprobability 1/" << TA << "=" << 1.0/TA << endl;};
        }
        else if (t>1) {
            if (Plot=="On") {outputLog << "Policy already initialized..." << endl;};
        };
        if (Plot=="On") {
            outputLog << "(To reduce output size, the PiTensorIndex vector indices and their corresponding TPiInd[] coordinates in pi-tensor space are not displayed)" << endl;
            
            
            /*************************************ABM II.3******************************************/
            outputLog << endl << "    ABM 2.3: Using current policy pi to select real action a=(Quant,Pinch) in present state s" << endl;
        }; // closes Plot condition
        // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 2.2 computed...";
        
        // We now use the current policy to perform an action $a$ in our state $s$. Note: There is no need to use \verb?PiStack?, all we need is when we have a vector of action each with a probability in $[0,1]$, we generate a random number from a uniform distribution: if that number is larger than the first probability it is subtracted from it, and if that new number is still larger than the second probability, we do likewise until we arrive to the probability corresponding to action $a$ to be picked up. This is because the sum of all these probabilities is $1$.
        TTenPiId[0]=Mu; TTenPiId[1]=Sig; TTenPiId[2]=RFALevel; TTenPiId[3]=RBALevel; TTenPiId[4]=Liquid; TTenPiId[5]=0; TTenPiId[6]=0;
        int Tid=0;
        Tid = get_tensor_index (7, TTenPiId, TPiDim); // Vector index of pi-tensor representation corresponding to present s and a=(0,0)
        if (Plot=="On") {
            outputLog << "TTenPiId[]={" ; for (int u=0; u<7; u++) {if (u<7-1) {outputLog << TTenPiId[u] << ", ";} else {outputLog << TTenPiId[u];};}; outputLog << "}" << endl;
            outputLog << "Tid=" << Tid << endl;
        }; // closes Plot condition
        
        double Tix=VRan[6];
        double Tix2=VRan[7];
        Act=0; Pos=Tid;
        if (Plot=="On") {outputLog << "Tix=" << Tix << ", Tix2=" << Tix2 << endl;};
        // outputDebug << ": a, ";
        if (Plot=="On") {outputLog << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 2.2" << endl;};
        
        if (Exploration==1) { // Epsilon-greedy: Choosing always the action with highest probability, but once in a while (with a small probability epsilon) choosing a (uniform) random one (which is a problem: the worst potential action can be taken as likely as the best, with in some cases devastating returns).
            // outputDebug << "b1, ";
            if (Plot=="On") {outputLog << "Epsilon-greedy method selected" << endl;};
            if (Tix<Epsilon) { // Exploration
                Pos=Tid + int(floor(TA*Tix2));
                Act=TPi[j][Pos];
                if (Plot=="On") {outputLog << "Exploration: Chose random action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
            }
            else { // Exploitation
                for (int n=Tid; n<Tid+TA; n++) { // Running over the actions a associated with that state s
                    if (Plot=="On") {outputLog << "For n=" << n << ", TPi[j][n]=" << floor(100*TPi[j][n]) << "%" << endl;};
                    if (TPi[j][n]>=Act) {Pos=n; Act=TPi[j][Pos];}; // We try and find the action with greatest probability (greedy move) PPP
                }; // closes for
                vector<int> PosVec; // Indices of equiprobable actions
                for (int n=Tid; n<Tid+TA; n++) {if (TPi[j][n]==Act) {PosVec.push_back(n);};}; // We record all the indices of equiprobable actions
                Pos=PosVec[int(floor(VRan[24]*PosVec.size()))]; // We choose a random indice among these!
                PosVec.erase(PosVec.begin(),PosVec.end());
                if (Plot=="On") {outputLog << "Exploitation epsilon-greedy: Chose greediest action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
            }; // closes else
        } // closes if
        
        else if (Exploration==2) { // Softmax: Same as above except the random choice in now not (uniform) equiprobable for all actions, but graded according to their estimated value. One method (using a Gibbs/Boltzmann distribution $e^{Q_t(a) / \tau} / \sum_{b=1}^n e^{Q_t(b) / \tau}$, with $\tau$ the temperature) allows to even shift this grading all the way back to epsilon greedy methods (when $\tau \rightarrow 0$).
            // outputDebug << "b2, ";
            if (Plot=="On") {outputLog << "Softmax method selected" << endl;};
            if (Tix<Epsilon) { // Exploration
                Pos=Tid + int(floor(TA*Tix2));
                Act=TPi[j][Pos];
                if (Plot=="On") {outputLog << "Exploration: Chose random action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
            }
            else { // Exploitation
                vector<int> N = Shuffle(TA);
                if (Plot=="On") {outputLog << "From n=" << Tid << " to " << Tid+TA-1 << ": Tix2=" << Tix2 << endl;};
                for (int n=Tid; n<Tid+TA; n++) { // Running over the actions a associated with that state s
                    if (Plot=="On") {outputLog << "For n=" << n << ", TPi[j][n]=" << floor(100*TPi[j][n]) << "%" << endl;};
                    int Newn = N[n-Tid]+Tid;
                    if (Plot=="On") {outputLog << "For n=" << Newn << ", TPi[j][Newn]=" << floor(100*TPi[j][Newn]) << "%: Tix2-=TPi[j][Newn]=" << Tix2 << endl;};
                    Tix2-=TPi[j][Newn];
                    if (Tix2<=0) {Pos=Newn; Act=TPi[j][Pos]; break;}; // We try and find a new action according to its graded probability
                }; // closes for
                if (Plot=="On") {outputLog << "Exploitation softmax: Chose action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
            }; // closes else
        } // closes else if
        
        else if (Exploration==3) { // Pursuit: The probability of selecting the greedy action $a_{t+1}=a_{t+1}^{\ast}$ at next time step is incremented a fraction $\beta$ of the way toward $1$, while the probabilities of selecting the other actions are decremented toward $0$, respectively shown by the following two expressions: $\pi_{t+1}(a_{t+1}^{\ast}) = \pi_{t}(a_{t+1}^{\ast}) + \beta \left[ 1- \pi_{t}(a_{t+1}^{\ast}) \right]$ and $\pi_{t+1}(a) = \pi_{t}(a) + \beta \left[ 0- \pi_{t}(a) \right] , \hspace{5mm} \forall a\neq a_{t+1}^{\ast}$.
            // outputDebug << "b3, ";
            if (Plot=="On") {outputLog << "Pursuit method selected" << endl;};
            if (Tix<Epsilon*(1-(double(t)/Time))) { // Exploration
                Pos=Tid + int(floor(TA*Tix2));
                Act=TPi[j][Pos];
                if (Plot=="On") {outputLog << "Exploration: Chose random action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
            }
            else { // Exploitation
                for (int n=Tid; n<Tid+TA; n++) { // Running over the actions a associated with that state s
                    if (Plot=="On") {outputLog << "For n=" << n << ", TPi[j][n]=" << floor(100*TPi[j][n]) << "%" << endl;};
                    if (TPi[j][n]>=Act) {Pos=n; Act=TPi[j][Pos];}; // We try and find the action with greatest probability (greedy move)
                }; // closes for
                vector<int> PosVec; // Indices of equiprobable actions
                for (int n=Tid; n<Tid+TA; n++) {if (TPi[j][n]==Act) {PosVec.push_back(n);};}; // We record all the indices of equiprobable actions
                Pos=PosVec[int(floor(VRan[24]*PosVec.size()))]; // We choose a random indice among these!
                PosVec.erase(PosVec.begin(),PosVec.end());
                if (Plot=="On") {outputLog << "Exploitation pursuit: Chose greediest action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
            }; // closes else
        } // closes else if
        
        else if (Exploration==4) { // This is for off-policy Watkin's Q-learning
            // outputDebug << "b4, ";
            if (Plot=="On") {outputLog << "Off-policy Watkin's Q-learning method selected" << endl;};
            for (int n=Tid; n<Tid+TA; n++) { // Running over the actions a associated with that state s
                if (Plot=="On") {outputLog << "For n=" << n << ", TPi[j][n]=" << floor(100*TPi[j][n]) << "%" << endl;};
                if (TPi[j][n]>=Act) {Pos=n; Act=TPi[j][Pos];}; // We try and find the action with greatest probability (greedy move)
            }; // closes for
            vector<int> PosVec; // Indices of equiprobable actions
            for (int n=Tid; n<Tid+TA; n++) {if (TPi[j][n]==Act) {PosVec.push_back(n);};}; // We record all the indices of equiprobable actions
            Pos=PosVec[int(floor(VRan[24]*PosVec.size()))]; // We choose a random indice among these!
            PosVec.erase(PosVec.begin(),PosVec.end());
            if (Plot=="On") {outputLog << "Exploitation Off-policy: Chose greediest action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
        }; // closes else if
        // outputDebug << "b, ";
        
        
        int Greenlight=1;
        if (Plot=="On") {outputLog << "Model-based RL switches Greenlight=1 if current state s is good wrt. previous states representations according to arg max_a Q(s,a)" << endl;};
        if (Plot=="On") {outputLog << "TradingWindowClock[j]/TradingWindow=" << TradingWindowClock[j] << "/" << TradingWindow << endl << "Greenlight=" << Greenlight << endl;};
        if ((TradingFrequencyCond=="On") && (t-Future>0)) {
            Greenlight=0; TradingWindowClock[j]+=1;
            if (TradingWindowClock[j]==1) {
                Qs[j].erase(Qs[j].begin(),Qs[j].end()); // Completely updating the model-based RL framework at every new trading window
                double Max=-999999;
                for (int u=0; u<TS*TA; u++) {
                    if (u%TA==0) {Max=-999999;};
                    if (TQ[j][u]>=Max) {Max=TQ[j][u];};
                    if ((u%TA==TA-1) && (Max!=0)) {Qs[j].push_back(Max);};
                    //if (u%TA==TA-1) {Qs[j].push_back(Max);};
                }; // closes u loop
                //string A = "TQ_"; string B=IntToString(t); string C=".txt"; A+=B+C; if (j==0) {PlotSTL(Qs[j], A.c_str(), "0");};
            }; // closes if
            double TSPercentile=0; double QsPresent=-999999; // QsPresent=arg max_a Q(s,a) for present state s
            for (int u=Tid; u<Tid+TA; u++) {if (TQ[j][u]>=QsPresent) {QsPresent=TQ[j][u];};}; // Defining QsPresent
            for (int u=0; u<int(Qs[j].size()); u++) {if (QsPresent>=Qs[j][u]) {TSPercentile+=1;};};
            TSPercentile*=100.0/int(Qs[j].size());
            if (Qs[j].size()==0) {TSPercentile=0;}; // Failsafe
            if (((TradingWindowClock[j]*100.0/TradingWindow)>=TSPercentile) || (Qs[j].size()==0) || (TradingWindowClock[j]==TradingWindow)) {Greenlight=1; TradingWindowClock[j]=0;};
            if (Plot=="On") {outputLog << "Qs[j].size()=" << Qs[j].size() << ", QsPresent=" << QsPresent << endl << "100*t/T=" << (TradingWindowClock[j]*100.0/TradingWindow) << " >= TSPercentile=" << TSPercentile << endl << "TradingWindowClock[j]/TradingWindow=" << TradingWindowClock[j] << "/" << TradingWindow << endl << "Greenlight=" << Greenlight << endl;};
        }; // closes if
        // outputDebug << "c, ";
        
        vector<int> TenPiCoo2 = get_tensor_coord (Pos, 7, TPiDim); // Getting the tensor coordinates from that index
        Quant=TenPiCoo2[5]; // Selection of first action from policy
        Pinch=TenPiCoo2[6]; // Selection of second action from policy
        if (Greenlight==0) {
            Quant=1; TenPiCoo2[5]=Quant;
            if (Plot=="On") {outputLog << "Greenlight is OFF..." << endl;}; // Holding position if the current state is not great enough according to BestStates
        }
        else {if (Plot=="On") {
            outputLog << "Greenlight is ON..." << endl;};
            ExitHorizons.push_back(t+Future); // JJJ4
        };
        
        
        if (NEBLossAversion>VRan[13]) {
            int FirstPinch=Pinch;
            if (Pinch==0) {Pinch=1;} // Loss aversion asking more to sell than to buy
            else {Pinch=2;};
            if (Plot=="On") {outputLog << "NEBLossAversion=" << NEBLossAversion << " : loss aversion => Pinch=" << FirstPinch << "->" << Pinch << endl;};
        }; // closes if
        
        if (NEBFear>VRan[13]) {
            int FirstQuant=Quant;
            if ((Sig==2) || (RFALevel==0) || (Liquid==0)) {Quant=0;}
            if (Plot=="On") {outputLog << "NEBFear=" << NEBFear << " : fear => Quant=" << FirstQuant << "->" << Quant << endl;};
        }; // closes if
        
        if (NEBFear2>VRan[13]) {
            int FirstQuant=Quant;
            if (Mu==0) {Quant=0;}
            if (Plot=="On") {outputLog << "NEBFear2=" << NEBFear2 << " : fear2 => Quant=" << FirstQuant << "->" << Quant << endl;};
        }; // closes if
        
        if (NEBGreed>VRan[13]) {
            int FirstQuant=Quant;
            if (Mu==2) {Quant=2;}
            if (Plot=="On") {outputLog << "NEBGreed=" << NEBGreed << " : fear => Quant=" << FirstQuant << "->" << Quant << endl;};
        }; // closes if
        
        
        // Cluster trading
        if ((AgentName<ClusterLimit) && (AgentName!=LeaderAgent) && (t>=LearningPhase)) {// WWW
            Quant=LeaderQuant;
            Pinch=LeaderPinch;
            cout << AgentName << ", ";
            if (Plot=="On") {outputLog << "i=" << AgentName << ", LeaderAgent=" << LeaderAgent << ", ClusterLimit=" << ClusterLimit << " => Quant=LeaderQuant=" << LeaderQuant << ", Pinch=LeaderPinch=" << LeaderPinch << " from LeaderAgent=" << LeaderAgent << endl;};
        };
        
        
        
        if (Plot=="On") {
            outputLog << "TenPiCoo2[]={" ; for (int u=0; u<int(TenPiCoo2.size()); u++) {if (u<int(TenPiCoo2.size())-1) {outputLog << TenPiCoo2[u] << ", ";} else {outputLog << TenPiCoo2[u];};}; outputLog << "}, ";
            outputLog << "Quant=TenPiCoo2[5]=" << TenPiCoo2[5] << ", Pinch=TenPiCoo2[6]=" << TenPiCoo2[6] << endl;
        }; // closes Plot condition
        // Defining the real action as tensor index to take according to above policy which selected Quant and Pinch
        if (t==1) {TAIndReal[0]=1; TAIndReal[1]=0;}
        else {TAIndReal[0]=Quant; TAIndReal[1]=Pinch;}; // Chosen in the policy above
        if (Plot=="On") {
            outputLog << "TADim[]={" ; for (int u=0; u<2; u++) {if (u<2-1) {outputLog << TADim[u] << ", ";} else {outputLog << TADim[u];};}; outputLog << "}" << endl;
            outputLog << "TAIndReal[]={" ; for (int u=0; u<2; u++) {if (u<2-1) {outputLog << TAIndReal[u] << ", ";} else {outputLog << TAIndReal[u];};}; outputLog << "}" << endl;
        }; // closes Plot condition
        // outputDebug << "d, ";
        
        ARealTensorIndex = get_tensor_index (2, TAIndReal, TADim); // Vector index in the tensor of all possible RL actions
        TAIndexReal[j].push_back(ARealTensorIndex); // Recording real action a_t (ABM#2 for Forecast)
        if (int(TAIndexReal[j].size()) - Future - 1 >0) {TAIndexReal[j].erase (TAIndexReal[j].begin());}; // OOO
        if (Plot=="On") {
            outputLog << "ARealTensorIndex=" << ARealTensorIndex << endl;
            outputLog << "TAIndexReal[j][]={" ; for (int u=0; u<int(TAIndexReal[j].size()); u++) {if (u<int(TAIndexReal[j].size())-1)  {outputLog << TAIndexReal[j][u] << ", ";} else {outputLog << TAIndexReal[j][u];};}; outputLog << "}, ";
            // Update the Bid and Ask that will be filed in the OB to take into account the Pinch
            outputLog << endl << "Updating the Bid and Ask that will be filed in the OB to take into account the Pinch" << endl;
        }; // closes Plot condition
        //double Spread=abs(Stocks[j].StockAskValue - Stocks[j].StockBidValue); // GGG
        double Spread=VSpread; // VSpread is the absolute value of average top bid minus average top ask
        if (Plot=="On") {
            outputLog << "ARealTensorIndex=" << ARealTensorIndex << endl;
            outputLog << "Stocks[j].StockBidValue=Mini=" << Stocks[j].StockBidValue << endl;
            outputLog << "Stocks[j].StockAskValue=Maxi=" << Stocks[j].StockAskValue << endl;
        }; // closes Plot condition
        // outputDebug << "e, ";
        
        // Gesture=0.2+0.8r
        // Stocks[j].StockBidValue = min(Pt,ˆPt)±Gesture*Spread
        // Stocks[j].StockAskValue = max(Pt,ˆPt)±Gesture*Spread
        
        double TransactionVirtualBid=ReflexiveVal_t;
        double TransactionVirtualAsk=ReflexiveVal_t;
        if (Exploration==4) {TransactionVirtualBid=Stocks[j].StockBidValue; TransactionVirtualAsk=Stocks[j].StockAskValue;}; // Backup to compute the virtual return (see below)
        if (Plot=="On") {outputLog << "Since Spread=" << Spread << "$, Pinch=" << Pinch << ", and Gesture=" << floor(Gesture*100-100) << "% of the Spread, we have:" << endl;};
        
        
        if (Pinch==0) {Stocks[j].StockBidValue+=Gesture*Spread; Stocks[j].StockAskValue-=Gesture*Spread;} // This corresponds to Pinch=-Gesture (more desperate)
        else if (Pinch==2) {Stocks[j].StockBidValue-=Gesture*Spread; Stocks[j].StockAskValue+=Gesture*Spread;}; // This corresponds to Pinch=+Gesture (less desperate)
        
        if (Trunk>-1) {
            if (Pinch==0) { // This corresponds to Pinch=-Gesture (more desperate)
                Stocks[j].StockBidValue=DigitTrunk (Mini+Gesture*Spread, Trunk, "Ceil");
                Stocks[j].StockAskValue=DigitTrunk (Maxi-Gesture*Spread, Trunk, "Ceil");
            }
            else if (Pinch==2) { // This corresponds to Pinch=+Gesture (less desperate)
                Stocks[j].StockBidValue=DigitTrunk (Mini-Gesture*Spread, Trunk, "Ceil");
                Stocks[j].StockAskValue=DigitTrunk (Maxi+Gesture*Spread, Trunk, "Ceil");
            };
            if (abs(Stocks[j].StockBidValue-ceil(Stocks[j].StockBidValue))>0) {cout << "ISSUE : Stocks[j].StockBidValue=" << Stocks[j].StockBidValue << endl;};
            if (abs(Stocks[j].StockAskValue-ceil(Stocks[j].StockAskValue))>0) {cout << "ISSUE : Stocks[j].StockAskValue=" << Stocks[j].StockAskValue << endl;};
        };
        
        if (Stocks[j].StockBidValue<0) {Stocks[j].StockBidValue=0;}; // Failsafe for the Pinch effect
        if (Stocks[j].StockAskValue<0) {Stocks[j].StockAskValue=0;}; //Failsafe for the Pinch effect
        if (Plot=="On") {
            outputLog << "Stocks[j].StockBidValue=" << Stocks[j].StockBidValue << endl;
            outputLog << "Stocks[j].StockAskValue=" << Stocks[j].StockAskValue << endl;
        }; // closes Plot condition
        double RealTransac=ReflexiveVal_t;
        int QuantityFraction=1;
        double AskVal=Stocks[j].StockQuantity/QuantityFraction;
        double BidVal=RFA/((NumberOfStocks+1)*Stocks[j].StockAskValue)/QuantityFraction; // Cannot buy more than its RFA divided by number of stocks
        //double BidVal=RFA/((NumberOfStocks)*Stocks[j].StockAskValue); // Cannot buy more than its RFA divided by number of stocks // GGG
        if ((RFA<0) || (Stocks[j].StockAskValue==0)) {BidVal=0;}; // Failsafe
        if (Quant==0) {QTradeReal=-int(AskVal); RealTransac=Stocks[j].StockAskValue;}
        else if (Quant==1) {QTradeReal=0; RealTransac=ReflexiveVal_t;}
        else if (Quant==2) {QTradeReal=int(BidVal); RealTransac=Stocks[j].StockBidValue;}
        QuantitiesReal[j].push_back(QTradeReal); // By default, but if there is an actual transaction then it is in OB's trading
        if (int(QuantitiesReal[j].size()) - Future - 1 >0) {QuantitiesReal[j].erase (QuantitiesReal[j].begin());}; // OOO
        TransactionPriceReal[j].push_back(RealTransac); // By default, but if there is an actual transaction then it is in OB's trading
        if (int(TransactionPriceReal[j].size()) - Future - 1 >0) {TransactionPriceReal[j].erase (TransactionPriceReal[j].begin());}; // OOO
        // outputDebug << "f, ";
        
        
        //ofstream outputTest("/Users/admin/Documents/GNT/SYMBA/ZULU.txt", ofstream::app); // ZZZ
        //if ((QTradeReal>0) && (QTradeReal>RFA/Stocks[j].StockBidValue)) {outputTest << "ErrorX Bid at i=" << AgentName << ", j=" << j << ", t=" << t << ": QTradeReal=" << QTradeReal << ", RFA/Stocks[j].StockBidValue=" << RFA << "/" << Stocks[j].StockBidValue << "=" << RFA/Stocks[j].StockBidValue << endl;};
        //if ((QTradeReal<0) && (abs(QTradeReal)>Stocks[j].StockQuantity)) {outputTest << "ErrorY Ask at i=" << AgentName << ", j=" << j << ", t=" << t << ": QTradeReal=" << QTradeReal << ", Stocks[j].StockQuantity=" << Stocks[j].StockQuantity << endl;};
        if (Plot=="On") {
            outputLog << "Since Quant=" << Quant << ", QTradeReal=" << QTradeReal << " (Stocks[j].StockQuantity=" << Stocks[j].StockQuantity << ", RFA=" << RFA << ", and Stocks[j].StockBidValue=" << Stocks[j].StockBidValue << ")";
            outputLog << endl << "QuantitiesReal[j][]={" ; for (int u=0; u<int(QuantitiesReal[j].size()); u++) {if (u<int(QuantitiesReal[j].size())-1)  {outputLog << QuantitiesReal[j][u] << ", ";} else {outputLog << QuantitiesReal[j][u];};}; outputLog << "}" << endl;
            outputLog << "RealTransac=" << RealTransac;
            outputLog << endl << "TransactionPriceReal[j][]={" ; for (int u=0; u<int(TransactionPriceReal[j].size()); u++) {if (u<int(TransactionPriceReal[j].size())-1)  {outputLog << TransactionPriceReal[j][u] << ", ";} else {outputLog << TransactionPriceReal[j][u];};}; outputLog << "}" << endl;
            QuantRealLog[j].push_back(Quant); // LLL
            PinchRealLog[j].push_back(Pinch); // LLL
        }; // closes Plot condition
        // outputDebug << "g" << endl;
        
        
        
        // Exploration5 retro-infering the right actions OptimalQuant and OptimalPinch at time t-Future by off-policy learning
        int OptimalQuant=0; int OptimalPinch=0; int TAid5 = 0;
        if (OffPolicy2==1) {
            // Finding OptimalQuant
            if (gsl_matrix_get (ReflexiveValues, j, t)>gsl_matrix_get (ReflexiveValues, j, t-Future)) {OptimalQuant=2; // Buy
                if (gsl_matrix_get (ReflexiveValues, j, t-Future+1)>gsl_matrix_get (ReflexiveValues, j, t-Future)) {OptimalPinch=0;} // Liberal gesture
                else if (gsl_matrix_get (ReflexiveValues, j, t-Future+1)==gsl_matrix_get (ReflexiveValues, j, t-Future)) {OptimalPinch=1;} // Neutral gesture
                else if (gsl_matrix_get (ReflexiveValues, j, t-Future+1)<gsl_matrix_get (ReflexiveValues, j, t-Future)) {OptimalPinch=2;}; // Tough gesture
            }
            
            else if (gsl_matrix_get (ReflexiveValues, j, t)==gsl_matrix_get (ReflexiveValues, j, t-Future)) {OptimalQuant=1; OptimalPinch=1;}// Hold, neutral gesture
            
            else if (gsl_matrix_get (ReflexiveValues, j, t)<gsl_matrix_get (ReflexiveValues, j, t-Future)) {OptimalQuant=0; // Sell
                if (gsl_matrix_get (ReflexiveValues, j, t-Future+1)>gsl_matrix_get (ReflexiveValues, j, t-Future)) {OptimalPinch=2;} // Tough gesture
                else if (gsl_matrix_get (ReflexiveValues, j, t-Future+1)==gsl_matrix_get (ReflexiveValues, j, t-Future)) {OptimalPinch=1;} // Neutral gesture
                else if (gsl_matrix_get (ReflexiveValues, j, t-Future+1)<gsl_matrix_get (ReflexiveValues, j, t-Future)) {OptimalPinch=0;}; // Liberal gesture
            };
            // Finding OptimalPinch (issue with hold because all Pinch should technically be updated equiprobably, but it's ok if we select one only)
            if (Plot=="On") {outputLog << "Since P(t)=" << gsl_matrix_get (ReflexiveValues, j, t) << ", P(t-T)=" << gsl_matrix_get (ReflexiveValues, j, t-Future) << ", P(t-T+1)=" << gsl_matrix_get (ReflexiveValues, j, t-Future+1) << ": OptimalQuant=" << OptimalQuant << " (0: short, 1: hold, 2: long) and OptimalPinch=" << OptimalPinch << " (0: liberal, 1: neutral, 2: tough gesture)" << endl;};
            TAInd5[0]=OptimalQuant; TAInd5[1]=OptimalPinch;
            if (Plot=="On") {outputLog << "TAInd5[0]=OptimalQuant=" << TAInd5[0] << ", TAInd5[1]=OptimalPinch=" << TAInd5[1] << endl;};
            TAid5 = get_tensor_index (2, TAInd5, TADim); // Vector index in the tensor of all possible RL actions
            //TAIndex5[j].push_back(TAid5); // Recording real action a_t (ABM#2 for Trading)
            //if (int(TAIndex5[j].size()) - Future - 1 >0) {TAIndex5[j].erase (TAIndex5[j].begin());};
        }; // closes if
        
        
        
        
        /*************************************ABM II.4******************************************/
        if (Plot=="On") {outputLog << endl << "    ABM 2.4: Randomly selecting a virtual action a=(VirtualQuant,VirtualPinch)" << endl;};
        // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 2.3 computed..." << endl;
        
        // Defining the virtual action as tensor index to take (not based on any exploration policy)
        int VirtualQuant, VirtualPinch;
        if (Exploration==4) {
            VirtualQuant=int(1000*VRan[8])%3; // 0, 1, 2, 3, 4, 5, or 6
            VirtualPinch=int(1000*VRan[9])%3; // 0, 1, or 2
            for (int i=0; i<100; i++) {int VQ=int(1000*VRan[10])%3; if (VQ!=Quant) {VirtualQuant=VQ; break;};}; // Making sure virtual action not like real one
            for (int i=0; i<100; i++) {int VP=int(1000*VRan[11])%3; if (VP!=Lag) {VirtualPinch=VP; break;};}; // Making sure virtual action not like real one
            if (Plot=="On") {outputLog << "VirtualQuant=" << VirtualQuant << ", VirtualPinch=" << VirtualPinch << endl;};
            double SpreadVirtual=abs(TransactionVirtualAsk - TransactionVirtualBid);
            if (VirtualPinch==0) {TransactionVirtualBid+=Gesture*SpreadVirtual; TransactionVirtualAsk-=Gesture*SpreadVirtual;} // This corresponds to VirtualPinch=-Gesture
            else if (VirtualPinch==2) {TransactionVirtualBid-=Gesture*SpreadVirtual; TransactionVirtualAsk+=Gesture*SpreadVirtual;}; // This corresponds to VirtualPinch=+Gesture
            if (TransactionVirtualBid<0) {TransactionVirtualBid=0;}; // Failsafe for the VirtualPinch effect
            if (TransactionVirtualAsk<0) {TransactionVirtualAsk=0;}; // Failsafe for the VirtualPinch effect
            double VirtualTransac=ReflexiveVal_t;
            if (VirtualQuant==0) {QTradeVirtual=-Stocks[j].StockQuantity; VirtualTransac=TransactionVirtualAsk;}
            else if (VirtualQuant==1) {QTradeVirtual=0; VirtualTransac=ReflexiveVal_t;}
            else if (VirtualQuant==2) {QTradeVirtual=int(floor(1*RFA/Stocks[j].StockBidValue)); VirtualTransac=TransactionVirtualBid;};
            if ((Stocks[j].StockBidValue==0) && (VirtualQuant==2)) {QTradeVirtual=0;} // For the case RFA/Stocks[j].StockBidValue=#IND
            if (Plot=="On") {outputLog << "Since VirtualQuant=" << VirtualQuant << ", QTradeVirtual=" << QTradeVirtual << " (Stocks[j].StockQuantity=" << Stocks[j].StockQuantity << ", RFA=" << RFA << ", and Stocks[j].StockBidValue=" << Stocks[j].StockBidValue << ")" << endl;};
            //if (t==1) {TAIndVirtual[0]=3; TAIndVirtual[1]=0;}
            TAIndVirtual[0]=VirtualQuant; TAIndVirtual[1]=VirtualPinch; // Chosen in the policy above
            if (Plot=="On") {outputLog << "TAIndVirtual[0]=VirtualQuant=" << TAIndVirtual[0] << ", TAIndVirtual[1]=VirtualPinch=" << TAIndVirtual[1] << endl;};
            int TAid = get_tensor_index (2, TAIndVirtual, TADim); // Vector index in the tensor of all possible RL actions
            TAIndexVirtual[j].push_back(TAid); // Recording real action a_t (ABM#2 for Trading)
            if (int(TAIndexVirtual[j].size()) - Future - 1 >0) {TAIndexVirtual[j].erase (TAIndexVirtual[j].begin());}; // OOO
            QuantitiesVirtual[j].push_back(QTradeVirtual); // This doest need to be in OB's trading
            if (int(QuantitiesVirtual[j].size()) - Future - 1 >0) {QuantitiesVirtual[j].erase (QuantitiesVirtual[j].begin());}; // OOO
            TransactionPriceVirtual[j].push_back(VirtualTransac); // By default, but if there is an actual transaction then it is in OB's trading
            if (int(TransactionPriceVirtual[j].size()) - Future - 1 >0) {TransactionPriceVirtual[j].erase (TransactionPriceVirtual[j].begin());}; // OOO
            if (Plot=="On") {
                outputLog << "TAid=" << TAid << endl;
                outputLog << "TAIndexVirtual[j][]={" ; for (int u=0; u<int(TAIndexVirtual[j].size()); u++) {if (u<int(TAIndexVirtual[j].size())-1)  {outputLog << TAIndexVirtual[j][u] << ", ";} else {outputLog << TAIndexVirtual[j][u];};}; outputLog << "}" << endl;
                outputLog << "QTradeVirtual=" << QTradeVirtual;
                outputLog << endl << "QuantitiesVirtual[j][]={" ; for (int u=0; u<int(QuantitiesVirtual[j].size()); u++) {if (u<int(QuantitiesVirtual[j].size())-1)  {outputLog << QuantitiesVirtual[j][u] << ", ";} else {outputLog << QuantitiesVirtual[j][u];};}; outputLog << "}" << endl;
                outputLog << "VirtualTransac=" << VirtualTransac << " (=ReflexiveVal_t)" << endl;
                outputLog << "TransactionPriceVirtual[j][]={" ; for (int u=0; u<int(TransactionPriceVirtual[j].size()); u++) {if (u<int(TransactionPriceVirtual[j].size())-1)  {outputLog << TransactionPriceVirtual[j][u] << ", ";} else {outputLog << TransactionPriceVirtual[j][u];};}; outputLog << "}" << endl;
                QuantVirtualLog[j].push_back(VirtualQuant); // LLL
                PinchVirtualLog[j].push_back(VirtualPinch); // LLL
            }; // closes Plot condition
        }; // closes Exploration==4 condition
        
        /*************************************ABM II.5******************************************/
        if (Plot=="On") {outputLog << endl << "    ABM 2.5: Calculating the return if t>Future as the cashflow difference if the action had been taken and if it had not" << endl;};
        
        // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 2.4 computed..." << endl;
        double DividendYield=2.276*0.01/Year;
        double BrokerFees=0.01*0.01; // 0.01% on the LSE
        if (t==1) {
            dTResultVec[j].push_back(0);
            for (int i=0; i<TS*TA; i++) {TQ[j].push_back(0); TNumberA[j].push_back(0);};
            if (Plot=="On") {outputLog << "Initialization of TQ[j][], TNumberA[j][], by Agent " << AgentName << ", t=" << t << endl;};
        };
        // We now calculate the return at time $t-T$. Let TransacQ(t-T) be the quantity of transaction at time $t-T$, TransacRFA(t-T) be the cash flow from RFA due to transaction at time $t-T$. Then the return is given by NAV(t) - NAV'(t), where NAV'(t) is the NAV(t) if the transaction at time $t-T$ had not taken place (Q(t-T) [P(t) - T(t-T)]), given by RFA(t) - TransacRFA(t-T) + (Q(t) - TransacQ(t-T))*MarketPrice(t). So the return is simply TransacRFA(t-T) + TransacQ(t-T)*MarketPrice(t)
        if (t-Future<=0) {if (Plot=="On") {outputLog << "We are not in the case t>Future (t=" << t << ", Future=" << Future << ")" << endl;}}
        else if (t-Future>0) {
            if (Plot=="On") {outputLog << "We are in the case t=" << t << ">" << Future << "=Future" << endl;};
            double dTResultReal=0; double dTResultVirtual=0;
            int Q=QuantitiesReal[j][max(0,int(QuantitiesReal[j].size()) - Future - 1)];
            double PastTransacPt=TransactionPriceReal[j][max(0,int(TransactionPriceReal[j].size()) - Future - 1)];
            dTResultReal=Q*(ReflexiveVal_t-PastTransacPt); // TransacRFA(t-T) + TransacQ(t-T)*MarketPrice(t). Note: the transaction price should be the one realized at OB (may be cancelled!) // JJJ
            dTResultReal+=Q*PastTransacPt*DividendYield*Future; // Dividends
            dTResultReal-=2*abs(Q)*PastTransacPt*BrokerFees; // Broker fees, one for in and one for out of position
            if (Exploration==4) {dTResultVirtual=QuantitiesVirtual[j][max(0,int(QuantitiesVirtual[j].size()) - Future - 1)] * (ReflexiveVal_t - TransactionPriceVirtual[j][max(0,int(TransactionPriceVirtual[j].size()) - Future - 1)]);};
            //TResultReal.push_back(dTResultReal);
            //if (int(TResultReal.size()) - Future - 1 >0) {TResultReal.erase (TResultReal.begin());}; // OOO
            //TResultVirtual.push_back(dTResultVirtual);
            //if (int(TResultVirtual.size()) - Future - 1 >0) {TResultVirtual.erase (TResultVirtual.begin());}; // OOO
            double dYield=100.0*(dTResultReal*NumberOfStocks*Year/Future)/(RFA+RBA); // Returns made by this action in % of agent NAV extended as adjusted return to all year
            double dMarketYield=(MarketPerformance-100)*Year/TimeSinceJan1st; // Or should we use Reflexive_t instead?
            double dTResult5=0;
            // Not really necessary since dTdisResult5 is always set at 4, regardless of dTResult5
            if (OffPolicy2==1) {
                int PseudoQuantity=ceil(RFAFirst/((NumberOfStocks+1)*ReflexiveVal_t)/QuantityFraction);
                dTResult5=PseudoQuantity * (ReflexiveVal_t - gsl_matrix_get (ReflexiveValues, j, t-Future)); // Approximation of Exploration5 by QuantitiesReal // JJJ
                dTResult5+=PseudoQuantity * gsl_matrix_get (ReflexiveValues, j, t-Future) * DividendYield * Future;
                dTResult5-=2*abs(PseudoQuantity) * gsl_matrix_get (ReflexiveValues, j, t-Future) * BrokerFees; // One for in and one for out of position
                if (Plot=="On") {outputLog << "PseudoQuantity=" << PseudoQuantity << ", P(t)=" << ReflexiveVal_t << ", P(t-T)=" << gsl_matrix_get (ReflexiveValues, j, t-Future) << " => dTResult5=" << dTResult5 << endl;};
            };
            if (Plot=="On") {
                outputLog << "At time t=" << t-Future << ", real placed order was $" << TransactionPriceReal[j][max(0,int(TransactionPriceReal[j].size()) - Future - 1)] << " (for a quantity " << QuantitiesReal[j][max(0,int(QuantitiesReal[j].size()) - Future - 1)] << " thereof)" << endl;
                if (Exploration==4) {outputLog << "At time t=" << t-Future << ", virtual placed order was $" << TransactionPriceVirtual[j][max(0,int(TransactionPriceVirtual[j].size()) - Future - 1)] << " (for a quantity " << QuantitiesVirtual[j][max(0,int(QuantitiesVirtual[j].size()) - Future - 1)] << " thereof)" << endl;};
                outputLog << "Finally today it is $" << ReflexiveVal_t << endl;
                outputLog << "So we get dTResultReal=$" << dTResultReal << ", dYield=" << dYield << "%, dMarketYield=" << dMarketYield << "%" << endl;
                if (Exploration==4) {outputLog << "And dTResultVirtual=$" << dTResultVirtual << endl;};
                //outputLog << "TResultReal[]={" ; for (int u=0; u<int(TResultReal.size()); u++) {if (u<int(TResultReal.size())-1)  {outputLog << TResultReal[u] << ", ";} else {outputLog << TResultReal[u];};}; outputLog << "}" << endl;
                //outputLog << "TResultVirtual[]={" ; for (int u=0; u<int(TResultVirtual.size()); u++) {if (u<int(TResultVirtual.size())-1)  {outputLog << TResultVirtual[u] << ", ";} else {outputLog << TResultVirtual[u];};}; outputLog << "}" << endl;
                
                
                
                
                /*************************************ABM II.6******************************************/
                outputLog << endl << "    ABM 2.6: Returns above now discretized according to past returns (cf. reinforcement comparison & actor-critic methods)" << endl;
            }; // closes Plot condition
            
            
            double dTResultPercentile=0; int dTResultVecSize=int(dTResultVec[j].size()); int dTResultDisReal=0;
            bool CondRecord=0;
            for (int i=0; i<int(ExitHorizons.size()); i++) {
                if (t==ExitHorizons[i]) {CondRecord=1;};
            };
            if (CondRecord==1) { // JJJ3
                for (int i=0; i<dTResultVecSize; i++) { // sorting the vector SmoothVec by ascending values
                    if (dTResultReal<dTResultVec[j][i]) {
                        dTResultVec[j].insert(dTResultVec[j].begin()+i, dTResultReal*(1 + 0.0001*VRan[27]));
                        dTResultPercentile=i*1.0/dTResultVecSize;
                        break;
                    };
                    if (i==dTResultVecSize-1) {
                        dTResultVec[j].push_back(dTResultReal*(1 + 0.0001*VRan[27]));
                        dTResultPercentile=1;
                        break;
                    };
                }; // closes for loop
                if ((dTResultReal>0) && (dTResultDisReal<0)) {dTResultDisReal=1;}; // Making sure a net profit is not accounted as a negative reward XYZ333
            }; // closes CondRecord condition
            
            if (dTResultPercentile<0.05) {dTResultDisReal=-4;}
            else if ((dTResultPercentile>=0.05) && (dTResultPercentile<0.25)) {dTResultDisReal=-2;}
            else if ((dTResultPercentile>=0.25) && (dTResultPercentile<0.5)) {dTResultDisReal=-1;}
            else if ((dTResultPercentile>=0.5) && (dTResultPercentile<0.75)) {dTResultDisReal=1;}
            else if ((dTResultPercentile>=0.75) && (dTResultPercentile<0.95)) {dTResultDisReal=2;}
            else if (dTResultPercentile>=0.95) {dTResultDisReal=4;};
            if (dTResultVecSize==1) {dTResultDisReal=1;}; // In the beginning it means nothing
            if (dTResultVecSize - Future - History >0) {dTResultVec[j].erase (dTResultVec[j].begin());}; // OOO
            
            
            if ((NEBPositivity>VRan[17]) && (dTResultDisReal>0)) {
                double FirstdTResultDisReal=dTResultDisReal;
                dTResultDisReal+=1; // Positivity bias: better recalling pleasant memories than unpleasant ones
                if (Plot=="On") {outputLog << "NEBPositivity=" << NEBPositivity << " : positivity bias switched on => dTResultDisReal=" << FirstdTResultDisReal << "->" << dTResultDisReal << endl;};
            };
            
            if ((NEBNegativity>VRan[17]) && (dTResultDisReal<0)) {
                double FirstdTResultDisReal=dTResultDisReal;
                dTResultDisReal+=1; // Positivity bias: better recalling pleasant memories than unpleasant ones
                if (Plot=="On") {outputLog << "NEBNegativity=" << NEBNegativity << " : negativity bias switched on => dTResultDisReal=" << FirstdTResultDisReal << "->" << dTResultDisReal << endl;};
            };
            
            if ((dTResultReal>0) && (dTResultDisReal<0)) {
                dTResultDisReal=1; // Making sure a net profit is not accounted as a negative reward
                if (Plot=="On") {outputLog << "Making sure a net profit is not accounted as a negative reward => dTResultDisReal=1" << endl;};
            };
            TResultDisReal[j].push_back(dTResultDisReal); // LLL
            
            //if (int(TResultDisReal.size()) - Future - 1 >0) {TResultDisReal.erase (TResultDisReal.begin());}; // OOO
            /*
             ofstream outputCheck("/Users/admin/Documents/GNT/SYMBA/Check.txt", ofstream::app);
             if (AgentName==5) {
             outputCheck << "At t=" << t << ": dTResultReal=" << dTResultReal*(1 + 0.0001*VRan[27]) << ", dTResultRealPercentile=" << ceil(100*dTResultRealPercentile) << "% => dTResultDisReal=" << dTResultDisReal << endl;
             outputCheck << "dTResultVec[]={" ; for (int u=0; u<int(dTResultVec.size()); u++) {if (u<int(dTResultVec.size())-1) {outputCheck << dTResultVec[u] << ", ";} else {outputCheck << dTResultVec[u];};}; outputCheck << "}" << endl;
             };
             */
            if (Plot=="On") {
                //outputLog << "TRMean=" << TRMean << ", TRVariance=" << TRVariance << ", TRStdDev=" << TRStdDev << endl;
                //outputLog << "BBMinus=" << BBMinus << ", BBPlus=" << BBPlus << endl;
                outputLog << "dTResultDisReal=" << dTResultDisReal << " (dTResultReal=" << dTResultReal << ")" << endl;
                outputLog << "TResultDisReal[j][]={" ; for (int u=0; u<int(TResultDisReal[j].size()); u++) {if (u<int(TResultDisReal[j].size())-1) {outputLog << TResultDisReal[j][u] << ", ";} else {outputLog << TResultDisReal[j][u];};}; outputLog << "}" << endl;
            }; // closes Plot condition
            // We do the same for the virtual results
            
            int dTResultDisVirtual=0;
            if (Exploration==4) {
                double dTResultPercentile=0; int dTResultVecSize=int(dTResultVec[j].size());
                for (int i=0; i<dTResultVecSize; i++) { // sorting the vector SmoothVec by ascending values
                    if (dTResultVirtual<dTResultVec[j][i]) {dTResultVec[j].insert(dTResultVec[j].begin()+i, dTResultVirtual*(1 + 0.0001*VRan[28])); dTResultPercentile=i*1.0/dTResultVecSize; break;};
                    if (i==dTResultVecSize-1) {dTResultVec[j].push_back(dTResultVirtual*(1 + 0.0001*VRan[28])); dTResultPercentile=1; break;};
                }; // closes for loop
                if (dTResultPercentile<0.05) {dTResultDisVirtual=-4;}
                else if ((dTResultPercentile>=0.05) && (dTResultPercentile<0.25)) {dTResultDisVirtual=-2;}
                else if ((dTResultPercentile>=0.25) && (dTResultPercentile<0.5)) {dTResultDisVirtual=-1;}
                else if ((dTResultPercentile>=0.5) && (dTResultPercentile<0.75)) {dTResultDisVirtual=1;}
                else if ((dTResultPercentile>=0.75) && (dTResultPercentile<0.95)) {dTResultDisVirtual=2;}
                else if (dTResultPercentile>=0.95) {dTResultDisVirtual=4;};
                if (dTResultVecSize==1) {dTResultDisVirtual=1;}; // In the beginning it means nothing
                if (dTResultVecSize - Future - History >0) {dTResultVec[j].erase (dTResultVec[j].begin());}; // OOO
                TResultDisVirtual[j].push_back(dTResultDisVirtual); // LLL
                if (Plot=="On") {
                    outputLog << "dTResultDisVirtual=" << dTResultDisVirtual << " (dTResultVirtual=" << dTResultVirtual << ")" << endl;
                    outputLog << "TResultDisVirtual[j][]={" ; for (int u=0; u<int(TResultDisVirtual[j].size()); u++) {if (u<int(TResultDisVirtual[j].size())-1) {outputLog << TResultDisVirtual[j][u] << ", ";} else {outputLog << TResultDisVirtual[j][u];};}; outputLog << "}" << endl;
                }; // closes plot condition
            }; // closes if
            
            double dTResultdis5=4;
            if (OffPolicy2==1) {
                if (Plot=="On") {outputLog << "dTResultdis5=" << dTResultdis5 << " (always)" << endl;};
            };
            
            
            /*************************************ABM II.7******************************************/
            if (Plot=="On") {outputLog << endl << "    ABM 2.7: Updating Q(s,a) based on these returns via SAM" << endl;};
            // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 2.6 computed..." << endl;
            
            // We now update the action-value function $Q(s,a)$ based on this return via the SAM or Sample Average Method
            vector<int> OldStates = get_tensor_coord (TSIndex[j][max(0,int(TSIndex[j].size()) - Future - 1)], 5, TSDim); // Coord in s-tensor of past state $s_{t-T}$
            vector<int> OldRealActions = get_tensor_coord (TAIndexReal[j][max(0,int(TAIndexReal[j].size()) - Future - 1)], 2, TADim); //Coord in a-tensor of past action $a_{t-T}$
            vector<int> OldVirtualActions, OldActions5;
            dTResultDisReal=TResultDisReal[j][max(0,int(TResultDisReal[j].size()) - Future - 1)]; // Past value t-Future time steps ago that will be used to update both Q() and pi() // TTT1
            
            if (Exploration==4) {
                OldVirtualActions = get_tensor_coord (TAIndexVirtual[j][max(0,int(TAIndexVirtual[j].size()) - Future - 1)], 2, TADim); //Coord in a-tensor of past action
                dTResultDisVirtual=TResultDisVirtual[j][max(0,int(TResultDisVirtual[j].size()) - Future - 1)]; // Past value t-Future time steps ago that will be used to update both Q() and pi() // TTT1
            };
            
            if (OffPolicy2==1) {OldActions5=get_tensor_coord (TAid5, 2, TADim);};
            
            if (Plot=="On") {
                outputLog << "OldStates[]={" ; for (int u=0; u<int(OldStates.size()); u++) {if (u<int(OldStates.size())-1) {outputLog << OldStates[u] << ", ";} else {outputLog << OldStates[u];};}; outputLog << "}" << endl;
                outputLog << "OldRealActions[]={" ; for (int u=0; u<int(OldRealActions.size()); u++) {if (u<int(OldRealActions.size())-1) {outputLog << OldRealActions[u] << ", ";} else {outputLog << OldRealActions[u];};}; outputLog << "}" << endl;
                if (Exploration==4) {outputLog << "OldVirtualActions[]={" ; for (int u=0; u<int(OldVirtualActions.size()); u++) {if (u<int(OldVirtualActions.size())-1) {outputLog << OldVirtualActions[u] << ", ";} else {outputLog << OldVirtualActions[u];};}; outputLog << "}" << endl;};
                outputLog << "We are at t=Future=" << t << endl;
            }; // closes Plot condition
            TQIndReal[0]=OldStates[0]; TQIndReal[1]=OldStates[1]; TQIndReal[2]=OldStates[2]; TQIndReal[3]=OldStates[3]; TQIndReal[4]=OldStates[4];//s
            TQIndReal[5]=OldRealActions[0]; TQIndReal[6]=OldRealActions[1]; // a
            if (Exploration==4) {TQIndVirtual[0]=OldStates[0]; TQIndVirtual[1]=OldStates[1]; TQIndVirtual[2]=OldStates[2]; TQIndVirtual[3]=OldStates[3]; TQIndVirtual[4]=OldStates[4]; // s
                TQIndVirtual[5]=OldVirtualActions[0]; TQIndVirtual[6]=OldVirtualActions[1];}; // a
            if (OffPolicy2==1) {TQInd5[0]=OldStates[0]; TQInd5[1]=OldStates[1]; TQInd5[2]=OldStates[2]; TQInd5[3]=OldStates[3]; TQInd5[4]=OldStates[4]; // s
                TQInd5[5]=OldActions5[0]; TQInd5[6]=OldActions5[1];}; // a
            int QRealTensorIndex = get_tensor_index (7, TQIndReal, TQDim); // Vector index of Q-tensor <=> $s_{t-T}$ x $a_{t-T}$
            int QVirtualTensorIndex=0; int QTensorIndex5=0;
            if (Exploration==4) {QVirtualTensorIndex = get_tensor_index (7, TQIndVirtual, TQDim);}; // Vector index of Q-tensor <=> $s_{t-T}$ x $a_{t-T}$
            if (OffPolicy2==1) {QTensorIndex5=get_tensor_index (7, TQInd5, TQDim);};
            if (Plot=="On") {
                outputLog << "TQIndReal[]={" ; for (int u=0; u<7; u++) {if (u<7-1) {outputLog << TQIndReal[u] << ", ";} else {outputLog << TQIndReal[u];};}; outputLog << "}" << endl;
                outputLog << "TQDim[]={" ; for (int u=0; u<7; u++) {if (u<7-1) {outputLog << TQDim[u] << ", ";} else {outputLog << TQDim[u];};}; outputLog << "}" << endl;
                outputLog << "QRealTensorIndex=" << QRealTensorIndex << endl;
                if (Exploration==4) {outputLog << "TQIndVirtual[]={" ; for (int u=0; u<7; u++) {if (u<7-1) {outputLog << TQIndVirtual[u] << ", ";} else {outputLog << TQIndVirtual[u];};}; outputLog << "}" << endl;
                    outputLog << "QVirtualTensorIndex=" << QVirtualTensorIndex << endl;};
                if (OffPolicy2==1) {outputLog << "TQInd5[]={" ; for (int u=0; u<7; u++) {if (u<7-1) {outputLog << TQInd5[u] << ", ";} else {outputLog << TQInd5[u];};}; outputLog << "}" << endl;
                    outputLog << "QTensorIndex5=" << QTensorIndex5 << endl;};
            }; // closes Plot condition
            
            if (CondRecord==1) {
                TNumberA[j][QRealTensorIndex]+=1; // Adding this action as taken in the ActionNumber tensor so as to count for SAM
                if (TNumberA[j][QRealTensorIndex]<=1) {TQ[j][QRealTensorIndex]=dTResultReal;} // No SAM if not yet populated (i.e =0)
                else {TQ[j][QRealTensorIndex]=(TQ[j][QRealTensorIndex]*(TNumberA[j][QRealTensorIndex]-1) + dTResultReal)/TNumberA[j][QRealTensorIndex];}; // SAM
                if (Plot=="On") {
                    outputLog << "Number of times this real action was previously taken in that former state: TNumberA[j][QRealTensorIndex]-1=" << TNumberA[j][QRealTensorIndex]-1 << ", so TQ[j][QRealTensorIndex]=" << TQ[j][QRealTensorIndex] << endl;
                }; // closes Plot condition
            };
            
            if (Exploration==4) {
                TNumberA[j][QVirtualTensorIndex]+=1; // Adding this action as taken in the ActionNumber tensor so as to count for SAM
                if (TNumberA[j][QVirtualTensorIndex]<=1) {TQ[j][QVirtualTensorIndex]=dTResultVirtual;} // No SAM if not yet populated
                else {TQ[j][QVirtualTensorIndex]=(TQ[j][QVirtualTensorIndex]*(TNumberA[j][QVirtualTensorIndex]-1) + dTResultVirtual)/TNumberA[j][QVirtualTensorIndex];}; // SAM
                if (Plot=="On") {
                    outputLog << "Number of times this virtual action was previously taken in that former state: TNumberA[j][QVirtualTensorIndex]-1=" << TNumberA[j][QVirtualTensorIndex]-1 << ", so TQ[j][QVirtualTensorIndex]=" << TQ[j][QVirtualTensorIndex] << endl;
                }; // closes Plot condition
            }; // closes if
            
            if (OffPolicy2==1) {
                TNumberA[j][QTensorIndex5]+=1; // Adding this action as taken in the ActionNumber tensor so as to count for SAM
                if (TNumberA[j][QTensorIndex5]<=1) {TQ[j][QTensorIndex5]=dTResult5;} // No SAM if not yet populated
                else {TQ[j][QTensorIndex5]=(TQ[j][QTensorIndex5]*(TNumberA[j][QTensorIndex5]-1) + dTResult5)/TNumberA[j][QTensorIndex5];}; // SAM
                if (Plot=="On") {
                    outputLog << "Number of times this off-policy action was previously taken in that former state: TNumberA[j][QTensorIndex5]-1=" << TNumberA[j][QTensorIndex5]-1 << ", so TQ[j][QTensorIndex5]=" << TQ[j][QTensorIndex5] << endl;
                }; // closes Plot condition
            };
            
            
            
            
            /*************************************ABM II.8******************************************/
            if (Plot=="On") {
                outputLog << endl << "    ABM 2.8: Updating existing policy pi(s,a) based on these discretized results" << endl;
                outputLog << "TA=" << TA << ", TS=" << TS << endl;
            };
            // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 2.7 computed..." << endl;
            
            // We now update existing policy pi(s,a) based on these discretized results
            TSInd[0]=OldStates[0]; TSInd[1]=OldStates[1]; TSInd[2]=OldStates[2]; TSInd[3]=OldStates[3]; TSInd[4]=OldStates[4];
            int OldSTensorIndex = get_tensor_index (5, TSInd, TSDim); // Vector index of the s-tensor associated with $s_{t-T}$
            if (Plot=="On") {
                outputLog << "State at time t-Future: TSInd[0]=OldStates[0]=" << TSInd[0] << ", TSInd[1]=OldStates[1]=" << TSInd[1] << ", TSInd[2]=OldStates[2]=" << TSInd[2] << ", TSInd[3]=OldStates[3]=" << TSInd[3] << ", TSInd[4]=OldStates[4]=" << TSInd[4] << endl;
                outputLog << "OldSTensorIndex=" << OldSTensorIndex << endl;
            }; // closes Plot condition
            int OldSCount=0;
            for (int i=0; i<int(TSIndex[j].size())-Future; i++) {if (TSIndex[j][i]==OldSTensorIndex) {OldSCount+=1;};}; // Nb of times $s_{t-T}$ seen
            TQIndReal[0]=OldStates[0]; TQIndReal[1]=OldStates[1]; TQIndReal[2]=OldStates[2]; TQIndReal[3]=OldStates[3]; TQIndReal[4]=OldStates[4]; TQIndReal[5]=0; TQIndReal[6]=0; // Using results from above on Q for s but with zero action
            if (Plot=="On") {
                outputLog << "Number of times state at t-Future encountered: OldSCount=" << OldSCount << endl;
                outputLog << "TQIndReal[0]=OldStates[0]=" << TQIndReal[0] << ", TQIndReal[1]=OldStates[1]=" << TQIndReal[1] << ", TQIndReal[2]=OldStates[2]=" << TQIndReal[2] << ", TQIndReal[3]=OldStates[3]=" << TQIndReal[3] << ", TQIndReal[4]=OldStates[4]=" << TQIndReal[4] << ", TQIndReal[5]=" << 0 << ", TQIndReal[6]=" << 0 << endl;
            }; // closes Plot condition
            int PiSTensorIndex = get_tensor_index (7, TQIndReal, TQDim); // Vector index of the Q-tensor associated with $s_{t-T}$
            if (Plot=="On") {
                outputLog << "PiSTensorIndex=" << PiSTensorIndex << " from which we span all actions to update their probability" << endl;
                //outputLog << "State s was encountered OldSCount=" << OldSCount << " times" << endl;
                outputLog << "Results where dTResultDisReal=" << dTResultDisReal << ", but with CondRecord=" << CondRecord << " => 0: not taken into account, 1: taken into account" << endl;
                if (Exploration==4) {outputLog << "And dTResultDisVirtual=" << dTResultDisVirtual << endl;};
                if (OffPolicy2==1) {outputLog << "And dTResultdis5=" << dTResultdis5 << endl;};
            }; // closes Plot condition
            double Sum=0;
            for (int i=PiSTensorIndex; i<PiSTensorIndex+TA; i++) {
                //if (Plot=="On") {outputLog << "Old TPi[j][" << i << "]=" << floor(100*TPi[j][i]) << "%" << endl;};
                Sum+=TPi[j][i];
            };
            if (Plot=="On") {
                outputLog << "Sum=" << 100*Sum << endl;
                if ((int(100*Sum)<99) || (int(100*Sum)>101)) {outputLog << "ErrorOriginalT[]" << endl;};
            };
            Sum=0;
            
            if (CondRecord==1) { // JJJ3
                // First for real
                for (int i=PiSTensorIndex; i<PiSTensorIndex+TA; i++) {
                    if (Plot=="On") {outputLog << "Old real TPi[j][" << i << "]=" << 100*TPi[j][i] << "% => ";};
                    if ((i==QRealTensorIndex) && (dTResultDisReal>=0)) {
                        for (int f=0; f<dTResultDisReal; f++) {TPi[j][i]=TPi[j][i]*(1-NEBLearningRate/TNumberA[j][QRealTensorIndex]) + NEBLearningRate/TNumberA[j][QRealTensorIndex];};
                    } // closes if
                    else if ((i!=QRealTensorIndex) && (dTResultDisReal>=0)) {
                        for (int f=0; f<dTResultDisReal; f++) {TPi[j][i]=TPi[j][i]*(1-NEBLearningRate/TNumberA[j][QRealTensorIndex]);};
                    } // closes else if
                    else if ((i==QRealTensorIndex) && (dTResultDisReal<0)) {
                        for (int f=0; f<abs(dTResultDisReal); f++) {TPi[j][i]=TPi[j][i]*(1-NEBLearningRate/TNumberA[j][QRealTensorIndex]);};
                    } // closes else if
                    else if ((i!=QRealTensorIndex) && (dTResultDisReal<0)) {
                        for (int f=0; f<abs(dTResultDisReal); f++) {TPi[j][i]=TPi[j][i]*(1-NEBLearningRate/TNumberA[j][QRealTensorIndex]) + NEBLearningRate/(TA-1)/TNumberA[j][QRealTensorIndex];};
                    }; // closes else if
                    if (Plot=="On") {outputLog << "New real TPi[j][" << i << "]=" << 100*TPi[j][i] << "%" << endl;};
                    Sum+=TPi[j][i];
                }; // closes i loop
                if (Plot=="On") {
                    outputLog << "(QRealTensorIndex=" << QRealTensorIndex << " was the real action update), TNumberA[j][QRealTensorIndex]=" << TNumberA[j][QRealTensorIndex] << endl;
                    outputLog << "Sum=" << 100*Sum << endl;
                    if ((int(100*Sum)<99) || (int(100*Sum)>101)) {outputLog << "ErrorRealT[]: For OldSCount=" << OldSCount << ", dTResultDisReal=" << dTResultDisReal << endl;};
                }; // closes Plot condition
                Sum=0;
            }; // closes if
            // Second for virtual
            //if ((Exploration==4) && (1-NEB8>VRan[19])) { // This is for off-policy Watkin's Q-learning
            if (Exploration==4) { // This is for off-policy Watkin's Q-learning
                for (int i=PiSTensorIndex; i<PiSTensorIndex+TA; i++) {
                    if (Plot=="On") {outputLog << "Old virtual TPi[j][" << i << "]=" << 100*TPi[j][i] << "% => ";};
                    if ((i==QVirtualTensorIndex) && (dTResultDisVirtual>=0)) {
                        for (int f=0; f<dTResultDisVirtual; f++) {TPi[j][i]=TPi[j][i]*(1-NEBLearningRate/TNumberA[j][QVirtualTensorIndex]) + NEBLearningRate/TNumberA[j][QVirtualTensorIndex];};
                    } // closes if
                    else if ((i!=QVirtualTensorIndex) && (dTResultDisVirtual>=0)) {
                        for (int f=0; f<dTResultDisVirtual; f++) {TPi[j][i]=TPi[j][i]*(1-NEBLearningRate/TNumberA[j][QVirtualTensorIndex]);};
                    } // closes else if
                    else if ((i==QVirtualTensorIndex) && (dTResultDisVirtual<0)) {
                        for (int f=0; f<abs(dTResultDisVirtual); f++) {TPi[j][i]=TPi[j][i]*(1-NEBLearningRate/TNumberA[j][QVirtualTensorIndex]);};
                    } // closes else if
                    else if ((i!=QVirtualTensorIndex) && (dTResultDisVirtual<0)) {
                        for (int f=0; f<abs(dTResultDisVirtual); f++) {TPi[j][i]=TPi[j][i]*(1-NEBLearningRate/TNumberA[j][QVirtualTensorIndex]) + NEBLearningRate/(TA-1)/TNumberA[j][QVirtualTensorIndex];};
                    }; // closes else if
                    if (Plot=="On") {outputLog << "New virtual TPi[j][" << i << "]=" << 100*TPi[j][i] << "%" << endl;};
                    //if (TPi[i]<0) {outputLog << "ErrorVirtual2a: at i=" << i << ", TPi[i]=" << TPi[i] << endl;};
                    Sum+=TPi[j][i];
                }; // closes i loop
                if (Plot=="On") {
                    outputLog << "(QVirtualTensorIndex=" << QVirtualTensorIndex << " was the virtual action update), TNumberA[j][QVirtualTensorIndex]=" << TNumberA[j][QVirtualTensorIndex] << endl;
                    outputLog << "Sum=" << 100*Sum << endl;
                    if ((int(100*Sum)<99) || (int(100*Sum)>101)) {outputLog << "ErrorVirtualT[]: For OldSCount=" << OldSCount << ", dTResultDisVirtual=" << dTResultDisVirtual << endl;};
                }; // closes Plot condition
                Sum=0;
            }; // closes if Exploration
            
            
            // Second for Exploration5
            if (OffPolicy2==1) {
                for (int i=PiSTensorIndex; i<PiSTensorIndex+TA; i++) {
                    if (Plot=="On") {outputLog << "Old Exploration5 TPi[j][" << i << "]=" << 100*TPi[j][i] << "% => ";};
                    if (i==QTensorIndex5) {
                        for (int f=0; f<dTResultdis5; f++) {TPi[j][i]=TPi[j][i]*(1-NEBLearningRate/TNumberA[j][QTensorIndex5]) + NEBLearningRate/TNumberA[j][QTensorIndex5];};
                    } // closes if
                    else if (i!=QTensorIndex5) {
                        for (int f=0; f<dTResultdis5; f++) {TPi[j][i]=TPi[j][i]*(1-NEBLearningRate/TNumberA[j][QTensorIndex5]);};
                    } // closes else if
                    if (Plot=="On") {outputLog << "New Exploration5 TPi[j][" << i << "]=" << 100*TPi[j][i] << "%" << endl;};
                    //if (TPi[i]<0) {outputLog << "ErrorVirtual2a: at i=" << i << ", TPi[i]=" << TPi[i] << endl;};
                    Sum+=TPi[j][i];
                }; // closes i loop
                if (Plot=="On") {
                    outputLog << "(QTensorIndex5=" << QTensorIndex5 << " was the off-policy action update), TNumberA[j][QTensorIndex5]=" << TNumberA[j][QTensorIndex5] << endl;
                    outputLog << "Sum=" << 100*Sum << endl;
                    if ((int(100*Sum)<99) || (int(100*Sum)>101)) {outputLog << "ErrorOffPolicyT[]: For OldSCount=" << OldSCount << ", dTResultdis5=" << dTResultdis5 << endl;};
                }; // closes Plot condition
                Sum=0;
            }; // closes if Exploration
            
            
            
            
        }; // closes if (t-Future>0)
        
        if (t==Future-1) { // Cancel all previous learnings (both Q(s,a) and Pi(s,a)) as they cannot be considered valid (either in state s which is represented according to Past*Future, or forecast which are for Lag*Future with Lag=1,2,3)
            for (int i=0; i<FS*FA; i++) {FQ[j][i]=0; FNumberA[j][i]=0; FPi[j][i]=1.0/FA;};
            for (int i=0; i<TS*TA; i++) {TQ[j][i]=0; TNumberA[j][i]=0; TPi[j][i]=1.0/TA;};
        };
        if (t<Future-1) {QTradeReal=0;}; // So agent doesn't trade before time
        if (Plot=="On") {outputLog << endl << endl;};
        // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 2.8 computed..." << endl;
        
        
        
    }; // closes method RL()
    
    
    
}; // closes Agent class






class Order { // Each line in the order book
public:
    bool Nature; // 0 for Bid and 1 for Ask
    int Agent; // Agent placing the order
    int Stock; // Stock to trade;
    int Q; // Quantity of the stock to trade
    double Price; // Price or less for bids, Price or more for asks
};

class StockOrderBook { // Order book for a given stock
public:
    vector<Order> Bids; // collection of all bids
    vector<Order> Asks; // collection of all asks
    int MetaorderLastAgent; // Agent placing the order
    int MetaorderNature; // 0 for Bid and 1 for Ask
    double Credit; // Extra RFA (for a buy) or stock quantity (for a sell) temporarily credited to the agent
    
    
    void Sort() { // Sorting the order book for the given stock, with Bids of decreasing prices and Asks of increasing prices (SS9)
        vector<Order> SortedBids, SortedAsks, SortedBidsNoId, SortedAsksNoId;
        vector<double> SortedListBids, SortedListAsks;
        // Saving the prices only of these orders Bids and Asks into two separate vectors SortedListBids and SortedListAsks
        for (int k=0; k<int(Bids.size()); k++) {SortedListBids.push_back(Bids[k].Price);}; // closes k loop
        for (int k=0; k<int(Asks.size()); k++) {SortedListAsks.push_back(Asks[k].Price);}; // closes k loop
        // Sorting them
        sort(SortedListBids.begin(), SortedListBids.begin() + Bids.size()); // sort() is by ascending order, bad for Bids
        reverse(SortedListBids.begin(), SortedListBids.begin() + Bids.size()); // reverse() puts it by descending order, good for Bids
        sort(SortedListAsks.begin(), SortedListAsks.begin() + Asks.size()); // sort() is by ascending order, good for Asks
        // Sorting the orders Bids and Asks according to these two sorted lists of prices SortedListBids and SortedListAsks, and saving them in new sorted orders SortedBids and SortedAsks
        for (int k=0; k<int(SortedListBids.size()); k++) {
            for (int m=0; m<int(Bids.size()); m++) {
                if (Bids[m].Price==SortedListBids[k]) {
                    SortedBids.push_back(Bids[m]);
                    Bids.erase(Bids.begin()+m); // To ensure this order is not taken again in the next count if another bid has the same price
                    m=int(SortedListBids.size());
                }; // closes if
            }; // closes m loop
        }; // closes k loop
        for (int k=0; k<int(SortedListAsks.size()); k++) {
            for (int m=0; m<int(Asks.size()); m++) {
                if (Asks[m].Price==SortedListAsks[k]) {
                    SortedAsks.push_back(Asks[m]);
                    Asks.erase(Asks.begin()+m); // To ensure this order is not taken again in the next count if another ask has the same price
                    m=int(SortedListAsks.size());
                }; // closes if
            }; // closes m loop
        }; // closes k loop
        Bids=SortedBids;
        Asks=SortedAsks;
    }; // closes method Sort()
    
    
    
    vector<double> Clear(vector<Agent> &Market, int t, int j, string Plot) { // Passing by ref to modify Market (SS11)
        double TempClearingPrice=0;
        double TotalStockQuantityTraded=0;
        int Count=0;
        double TempClearingPriceLast=0;
        double TotalStockQuantityTradedLast=0;
        int CountLast=0;
        int LevelCount=0;
        double BrokerFees=0.01*0.01; // 0.01% on the LSE
        //double AverageBids=0;
        //double AverageAsks=0;
        //ofstream outputTest("/Users/admin/Documents/GNT/SYMBA/ZULU.txt", ofstream::app); int NumberOfAgents=30; int NumberOfStocks=2;
        ofstream outputLog(Machine+"SimLog.txt", ofstream::app);
        if (Plot=="On") {outputLog << endl << endl << "***************CLEAR() AT TIME t=" << t << "**************" << endl;};
        for (int k=-1; k<int(min(Bids.size(), Asks.size())); k++) {
            if (int(min(Bids.size(), Asks.size()))==0) {/*outputLog << "Function Clear() not activated" << endl;*/ break;}
            if (k==-1) {k=0;};
            //if (Plot=="On") {outputLog << "New loop: k=" << k << ", Bids.size()=" << Bids.size() << ", Asks.size()=" << Asks.size() << ", min(Bids.size(), Asks.size())=" << min(Bids.size(), Asks.size()) << endl;};
            int TempQuantity=int(min (Bids[k].Q, Asks[k].Q)); // KKK
            //int TempQuantity=max(0, int(min (Bids[k].Q, Asks[k].Q))); // KKK
            TempClearingPrice = (Bids[k].Price + Asks[k].Price)/2; // Transaction taken at midprice, a neutral and common solution
            bool C1 = Bids[k].Price>=Asks[k].Price; // Bid is larger than Ask at a given level k of the order book
            bool C2 = Market[Asks[k].Agent].Stocks[Asks[k].Stock].StockQuantity > 0; // Seller has non-zero quantity of this stock
            bool C3 = Market[Bids[k].Agent].RFA >= TempClearingPrice * TempQuantity; // Buyer has enough RFA to buy
            bool C4 = Market[Asks[k].Agent].Stocks[Asks[k].Stock].StockQuantity >= TempQuantity; // Seller has enough quantities to sell
            if (C1 && C2 && C3 && C4) {
                //int QBidsBefore=Market[Bids[k].Agent].Stocks[Bids[k].Stock].StockQuantity;
                //int QAsksBefore=Market[Asks[k].Agent].Stocks[Asks[k].Stock].StockQuantity;
                // This is the OB clearing per se
                Market[Bids[k].Agent].RFA -= TempClearingPrice * TempQuantity;
                Market[Asks[k].Agent].RFA += TempClearingPrice * TempQuantity;
                Market[Bids[k].Agent].RFA*=1-BrokerFees; // Broker fees
                Market[Asks[k].Agent].RFA*=1-BrokerFees; // Broker fees
                Market[Bids[k].Agent].Stocks[Bids[k].Stock].StockQuantity += TempQuantity;
                Market[Asks[k].Agent].Stocks[Asks[k].Stock].StockQuantity -= TempQuantity;
                Market[Bids[k].Agent].TradingWindowClock.push_back(Market[Bids[k].Agent].Future+t); // JJJ3
                Market[Asks[k].Agent].TradingWindowClock.push_back(Market[Asks[k].Agent].Future+t); // JJJ3
                TotalStockQuantityTraded+=TempQuantity;
                Count+=1;
                //for (int i=0; i<NumberOfAgents ; i++) {for (int j=0; j<NumberOfStocks; j++) {if (Market[i].Stocks[j].StockQuantity <0) {outputTest << "Error Z Clear(): TempQuantity=min (Bids[k].Q, Asks[k].Q)=" << TempQuantity << ", Bids[k].Q=" << Bids[k].Q << ", Asks[k].Q=" << Asks[k].Q << ", Market[Bids[k].Agent].Stocks[Bids[k].Stock].StockQuantity=" << QBidsBefore << "->" << Market[Bids[k].Agent].Stocks[Bids[k].Stock].StockQuantity << ", Market[Asks[k].Agent].Stocks[Asks[k].Stock].StockQuantity=" << QAsksBefore << "->" << Market[Asks[k].Agent].Stocks[Asks[k].Stock].StockQuantity << endl;};};};
                // This is for method RL()
                int TransactionPriceRealSize = int(Market[Bids[k].Agent].TransactionPriceReal[j].size());
                int TransactionPriceRealSize2 = int(Market[Asks[k].Agent].TransactionPriceReal[j].size());
                int QuantitiesRealSize = int(Market[Bids[k].Agent].QuantitiesReal[j].size());
                int QuantitiesRealSize2 = int(Market[Asks[k].Agent].QuantitiesReal[j].size());
                Market[Bids[k].Agent].TransactionPriceReal[j][TransactionPriceRealSize-1]=TempClearingPrice;
                Market[Asks[k].Agent].TransactionPriceReal[j][TransactionPriceRealSize2-1]=TempClearingPrice;
                Market[Bids[k].Agent].QuantitiesReal[j][QuantitiesRealSize-1]=TempQuantity;
                Market[Asks[k].Agent].QuantitiesReal[j][QuantitiesRealSize2-1]=-TempQuantity;
                
                // This computes the average of all cleared Bids and all cleared Asks
                //AverageBids+=Bids[k].Price;
                //AverageAsks+=Asks[k].Price;
                
                TempClearingPriceLast=TempClearingPrice;
                TotalStockQuantityTradedLast=TotalStockQuantityTraded;
                CountLast=Count;
                
                // Output
                if (Plot=="On") {
                    outputLog << "Level " << LevelCount << ": Agent " << Bids[k].Agent << " with Q=" << Market[Bids[k].Agent].Stocks[Bids[k].Stock].StockQuantity - TempQuantity << "->" << Market[Bids[k].Agent].Stocks[Bids[k].Stock].StockQuantity << " longs " << TempQuantity << " of Stock " << Bids[k].Stock << " from Agent " << Asks[k].Agent << " with Q=" << Market[Asks[k].Agent].Stocks[Asks[k].Stock].StockQuantity+TempQuantity << "->" << Market[Asks[k].Agent].Stocks[Asks[k].Stock].StockQuantity << " for $" << TempClearingPrice << " (midprice from " << Asks[k].Price << "~" << Bids[k].Price << ")" << endl;
                }; // closes if
                
                // Proper draft of all orders
                if (Bids[k].Q>TempQuantity) {Bids[k].Q-=TempQuantity; Asks.erase(Asks.begin()+k); k-=1;} // Buyer stays in business if still eager to buy more from next seller
                else if (Asks[k].Q>TempQuantity) {Asks[k].Q-=TempQuantity; Bids.erase(Bids.begin()+k); k-=1;}; // Seller stays in business if still eager to sell more to next buyer
                //if (Plot=="On") {outputLog << "Now k=" << k << ", Bids.size()=" << Bids.size() << ", Asks.size()=" << Asks.size() << ", min(Bids.size(), Asks.size())=" << min(Bids.size(), Asks.size()) << endl;};
            } // closes if
            
            else if (!C1) {if (Plot=="On") {outputLog << "Level " << LevelCount << ": not cleared since bid=" << Bids[k].Price << " < ask=" << Asks[k].Price << endl;};}
            else if (!C2) {if (Plot=="On") {outputLog << "Level " << LevelCount << ": not cleared since seller out of stock" << endl;};}
            else if (!C3) {if (Plot=="On") {outputLog << "Level " << LevelCount << ": not cleared since buyer hasn't enough RFA to buy ($" << Market[Bids[k].Agent].RFA << ")" << endl;};}
            else if (!C4) {if (Plot=="On") {outputLog << "Level " << LevelCount << ": not cleared since seller hasn't enough quantities to sell (" << Market[Asks[k].Agent].Stocks[Asks[k].Stock].StockQuantity << ")" << endl;};};
            
            /*
             if ((k<0) || (OBSize<=0)) {
             //if (Plot=="On") {outputLog << "We are in a case (k<0) || (OBSize<=0)" << endl;};
             k=0; OBSize=0;
             //if (Plot=="On") {outputLog << "So we change k=" << k << ", OBSize=" << OBSize << endl;};
             }; // closes if
             */
            LevelCount+=1;
        }; // closes k loop
        
        // Outputting clearing price
        //if ((Count!=0) && ((Plot=="On") || (Plot=="Mini"))) {outputLog << endl << "CLEARING PRICE FOR STOCK " << OBStock << " AT TIME t=" << t << ": $" << TempClearingPrice << endl<< endl << endl;};
        
        // Returning result
        //if (Count>0) {AverageBids/=Count; AverageAsks/=Count;};
        vector<double> Res;
        Res.push_back(TempClearingPriceLast); // Market price given by the latest (lowest) clearing price
        Res.push_back(TotalStockQuantityTradedLast); // Total stock quantity that was traded at this time step t
        Res.push_back(CountLast); // Indicator of emptyness of OB
        //ofstream outputOBLevels(Machine+"OBlevels.csv", ofstream::app);
        //outputOBLevels << CountLast << endl;
        //Res.push_back(AverageBids); // Average of all Bids
        //Res.push_back(AverageAsks); // Average of all Asks
        if (Plot=="On") {outputLog << endl << "TempClearingPriceLast=" << Res[0] << ", TotalStockQuantityTradedLast=" << Res[1] << ", CountLast=" << Res[2] << endl;};
        outputLog.close();
        //outputOBLevels.close();
        return Res;
    }; // closes method Clear()
    
    
    // Selecting a random level in the OB and shifting it to a metaorder (if at least one business order in the OB)
    int MetaorderInjection(vector<Agent> &Market, int SharesOutstanding, int MetaorderImpact, int OBLevelSize) {
        ofstream outputMetalog(Machine+"MetaorderLog.txt", ofstream::app);
        int Res=-1;
        if (OBLevelSize>-1) {
            vector<int> J = Shuffle2(OBLevelSize);
            int MetaorderLevel=J[0]; // OB level where the metaorder will be injected (whether in the ask or the bid order)
            int BidOrAskMetaorder=int(time(0))%2; // If 0 then bid, if 1 then ask
            if ((BidOrAskMetaorder==0) && (OBLevelSize>-1)) {
                Bids[MetaorderLevel].Q=int(MetaorderImpact*SharesOutstanding/100);
                Market[Bids[MetaorderLevel].Agent].RFA+=ceil(Bids[MetaorderLevel].Q*(Bids[MetaorderLevel].Price+Asks[MetaorderLevel].Price)/2); // This agent is made richer in order to buy
                Res=MetaorderLevel;
                outputMetalog << "Bid metaorder injected at OB level " << MetaorderLevel << ": agent " << Bids[MetaorderLevel].Agent << " buys " << Bids[MetaorderLevel].Q << " for £" << Bids[MetaorderLevel].Price << endl;
                MetaorderLastAgent=Bids[MetaorderLevel].Agent; // Agent placing the order
                MetaorderNature=0; // 0 for Bid and 1 for Ask
                Credit=ceil(Bids[MetaorderLevel].Q*(Bids[MetaorderLevel].Price+Asks[MetaorderLevel].Price)/2); // Extra RFA (for a buy) or stock quantity (for a sell) temporarily credited to the agent
            }
            else if ((BidOrAskMetaorder==1) && (OBLevelSize>-1)) {
                Asks[MetaorderLevel].Q=int(MetaorderImpact*SharesOutstanding/100);
                Market[Asks[MetaorderLevel].Agent].Stocks[0].StockQuantity+=Asks[MetaorderLevel].Q+1; // This agent is made richer in order to sell
                Res=MetaorderLevel;
                outputMetalog << "Ask metaorder injected at OB level " << MetaorderLevel << ": agent " << Asks[MetaorderLevel].Agent << " sells " << Asks[MetaorderLevel].Q << " for £" << Asks[MetaorderLevel].Price << endl;
                MetaorderLastAgent=Asks[MetaorderLevel].Agent; // Agent placing the order
                MetaorderNature=1; // 0 for Bid and 1 for Ask
                Credit=double(Asks[MetaorderLevel].Q+1); // Extra RFA (for a buy) or stock quantity (for a sell) temporarily credited to the agent
            };
        };
        outputMetalog.close();
        return Res;
    }; // closes method MetaorderInjection()
    
    
    
    
    
    
}; // closes class StockOrderBook









// This plots the distribution of an STL vector taken as parameter in Mathematica. Setting Xmin=Xmax=0 finds these bounds automatically
vector<vector<double> > Distribution (vector<double> X, int Precision, double Xmin, double Xmax, string Name, string Plot) {
    vector<vector<double> > DistRes;
    // Randomization of output file name
    string A2 = Machine;
    string B2=Name;
    if (Name=="000") {int B = int(1000*((STLRandom(5, "Uniform", "NoPlot"))[2])); B2=IntToString(B);};
    A2+=B2;
    const char * Title = A2.c_str();
    ofstream output5(Title);
    
    int SeriesSize = int(X.size());
    double res;
    vector<double> X1, X2, Xbin, Xgrid; X1=X; X2=X;
    sort (X1.begin(), X1.begin() + SeriesSize); sort (X2.begin(), X2.begin() + SeriesSize);
    // After sorting the series, we take the first and last element as min and max resp.
    if ((Xmin==0) && (Xmax==0)) {Xmin = X1[0]; Xmax = X1[SeriesSize-1];};
    
    // Now generating the distribution itself in Excel
    for (int i=0; i<Precision; i+=1) {
        res=0;
        double BinStep = (Xmax-Xmin)/Precision;
        Xgrid.push_back(Xmin + i*BinStep/Precision);
        if (Plot=="On") {output5 << Xmin + i*BinStep << "\t";};
        for (int j=0; j<SeriesSize; j+=1) {
            if((X2[j] >= Xmin + i*BinStep) && (X2[j] < Xmin + (i+1)*BinStep)) {res=res+1;};
        }
        Xbin.push_back(res);
        if (Plot=="On") {output5 << Xbin[i] << endl;}
    }
    if (Plot=="On") {output5.close();}
    
    DistRes.push_back(Xgrid);
    DistRes.push_back(Xbin);
    return DistRes;
}




// This plots the distribution of an STL vector taken as parameter in Mathematica. Setting Xmin=Xmax=0 finds these bounds automatically
gsl_matrix * GSLDistribution (gsl_matrix * M, int Precision) {
    int SeriesSize1 = int(M->size1);
    int SeriesSize2 = int(M->size2);
    gsl_matrix * D = gsl_matrix_calloc (SeriesSize1*2, Precision);
    double Xmin, Xmax;
    // Determining the min and max of each time series j from gsl_matrix M
    vector<double> X;
    for (int j=0; j<SeriesSize1; j++) {
        for (int t=0; t<SeriesSize2; t++) {X.push_back(gsl_matrix_get(M,j,t));};
        sort (X.begin(), X.begin() + SeriesSize2);
        Xmin = X[0]; Xmax = X[SeriesSize2-1];
        double BinStep = (Xmax-Xmin+1)/Precision; // The +1 term is what fixed the whole range issue
        // Generating the distribution
        for (int k=0; k<Precision; k+=1) {
            double res=0;
            gsl_matrix_set(D, 2*j, k, Xmin+k*BinStep);
            for (int t=0; t<SeriesSize2; t+=1) {
                if ((X[t] >= Xmin + k*BinStep) && (X[t] < Xmin + (k+1)*BinStep)) {res+=1;};
            }; // closes t loop
            gsl_matrix_set(D, 2*j+1, k, res);
        }; // closes k loop
        //gsl_matrix_set(D, 2*j+1, Precision-1, gsl_matrix_get(D, 2*j+1, Precision-1)+1); // Adding one bin count at Xmax
        X.clear();
    }; // closes j loop
    return D;
}













// This function returns a lognormal random walk simulating an asset: dS/S= mu.dt + sig.N(0,1).sqrt(dt) + \mathcal P(k)
vector<double> RandomWalk(double InitialValue, double Drift, double Volatility, int TimeSteps, double Seed) {
    double Skewness=0.00;
    double FundamentalValue=InitialValue/5; // Asset cannot go below 20% of its fundamental value
    Drift=0.00;
    vector<double> S,S2, StandardNormal, Jumps;
    //gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937); // or gsl_rng_default instead of gsl_rng_mt19937
    gsl_rng * r = make_rng();
    gsl_rng_set(r, static_cast<unsigned long int>(Seed)); // gsl_rng_set(const gsl_rng * r, unsigned long int s)
    for (int t=0; t<TimeSteps; t++) {Jumps.push_back(gsl_ran_poisson (r, 10));};
    for (int t=0; t<TimeSteps; t++) {StandardNormal.push_back(Skewness+gsl_ran_gaussian (r, 1));};
    S.push_back(InitialValue); // The +20 is there to ensure we have an entry price that is not negligeable
    S2.push_back(S[0]+FundamentalValue);
    double dt=1;
    int Condition=0;
    for (int t=1; t<TimeSteps; t++) {
        if ((Jumps[t]) > 15) {Condition=1;};
        if ((StandardNormal[t]) >= 0) {Condition *= 1;};
        if ((StandardNormal[t]) < 0) {Condition *= (-1);};
        //S.push_back(floor(S[t-1] + Volatility*(S[t-1])*(StandardNormal[t]) + (S[t-1])*Condition*(Jumps[t])/100));
        S.push_back(ceil((S[t-1])*(1 + Drift*dt + Volatility*(StandardNormal[t])*sqrt(dt) + Condition*(Jumps[t])/100)));
        Condition=0;
        S2.push_back(S[t]+FundamentalValue);
    };
    return S2;
};



// Mu is the poisson-length between each jump, and Sig the intensity of the jump
vector<double> PoissonRandomWalk (double S0, int Mu, int Sig, int Time, double Seed) {
    vector<double> S, VSig; vector<int> VMu;
    //gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937); // or gsl_rng_default instead of gsl_rng_mt19937
    gsl_rng * r = make_rng();
    gsl_rng_set(r, static_cast<unsigned long int>(Seed)); // gsl_rng_set(const gsl_rng * r, unsigned long int s)
    for (int t=0; t<Time; t++) {
        S.push_back(S0);
        VMu.push_back(1+gsl_ran_poisson (r, Mu));
        VSig.push_back((1+gsl_ran_poisson (r, Sig)) * (gsl_ran_gaussian (r, 1)));
    };
    //PlotSTLInt(VMu, "VMu.txt"); PlotSTL(VSig, "VSig.txt", "NoXL");
    int P=VMu[0];
    for (int t=1; t<Time; t++) {
        if (P<=0) {P=VMu[t]; S[t] = S[t-1] + VSig[t];}
        else {P-=1; S[t] = S[t-1];};
        if (S[t]<=0.2*S[0]) {S[t] = S[t-1];}; // Stock cannot go below treshold of 20% of initial value
    };
    return S;
};





// This function returns a lognormal random walk that is cointegrated to another
vector<double> CointegratedWalk(vector<double> Master, double Leash, double LeashVolatility, double Accuracy) {
    double Val;
    vector<double> S, L;
    vector<double> StandardNormal = STLRandom((int(Master.size())), "StandardNormal2", "NoPlot"); // Time series of size Master.size() made of values gsl_ran_gaussian (r, 1)
    S.push_back(Master[0]);
    L.push_back(Leash);
    //gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937); // or gsl_rng_default instead of gsl_rng_mt19937
    gsl_rng * r = make_rng();
    gsl_rng_set(r, static_cast<unsigned long int>(time(0))); //gsl_rng_set(const gsl_rng * r, unsigned long int s)
    for (int t=1; t<(int(Master.size())); t++) {
        L.push_back((L[t-1])*(1 + LeashVolatility*(gsl_ran_gaussian (r, 1)))); // The leash is elastic with its own vol
        Val = S[t-1] + LeashVolatility*(S[t-1])*(StandardNormal[t]);
        // Asset comes back if too far
        if ((Master[t] - Val) > ((L[t])*(Master[t]))) {Val=S[t-1] + LeashVolatility*(S[t-1])*(abs(StandardNormal[t]));};
        if ((Val - Master[t]) > ((L[t])*(Master[t]))) {Val=S[t-1] + LeashVolatility*(S[t-1])*(-abs(StandardNormal[t]));};
        if ((gsl_ran_poisson (r, 1))*(gsl_ran_poisson (r, 1)) >= Accuracy) {Val=Master[t]*(1 + 0.10*StandardNormal[t]);}; //This poisson event is computed as happening 0.2, 0.75, 1, 2 per year for Accuracy=13,11,10,9 on average and brings back the cointegrated walk to the Master time series (after a news for example that reveals TruePresentValue)
        S.push_back(Val);
    };
    return S;
};






double Compare (vector<double> X1, vector<double> X2) {
    int Size = int(min(X1.size(),X2.size()));
    double Result=0;
    for (int i=0; i<Size; i++) {
        Result += abs(X1[i]-X2[i]);
    }; // closes i loop
    return 100*Result/Size;
}; // closes Compare()






// This function returns the Future parameters for a NumberOfAgents and a simulation of time Time, according to a quasi-hyperbolic discount function F with parameters Beta and Delta: F(0)=1, F(t)=Beta*pow(Delta, t) with 0 < Beta,Delta < 1
vector<int> Quasihyperbolic(int NumberOfAgents, int Time, double Beta, double Delta, double Seed) {
    vector<int> Res;
    //gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937); // or gsl_rng_default instead of gsl_rng_mt19937
    gsl_rng * r = make_rng();
    gsl_rng_set(r, static_cast<unsigned long int>(Seed)); // gsl_rng_set(const gsl_rng * r, unsigned long int s)
    for (int i=0; i<NumberOfAgents+10; i++) {
        while (i>0) {
            double x = gsl_rng_uniform (r);
            double y = gsl_rng_uniform (r);
            if (y<=Beta*pow(Delta, floor(10*x))) {Res.push_back(x*Time); break;};
        };
    };
    return Res;
};









// Outputs the agent parameters for the stitch into one file: for each agent, the first line are integers, the second doubles, the third vector<int> Accuracy, the fourth vector<double> Leash, the fifth vector<double> LeashVol, the sixth vector<vector<int>> FPi[j][i], and the seventh vector<vector<int>> TPi[j][i]
void AgentParametersOutput (vector<Agent> &Market, int NumberOfAgents, int NumberOfStocks, string OutputName) {
    string RootName = Machine; RootName+=OutputName;
    ofstream outputParameters(RootName.c_str(), ofstream::app);
    for (int i=0; i<NumberOfAgents ; i++) {
        // int Future, History, TradingWindow, LiquidationFloor, Exploration, NEB, NEB1, NEBLossAversion, NEB3, NEB4, NEBPositivity, NEB6
        outputParameters << Market[i].Future << ", " << Market[i].History << ", " << Market[i].TradingWindow << ", " << Market[i].LiquidationFloor << ", " << Market[i].Exploration << ", " << Market[i].NEB << ", " << Market[i].NEBLossAversion << ", " << Market[i].NEBPositivity << endl;
        
        // double Reflexivity, Versatility, Gesture, Epsilon, Human, NEB1p, NEB3p, NEB4p, NEBPositivity, NEBLearningRate, NEB8, NEB9
        outputParameters << Market[i].Reflexivity << ", " << Market[i].Versatility << ", " << Market[i].Gesture << ", " << Market[i].Epsilon << ", " << Market[i].Human << ", " << Market[i].NEBPositivity << ", " << Market[i].NEBLearningRate << endl;
        
        // vector<int> Accuracy
        for (int j=0; j<NumberOfStocks ; j++) {
            if (j<NumberOfStocks-1) {outputParameters << Market[i].Accuracy[j] << ", ";}
            else {outputParameters << Market[i].Accuracy[j] << endl;};
        }; // closes j loop
        
        // vector<double> Leash
        for (int j=0; j<NumberOfStocks ; j++) {
            if (j<NumberOfStocks-1) {outputParameters << Market[i].Leash[j] << ", ";}
            else {outputParameters << Market[i].Leash[j] << endl;};
        }; // closes j loop
        
        // vector<double> LeashVol
        for (int j=0; j<NumberOfStocks ; j++) {
            if (j<NumberOfStocks-1) {outputParameters << Market[i].LeashVol[j] << ", ";}
            else {outputParameters << Market[i].LeashVol[j] << endl;};
        }; // closes j loop
        
        // vector<vector<int>> FPi[j][i]
        for (int j=0; j<NumberOfStocks ; j++) {
            for (int k=0; k < Market[0].FS * Market[0].FA; k++) {
                if (j<Market[0].FS * Market[0].FA * NumberOfStocks - 1) {outputParameters << Market[i].FPi[j][k] << ", ";}
                else {outputParameters << Market[i].FPi[j][k] << endl;};
            }; // closes k loop
        }; // closes j loop
        
        // vector<vector<int>> TPi[j][i]
        for (int j=0; j<NumberOfStocks ; j++) {
            for (int k=0; k < Market[0].TS * Market[0].TA; k++) {
                if (j<Market[0].TS * Market[0].TA * NumberOfStocks - 1) {outputParameters << Market[i].TPi[j][k] << ", ";}
                else {outputParameters << Market[i].TPi[j][k] << endl;};
            }; // closes k loop
        }; // closes j loop
        
    }; // closes i loop
    
    
    outputParameters.close();
}; // closes OutputParameters()










// Outputs the agent parameters of agents in a line
void AgentParameters (vector<Agent> &Market, int NumberOfAgents, int NumberOfAgentsFull, int NumberOfStocks, int Time) {
    ofstream outputParameters(Machine+"AgentParameters.txt", ofstream::app);
    vector<int> J = Shuffle(NumberOfAgentsFull);
    outputParameters << "I=" << NumberOfAgentsFull << ", J=" << NumberOfStocks << ", T=" << Time << endl;
    outputParameters << "***********************" << endl;
    outputParameters << "Selected agents: i={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << J[i] << ", ";} else {outputParameters << J[i] << "}" << endl;};};
    outputParameters << "Future={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Future << ", ";} else {outputParameters << Market[J[i]].Future << "}" << endl;};};
    outputParameters << "History={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].History << ", ";} else {outputParameters << Market[J[i]].History << "}" << endl;};};
    outputParameters << "TradingWindow={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].TradingWindow << ", ";} else {outputParameters << Market[J[i]].TradingWindow << "}" << endl;};};
    outputParameters << "Exploration={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Exploration << ", ";} else {outputParameters << Market[J[i]].Exploration << "}" << endl;};};
    outputParameters << "Epsilon={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Epsilon << ", ";} else {outputParameters << Market[J[i]].Epsilon << "}" << endl;};};
    outputParameters << "Gesture={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Gesture << ", ";} else {outputParameters << Market[J[i]].Gesture << "}" << endl;};};
    outputParameters << "Versatility={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Versatility << ", ";} else {outputParameters << Market[J[i]].Versatility << "}" << endl;};};
    outputParameters << "Reflexivity={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Reflexivity << ", ";} else {outputParameters << Market[J[i]].Reflexivity << "}" << endl;};};
    //outputParameters << "DiscountFactor={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].DiscountFactor << ", ";} else {outputParameters << Market[J[i]].DiscountFactor << "}" << endl;};};
    outputParameters << "RFA={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].RFAFirst << ", ";} else {outputParameters << Market[J[i]].RFAFirst << "}" << endl;};};
    outputParameters << "RBA={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].RBAFirst << ", ";} else {outputParameters << Market[J[i]].RBAFirst << "}" << endl;};};
    outputParameters << "Bankruptcy={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Bankruptcy << ", ";} else {outputParameters << Market[J[i]].Bankruptcy << "}" << endl;};};
    outputParameters << "LiquidationFloor={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].LiquidationFloor << ", ";} else {outputParameters << Market[J[i]].LiquidationFloor << "}" << endl;};};
    outputParameters << "Human={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Human << ", ";} else {outputParameters << Market[J[i]].Human << "}" << endl;};};
    outputParameters << "NEBLossAversion={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].NEBLossAversion << ", ";} else {outputParameters << Market[J[i]].NEBLossAversion << "}" << endl;};};
    outputParameters << "NEBPositivity={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].NEBPositivity << ", ";} else {outputParameters << Market[J[i]].NEBPositivity << "}" << endl;};};
    outputParameters << "NEBLearningRate={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].NEBLearningRate << ", ";} else {outputParameters << Market[J[i]].NEBLearningRate << "}" << endl;};};
    outputParameters << "NEBPositivity={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].NEBPositivity << ", ";} else {outputParameters << Market[J[i]].NEBPositivity << "}" << endl;};};
    outputParameters << "***********************" << endl;
    outputParameters.close();
}; // closes OutputParameters()






// Outputs the RL variables (mostly vectors) of agents in a line, for a given stock j
void AgentVariables (vector<Agent> &Market, int NumberOfAgents, int NumberOfAgentsFull, int NumberOfStocks, int Time, gsl_matrix * ReflexiveValues, int j) {
    vector<int> J = Shuffle(NumberOfAgentsFull);
    string A = Machine+"AgentVariables_j";
    string B=IntToString(j);
    string C= ".txt";
    A+=B+C;
    const char * Title=A.c_str();
    ofstream outputParameters(Title, ofstream::app);
    outputParameters << "***********************" << "j=" << j << "***********************" << endl;
    outputParameters << "I=" << NumberOfAgentsFull << ", J=" << NumberOfStocks << ", T=" << Time << endl;
    outputParameters << "Selected agents: i={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << J[i] << ", ";} else {outputParameters << J[i] << "}" << endl;};};
    outputParameters << "Q={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Stocks[j].StockQuantity << ", ";} else {outputParameters << Market[J[i]].Stocks[j].StockQuantity << "}" << endl;};};
    outputParameters << "Bid={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Stocks[j].StockBidValue << ", ";} else {outputParameters << Market[J[i]].Stocks[j].StockBidValue << "}" << endl;};};
    outputParameters << "Ask={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Stocks[j].StockAskValue << ", ";} else {outputParameters << Market[J[i]].Stocks[j].StockAskValue << "}" << endl;};};
    outputParameters << "Leash={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Leash[j] << ", ";} else {outputParameters << Market[J[i]].Leash[j] << "}" << endl;};};
    outputParameters << "LeashVol={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].LeashVol[j] << ", ";} else {outputParameters << Market[J[i]].LeashVol[j] << "}" << endl;};};
    outputParameters << "Accuracy={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Accuracy[j] << ", ";} else {outputParameters << Market[J[i]].Accuracy[j] << "}" << endl;};};
    outputParameters << "***********************" << endl;
    outputParameters << "ReflexiveValues={"; for (int t=0; t<Time ; t++) {if (t!=Time-1)  {outputParameters << gsl_matrix_get (ReflexiveValues, j, t) << ", ";} else {outputParameters << gsl_matrix_get (ReflexiveValues, j, t) << "}" << endl;};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "FSIndex[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].FSIndex[j].size()) ; u++) {if (u!=int(Market[J[i]].FSIndex[j].size())-1)  {outputParameters << Market[J[i]].FSIndex[j][u] << ", ";} else {outputParameters << Market[J[i]].FSIndex[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "dRoughVec[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].dRoughVec[j].size()) ; u++) {if (u!=int(Market[J[i]].dRoughVec[j].size())-1)  {outputParameters << Market[J[i]].dRoughVec[j][u] << ", ";} else {outputParameters << Market[J[i]].dRoughVec[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "dSmoothVec[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].dSmoothVec[j].size()) ; u++) {if (u!=int(Market[J[i]].dSmoothVec[j].size())-1)  {outputParameters << Market[J[i]].dSmoothVec[j][u] << ", ";} else {outputParameters << Market[J[i]].dSmoothVec[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "Lvol[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].LvolLog[j].size()) ; u++) {if (u!=int(Market[J[i]].LvolLog[j].size())-1)  {outputParameters << Market[J[i]].LvolLog[j][u] << ", ";} else {outputParameters << Market[J[i]].LvolLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "Svol[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].SvolLog[j].size()) ; u++) {if (u!=int(Market[J[i]].SvolLog[j].size())-1)  {outputParameters << Market[J[i]].SvolLog[j][u] << ", ";} else {outputParameters << Market[J[i]].SvolLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "Rough[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].RoughLog[j].size()) ; u++) {if (u!=int(Market[J[i]].RoughLog[j].size())-1)  {outputParameters << Market[J[i]].RoughLog[j][u] << ", ";} else {outputParameters << Market[J[i]].RoughLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "Smooth[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].SmoothLog[j].size()) ; u++) {if (u!=int(Market[J[i]].SmoothLog[j].size())-1)  {outputParameters << Market[J[i]].SmoothLog[j][u] << ", ";} else {outputParameters << Market[J[i]].SmoothLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "Reflexive[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].ReflexiveLog[j].size()) ; u++) {if (u!=int(Market[J[i]].ReflexiveLog[j].size())-1)  {outputParameters << Market[J[i]].ReflexiveLog[j][u] << ", ";} else {outputParameters << Market[J[i]].ReflexiveLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "FAIndexReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].FAIndexReal[j].size()) ; u++) {if (u!=int(Market[J[i]].FAIndexReal[j].size())-1)  {outputParameters << Market[J[i]].FAIndexReal[j][u] << ", ";} else {outputParameters << Market[J[i]].FAIndexReal[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "ToolReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].ToolRealLog[j].size()) ; u++) {if (u!=int(Market[J[i]].ToolRealLog[j].size())-1)  {outputParameters << Market[J[i]].ToolRealLog[j][u] << ", ";} else {outputParameters << Market[J[i]].ToolRealLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "LagReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].LagRealLog[j].size()) ; u++) {if (u!=int(Market[J[i]].LagRealLog[j].size())-1)  {outputParameters << Market[J[i]].LagRealLog[j][u] << ", ";} else {outputParameters << Market[J[i]].LagRealLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "WeightReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].WeightRealLog[j].size()) ; u++) {if (u!=int(Market[J[i]].WeightRealLog[j].size())-1)  {outputParameters << Market[J[i]].WeightRealLog[j][u] << ", ";} else {outputParameters << Market[J[i]].WeightRealLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "FAIndexVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].FAIndexVirtual[j].size()) ; u++) {if (u!=int(Market[J[i]].FAIndexVirtual[j].size())-1)  {outputParameters << Market[J[i]].FAIndexVirtual[j][u] << ", ";} else {outputParameters << Market[J[i]].FAIndexVirtual[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "ToolVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].ToolVirtualLog[j].size()) ; u++) {if (u!=int(Market[J[i]].ToolVirtualLog[j].size())-1)  {outputParameters << Market[J[i]].ToolVirtualLog[j][u] << ", ";} else {outputParameters << Market[J[i]].ToolVirtualLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "LagVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].LagVirtualLog[j].size()) ; u++) {if (u!=int(Market[J[i]].LagVirtualLog[j].size())-1)  {outputParameters << Market[J[i]].LagVirtualLog[j][u] << ", ";} else {outputParameters << Market[J[i]].LagVirtualLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "WeightVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].WeightVirtualLog[j].size()) ; u++) {if (u!=int(Market[J[i]].WeightVirtualLog[j].size())-1)  {outputParameters << Market[J[i]].WeightVirtualLog[j][u] << ", ";} else {outputParameters << Market[J[i]].WeightVirtualLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "ForecastReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].ForecastReal[j].size()) ; u++) {if (u!=int(Market[J[i]].ForecastReal[j].size())-1)  {outputParameters << Market[J[i]].ForecastReal[j][u] << ", ";} else {outputParameters << Market[J[i]].ForecastReal[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "dFResultVec[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].dFResultVec[j].size()) ; u++) {if (u!=int(Market[J[i]].dFResultVec[j].size())-1)  {outputParameters << Market[J[i]].dFResultVec[j][u] << ", ";} else {outputParameters << Market[J[i]].dFResultVec[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "FResultDisReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].FResultDisReal[j].size()) ; u++) {if (u!=int(Market[J[i]].FResultDisReal[j].size())-1)  {outputParameters << Market[J[i]].FResultDisReal[j][u] << ", ";} else {outputParameters << Market[J[i]].FResultDisReal[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "ForecastVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].ForecastVirtual[j].size()) ; u++) {if (u!=int(Market[J[i]].ForecastVirtual[j].size())-1)  {outputParameters << Market[J[i]].ForecastVirtual[j][u] << ", ";} else {outputParameters << Market[J[i]].ForecastVirtual[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "FResultDisVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].FResultDisVirtual[j].size()) ; u++) {if (u!=int(Market[J[i]].FResultDisVirtual[j].size())-1)  {outputParameters << Market[J[i]].FResultDisVirtual[j][u] << ", ";} else {outputParameters << Market[J[i]].FResultDisVirtual[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "TSIndex[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].TSIndex[j].size()) ; u++) {if (u!=int(Market[J[i]].TSIndex[j].size())-1)  {outputParameters << Market[J[i]].TSIndex[j][u] << ", ";} else {outputParameters << Market[J[i]].TSIndex[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "dMuPosVec[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].dMuPosVec[j].size()) ; u++) {if (u!=int(Market[J[i]].dMuPosVec[j].size())-1)  {outputParameters << Market[J[i]].dMuPosVec[j][u] << ", ";} else {outputParameters << Market[J[i]].dMuPosVec[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "dMuNegVec[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].dMuNegVec[j].size()) ; u++) {if (u!=int(Market[J[i]].dMuNegVec[j].size())-1)  {outputParameters << Market[J[i]].dMuNegVec[j][u] << ", ";} else {outputParameters << Market[J[i]].dMuNegVec[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "dSigVec[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].dSigVec[j].size()) ; u++) {if (u!=int(Market[J[i]].dSigVec[j].size())-1)  {outputParameters << Market[J[i]].dSigVec[j][u] << ", ";} else {outputParameters << Market[J[i]].dSigVec[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "LiquidPercentVec[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].LiquidPercentVec[j].size()) ; u++) {if (u!=int(Market[J[i]].LiquidPercentVec[j].size())-1)  {outputParameters << Market[J[i]].LiquidPercentVec[j][u] << ", ";} else {outputParameters << Market[J[i]].LiquidPercentVec[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "Mu[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].MuLog[j].size()) ; u++) {if (u!=int(Market[J[i]].MuLog[j].size())-1)  {outputParameters << Market[J[i]].MuLog[j][u] << ", ";} else {outputParameters << Market[J[i]].MuLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "Sig[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].SigLog[j].size()) ; u++) {if (u!=int(Market[J[i]].SigLog[j].size())-1)  {outputParameters << Market[J[i]].SigLog[j][u] << ", ";} else {outputParameters << Market[J[i]].SigLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "RFALevel[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].RFALog[j].size()) ; u++) {if (u!=int(Market[J[i]].RFALog[j].size())-1)  {outputParameters << Market[J[i]].RFALog[j][u] << ", ";} else {outputParameters << Market[J[i]].RFALog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "RBALevel[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].RBALog[j].size()) ; u++) {if (u!=int(Market[J[i]].RBALog[j].size())-1)  {outputParameters << Market[J[i]].RBALog[j][u] << ", ";} else {outputParameters << Market[J[i]].RBALog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "Liquidity[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].LiquidityLog[j].size()) ; u++) {if (u!=int(Market[J[i]].LiquidityLog[j].size())-1)  {outputParameters << Market[J[i]].LiquidityLog[j][u] << ", ";} else {outputParameters << Market[J[i]].LiquidityLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "TAIndexReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].TAIndexReal[j].size()) ; u++) {if (u!=int(Market[J[i]].TAIndexReal[j].size())-1)  {outputParameters << Market[J[i]].TAIndexReal[j][u] << ", ";} else {outputParameters << Market[J[i]].TAIndexReal[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "QuantReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].QuantRealLog[j].size()) ; u++) {if (u!=int(Market[J[i]].QuantRealLog[j].size())-1)  {outputParameters << Market[J[i]].QuantRealLog[j][u] << ", ";} else {outputParameters << Market[J[i]].QuantRealLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "PinchReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].PinchRealLog[j].size()) ; u++) {if (u!=int(Market[J[i]].PinchRealLog[j].size())-1)  {outputParameters << Market[J[i]].PinchRealLog[j][u] << ", ";} else {outputParameters << Market[J[i]].PinchRealLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "TAIndexVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].TAIndexVirtual[j].size()) ; u++) {if (u!=int(Market[J[i]].TAIndexVirtual[j].size())-1)  {outputParameters << Market[J[i]].TAIndexVirtual[j][u] << ", ";} else {outputParameters << Market[J[i]].TAIndexVirtual[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "QuantVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].QuantVirtualLog[j].size()) ; u++) {if (u!=int(Market[J[i]].QuantVirtualLog[j].size())-1)  {outputParameters << Market[J[i]].QuantVirtualLog[j][u] << ", ";} else {outputParameters << Market[J[i]].QuantVirtualLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "PinchVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].PinchVirtualLog[j].size()) ; u++) {if (u!=int(Market[J[i]].PinchVirtualLog[j].size())-1)  {outputParameters << Market[J[i]].PinchVirtualLog[j][u] << ", ";} else {outputParameters << Market[J[i]].PinchVirtualLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "QuantitiesReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].QuantitiesReal[j].size()) ; u++) {if (u!=int(Market[J[i]].QuantitiesReal[j].size())-1)  {outputParameters << Market[J[i]].QuantitiesReal[j][u] << ", ";} else {outputParameters << Market[J[i]].QuantitiesReal[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "QuantitiesVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].QuantitiesVirtual[j].size()) ; u++) {if (u!=int(Market[J[i]].QuantitiesVirtual[j].size())-1)  {outputParameters << Market[J[i]].QuantitiesVirtual[j][u] << ", ";} else {outputParameters << Market[J[i]].QuantitiesVirtual[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "TransactionPriceReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].TransactionPriceReal[j].size()) ; u++) {if (u!=int(Market[J[i]].TransactionPriceReal[j].size())-1)  {outputParameters << Market[J[i]].TransactionPriceReal[j][u] << ", ";} else {outputParameters << Market[J[i]].TransactionPriceReal[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "TransactionPriceVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].TransactionPriceVirtual[j].size()) ; u++) {if (u!=int(Market[J[i]].TransactionPriceVirtual[j].size())-1)  {outputParameters << Market[J[i]].TransactionPriceVirtual[j][u] << ", ";} else {outputParameters << Market[J[i]].TransactionPriceVirtual[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    double dTRmean=0;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "dTResultVec[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].dTResultVec[j].size()) ; u++) {if (u!=int(Market[J[i]].dTResultVec[j].size())-1)  {outputParameters << Market[J[i]].dTResultVec[j][u] << ", ";} else {outputParameters << Market[J[i]].dTResultVec[j][u] << "}";}; dTRmean+=Market[J[i]].dTResultVec[j][u]/int(Market[J[i]].dTResultVec[j].size());}; outputParameters << " => <" << dTRmean << ">" << endl;}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "TResultDisReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].TResultDisReal[j].size()) ; u++) {if (u!=int(Market[J[i]].TResultDisReal[j].size())-1)  {outputParameters << Market[J[i]].TResultDisReal[j][u] << ", ";} else {outputParameters << Market[J[i]].TResultDisReal[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "TResultDisVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].TResultDisVirtual[j].size()) ; u++) {if (u!=int(Market[J[i]].TResultDisVirtual[j].size())-1)  {outputParameters << Market[J[i]].TResultDisVirtual[j][u] << ", ";} else {outputParameters << Market[J[i]].TResultDisVirtual[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;     outputParameters.close();
    
}; // closes OutputParameters()



gsl_matrix *  StockMoments (gsl_matrix * ReflexiveValues, gsl_matrix * Spread, gsl_matrix * TotalStockQuantityTraded, int Time, int Lag, int Delta, int CorrLag, int j) {
    gsl_matrix * StockIndexMoments = gsl_matrix_calloc (12, Time);
    gsl_matrix * StocksDistributions = gsl_matrix_calloc (24, 1000);
    gsl_matrix * AutoCorrLagP = gsl_matrix_alloc (24, Time);
    vector<double> RealizedVolatility, SpreadV, Volume, Stock;
    for (int t=0; t<Time; t++) {Stock.push_back(gsl_matrix_get(ReflexiveValues, j, t));};
    vector<vector<double> > MomentsStockIndex = DynamicMoments (Stock, 0, int(Stock.size()), Lag, Delta);
    for (int t=0; t<Time; t++) {
        gsl_matrix_set(StockIndexMoments, 0, t, Stock[t]);                              // Raw Stock
        gsl_matrix_set(StockIndexMoments, 1, t, MomentsStockIndex[0][t]); // Mean of Stock
        gsl_matrix_set(StockIndexMoments, 2, t, MomentsStockIndex[4][t]); // StdDev of Stock
        gsl_matrix_set(StockIndexMoments, 3, t, MomentsStockIndex[2][t]); // Skewness of Stock
        gsl_matrix_set(StockIndexMoments, 4, t, MomentsStockIndex[3][t]); // Kurtosis of Stock
        gsl_matrix_set(StockIndexMoments, 5, t, MomentsStockIndex[5][t]); // Stock + 3 Sigma
        gsl_matrix_set(StockIndexMoments, 6, t, MomentsStockIndex[6][t]); // Stock - 3 Sigma
        gsl_matrix_set(StockIndexMoments, 7, t, MomentsStockIndex[7][t]); // Realized volatility of Stock
        gsl_matrix_set(StockIndexMoments, 8, t, MomentsStockIndex[8][t]); // Autocorrelations of lag Delta in bsp for Stock
        gsl_matrix_set(StockIndexMoments, 9, t, MomentsStockIndex[9][t]); // Squarred of daily log returns
        gsl_matrix_set(StockIndexMoments, 10, t, MomentsStockIndex[10][t]); // Squarred of weekly log returns
        gsl_matrix_set(StockIndexMoments, 11, t, MomentsStockIndex[11][t]); // Squarred of bi-weekly log returns
        RealizedVolatility.push_back(gsl_matrix_get(StockIndexMoments, 7, t));
        SpreadV.push_back(gsl_matrix_get(Spread, j, t));
        Volume.push_back(gsl_matrix_get(TotalStockQuantityTraded, j, t));
    };
    StocksDistributions = GSLDistribution(StockIndexMoments, 1000);
    AutoCorrLagP = AutocorrelationLagP (StockIndexMoments);
    string B=IntToString(j); string C= ".xls";
    string A = "Moments_j"; A+=B+C;
    PlotGSLMatrix(StockIndexMoments, A.c_str(), 1);
    A = "Distributions_j"; A+=B+C;
    PlotGSLMatrix(StocksDistributions, A.c_str(), 1);
    A = "AutoCorrLagP_j"; A+=B+C;
    PlotGSLMatrix(AutoCorrLagP, A.c_str(), 1);
    
    // Outputing correlations
    vector<double> CorrRealizedVolatilityVolume = DynamicCorrelation (RealizedVolatility, Volume, CorrLag);
    vector<double> CorrRealizedVolatilitySpread = DynamicCorrelation (RealizedVolatility, SpreadV, CorrLag);
    vector<double> CorrSpreadVolume = DynamicCorrelation (SpreadV, Volume, CorrLag);
    int Minsize = min(int(CorrRealizedVolatilityVolume.size()), min(int(CorrRealizedVolatilitySpread.size()), int(CorrSpreadVolume.size())));
    gsl_matrix * Correlations = gsl_matrix_alloc (3, Minsize);
    for (int t=0; t<Minsize; t++) {
        gsl_matrix_set (Correlations, 0, t, CorrRealizedVolatilityVolume[t]);
        gsl_matrix_set (Correlations, 1, t, CorrRealizedVolatilitySpread[t]);
        gsl_matrix_set (Correlations, 2, t, CorrSpreadVolume[t]);
    };
    A = "Correlations_j"; A+=B+C;
    PlotGSLMatrix(Correlations, A.c_str(), 1);
    
    return StockIndexMoments;
}; // closes function StockMoments()










// Policy distances (NON-MULTIVARIATE!)
// The simulation can output matrices of dimension IxI of the average of the absolute differences between each value P(s,a) of two agents. This is what we call Policy Distance (PD): FPolicyDistances, TPolicyDistances, MostSuccessfulFPolicyDistances, MostSuccessfulTPolicyDistances. But the outputs are 3d plots that are hard to read. So we should create a measure of policy distance not between an agent i1 and an agent i2, but between an agent i1 and many other agents (either all others or just a subgroup, such as the most successful percentile) as the average all its policy distances with all other agents, and then we should rank these agents by their increasing Average Policy Distances (APD). Then we know that an agent with a large APD compared to the market (or subset thereof) Average of other APD's (AAPD) has a largely different policy than the rest of the market or a subset thereof. Then we can extract features about the parameters of these most different or most similar agents to market or submarket.
gsl_matrix * PolicyDistances (vector<Agent> Market, int NumberOfAgents, int NumberOfStocks, int Time, int Percentile, int t) {
    // Most successful and less successful agents
    vector<double> CapitalEvolution;
    for (int i=0; i<NumberOfAgents ; i++) {
        double x = Market[i].Capital();
        double y = Market[i].RFAFirst + Market[i].RBAFirst;
        CapitalEvolution.push_back(100*x/y);
    }; // closes i loop
    vector<int> DescendingRank;
    for (int k=0; k<NumberOfAgents; k++) {
        double Bar=0; int BarX=0;
        for (int i=0; i<NumberOfAgents; i++) {
            if ((CapitalEvolution[i]>=Bar) && (CapitalEvolution[i]>=0)) {Bar=CapitalEvolution[i]; BarX=i;};
        }; // closes i loop
        CapitalEvolution[BarX]=-1;
        DescendingRank.push_back(BarX);
    }; // closes k loop
    int AgentNumberPercentile=Percentile*NumberOfAgents/100;
    if (AgentNumberPercentile<=5) {AgentNumberPercentile=5;};
    
    gsl_matrix * Res = gsl_matrix_calloc (10, AgentNumberPercentile);
    gsl_matrix * FPolicyDistancesMostSuccessfulToMarket = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix * TPolicyDistancesMostSuccessfulToMarket = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix * FPolicyDistancesLessSuccessfulToMarket = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix * TPolicyDistancesLessSuccessfulToMarket = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix * FPolicyDistancesMostSuccessfulToMostSuccessful = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix * TPolicyDistancesMostSuccessfulToMostSuccessful = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix * FPolicyDistancesLessSuccessfulToLessSuccessful = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix * TPolicyDistancesLessSuccessfulToLessSuccessful = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix * FPolicyDistancesMostSuccessfulToLessSuccessful = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix * TPolicyDistancesMostSuccessfulToLessSuccessful = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix * FPolicyDistancesLessSuccessfulToMostSuccessful = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix * TPolicyDistancesLessSuccessfulToMostSuccessful = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix * FPolicyDistancesM = gsl_matrix_calloc (NumberOfAgents, NumberOfAgents);
    gsl_matrix * TPolicyDistancesM = gsl_matrix_calloc (NumberOfAgents, NumberOfAgents);
    
    for (int j=0; j<NumberOfStocks; j++) {
        // Computing FPolicyDistancesM and TPolicyDistancesM as the IxI matrices of all PD of each agent with each agent
        int SizeFPi = int(Market[0].FPi[j].size());
        int SizeTPi = int(Market[0].TPi[j].size());
        for (int i1=0; i1<NumberOfAgents; i1++) {
            for (int i2=i1+1; i2<NumberOfAgents; i2++) {
                cout << "****************************" << endl;
                cout << "** STOCK MARKET SIMULATOR **" << endl;
                cout << "****************************" << endl;
                cout << "I=" << NumberOfAgents << ", J=" << NumberOfStocks << ", T=" << Time << endl  << endl;
                cout.setf(ios::fixed);
                cout << "PD step             : " << i1*NumberOfAgents + i2 << "/" << NumberOfAgents*NumberOfAgents << " (" << 100.0*(i1*NumberOfAgents + i2)/(NumberOfAgents*NumberOfAgents) << "%)" << endl;
                
                double ResultF=0; double ResultT=0;
                for (int k=0; k<SizeFPi; k++) {ResultF += abs(Market[i1].FPi[j][k] - Market[i2].FPi[j][k]);};
                gsl_matrix_set (FPolicyDistancesM, i1, i2, 100*ResultF/SizeFPi);
                gsl_matrix_set (FPolicyDistancesM, i2, i1, 100*ResultF/SizeFPi);
                for (int k=0; k<SizeTPi; k++) {ResultT += abs(Market[i1].TPi[j][k] - Market[i2].TPi[j][k]);};
                gsl_matrix_set (TPolicyDistancesM, i1, i2, 100*ResultT/SizeTPi);
                gsl_matrix_set (TPolicyDistancesM, i2, i1, 100*ResultT/SizeTPi);
            }; // closes i2 loop
        }; // closes i1 loop
        
        // Computing FPolicyDistancesMostSuccessfulToMarket and TPolicyDistancesMostSuccessfulToMarket as matrices of j APD of most successful agents with all agents
        for (int i1=0; i1<AgentNumberPercentile; i1++) {
            double ResultF=0; double ResultT=0;
            for (int i2=0; i2<NumberOfAgents; i2++) {
                ResultF += gsl_matrix_get (FPolicyDistancesM, Market[DescendingRank[i1]].AgentName, i2);
                ResultT += gsl_matrix_get (TPolicyDistancesM, Market[DescendingRank[i1]].AgentName, i2);
            }; // closes i2 loop
            gsl_matrix_set (FPolicyDistancesMostSuccessfulToMarket, i1, j, ResultF/(NumberOfAgents-1));
            gsl_matrix_set (TPolicyDistancesMostSuccessfulToMarket, i1, j, ResultT/(NumberOfAgents-1));
        }; // closes i1 loop
        
        // Computing FPolicyDistancesLessSuccessfulToMarket and TPolicyDistancesLessSuccessfulToMarket as matrices of j APD of less successful agents with all agents
        for (int i1=0; i1<AgentNumberPercentile; i1++) {
            double ResultF=0; double ResultT=0;
            for (int i2=0; i2<NumberOfAgents; i2++) {
                ResultF += gsl_matrix_get (FPolicyDistancesM, Market[DescendingRank[i1+NumberOfAgents-AgentNumberPercentile]].AgentName, i2);
                ResultT += gsl_matrix_get (TPolicyDistancesM, Market[DescendingRank[i1+NumberOfAgents-AgentNumberPercentile]].AgentName, i2);
            }; // closes i2 loop
            gsl_matrix_set (FPolicyDistancesLessSuccessfulToMarket, i1, j, ResultF/(NumberOfAgents-1));
            gsl_matrix_set (TPolicyDistancesLessSuccessfulToMarket, i1, j, ResultT/(NumberOfAgents-1));
        }; // closes i1 loop
        
        // Computing FPolicyDistancesLessSuccessfulToMarket and TPolicyDistancesLessSuccessfulToMarket as matrices of j APD of most successful agents with most successful agents
        for (int i1=0; i1<AgentNumberPercentile; i1++) {
            double ResultF=0; double ResultT=0;
            for (int i2=0; i2<AgentNumberPercentile; i2++) {
                ResultF += gsl_matrix_get (FPolicyDistancesM, Market[DescendingRank[i1]].AgentName, Market[DescendingRank[i2]].AgentName);
                ResultT += gsl_matrix_get (TPolicyDistancesM, Market[DescendingRank[i1]].AgentName, Market[DescendingRank[i2]].AgentName);
            }; // closes i2 loop
            gsl_matrix_set (FPolicyDistancesMostSuccessfulToMostSuccessful, i1, j, ResultF/(AgentNumberPercentile-1));
            gsl_matrix_set (TPolicyDistancesMostSuccessfulToMostSuccessful, i1, j, ResultT/(AgentNumberPercentile-1));
        }; // closes i1 loop
        
        // Computing FPolicyDistancesLessSuccessfulToMarket and TPolicyDistancesLessSuccessfulToMarket as matrices of j APD of less successful agents with less successful agents
        for (int i1=0; i1<AgentNumberPercentile; i1++) {
            double ResultF=0; double ResultT=0;
            for (int i2=0; i2<AgentNumberPercentile; i2++) {
                ResultF += gsl_matrix_get (FPolicyDistancesM, Market[DescendingRank[i1+NumberOfAgents-AgentNumberPercentile]].AgentName, Market[DescendingRank[i2+NumberOfAgents-AgentNumberPercentile]].AgentName);
                ResultT += gsl_matrix_get (TPolicyDistancesM, Market[DescendingRank[i1+NumberOfAgents-AgentNumberPercentile]].AgentName, Market[DescendingRank[i2+NumberOfAgents-AgentNumberPercentile]].AgentName);
            }; // closes i2 loop
            gsl_matrix_set (FPolicyDistancesLessSuccessfulToLessSuccessful, i1, j, ResultF/(AgentNumberPercentile-1));
            gsl_matrix_set (TPolicyDistancesLessSuccessfulToLessSuccessful, i1, j, ResultT/(AgentNumberPercentile-1));
        }; // closes i1 loop
        
        // Computing FPolicyDistancesLessSuccessfulToMarket and TPolicyDistancesLessSuccessfulToMarket as matrices of j APD of most successful agents with less successful agents
        for (int i1=0; i1<AgentNumberPercentile; i1++) {
            double ResultF=0; double ResultT=0;
            for (int i2=0; i2<AgentNumberPercentile; i2++) {
                ResultF += gsl_matrix_get (FPolicyDistancesM, Market[DescendingRank[i1]].AgentName, Market[DescendingRank[i2+NumberOfAgents-AgentNumberPercentile]].AgentName);
                ResultT += gsl_matrix_get (TPolicyDistancesM, Market[DescendingRank[i1]].AgentName, Market[DescendingRank[i2+NumberOfAgents-AgentNumberPercentile]].AgentName);
            }; // closes i2 loop
            gsl_matrix_set (FPolicyDistancesMostSuccessfulToLessSuccessful, i1, j, ResultF/(AgentNumberPercentile-1));
            gsl_matrix_set (TPolicyDistancesMostSuccessfulToLessSuccessful, i1, j, ResultT/(AgentNumberPercentile-1));
        }; // closes i1 loop
        
        // Computing FPolicyDistancesLessSuccessfulToMarket and TPolicyDistancesLessSuccessfulToMarket as matrices of j APD of less successful agents with most successful agents
        for (int i1=0; i1<AgentNumberPercentile; i1++) {
            double ResultF=0; double ResultT=0;
            for (int i2=0; i2<AgentNumberPercentile; i2++) {
                ResultF += gsl_matrix_get (FPolicyDistancesM, Market[DescendingRank[i1+NumberOfAgents-AgentNumberPercentile]].AgentName, Market[DescendingRank[i2]].AgentName);
                ResultT += gsl_matrix_get (TPolicyDistancesM, Market[DescendingRank[i1+NumberOfAgents-AgentNumberPercentile]].AgentName, Market[DescendingRank[i2]].AgentName);
            }; // closes i2 loop
            gsl_matrix_set (FPolicyDistancesLessSuccessfulToMostSuccessful, i1, j, ResultF/(AgentNumberPercentile-1));
            gsl_matrix_set (TPolicyDistancesLessSuccessfulToMostSuccessful, i1, j, ResultT/(AgentNumberPercentile-1));
        }; // closes i1 loop
        
    }; // closes j loop
    
    for (int k=0; k<AgentNumberPercentile; k++) {
        gsl_matrix_set (Res, 0, k, gsl_matrix_get (FPolicyDistancesMostSuccessfulToMostSuccessful, k, 0));
        gsl_matrix_set (Res, 1, k, gsl_matrix_get (FPolicyDistancesMostSuccessfulToMarket, k, 0));
        gsl_matrix_set (Res, 2, k, gsl_matrix_get (FPolicyDistancesMostSuccessfulToLessSuccessful, k, 0));
        gsl_matrix_set (Res, 3, k, gsl_matrix_get (FPolicyDistancesLessSuccessfulToMarket, k, 0));
        gsl_matrix_set (Res, 4, k, gsl_matrix_get (FPolicyDistancesLessSuccessfulToLessSuccessful, k, 0));
        gsl_matrix_set (Res, 5, k, gsl_matrix_get (TPolicyDistancesMostSuccessfulToMostSuccessful, k, 0));
        gsl_matrix_set (Res, 6, k, gsl_matrix_get (TPolicyDistancesMostSuccessfulToMarket, k, 0));
        gsl_matrix_set (Res, 7, k, gsl_matrix_get (TPolicyDistancesMostSuccessfulToLessSuccessful, k, 0));
        gsl_matrix_set (Res, 8, k, gsl_matrix_get (TPolicyDistancesLessSuccessfulToMarket, k, 0));
        gsl_matrix_set (Res, 9, k, gsl_matrix_get (TPolicyDistancesLessSuccessfulToLessSuccessful, k, 0));
    }; // closes k loop
    PlotGSLMatrix(Res, "PDistances.xls", 1);
    
    // Matrix memory freeing
    gsl_matrix_free(FPolicyDistancesM);
    gsl_matrix_free(FPolicyDistancesMostSuccessfulToMostSuccessful);
    gsl_matrix_free(FPolicyDistancesMostSuccessfulToMarket);
    gsl_matrix_free(FPolicyDistancesMostSuccessfulToLessSuccessful);
    gsl_matrix_free(FPolicyDistancesLessSuccessfulToMarket);
    gsl_matrix_free(FPolicyDistancesLessSuccessfulToLessSuccessful);
    gsl_matrix_free(TPolicyDistancesM);
    gsl_matrix_free(TPolicyDistancesMostSuccessfulToMostSuccessful);
    gsl_matrix_free(TPolicyDistancesMostSuccessfulToMarket);
    gsl_matrix_free(TPolicyDistancesMostSuccessfulToLessSuccessful);
    gsl_matrix_free(TPolicyDistancesLessSuccessfulToMarket);
    gsl_matrix_free(TPolicyDistancesLessSuccessfulToLessSuccessful);
    
    return Res;
    
};




void PoliciesAverages (vector<Agent> Market, int NumberOfAgents, int NumberOfStocks, int Percentile) {
    vector<gsl_matrix *> Res;
    // Most successful and less successful agents
    vector<double> CapitalEvolution;
    for (int i=0; i<NumberOfAgents ; i++) {
        double x = Market[i].Capital();
        double y = Market[i].RFAFirst + Market[i].RBAFirst;
        CapitalEvolution.push_back(100*x/y);
    }; // closes i loop
    vector<int> DescendingRank;
    for (int k=0; k<NumberOfAgents; k++) {
        double Bar=0; int BarX=0;
        for (int i=0; i<NumberOfAgents; i++) {
            if ((CapitalEvolution[i]>=Bar) && (CapitalEvolution[i]>=0)) {Bar=CapitalEvolution[i]; BarX=i;};
        }; // closes i loop
        CapitalEvolution[BarX]=-1;
        DescendingRank.push_back(BarX);
    }; // closes k loop
    int AgentNumberPercentile=Percentile*NumberOfAgents/100;
    if (AgentNumberPercentile<=5) {AgentNumberPercentile=5;};
    
    gsl_matrix * FPoliciesAverageBest = gsl_matrix_calloc (Market[0].FA, Market[0].FS); // Heat maps
    gsl_matrix * FPoliciesAverageWorst = gsl_matrix_calloc (Market[0].FA, Market[0].FS); // Heat maps
    gsl_matrix * TPoliciesAverageBest = gsl_matrix_calloc (Market[0].TA, Market[0].TS); // Heat maps
    gsl_matrix * TPoliciesAverageWorst = gsl_matrix_calloc (Market[0].TA, Market[0].TS); // Heat maps
    
    for (int j=0; j<NumberOfStocks; j++) {
        // Computing the F() policies averages as heat maps for Best and Worst populations
        for (int n=0; n<AgentNumberPercentile; n++) {
            for (int i=0; i<(Market[0].FS)*(Market[0].FA); i++) {
                gsl_matrix_set (FPoliciesAverageBest, i%Market[0].FA, i/Market[0].FA, Market[DescendingRank[n]].FPi[j][i]/AgentNumberPercentile + gsl_matrix_get (FPoliciesAverageBest, i%Market[0].FA, i/Market[0].FA));
            }; // closes i loop
            for (int i=0; i<(Market[0].FS)*(Market[0].FA); i++) {
                gsl_matrix_set (FPoliciesAverageWorst, i%Market[0].FA, i/Market[0].FA, Market[DescendingRank[n+NumberOfAgents-AgentNumberPercentile]].FPi[j][i]/AgentNumberPercentile + gsl_matrix_get (FPoliciesAverageWorst, i%Market[0].FA, i/Market[0].FA));
            }; // closes i loop
        }; // closes n loop
        
        // Computing the T() policies averages as heat maps for Best and Worst populations
        for (int n=0; n<AgentNumberPercentile; n++) {
            for (int i=0; i<(Market[0].TS)*(Market[0].TA); i++) {
                gsl_matrix_set (TPoliciesAverageBest, i%Market[0].TA, i/Market[0].TA, Market[DescendingRank[n]].TPi[j][i]/AgentNumberPercentile + gsl_matrix_get (TPoliciesAverageBest, i%Market[0].TA, i/Market[0].TA));
            }; // closes i loop
            for (int i=0; i<(Market[0].TS)*(Market[0].TA); i++) {
                gsl_matrix_set (TPoliciesAverageWorst, i%Market[0].TA, i/Market[0].TA, Market[DescendingRank[n+NumberOfAgents-AgentNumberPercentile]].TPi[j][i]/AgentNumberPercentile + gsl_matrix_get (TPoliciesAverageWorst, i%Market[0].TA, i/Market[0].TA));
            }; // closes i loop
        }; // closes n loop
        
    }; // closes j loop
    
    PlotGSLMatrix(FPoliciesAverageBest, "FPoliciesAverageBest.xls", 1);
    PlotGSLMatrix(FPoliciesAverageWorst, "FPoliciesAverageWorst", 1);
    PlotGSLMatrix(TPoliciesAverageBest, "TPoliciesAverageBest.xls", 1);
    PlotGSLMatrix(TPoliciesAverageWorst, "TPoliciesAverageWorst", 1);
    
    // Matrix memory freeing
    gsl_matrix_free(FPoliciesAverageBest);
    gsl_matrix_free(FPoliciesAverageWorst);
    gsl_matrix_free(TPoliciesAverageBest);
    gsl_matrix_free(TPoliciesAverageWorst);
    
};



// Function to Simulate a given Market over a given number of time steps
//pNEBLossAversion24 => pNEBLossAversion, pNEBDelayDiscounting56=>pNEBPositivity, pNEB89=>pNEBDelayDiscounting
vector<gsl_matrix *> MarketSimulator (int NumberOfAgents, int NumberOfStocks, int Time, double Rate, string Plot, string PDCondition, string TypeNEB, int HPGesture, double HPTrueMu, int HPAccuracy, int LiquidationFloor, string LeaderType, int ClusterLimit, int pNEB, double LearningRateScale, int s, int Trunk, int MetaorderImpact) {
    
    int pBots=100; int pNEBLossAversion=0; int pNEBPositivity=0; int pNEBNegativity=0; int pNEBDelayDiscounting=0; int pNEBFear=0; int pNEBFear2=0; int pNEBGreed=0; int pNEBLearningRate=0;
    if (TypeNEB=="Classic") {pBots=80; pNEBLossAversion=3; pNEBPositivity=3; pNEBNegativity=3; pNEBDelayDiscounting=3; pNEBFear=3; pNEBGreed=3; pNEBLearningRate=2;}
    else if (TypeNEB=="Algorithmic") {pBots=100; pNEBLossAversion=0; pNEBPositivity=0; pNEBNegativity=0; pNEBDelayDiscounting=0; pNEBFear=0; pNEBGreed=0; pNEBLearningRate=0;}
    else if (TypeNEB=="Human") {pBots=100-pNEB; pNEBLossAversion=pNEB-(6*pNEB/7); pNEBPositivity=pNEB/7; pNEBNegativity=pNEB/7; pNEBDelayDiscounting=pNEB/7; pNEBFear=pNEB/7; pNEBGreed=pNEB/7; pNEBLearningRate=pNEB/7;}
    else if (TypeNEB=="LossAversion") {pBots=100-pNEB; pNEBLossAversion=pNEB; pNEBPositivity=0; pNEBNegativity=0; pNEBDelayDiscounting=0; pNEBFear=0; pNEBGreed=0; pNEBLearningRate=0;}
    else if (TypeNEB=="Positivity") {pBots=100-pNEB; pNEBLossAversion=0; pNEBPositivity=pNEB; pNEBNegativity=0; pNEBDelayDiscounting=0; pNEBFear=0; pNEBGreed=0; pNEBLearningRate=0;}
    else if (TypeNEB=="Negativity") {pBots=100-pNEB; pNEBLossAversion=0; pNEBPositivity=0; pNEBNegativity=pNEB; pNEBDelayDiscounting=0; pNEBFear=0; pNEBGreed=0; pNEBLearningRate=0;}
    else if (TypeNEB=="DelayDiscounting") {pBots=100-pNEB; pNEBLossAversion=0; pNEBPositivity=0; pNEBNegativity=0; pNEBDelayDiscounting=pNEB; pNEBFear=0; pNEBGreed=0; pNEBLearningRate=0;}
    else if (TypeNEB=="Fear") {pBots=100-pNEB; pNEBLossAversion=0; pNEBPositivity=0; pNEBNegativity=0; pNEBDelayDiscounting=0; pNEBFear=pNEB; pNEBGreed=0; pNEBLearningRate=0;}
    else if (TypeNEB=="Fear2") {pBots=100-pNEB; pNEBLossAversion=0; pNEBPositivity=0; pNEBNegativity=0; pNEBDelayDiscounting=0; pNEBFear=0; pNEBFear2=pNEB; pNEBGreed=0; pNEBLearningRate=0;}
    else if (TypeNEB=="Greed") {pBots=100-pNEB; pNEBLossAversion=0; pNEBPositivity=0; pNEBNegativity=0; pNEBDelayDiscounting=0; pNEBFear=0; pNEBGreed=pNEB; pNEBLearningRate=0;}
    else if (TypeNEB=="LearningRate") {pBots=100-pNEB; pNEBLossAversion=0; pNEBPositivity=0; pNEBNegativity=0; pNEBDelayDiscounting=0; pNEBFear=0; pNEBGreed=0; pNEBLearningRate=pNEB;};
    
    
    
    string VersatilityCondition="Off"; string TradingFrequencyCond="On"; int Percentile=10; string OutputName="OutputCondOff";
    // ofstream outputDebug("/Users/admin/Documents/GNT/SYMBA/Debug.txt", ofstream::app);
    ofstream outputLog(Machine+"SimLog.txt", ofstream::app);
    ofstream outputMetalog(Machine+"MetaorderLog.txt", ofstream::app);
    ofstream outputClustered(Machine+"ClusteredBestWorstAgent.csv", ofstream::app);
    if (Plot=="On") {
        outputLog << "**********SUM UP OF SIMULATION**********" << endl;
        outputLog << "SS1: Initialization at t=0" << endl;
        outputLog << "SS2: Start of t-loop t>=1" << endl;
        outputLog << "SS3: Computing ReflexiveValues" << endl;
        //outputLog << "SS4: Computing ReferenceValues" << endl;
        outputLog << "SS5: Computing interest rates" << endl;
        outputLog << "SS6: RL()" << endl;
        outputLog << "SS7: MarketEvolution.push_back(Market)" << endl;
        outputLog << "SS8: Filling of OB" << endl;
        outputLog << "SS9: Sort()" << endl;
        outputLog << "SS10: Output of OB" << endl;
        outputLog << "SS11: Clear()" << endl;
        outputLog << "SS12: Market indicators" << endl;
        outputLog << "SS13: End of t-loop" << endl;
        outputLog << "*********************************************" << endl << endl;
    };
    
    /*
     int t=0;
     string N1="/Users/admin/Documents/GNT/SYMBA/SimLog";
     string N2=IntToString(t);
     string N3=".txt";
     N1+=N2+N3;
     const char * Title = N1.c_str();
     ofstream outputLog(Title, ofstream::app);
     */
    // INITIALIZING CLOCK (SS1)
    int TickStart=int(time(0)); // Starting the clock to measure time of computation
    int TickEnd=TickStart;
    int CountBrankrupcies=0;
    double DividendYield=2.276*0.01; // According to http://indexarb.com/dividendYieldSorteddj.html
    double ClusteredBestAgent=0; // Number of  best agents that are within the ClusterLimit population
    double ClusteredWorstAgent=0; // Number of  best agents that are within the ClusterLimit population
    int SharesOutstanding=0; // For stock j=0 only (study of metaorder impact)
    int LastMetaorder=0; // Last time step a metaorder was injected in the OB
    vector<int> MetaorderInjectionTimes; // Times at which metaorders are injected
    //vector<int> MetaorderInjectionAmplitudes; // Effective volume transacted at metaorder injection
    
    //system("CLS");
    cout << "INITIALIZED CLOCK..." << endl;
    
    //gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937); // or gsl_rng_default instead of gsl_rng_mt19937
    gsl_rng * r = make_rng();
    gsl_rng_set(r, static_cast<unsigned long int>(time(0))); // gsl_rng_set(const gsl_rng * r, unsigned long int s)
    vector<double> StockVlm, StockVlmInit, StockIndex, StockLOD, TotalStockQuantities;
    //string TimeScope="NoTimeScope";
    
    // INITIALIZING EACH AGENT
    gsl_matrix * ReflexiveValues = gsl_matrix_calloc (NumberOfStocks, Time); // Computed as the average of all bids
    gsl_matrix * AverageBids = gsl_matrix_calloc (NumberOfStocks, Time); // Computed as the average of all bids
    gsl_matrix * AverageAsks = gsl_matrix_calloc (NumberOfStocks, Time); // Computed as the average of all asks
    gsl_matrix * AverageTopBids = gsl_matrix_calloc (NumberOfStocks, Time); // Computed as the average of all bids in buisness
    gsl_matrix * AverageTopAsks = gsl_matrix_calloc (NumberOfStocks, Time); // Computed as the average of all asks in buisness
    gsl_matrix * MedianTopBids = gsl_matrix_calloc (NumberOfStocks, Time); // Computed as the median of all bids in buisness
    gsl_matrix * MedianTopAsks = gsl_matrix_calloc (NumberOfStocks, Time); // Computed as the median of all asks in buisness
    gsl_matrix * HighestBid = gsl_matrix_calloc (NumberOfStocks, Time);
    gsl_matrix * LowestAsk = gsl_matrix_calloc (NumberOfStocks, Time);
    gsl_matrix * GSLTrueValues = gsl_matrix_calloc (NumberOfStocks, Time);
    gsl_matrix * Spread = gsl_matrix_calloc (NumberOfStocks, Time);
    gsl_matrix * TotalStockQuantityTraded = gsl_matrix_calloc (NumberOfStocks, Time);
    gsl_matrix * VolumesDistributions = gsl_matrix_alloc (2*NumberOfStocks, Time);
    gsl_matrix * ReturnsDistributions = gsl_matrix_alloc (2*NumberOfStocks, Time);
    gsl_matrix * AutoCorrLagPVolumes = gsl_matrix_alloc (2*NumberOfStocks, Time);
    gsl_matrix * CapitalEvolution = gsl_matrix_calloc (NumberOfAgents, Time);
    gsl_matrix * Bankruptcies = gsl_matrix_calloc (2, Time); // Bankruptcies
    gsl_matrix_set(Bankruptcies, 1, 0, 100);
    gsl_matrix * BiasedDelta = gsl_matrix_alloc (NumberOfStocks, Time);
    gsl_matrix * TrueDelta = gsl_matrix_alloc (NumberOfStocks, Time);
    vector<gsl_matrix *> PDs;
    vector<gsl_matrix *> PAverages;
    
    // Generating the Future parameter of each agent according to quasi-hyperbolic discounting
    //vector<double> Gen = STLRandom (NumberOfStocks, "Uniform", "NoPlot");
    int LearningPhase=1000;
    /*
     int Tau=Time;
     if (Time>2*Year) {Tau=8*Month;}
     else if ((Time<=2*Year) && (Time>1*Year)) {Tau=4*Month;}
     else if ((Time<=1*Year) && (Time>6*Month)) {Tau=3*Month;}
     else if ((Time<=6*Month) && (Time>1*Month)) {Tau=2*Month;}
     else if (Time<=1*Month) {Tau=Week;}
     //Tau/=2; // AAA
     vector<int> QH = Quasihyperbolic(NumberOfAgents, Tau, 1, 0.8, 1000*(gsl_rng_uniform (r)+0.0001)); // Beta=1, Delta=0.8 for the Hyperbolic discounting model
     vector<int> Reflex = Quasihyperbolic(NumberOfAgents, 100, 1, 0.75, 1000*(gsl_rng_uniform (r)+0.0001)); // e.g Beta=1, Delta=0.75 yields ~75% of agents with reflexivity below 0.5 // GGG
     //gsl_matrix * QHPlot = gsl_matrix_calloc (1, int(QH.size())); for (int i=0; i<int(QH.size()); i++) {gsl_matrix_set(QHPlot, 0, i, QH[i]);}; PlotGSLMatrix(QHPlot, "QHPlot.xls", 1); exit(0);
     //vector<int> Quasihyperbolic(int NumberOfAgents, int Time, double Beta, double Delta, double Seed) {
     */
    
    vector<Agent> Market;
    int Counter=1;
    for (int j=0; j<NumberOfStocks; j++) {
        StockVlmInit.push_back(0);
    }; // closes j loop
    for (int i=0; i<NumberOfAgents; i++) {
        //outputLog << "Agent " << i << ", ";
        Agent TempAgent;
        vector<Stock> TempStocks;
        TempAgent.AgentName=i;
        TempAgent.RFA=1000*abs((gsl_ran_gaussian (r, 10)));
        TempAgent.RFAFirst=TempAgent.RFA; // RFA at time t=0
        //TempAgent.DiscountFactor = gsl_rng_uniform (r);
        //TempAgent.Future=QH[i]+Week;
        TempAgent.Future=ceil(6*Month*gsl_rng_uniform (r))+Week;
        //if (Time<=1*Month) {TempAgent.Future=int(ceil(Week*gsl_rng_uniform (r)));}; // For very small Time, we go for a basic uniform distribution
        //if ((TempAgent.Future<2) || (TempAgent.Future>Tau)) { // Failsafe
        //outputLog << "Agent " << i << ": Future=" << TempAgent.Future << " (QH[" << i << "]=" << QH[i] << ") replaced to Week" << endl;
        //TempAgent.Future=Week;
        //};
        TempAgent.History=abs(Week+int(gsl_rng_uniform (r)*(Time-TempAgent.Future-1-(2*Week)))); //Week to Time-Future-1
        TempAgent.Reflexivity = gsl_rng_uniform (r);
        
        //TempAgent.Reflexivity*=0.75; // GGG
        //TempAgent.Reflexivity=0.01*Reflex[i]; // GGG
        //TempAgent.TradingWindow=int(0.25*TempAgent.Future + (gsl_rng_uniform(r) * ((1.5*TempAgent.Future)-(0.25*TempAgent.Future))));
        //TempAgent.TradingWindow=int(TempAgent.Future*(0.33 + 0.66*gsl_rng_uniform(r))); // EEE
        //if (TempAgent.TradingWindow<=2) {TempAgent.TradingWindow=3;}; // Lower bound
        //if (TempAgent.TradingWindow>TempAgent.Future) {TempAgent.TradingWindow=TempAgent.Future;}; // Higher bound
        TempAgent.TradingWindow=5+gsl_rng_uniform(r)*TempAgent.Future;
        TempAgent.Bankruptcy = 1; // Bankruptcy = 0 is effective bankruptcy
        //TempAgent.LiquidationFloor = 80 + 10*gsl_rng_uniform (r);
        TempAgent.LiquidationFloor = LiquidationFloor + 10*gsl_rng_uniform (r);
        TempAgent.BiasedValues = gsl_matrix_calloc (NumberOfStocks, Time);
        TempAgent.Versatility=1+floor(10*abs(gsl_rng_uniform (r)));
        //TempAgent.Gesture=floor(100*(4*gsl_rng_uniform (r)/5 + 0.2))/100; // Commercial gesture of the agent in percent of the spread
        TempAgent.Gesture=0.2+0.8*gsl_rng_uniform (r);
        //TempAgent.Gesture*=2;
        if (HPGesture==0) {TempAgent.Gesture*=1;}
        else if (HPGesture==1) {TempAgent.Gesture*=1.5;}
        else if (HPGesture==2) {TempAgent.Gesture*=2;}
        else if (HPGesture==3) {TempAgent.Gesture*=3;};
        //TempAgent.Exploration=int((1000*(gsl_rng_uniform (r))))%4 + 1;
        TempAgent.Exploration=1;
        TempAgent.OffPolicy=1;
        //TempAgent.Epsilon=0.1+0.15*gsl_rng_uniform (r); // 0.3r: The epsilon factor must be small
        TempAgent.Epsilon=0;
        TempAgent.Quant=0;
        TempAgent.Pinch=0;
        TempAgent.FS=27; // Number of states for Forecast()
        TempAgent.FA=27; // Number of actions for Forecast()
        TempAgent.TS=108; // Number of states for Trade()
        TempAgent.TA=9; // Number of actions for Trade()
        TempAgent.FSDim[0]=3; TempAgent.FSDim[1]=3; TempAgent.FSDim[2]=3;
        TempAgent.FADim[0]=3; TempAgent.FADim[1]=3; TempAgent.FADim[2]=3;
        TempAgent.FQDim[0]=3; TempAgent.FQDim[1]=3; TempAgent.FQDim[2]=3; TempAgent.FQDim[3]=3; TempAgent.FQDim[4]=3; TempAgent.FQDim[5]=3;
        TempAgent.FPiDim[0]=3; TempAgent.FPiDim[1]=3; TempAgent.FPiDim[2]=3; TempAgent.FPiDim[3]=3; TempAgent.FPiDim[4]=3; TempAgent.FPiDim[5]=3;
        //TempAgent.FPiDim[0]=4; TempAgent.FPiDim[1]=4; TempAgent.FPiDim[2]=3; TempAgent.FPiDim[3]=3; // DDD2
        TempAgent.TSDim[0]=3; TempAgent.TSDim[1]=3; TempAgent.TSDim[2]=2; TempAgent.TSDim[3]=2; TempAgent.TSDim[4]=3;
        TempAgent.TADim[0]=3; TempAgent.TADim[1]=3;
        TempAgent.TQDim[0]=3; TempAgent.TQDim[1]=3; TempAgent.TQDim[2]=2; TempAgent.TQDim[3]=2; TempAgent.TQDim[4]=3; TempAgent.TQDim[5]=3; TempAgent.TQDim[6]=3;
        TempAgent.TPiDim[0]=3; TempAgent.TPiDim[1]=3; TempAgent.TPiDim[2]=2; TempAgent.TPiDim[3]=2; TempAgent.TPiDim[4]=3; TempAgent.TPiDim[5]=3; TempAgent.TPiDim[6]=3;
        // RL() and LOG of vectors of vectors
        for (int j=0; j<NumberOfStocks; j++) {
            TempAgent.TradingWindowClock.push_back(0);
            vector<double> V; TempAgent.dRoughVec.push_back(V); TempAgent.dSmoothVec.push_back(V); TempAgent.dMuPosVec.push_back(V); TempAgent.dMuNegVec.push_back(V); TempAgent.dSigVec.push_back(V); TempAgent.LiquidPercentVec.push_back(V); TempAgent.dFResultVec.push_back(V); TempAgent.dTResultVec.push_back(V); TempAgent.FPi.push_back(V); TempAgent.TPi.push_back(V); TempAgent.FQ.push_back(V); TempAgent.TQ.push_back(V); TempAgent.ForecastReal.push_back(V); TempAgent.ForecastVirtual.push_back(V); TempAgent.Forecast5.push_back(V); TempAgent.TransactionPriceReal.push_back(V); TempAgent.TransactionPriceVirtual.push_back(V); TempAgent.Qs.push_back(V); // RL
            for (int u=0; u<TempAgent.TS; u++) {TempAgent.Qs[j].push_back(0);}; // Initializatioon of model-based RL framework
            TempAgent.LvolLog.push_back(V); TempAgent.SvolLog.push_back(V); // LOG
            vector<int> U; TempAgent.FSIndex.push_back(U); TempAgent.TSIndex.push_back(U); TempAgent.FAIndexReal.push_back(U); TempAgent.FAIndexVirtual.push_back(U); TempAgent.TAIndexReal.push_back(U); TempAgent.TAIndexVirtual.push_back(U); TempAgent.QuantitiesReal.push_back(U); TempAgent.QuantitiesVirtual.push_back(U); TempAgent.FNumberA.push_back(U); TempAgent.TNumberA.push_back(U); // RL
            TempAgent.RoughLog.push_back(U); TempAgent.SmoothLog.push_back(U); TempAgent.ReflexiveLog.push_back(U); TempAgent.ToolRealLog.push_back(U); TempAgent.LagRealLog.push_back(U); TempAgent.WeightRealLog.push_back(U); TempAgent.ToolVirtualLog.push_back(U); TempAgent.LagVirtualLog.push_back(U); TempAgent.WeightVirtualLog.push_back(U); TempAgent.MuLog.push_back(U); TempAgent.SigLog.push_back(U); TempAgent.RFALog.push_back(U); TempAgent.RBALog.push_back(U); TempAgent.LiquidityLog.push_back(U); TempAgent.PinchRealLog.push_back(U); TempAgent.QuantRealLog.push_back(U); TempAgent.PinchVirtualLog.push_back(U); TempAgent.QuantVirtualLog.push_back(U); TempAgent.FResultDisReal.push_back(U); TempAgent.FResultDisVirtual.push_back(U); TempAgent.TResultDisReal.push_back(U); TempAgent.TResultDisVirtual.push_back(U); // LOG
        }; // closes j loop
        
        
        // Attribution of biases (for all: if NEBp>r then NEB is activated)
        TempAgent.Human=100*gsl_rng_uniform (r); // If it is smaller than pBots then it is a bot
        TempAgent.NEBLossAversion=0.5*gsl_rng_uniform (r);
        TempAgent.NEBPositivity=0.5*gsl_rng_uniform (r);
        TempAgent.NEBNegativity=0.5*gsl_rng_uniform (r);
        // DelayDiscounting is set below
        TempAgent.NEBLearningRate=0.05+0.2*gsl_rng_uniform (r);
        TempAgent.NEBFear=0.2*gsl_rng_uniform (r);
        TempAgent.NEBGreed=0.2*gsl_rng_uniform (r);
        
        
        // 50% have nothing, 10% PROFILE1 (NEB124), 10% PROFILE2 (NEB356), 5% PROFILE3 (NEB89), 5% PROFILE12 (NEB123456), 5% PROFILE13 (NEB12489), 5% PROFILE23 (NEB35689), 10% everything
        //pNEBLossAversion24 => pNEBLossAversion, pNEBDelayDiscounting56=>pNEBPositivity, pNEB89=>pNEBDelayDiscounting
        if (TempAgent.Human<pBots) {
            TempAgent.NEB=0; // Bot
            TempAgent.NEBLossAversion=-1;
            TempAgent.NEBPositivity=-1;
            TempAgent.NEBNegativity=-1;
            TempAgent.NEBFear=-1;
            TempAgent.NEBGreed=-1;
        }
        else if ((TempAgent.Human>=pBots) && (TempAgent.Human<pBots + pNEBLossAversion)) {
            TempAgent.NEB=1; // Loss aversion
            TempAgent.NEBPositivity=-1;
            TempAgent.NEBNegativity=-1;
            TempAgent.NEBFear=-1;
            TempAgent.NEBGreed=-1;
        }
        else if ((TempAgent.Human>=pBots + pNEBLossAversion) && (TempAgent.Human<pBots + pNEBLossAversion + pNEBPositivity)) {
            TempAgent.NEB=2; // Positivity
            TempAgent.NEBLossAversion=-1;
            TempAgent.NEBNegativity=-1;
            TempAgent.NEBFear=-1;
            TempAgent.NEBGreed=-1;
        }
        else if ((TempAgent.Human>=pBots + pNEBLossAversion + pNEBPositivity) && (TempAgent.Human<pBots + pNEBLossAversion + pNEBPositivity + pNEBNegativity)) {
            TempAgent.NEB=3; // Negativity
            TempAgent.NEBLossAversion=-1;
            TempAgent.NEBPositivity=-1;
            TempAgent.NEBFear=-1;
            TempAgent.NEBGreed=-1;
        }
        else if ((TempAgent.Human>=pBots + pNEBLossAversion + pNEBPositivity + pNEBNegativity) && (TempAgent.Human<pBots + pNEBLossAversion + pNEBPositivity + pNEBNegativity + pNEBDelayDiscounting)) {
            TempAgent.NEB=4; // Delay discounting
            TempAgent.NEBLossAversion=-1;
            TempAgent.NEBPositivity=-1;
            TempAgent.NEBNegativity=-1;
            TempAgent.NEBFear=-1;
            TempAgent.NEBGreed=-1;
            TempAgent.Future=ceil(Week*gsl_rng_uniform (r))+Week;
            TempAgent.History=abs(Week+int(gsl_rng_uniform (r)*(Time-TempAgent.Future-1-(2*Week)))); //2*Week to Time-Future-1
            TempAgent.TradingWindow=5+gsl_rng_uniform(r)*TempAgent.Future;
        }
        else if ((TempAgent.Human>=pBots + pNEBLossAversion + pNEBPositivity + pNEBNegativity + pNEBDelayDiscounting) && (TempAgent.Human<pBots + pNEBLossAversion + pNEBPositivity + pNEBNegativity + pNEBDelayDiscounting + pNEBFear)) {
            TempAgent.NEB=5; // Fear
            TempAgent.NEBLossAversion=-1;
            TempAgent.NEBPositivity=-1;
            TempAgent.NEBGreed=-1;
        }
        else if ((TempAgent.Human>=pBots + pNEBLossAversion + pNEBPositivity + pNEBNegativity + pNEBDelayDiscounting + pNEBFear) && (TempAgent.Human<pBots + pNEBLossAversion + pNEBPositivity + pNEBNegativity + pNEBDelayDiscounting + pNEBFear + pNEBGreed)) {
            TempAgent.NEB=6; // Greed
            TempAgent.NEBLossAversion=-1;
            TempAgent.NEBPositivity=-1;
            TempAgent.NEBFear=-1;
        }
        else if ((TempAgent.Human>=pBots + pNEBLossAversion + pNEBPositivity + pNEBNegativity + pNEBDelayDiscounting + pNEBFear + pNEBGreed) && (TempAgent.Human<pBots + pNEBLossAversion + pNEBPositivity + pNEBNegativity + pNEBDelayDiscounting + pNEBFear + pNEBGreed + pNEBLearningRate)) {
            TempAgent.NEB=7; // Learning
            TempAgent.NEBLossAversion=-1;
            TempAgent.NEBPositivity=-1;
            TempAgent.NEBFear=-1;
            TempAgent.NEBGreed=-1;
            TempAgent.NEBLearningRate*=LearningRateScale;
        }
        else if ((TempAgent.Human>=pBots + pNEBLossAversion + pNEBPositivity + pNEBNegativity + pNEBDelayDiscounting + pNEBFear2) && (TempAgent.Human<pBots + pNEBLossAversion + pNEBPositivity + pNEBNegativity + pNEBDelayDiscounting + pNEBFear2 + pNEBGreed)) {
            TempAgent.NEB=8; // Fear2
            TempAgent.NEBLossAversion=-1;
            TempAgent.NEBPositivity=-1;
            TempAgent.NEBFear=-1;
            TempAgent.NEBGreed=-1;
        };
        
        for (int j=0; j<NumberOfStocks; j++) {
            TempAgent.Leash.push_back(0.001*gsl_rng_uniform (r));
            TempAgent.LeashVol.push_back(5*0.01*gsl_rng_uniform (r)); // TTT4
            //TempAgent.Accuracy.push_back(7+int((1000*(gsl_rng_uniform (r))))%5); // Prob. of bubble burst, on average 0.2,0.75,1,2 per year for Accuracy=13,11,10,9 resp.
            //TempAgent.Accuracy.push_back(10+gsl_rng_uniform (r)); // TTT3
            TempAgent.Accuracy.push_back(HPAccuracy+gsl_rng_uniform (r)); // TTT3
            Stock TempStock;
            TempStock.StockName = j;
            TempStock.StockQuantity = abs(int(gsl_ran_gaussian (r, 100)));
            TempStock.StockQuantityFirst = TempStock.StockQuantity;
            TempStocks.push_back(TempStock);
            Counter += 1;
        }; // closes j loop
        TempAgent.NAVOnJanFirst=TempAgent.Capital(); // NAV on Jan 1st
        TempAgent.Stocks = TempStocks;
        //TempStocks.erase(TempStocks.begin(),TempStocks.end());
        Market.push_back(TempAgent);
        SharesOutstanding+=TempAgent.Stocks[0].StockQuantity;
    }; // closes i loop
    //outputLog << "INITIALIZED AGENTS..." << endl;
    cout << "INITIALIZED AGENTS..." << endl;
    
    
    // GENERATING THE STOCKS TRUE VALUES
    vector<double> Gen2 = STLRandom (NumberOfStocks, "Uniform", "NoPlot");
    //double HPTrueMu=0.1;
    for (int j=0; j<NumberOfStocks; j++) {
        //vector<double> S = PoissonRandomWalk (100, Week+3*Month*Gen2[j], int(0.1*Time/250+0.9*Gen2[j]*Time/250), Time, 1000*Gen2[j]);
        vector<double> S = PoissonRandomWalk (100, Week+HPTrueMu*3*Month*Gen2[j], int(0.1*Time/250+0.9*Gen2[j]*Time/250), Time, 1000*Gen2[j]);
        for (int t=0; t<Time; t++) {gsl_matrix_set (GSLTrueValues, j, t, S[t]);};
    };
    vector<vector<double> > TrueValues = GSLMatrixToSTLMatrix(GSLTrueValues);
    cout << "INITIALIZED STOCKS TRUE VALUES..." << endl;
    
    
    // INITIALIZING THE AGENTS BIASED VALUES
    for (int i=0; i<NumberOfAgents; i++) {
        for (int j=0; j<NumberOfStocks; j++) {
            vector<double> C = CointegratedWalk(TrueValues[j], Market[i].Leash[j], Market[i].LeashVol[j], Market[i].Accuracy[j]); // Master, Leash, LeashVolatility, Accuracy
            // vector<double> C = CointegratedWalk(TrueValues[j], 0.0005*2*gsl_rng_uniform (r), 0.005*2*gsl_rng_uniform (r), Accuracy); // Master, Leash, LeashVolatility, Accuracy
            //if ((j==0) || (j== NumberOfStocks-1)) {PlotSTL(C, "000");};
            for (int t=0; t<Time; t++) {
                gsl_matrix_set(Market[i].BiasedValues,j,t, C[t]);
            }; // closes t loop
            Market[i].Stocks[j].StockBidValue = gsl_matrix_get(Market[i].BiasedValues,j,0) + abs(gsl_ran_gaussian (r, 1));
            //Market[i].Stocks[j].StockAskValue = (Market[i].Stocks[j].StockBidValue)*(1 + Market[i].ProfitMargin + Market[i].RiskPremium);
            Market[i].Stocks[j].StockAskValue = (Market[i].Stocks[j].StockBidValue)*(1 + 0.02*gsl_rng_uniform (r));
            C.clear();
        }; // closes j loop
        Market[i].RBAFirst=Market[i].StockHoldings(); // RBA at time t=0
    }; // clloses i loop
    cout << "INITIALIZED AGENTS BIASED VALUES..." << endl;
    
    
    
    // Checking on distribution of agents Future parameter
    //gsl_matrix * M=gsl_matrix_alloc (1, NumberOfAgents);
    //for (int i=0; i<NumberOfAgents; i++) {gsl_matrix_set (M, 0, i, Market[i].Future);};
    //gsl_matrix * D = gsl_matrix_alloc (2, NumberOfAgents);
    //D = GSLDistribution(M, 10);
    //PlotGSLMatrix(M, "MM.xls", 1);
    //PlotGSLMatrix(D, "DD.xls", 1);
    
    
    
    // INITILAZING THE MARKET SIMULATION
    // Initializing reflexive values as averages at t=0, and reference values = f(biasedtrue,reflexive)
    vector<double> VSpread;
    for (int j=0; j<NumberOfStocks; j++) {
        double InitValBid=0;
        for (int i=0; i<NumberOfAgents; i++) {
            InitValBid+=Market[i].Stocks[j].StockBidValue;
        }; // closes i loop
        gsl_matrix_set(ReflexiveValues, j, 0, InitValBid/NumberOfAgents);
        gsl_matrix_set(AverageBids, j, 0, InitValBid/NumberOfAgents);
        gsl_matrix_set(AverageAsks, j, 0, InitValBid/NumberOfAgents);
        gsl_matrix_set(AverageTopBids, j, 0, InitValBid/NumberOfAgents);
        gsl_matrix_set(AverageTopAsks, j, 0, InitValBid/NumberOfAgents);
        gsl_matrix_set(MedianTopBids, j, 0, InitValBid/NumberOfAgents);
        gsl_matrix_set(MedianTopAsks, j, 0, InitValBid/NumberOfAgents);
        gsl_matrix_set(HighestBid, j, 0, InitValBid/NumberOfAgents);
        gsl_matrix_set(LowestAsk, j, 0, InitValBid/NumberOfAgents);
        VSpread.push_back(0);
    }; // closes j loop
    
    // Initializing TotalStockQuantities, the vector of total quantity of each stock
    for (int j=0; j<NumberOfStocks; j++) {
        double TempStockQuantity=0;
        for (int i=0; i<NumberOfAgents; i++) {
            TempStockQuantity+=Market[i].Stocks[j].StockQuantity;
        }; // closes i loop
        TotalStockQuantities.push_back(TempStockQuantity);
    }; // closes j loop
    // Initializing Market Indicators
    for (int i=0; i<NumberOfAgents; i++) {gsl_matrix_set(CapitalEvolution, i, 0, 100);};
    cout << "INITIALIZED MARKET INDICATORS..." << endl;
    
    // INITIALIZATION OF MARKETEVOLUTION
    vector<vector<Agent> > MarketEvolution;
    MarketEvolution.push_back(Market);
    cout << "INITIALIZED MARKET EVOLUTION..." << endl;
    
    
    // INITIALIZATION OF THE ORDER BOOK
    vector<vector<StockOrderBook> > FullOrderBook; // This is our final order book with many pages t
    vector<StockOrderBook> TempOrderBook; // This is the order book at time t, composed of many order books (one for each stock)
    for (int j=0; j<NumberOfStocks; j++) {
        StockOrderBook TempStockOrderBook;
        vector<int> J = Shuffle(NumberOfAgents);
        for (int i=0; i<NumberOfAgents ; i++) {
            Order TempOrderBid;
            TempOrderBid.Nature=0; // 0 for Bid, 1 for Ask
            TempOrderBid.Agent=J[i];
            TempOrderBid.Stock=j;
            TempOrderBid.Q=MarketEvolution[0][J[i]].Stocks[j].StockQuantity; // CHECK!!! For now the buyer buys all the seller has !!!
            TempOrderBid.Price=MarketEvolution[0][J[i]].Stocks[j].StockBidValue;
            Order TempOrderAsk;
            TempOrderAsk.Nature=1; // 0 for Bid, 1 for Ask
            TempOrderAsk.Agent=J[i];
            TempOrderAsk.Stock=j;
            TempOrderAsk.Q=MarketEvolution[0][J[i]].Stocks[j].StockQuantity; // CHECK!!! For now the seller sells all he has !!!
            TempOrderAsk.Price=MarketEvolution[0][J[i]].Stocks[j].StockAskValue;
            // Now filling the final order book of stock j
            TempStockOrderBook.Bids.push_back(TempOrderBid);
            TempStockOrderBook.Asks.push_back(TempOrderAsk);
        }; // closes i loop
        TempOrderBook.push_back(TempStockOrderBook);
        TempStockOrderBook.Bids.clear();
        TempStockOrderBook.Asks.clear();
    }; // closes j loop
    FullOrderBook.push_back(TempOrderBook); // FullOrderBook is filled on first page at t=0
    cout << "INITIALIZED ORDER BOOK..." << endl;
    
    
    // INITIALIZATION t=0 FINISHES HERE
    int TickMid=(int(time(0))-TickStart); // Checking the clock for time of computation of initialization
    //ofstream outputT2("/Users/admin/Documents/GNT/SYMBA/Details.txt", ofstream::app);
    system("CLS");
    cout << "****************************" << endl;
    cout << "** STOCK MARKET SIMULATOR **" << endl;
    cout << "****************************" << endl;
    cout << "I=" << NumberOfAgents << ", J=" << NumberOfStocks << ", T=" << Time << endl << endl;
    cout.setf(ios::fixed);
    cout << "Initialization : " << TickMid << "s" << endl;
    if (Plot=="On") {outputLog << "INITIALIZATION t=0 FINISHED..." << endl << endl;};
    
    
    
    
    // ***** MARKET SIMULATION STARTS HERE WITH t-LOOP ***** (SS2)
    if (Plot=="On") {outputLog << "***** MARKET SIMULATION WITH t-LOOP STARTS HERE *****" << endl;};
    for (int t=1; t<Time; t++) {
        for (int j=0; j<NumberOfStocks; j++) {
            // REFLEXIVE VALUES = CLEARING PRICE OF ORDER BOOK (SS3)
            gsl_matrix_set(ReflexiveValues, j, t, gsl_matrix_get(ReflexiveValues, j, t-1));
            if (Plot=="On") {outputLog << "t=" << t << ", j=" << j << ": ReflexiveValues computed..." << endl;};
        }; // closes j loop
        // outputDebug << "t=" << t << ": ReflexiveValues computed..." << endl;
        
        
        // INTEREST RATES (SS5)
        for (int i=0; i<NumberOfAgents; i++) {
            Market[i].RFA*=1+(Rate/Year); // Saving account banking rate, modeled w/o inflation
            Market[i].RFA+=Market[i].StockHoldings()*DividendYield/Year; // Modeled as compounded dividend yield
        }; // closes i loop
        if (Plot=="On") {outputLog << "t=" << t << ", Interest rates computed..." << endl << endl;};
        // outputDebug << "t=" << t << ": Interest rates computed..." << endl;
        
        // REINITIALIZATION OF THE AGENTS' PORTFOLIO AFTER LEARNING PHASE
        for (int i=0; i<NumberOfAgents ; i++) {
            if (t==LearningPhase) {
                Market[i].Bankruptcy=1;
                Market[i].RFA=Market[i].RFAFirst;
                for (int j=0; j<NumberOfStocks; j++) {
                    Market[i].Stocks[j].StockQuantity=Market[i].Stocks[j].StockQuantityFirst;
                };
            };
        };
        
        // RECORD OF LEADER ACTIONS FOR CLUSTER TRADING
        int LeaderAgent=0; // "Random" agent is the leader for cluster trading
        double MinCap=9999999; double MaxCap=-9999999; int WorstAgent=0; int BestAgent=0;
        for (int i=0; i<NumberOfAgents; i++) {
            if (Market[i].Capital()<MinCap) {MinCap=Market[i].Capital(); WorstAgent=i;};
            if (Market[i].Capital()>MaxCap) {MaxCap=Market[i].Capital(); BestAgent=i;};
        };
        if (LeaderType=="Worst") {LeaderAgent=WorstAgent;};
        if (LeaderType=="Best") {LeaderAgent=BestAgent;};
        if (LeaderType=="Static") {LeaderAgent=0;};
        if ((BestAgent<ClusterLimit) && (t>=LearningPhase)) {ClusteredBestAgent+=1;}; // Nb of LeaderAgent that are part of the ClusterLimit
        if ((WorstAgent<ClusterLimit) && (t>=LearningPhase)) {ClusteredWorstAgent+=1;}; // Nb of LeaderAgent that are part of the ClusterLimit
        
        if ((Market[LeaderAgent].Bankruptcy==0) && (ClusterLimit>1)) { // Agent leader cannot be bankrupt
            Market[LeaderAgent].RFA=Market[LeaderAgent].RFAFirst;
            Market[LeaderAgent].Bankruptcy=1;
            for (int j=0; j<NumberOfStocks; j++) {Market[LeaderAgent].Stocks[j].StockQuantity=Market[LeaderAgent].Stocks[j].StockQuantityFirst;};
        };
        int LeaderQuant=Market[LeaderAgent].Quant;
        int LeaderPinch=Market[LeaderAgent].Pinch;
        if (t<=2) {LeaderQuant=0; LeaderPinch=0;};
        
        
        
        
        // REINFORCEMENT LEARNING RL() (SS6)
        gsl_matrix_set(Bankruptcies, 1, t, gsl_matrix_get(Bankruptcies, 1, t-1));
        for (int j=0; j<NumberOfStocks; j++) {
            for (int i=0; i<NumberOfAgents; i++) {
                string NewPlot=Plot;
                if (LeaderType=="Noise") {
                    LeaderQuant=int((1000*(gsl_rng_uniform (r))))%3; //0, 1, 2
                    LeaderPinch=int((1000*(gsl_rng_uniform (r))))%3; //0, 1, 2
                };
                if (Market[i].Bankruptcy!=0) {Market[i].RL(j, t, Rate, ReflexiveValues, VSpread[j], gsl_matrix_get(TotalStockQuantityTraded, j, t-1), Time, NumberOfStocks, TradingFrequencyCond, NewPlot, VersatilityCondition, gsl_matrix_get(Bankruptcies, 1, t), t-(t-LearningPhase)%Year, LearningPhase, LeaderAgent, LeaderQuant, LeaderPinch, ClusterLimit, Trunk);}; // XYZ1
                
            }; // closes i loop
        }; // closes j loop
        if (Plot=="On") {outputLog << "t=" << t << ", RL() method ran in all agents..." << endl << endl;};
        // outputDebug << "t=" << t << ": RL() computed..." << endl;
        
        
        // FILLING THE ORDER BOOK (SS8)
        TempOrderBook.clear();
        for (int j=0; j<NumberOfStocks; j++) {
            StockOrderBook TempStockOrderBook;
            vector<int> J = Shuffle(NumberOfAgents);
            for (int i=0; i<NumberOfAgents ; i++) {
                bool B1 = Market[(J[i])].Bankruptcy != 1; // Agent is bankrupt
                bool B2 = Market[(J[i])].QuantitiesReal[j][int(Market[(J[i])].QuantitiesReal[j].size())-1]==0; // Agent is not trading at this time step
                if (Plot=="On") {
                    outputLog << "Filling TempOrderBook of Agent " << J[i] << ", Stock " << j << ", at t=" << t << ": ";
                    outputLog << "(Bankruptcy=" << B1 << ", Waiting=" << B2 << "): ";
                };
                if (B1 || B2) {if (Plot=="On") {outputLog << "order cancelled" << endl;}; continue;}
                else if ((Market[(J[i])].QuantitiesReal[j][int(Market[(J[i])].QuantitiesReal[j].size()) - 1]>0) && !B1 && !B2) {
                    Order TempOrderBid;
                    TempOrderBid.Nature=0; // 0 for Bid, 1 for Ask
                    TempOrderBid.Agent=Market[(J[i])].AgentName;
                    TempOrderBid.Stock=j;
                    TempOrderBid.Price=Market[(J[i])].Stocks[j].StockBidValue;
                    //if (Trunk>-1) {TempOrderBid.Price=DigitTrunk (TempOrderBid.Price, Trunk, "Ceil");}; // Trunk according to the number of significant digits
                    TempOrderBid.Q=Market[(J[i])].QuantitiesReal[j][int(Market[(J[i])].QuantitiesReal[j].size()) - 1];
                    TempStockOrderBook.Bids.push_back(TempOrderBid);
                    if (Plot=="On") {outputLog << "Bid Q=" << TempOrderBid.Q << " at $" << TempOrderBid.Price << " (RFA=$" << Market[(J[i])].RFA << ")" << endl;};
                } // closes else if
                else if ((Market[(J[i])].QuantitiesReal[j][int(Market[(J[i])].QuantitiesReal[j].size()) - 1]<0) && !B1 && !B2) {
                    Order TempOrderAsk;
                    TempOrderAsk.Nature=1; // 0 for Bid, 1 for Ask
                    TempOrderAsk.Agent=Market[(J[i])].AgentName;
                    TempOrderAsk.Stock=j;
                    TempOrderAsk.Price=Market[(J[i])].Stocks[j].StockAskValue;
                    //if (Trunk>-1) {TempOrderAsk.Price=DigitTrunk (TempOrderAsk.Price, Trunk, "Ceil");}; // Trunk according to the number of significant digits
                    TempOrderAsk.Q=abs(Market[(J[i])].QuantitiesReal[j][int(Market[(J[i])].QuantitiesReal[j].size()) - 1]);
                    TempStockOrderBook.Asks.push_back(TempOrderAsk);
                    if (Plot=="On") {outputLog << "Ask Q=" << TempOrderAsk.Q << " at $" << TempOrderAsk.Price << " (QuantitiesReal[j]=" << Market[(J[i])].QuantitiesReal[j][int(Market[(J[i])].QuantitiesReal[j].size()) - 1] << ")" << endl;};
                }; // closes else if
            }; // closes i loop
            TempStockOrderBook.Sort(); // (SS9)
            if (Plot=="On") {outputLog << "t=" << t << ", Sorted OB[j=" << j << "]" << endl << endl;};
            
            TempOrderBook.push_back(TempStockOrderBook);
            TempStockOrderBook.Bids.clear();
            TempStockOrderBook.Asks.clear();
        }; // closes j loop
        FullOrderBook.push_back(TempOrderBook); // FullOrderBook is filled at page t
        if (Plot=="On") {outputLog << "t=" << t << ", Filled the OB..." << endl;};
        // outputDebug << "t=" << t << ": OB filled..." << endl;
        
        
        // RECORDING AVERAGE BIDS AND ASKS
        // The market price is taken as the lowest of clearing prices (not the median), and the bid and ask are taken as the average (not median) of all bid and ask prices of orders that will be transacted: hence the bid and ask being sometimes not below and above of market price
        for (int j=0; j<NumberOfStocks; j++) {
            double MeanBid=0; double MeanAsk=0; double MeanTopBid=0; double MeanTopAsk=0; int Level=-1;
            int OBSize=int(min(FullOrderBook[t][j].Asks.size(), FullOrderBook[t][j].Bids.size()));
            if (OBSize==0) {
                gsl_matrix_set(AverageTopBids, j, t, gsl_matrix_get (AverageTopBids, j, t-1));
                gsl_matrix_set(AverageTopAsks, j, t, gsl_matrix_get (AverageTopAsks, j, t-1));
                gsl_matrix_set(MedianTopBids, j, t, gsl_matrix_get (MedianTopBids, j, t-1));
                gsl_matrix_set(MedianTopAsks, j, t, gsl_matrix_get (MedianTopAsks, j, t-1));
                gsl_matrix_set(AverageBids, j, t, gsl_matrix_get (AverageBids, j, t-1));
                gsl_matrix_set(AverageAsks, j, t, gsl_matrix_get (AverageAsks, j, t-1));
                gsl_matrix_set(HighestBid, j, t, gsl_matrix_get (HighestBid, j, t-1));
                gsl_matrix_set(LowestAsk, j, t, gsl_matrix_get (LowestAsk, j, t-1));
            }
            else {
                for (int k=0; k<OBSize; k++) {if (FullOrderBook[t][j].Bids[k].Price>=FullOrderBook[t][j].Asks[k].Price) {Level+=1;}}; // Level
                gsl_matrix_set(AverageTopBids, j, t, gsl_matrix_get (AverageTopBids, j, t-1));
                gsl_matrix_set(AverageTopAsks, j, t, gsl_matrix_get (AverageTopAsks, j, t-1));
                gsl_matrix_set(HighestBid, j, t, gsl_matrix_get (HighestBid, j, t-1));
                gsl_matrix_set(LowestAsk, j, t, gsl_matrix_get (LowestAsk, j, t-1));
                if (Level>-1) {
                    for (int k=0; k<Level+1; k++) {
                        MeanTopBid+=FullOrderBook[t][j].Bids[k].Price/(Level+1);
                        MeanTopAsk+=FullOrderBook[t][j].Asks[k].Price/(Level+1);
                    }; // closes k loop
                    gsl_matrix_set(AverageTopBids, j, t, MeanTopBid);
                    gsl_matrix_set(AverageTopAsks, j, t, MeanTopAsk);
                    gsl_matrix_set(MedianTopBids, j, t, FullOrderBook[t][j].Bids[(Level+1)/2].Price);
                    gsl_matrix_set(MedianTopAsks, j, t, FullOrderBook[t][j].Asks[(Level+1)/2].Price);
                    gsl_matrix_set(HighestBid, j, t, FullOrderBook[t][j].Bids[0].Price);
                    gsl_matrix_set(LowestAsk, j, t, FullOrderBook[t][j].Asks[0].Price);
                }; // closes if
                for (int k=0; k<OBSize; k++) {
                    MeanBid+=FullOrderBook[t][j].Bids[k].Price/OBSize;
                    MeanAsk+=FullOrderBook[t][j].Asks[k].Price/OBSize;
                }; // closes k loop
                gsl_matrix_set(AverageBids, j, t, MeanBid);
                gsl_matrix_set(AverageAsks, j, t, MeanAsk);
            }; // closes else
            VSpread.at(j)=abs(gsl_matrix_get(AverageBids, j, t) - gsl_matrix_get(AverageAsks, j, t));
            /*
             if (HPSpread==0) {VSpread.at(j)=abs(gsl_matrix_get(AverageBids, j, t) - gsl_matrix_get(AverageAsks, j, t));} // whole average definition of spread XXX3 // Small volatility
             else if (HPSpread==1) {VSpread.at(j)=abs(gsl_matrix_get(AverageTopBids, j, t) - gsl_matrix_get(AverageTopAsks, j, t));} // average definition of spread XXX1 // Too smooth
             else if (HPSpread==2) {VSpread.at(j)=abs(gsl_matrix_get(MedianTopBids, j, t) - gsl_matrix_get(MedianTopAsks, j, t));} // median definition of spread XXX2 // Too smooth
             else {VSpread.at(j)=gsl_matrix_get(HighestBid, j, t) - gsl_matrix_get(LowestAsk, j, t); if (VSpread[j]<0) {VSpread.at(j)=0;}}; // Formal definition (difference between the highest bid and lowest ask) XXX4 // Too volatile
             */
            double FormalDefinitionSpread=gsl_matrix_get(HighestBid, j, t) - gsl_matrix_get(LowestAsk, j, t); if (FormalDefinitionSpread<0) {FormalDefinitionSpread=0;};
            if (Plot=="On") {
                outputLog << "gsl_matrix_get(AverageTopBids, j, t)=" << gsl_matrix_get(AverageTopBids, j, t) << ", gsl_matrix_get(AverageTopAsks, j, t)=" << gsl_matrix_get(AverageTopAsks, j, t) << ", Type-I VSpread[j]=" << abs(gsl_matrix_get(AverageTopBids, j, t) - gsl_matrix_get(AverageTopAsks, j, t)) << endl;
                outputLog << "gsl_matrix_get(MedianTopBids, j, t)=" << gsl_matrix_get(MedianTopBids, j, t) << ", gsl_matrix_get(MedianTopAsks, j, t)=" << gsl_matrix_get(MedianTopAsks, j, t) << ", Type-II VSpread[j]=" << abs(gsl_matrix_get(MedianTopBids, j, t) - gsl_matrix_get(MedianTopAsks, j, t)) << endl;
                outputLog << "gsl_matrix_get(AverageBids, j, t)=" << gsl_matrix_get(AverageBids, j, t) << ", gsl_matrix_get(AverageAsks, j, t)=" << gsl_matrix_get(AverageAsks, j, t) << ", Type-III VSpread[j]=" << abs(gsl_matrix_get(AverageBids, j, t) - gsl_matrix_get(AverageAsks, j, t)) << endl;
                outputLog << "gsl_matrix_get(HighestBid, j, t)=" << gsl_matrix_get(HighestBid, j, t) << ", gsl_matrix_get(LowestAsk, j, t)=" << gsl_matrix_get(LowestAsk, j, t) << ", Type-IV VSpread[j]=" << FormalDefinitionSpread << endl;
                outputLog << "VSpread[j]=" << VSpread[j] << endl;
            };
        }; // closes j loop
        // outputDebug << "t=" << t << ": Average bids and asks computed..." << endl;
        
        
        // LOG OF CREDIT EVENTS
        // Bankruptcy = 0 is effective bankruptcy
        int LearningPhase2=10;
        double NAVOfNonBankruptAgents=0; int NbOfNonBankruptAgents=0; double MarketPerformanceJan1st=0;
        if ((t>=LearningPhase2) && ((t-LearningPhase2)%Year2==0)) {
            for (int i=0; i<NumberOfAgents ; i++) {
                Market[i].NAVOnJanFirst=Market[i].Capital(); // Recording all Jan 1st NAVs
            };
            MarketPerformanceJan1st=gsl_matrix_get(ReflexiveValues, 0, t);
        };
        if ((t>LearningPhase2) && ((t-LearningPhase2)%Year2>Week)) {
            for (int i=0; i<NumberOfAgents ; i++) {
                if (Market[i].Bankruptcy!=0) {NAVOfNonBankruptAgents+=100*Market[i].Capital()/Market[i].NAVOnJanFirst; NbOfNonBankruptAgents+=1;};
            };
            NAVOfNonBankruptAgents/=NbOfNonBankruptAgents; // Mean NAV ratio (%) since Jan 1st of all non-bankrupt agents
            if (NbOfNonBankruptAgents==0) {NAVOfNonBankruptAgents=100;}; // Failsafe
            gsl_matrix_set(Bankruptcies, 1, t, NAVOfNonBankruptAgents);
            //outputLog << "At t=" << t << ", NAVOfNonBankruptAgents=" << NAVOfNonBankruptAgents << ", NbOfNonBankruptAgents=" << NbOfNonBankruptAgents << ", Pt(Jan1st)=" << MarketPerformanceJan1st << ", Pt=" << gsl_matrix_get(ReflexiveValues, 0, t) << endl;
        };
        
        int Discount=max(0, 100-int(NAVOfNonBankruptAgents));
        for (int i=0; i<NumberOfAgents ; i++) {
            if ((Market[i].Bankruptcy!=0) && (t>LearningPhase2) && ((t-LearningPhase2)%Year2>Week) && (100*Market[i].Capital()/Market[i].NAVOnJanFirst+Discount<=Market[i].LiquidationFloor)) {
                CountBrankrupcies+=1;
                //outputLog << "At t=" << t << ", NAV of agent " << i << " reached " << 100*Market[i].Capital()/Market[i].NAVOnJanFirst << "% (while average agent performance was " << NAVOfNonBankruptAgents << "% since Jan 1st at t=" << t-(t-LearningPhase2)%Year2 << "), and got bankrupt (from NAV=$" << Market[i].Capital() << ") due to agent LiquidationFloor=" << Market[i].LiquidationFloor << "%" << endl;
                Market[i].Liquidation(i, Market);
            };
        };
        gsl_matrix_set(Bankruptcies, 0, t, 100.0*CountBrankrupcies/NumberOfAgents);
        
        
        
        if (Plot=="On") {outputLog << "t=" << t << ", Log of credit events completed..." << endl;};
        // outputDebug << "t=" << t << ": Credit events computed..." << endl;
        
        
        
        
        // METAORDERS
        // Counting the OB business levels
        //if (t>1000) {Plot="On";};
        int OBLevelSize=0;
        int CountOB=0; int OBLevel=0;
        for (int k=0; k<int(min(FullOrderBook[t][0].Bids.size(), FullOrderBook[t][0].Asks.size())); k++) {
            if ((FullOrderBook[t][0].Bids[k].Price < FullOrderBook[t][0].Asks[k].Price) && (OBLevel==0)) {OBLevel=1; OBLevelSize=k-1;};
            CountOB=1;
        }; // closes k loop
        //ofstream outputOBLevels(Machine+"OBlevels.log", ofstream::app); outputOBLevels << " At t=" << t << " OB size =" << OBLevelSize << endl;;
        // Meta order injection
        if ((t>500) && (t-LastMetaorder>=6*Month) && (MetaorderImpact>0)) {
            int Res=FullOrderBook[t][0].MetaorderInjection(Market, SharesOutstanding, MetaorderImpact, OBLevelSize);
            if (Res>-1) {
                LastMetaorder=t;
                outputMetalog << "t=" << t << ", metaorder successfully injected at OB level " << Res << ", last metaorder was at t=" << LastMetaorder << endl;
                if ((t-LearningPhase>1000) && (Time-t>7*Month)) {
                    MetaorderInjectionTimes.push_back(LastMetaorder-LearningPhase);
                    outputMetalog << "And metaorder recorded in MetaorderInjectionTimes[]..." << endl;
                };
            };
        };
        
        
        // OUTPUTING THE ORDER BOOK AT TIME t (SS10)
        if (Plot=="On") {
            for (int j=0; j<NumberOfStocks; j++) {
                int CountOB=0;
                int OBLevel=0;
                outputLog << endl << "ORDER BOOK OF STOCK " << j << " AT TIME t=" << t << endl;
                if (int(min(FullOrderBook[t][j].Bids.size(), FullOrderBook[t][j].Asks.size()))==0) {outputLog << "t=" << t << ", the Order Book has zero size..." << endl;};
                for (int k=0; k<int(min(FullOrderBook[t][j].Bids.size(), FullOrderBook[t][j].Asks.size())); k++) { // ABCD
                    if ((FullOrderBook[t][j].Bids[k].Price < FullOrderBook[t][j].Asks[k].Price) && (OBLevel==0)) {
                        outputLog << "============================================================================================" << endl;
                        OBLevel=1;
                    };
                    CountOB=1;
                    outputLog << "Level " << k << ": Agent " << FullOrderBook[t][j].Bids[k].Agent << " bids " << FullOrderBook[t][j].Bids[k].Q << " Stocks " << FullOrderBook[t][j].Bids[k].Stock << " at $" << FullOrderBook[t][j].Bids[k].Price;
                    outputLog << " || ";
                    outputLog << "Agent " << FullOrderBook[t][j].Asks[k].Agent << " asks " << FullOrderBook[t][j].Asks[k].Q << " Stocks " << FullOrderBook[t][j].Asks[k].Stock << " at $" << FullOrderBook[t][j].Asks[k].Price << endl;
                    //}; // closes if
                }; // closes k loop
                if (CountOB==0) {outputLog << "Order book is empty!" << endl;};
                outputLog << endl;
            }; // closes j loop
            outputLog << "t=" << t << ", Outputed the Order Book at time t" << endl;
        }; // closes Plot condition
        // outputDebug << "t=" << t << ": OB outputted..." << endl;
        
        
        
        
        
        
        // ACTUAL TRADING VIA CLEAR() METHOD OF THE ORDER BOOK (SS11)
        for (int j=0; j<NumberOfStocks; j++) {
            // REMOVING METAORDERS EFFECTS (CREDITING THE METAORDER AGENT WITH ITS PREVIOUS Q OR RFA) WHICH MAY BRING THIS AGENT BANKRUPT
            double CreditBeforeTrading=0; double CreditAfterTrading=0;
            // Before trading
            if ((t==LastMetaorder) && (MetaorderImpact>0) && (FullOrderBook[t][0].MetaorderNature==0)) { // 0 for Bid and 1 for Ask
                CreditBeforeTrading=Market[FullOrderBook[t][0].MetaorderLastAgent].RFA;
                outputMetalog << "Before metaorder transaction, agent " << FullOrderBook[t][0].MetaorderLastAgent << " had a (meta) RFA of £" << CreditBeforeTrading << endl;
            }
            else if ((t==LastMetaorder) && (MetaorderImpact>0) && (FullOrderBook[t][0].MetaorderNature==1)) { // 0 for Bid and 1 for Ask
                CreditBeforeTrading=double(Market[FullOrderBook[t][0].MetaorderLastAgent].Stocks[0].StockQuantity);
                outputMetalog << "Before metaorder transaction, agent " << FullOrderBook[t][0].MetaorderLastAgent << " had a (meta) number of " << int(CreditBeforeTrading) << " stocks" << endl;
            };
            // ACTUAL TRADING
            vector<double> OBClearing = FullOrderBook[t][j].Clear(Market, t, j, Plot); // Actual transations
            
            // After trading
            if ((t==LastMetaorder) && (MetaorderImpact>0) && (FullOrderBook[t][0].MetaorderNature==0)) { // 0 for Bid and 1 for Ask
                CreditAfterTrading=Market[FullOrderBook[t][0].MetaorderLastAgent].RFA;
                Market[FullOrderBook[t][0].MetaorderLastAgent].RFA-=FullOrderBook[t][0].Credit-(CreditBeforeTrading-CreditAfterTrading);
                outputMetalog << "After metaorder transaction, agent " << FullOrderBook[t][0].MetaorderLastAgent << " had its unused cash of £" << FullOrderBook[t][0].Credit-(CreditBeforeTrading-CreditAfterTrading) << " removed and thus RFA lowered to £" << Market[FullOrderBook[t][0].MetaorderLastAgent].RFA << endl;
            }
            else if ((t==LastMetaorder) && (MetaorderImpact>0) && (FullOrderBook[t][0].MetaorderNature==1)) { // 0 for Bid and 1 for Ask
                CreditAfterTrading=double(Market[FullOrderBook[t][0].MetaorderLastAgent].Stocks[0].StockQuantity);
                Market[FullOrderBook[t][0].MetaorderLastAgent].Stocks[0].StockQuantity-=int(FullOrderBook[t][0].Credit)-(int(CreditBeforeTrading)-int(CreditAfterTrading));
                outputMetalog << "After metaorder transaction, agent " << FullOrderBook[t][0].MetaorderLastAgent << " had its unused number of " << int(FullOrderBook[t][0].Credit)-(int(CreditBeforeTrading)-int(CreditAfterTrading)) << " stocks removed and thus lowered to " << Market[FullOrderBook[t][0].MetaorderLastAgent].Stocks[0].StockQuantity << endl;
            };
            
            if (OBClearing[2]==0) { // If OB is empty
                gsl_matrix_set(ReflexiveValues, j, t, gsl_matrix_get (ReflexiveValues, j, t-1));
                gsl_matrix_set(AverageBids, j, t, gsl_matrix_get (AverageBids, j, t-1));
                gsl_matrix_set(AverageAsks, j, t, gsl_matrix_get (AverageAsks, j, t-1));
                gsl_matrix_set(TotalStockQuantityTraded, j, t, 0);
            } // closes if
            else {  // If OB is non-empty
                gsl_matrix_set(ReflexiveValues, j, t, OBClearing[0]); // Lowest of clearing prices is the market price, not the median
                // Calculating the average bids and asks at the top levels of the OB
                double MeanBid=0; int CountBid=int(FullOrderBook[t][j].Bids.size());
                for (int k=0; k<CountBid; k++) {MeanBid+=FullOrderBook[t][j].Bids[k].Price/CountBid;};
                double MeanAsk=0; int CountAsk=int(FullOrderBook[t][j].Asks.size());
                for (int k=0; k<CountAsk; k++) {MeanAsk+=FullOrderBook[t][j].Asks[k].Price/CountAsk;};
                if (CountBid==0) {MeanBid=gsl_matrix_get(AverageBids, j, t-1);};
                if (CountAsk==0) {MeanAsk=gsl_matrix_get(AverageAsks, j, t-1);};
                gsl_matrix_set(AverageBids, j, t, MeanBid);
                gsl_matrix_set(AverageAsks, j, t, MeanAsk);
                //gsl_matrix_set(TotalStockQuantityTraded, j, t, 10000*OBClearing[1]/TotalStockQuantities[j]); // Quantity of stock traded at time t in basis points of total quantity of stock RPR
                gsl_matrix_set(TotalStockQuantityTraded, j, t, OBClearing[1]); // Quantity of stock traded at time t in nominal
                if ((Plot=="On") && (Plot=="Mini")) {outputLog << "t=" << t << ", Actual trading via StockOrderBook.Clear()" << endl;};
            }; // closes else
            // FAILSAFE
            if (gsl_matrix_get(ReflexiveValues, j, t)<0.000001) {
                gsl_matrix_set(ReflexiveValues, j, t, 0.5*gsl_matrix_get(ReflexiveValues, j, t-1));
                //outputLog << "At t=" << t << " P(t) was zero => now P(t)=P(t-1)/2" << endl;
            };
        }; // closes j loop
        if (Plot=="On") {outputLog << "t=" << t << ", Actual trading via StockOrderBook.Clear()" << endl;};
        // outputDebug << "t=" << t << ": Clear() computed..." << endl;
        
        
        
        
        
        //MetaorderLastAgent=Bids[MetaorderLevel].Agent; // Agent placing the order
        //MetaorderNature=0; // 0 for Bid and 1 for Ask
        //Credit=ceil(Bids[MetaorderLevel].Q*(Bids[MetaorderLevel].Price+Asks[MetaorderLevel].Price)/2); // Extra RFA (for a buy) or stock quantity (for a sell) temporarily credited to the agent
        //Credit=double(Asks[MetaorderLevel].Q+1);
        
        
        
        // MARKET INDICATORS (SS12)
        // Capital evolution
        for (int i=0; i<NumberOfAgents ; i++) {
            double x = Market[i].Capital();
            double y = MarketEvolution[0][i].Capital();
            gsl_matrix_set(CapitalEvolution, i, t, (Market[i].Bankruptcy)*100*x/y);
        }; // closes i loop
        if (Plot=="On") {outputLog << "t=" << t << ", Indicator: CapitalEvolution" << endl;};
        
        // RELATIVE MAX DRAWDOWN
        // The (relative) max drawdown is the largest (adjusted to market) cumulated loss over a given period. We consider each year (after the second year) the lowest value of gsl_matrix * CapitalEvolution, and if this measure is below 100-Drawdown the agent reaches bankruptcy.
        
        
        // Computing the average of the biased values of all agents
        for (int j=0; j<NumberOfStocks; j++) {
            double Mean=0;
            for (int i=0; i<NumberOfAgents; i++) {
                Mean+=abs(-100+100*gsl_matrix_get (Market[i].BiasedValues, j, t)/gsl_matrix_get (ReflexiveValues, j, t)); // Absolute percentage of difference between agent biased and market price
            }; // closes i loop
            gsl_matrix_set (BiasedDelta, j, t, Mean/NumberOfAgents);
            gsl_matrix_set (TrueDelta, j, t, abs(-100+100*gsl_matrix_get (GSLTrueValues, j, t)/gsl_matrix_get (ReflexiveValues, j, t))); // Absolute percentage of difference between true and market price
        }; // closes j loop
        
        
        
        // CONSOLE INFO SCREEN
        TickEnd=(int(time(0))-TickStart); // Closing the clock for time of computation
        //system("CLS");
        cout << "****************************" << endl;
        cout << "** STOCK MARKET SIMULATOR **" << endl;
        cout << "****************************" << endl;
        cout << "I=" << NumberOfAgents << ", J=" << NumberOfStocks << ", T=" << Time << ", s=" << s << endl  << endl;
        cout.setf(ios::fixed);
        cout << "Initialization : " << TickMid << "s" << endl;
        cout << "Timestep       : " << t+1 << "/" << Time << endl;
        cout << "Elapsed time   : " << TickEnd/3600 << " h " << TickEnd/60 << " min " << TickEnd%60 << " s" << endl;
        cout << "Remaining time : " << (Time*TickEnd/t - TickEnd)/3600 << " h " << ((Time*TickEnd/t - TickEnd)%3600)/60 << " min " << (Time*TickEnd/t - TickEnd)%60 << " s" << endl << endl;
        cout << "Progress : " << 100*double(t)/(double(Time)) << " %" << endl;
        if (TickEnd%4==0) {cout << "                                    |" << endl;}
        else if (TickEnd%4==1) {cout << "                                    /" << endl;}
        else if (TickEnd%4==2) {cout << "                                    --" << endl;}
        else if (TickEnd%4==3) {cout << "                                    /" << endl;};
        if ((PDCondition=="PDOn") && (t>=1000) && ((t-1000)%Year2==0)) {outputLog << "Computing PDistances..." << endl;};
        
        // END OF t-LOOP (SS13)
        if (Plot=="On") {
            outputLog << "t=" << t << ", End of t-loop" << endl;
            outputLog << "***************************************************" << endl << endl << endl;
            // outputDebug << "t=" << t << ": End of t-loop" << endl;
        }; // closes Plot condition
        
        
        // Policy distances
        if ((PDCondition=="PDOn") && (t>=1000) && ((t-1000)%Year2==0)) {
            PDs.push_back(PolicyDistances (Market, NumberOfAgents, NumberOfStocks, Time, Percentile, t));
        };
        
        
    }; // closes t loop
    cout << "Closed t-loop..." << endl;
    outputMetalog.close();
    
    
    // Policy heat maps for best and worst populations
    if (PDCondition=="PDOn") {PoliciesAverages (Market, NumberOfAgents, NumberOfStocks, Percentile);};
    
    // Outputing the agents parameters and especially policies for eventual re-use
    if (OutputName!="OutputCondOff") {AgentParametersOutput (Market, NumberOfAgents, NumberOfStocks, OutputName);};
    
    // Outputting the percentage during the simulation of best and worst agents that belonged to the agent population within ClusterLimit
    outputClustered << 100*ClusteredBestAgent/(Time-LearningPhase) << '\t' << 100*ClusteredWorstAgent/(Time-LearningPhase) << endl;
    
    
    
    
    // MarketBidAskTopBidTopAskTrue of j stocks, and Spread of j stocks
    gsl_matrix * MarketBidAskTrue = gsl_matrix_alloc (12*NumberOfStocks, Time);
    gsl_matrix * MarketBidAskTrue2 = gsl_matrix_alloc (12*NumberOfStocks, Time-LearningPhase); // Without LearningPhase
    for (int j=0; j<NumberOfStocks; j++) {
        int Modulo=12*j;
        for (int t=0; t<Time; t++) {
            gsl_matrix_set (MarketBidAskTrue, 0+Modulo, t, gsl_matrix_get (ReflexiveValues, j, t));
            //gsl_matrix_set (MarketBidAskTrue, 1+Modulo, t, gsl_matrix_get (AverageBids, j, t));
            //gsl_matrix_set (MarketBidAskTrue, 2+Modulo, t, gsl_matrix_get (AverageAsks, j, t));
            //gsl_matrix_set (MarketBidAskTrue, 1+Modulo, t, gsl_matrix_get (AverageTopBids, j, t));
            //gsl_matrix_set (MarketBidAskTrue, 2+Modulo, t, gsl_matrix_get (AverageTopAsks, j, t));
            gsl_matrix_set (MarketBidAskTrue, 1+Modulo, t, gsl_matrix_get (GSLTrueValues, j, t));
            gsl_matrix_set (MarketBidAskTrue, 2+Modulo, t, gsl_matrix_get (HighestBid, j, t));
            gsl_matrix_set (MarketBidAskTrue, 3+Modulo, t, gsl_matrix_get (LowestAsk, j, t));
            if (t==0) {gsl_matrix_set (MarketBidAskTrue, 4+Modulo, t, 1);}
            else {gsl_matrix_set (MarketBidAskTrue, 4+Modulo, t, gsl_matrix_get (ReflexiveValues, j, t)/gsl_matrix_get (ReflexiveValues, j, t-1));}; // Returns P(t)/P(t-1)
            gsl_matrix_set (Spread, j, t, 100*VSpread[j]/gsl_matrix_get (ReflexiveValues, j, t));
            gsl_matrix_set (MarketBidAskTrue, 5+Modulo, t, 100*VSpread[j]/gsl_matrix_get (ReflexiveValues, j, t)); // Spread pct
            //gsl_matrix_set (MarketBidAskTrue, 6+Modulo, t, gsl_matrix_get (TotalStockQuantityTraded, j, t)); // Quantity of stock traded at time t in basis points of total quantity of stock RPR
            gsl_matrix_set (MarketBidAskTrue, 6+Modulo, t, gsl_matrix_get (TotalStockQuantityTraded, j, t)); // Quantity of stock traded at time t in nominal
            gsl_matrix_set (MarketBidAskTrue, 7+Modulo, t, gsl_matrix_get(Bankruptcies, 0, t)); // Bankruptcy pct
            gsl_matrix_set (MarketBidAskTrue, 8+Modulo, t, log(gsl_matrix_get(MarketBidAskTrue, 4+Modulo, t))); // Log-returns log[P(t)/P(t-1)]
            gsl_matrix_set (MarketBidAskTrue, 9+Modulo, t, abs(log(gsl_matrix_get(MarketBidAskTrue, 4+Modulo, t)))); // Absolute log-returns abs(log[P(t)/P(t-1)])
            //if ((t==0) || (gsl_matrix_get (TotalStockQuantityTraded, j, t-1)==0)) {gsl_matrix_set (MarketBidAskTrue, 10+Modulo, t, 1);}
            //else {gsl_matrix_set (MarketBidAskTrue, 10+Modulo, t, gsl_matrix_get (TotalStockQuantityTraded, j, t)/gsl_matrix_get (TotalStockQuantityTraded, j, t-1));}; // Volume returns
            gsl_matrix_set (MarketBidAskTrue, 10+Modulo, t, gsl_matrix_get (TotalStockQuantityTraded, j, t)); // Volume
            gsl_matrix_set (MarketBidAskTrue, 11+Modulo, t, gsl_matrix_get(Bankruptcies, 1, t)); // Market performance in % since Jan 1st
        }; // closes t loop
    }; // closes j loop
    
    // Populating MarketBidAskTrue2, which is like MarketBidAskTrue but without the first LearningPhase time steps
    for (int j=0; j<12*NumberOfStocks; j++) {
        for (int t=0; t<Time-LearningPhase; t++) {
            gsl_matrix_set(MarketBidAskTrue2, j, t, gsl_matrix_get(MarketBidAskTrue, j, t+LearningPhase));
        };
    };
    
    
    
    // Metrics: log-returns (& AC), absolute log-returns (& AC), volatilities (& AC), volumes (& AC)
    // Calibration: we can make distributions of stacked real/simulated log-returns, absolute log-returns, volatilities, volumes, but not AC's thereof.
    gsl_matrix * Moments = gsl_matrix_calloc (54*NumberOfStocks, Time);
    gsl_matrix * Moments2 = gsl_matrix_calloc (54*NumberOfStocks, Time-LearningPhase);
    for (int j=0; j<NumberOfStocks; j++) {
        int Modulo=12*j; // Size1 of MarketBidAskTrue above
        int ModuloMoments=54*j; // Size1 of Moments above
        for (int t=0; t<Time; t++) {
            gsl_matrix_set (Moments, 0+ModuloMoments, t, gsl_matrix_get (MarketBidAskTrue, 8+Modulo, t)); // Log-returns log[P(t)/P(t-1)]
            gsl_matrix_set (Moments, 1+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, Week, 8+Modulo)); // AC of log-returns at lag of Week
            gsl_matrix_set (Moments, 2+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, 2*Week, 8+Modulo)); // AC of log-returns at lag of 2*Week
            gsl_matrix_set (Moments, 3+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, Month, 8+Modulo)); // AC of log-returns at lag of Month
            gsl_matrix_set (Moments, 4+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, 3*Month, 8+Modulo)); // AC of log-returns at lag of 3*Month
            gsl_matrix_set (Moments, 5+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, 6*Month, 8+Modulo)); // AC of log-returns at lag of 6*Month
            gsl_matrix_set (Moments, 6+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, Year, 8+Modulo)); // AC of log-returns at lag of Year
            gsl_matrix_set (Moments, 7+ModuloMoments, t, gsl_matrix_get (MarketBidAskTrue, 9+Modulo, t)); // Absolute log-returns abs(log[P(t)/P(t-1)])
            gsl_matrix_set (Moments, 8+ModuloMoments, t, abs(MALAutoCorrelation (MarketBidAskTrue, t, Week, 9+Modulo))); // AC of absolute log-returns at lag of Week
            gsl_matrix_set (Moments, 9+ModuloMoments, t, abs(MALAutoCorrelation (MarketBidAskTrue, t, 2*Week, 9+Modulo))); // AC of absolute log-returns at lag of 2*Week
            gsl_matrix_set (Moments, 10+ModuloMoments, t, abs(MALAutoCorrelation (MarketBidAskTrue, t, Month, 9+Modulo))); // AC of absolute log-returns at lag of Month
            gsl_matrix_set (Moments, 11+ModuloMoments, t, abs(MALAutoCorrelation (MarketBidAskTrue, t, 3*Month, 9+Modulo))); // AC of absolute log-returns at lag of 3*Month
            gsl_matrix_set (Moments, 12+ModuloMoments, t, abs(MALAutoCorrelation (MarketBidAskTrue, t, 6*Month, 9+Modulo))); // AC of absolute log-returns at lag of 6*Month
            gsl_matrix_set (Moments, 13+ModuloMoments, t, abs(MALAutoCorrelation (MarketBidAskTrue, t, Year, 9+Modulo))); // AC of absolute log-returns at lag of Year
            gsl_matrix_set (Moments, 14+ModuloMoments, t, MALVolatility (MarketBidAskTrue, t, Week, 0+Modulo)); // Volatility at lag of Week
            gsl_matrix_set (Moments, 16+ModuloMoments, t, MALVolatility (MarketBidAskTrue, t, 2*Week, 0+Modulo)); // Volatility at lag of 2*Week
            gsl_matrix_set (Moments, 18+ModuloMoments, t, MALVolatility (MarketBidAskTrue, t, Month, 0+Modulo)); // Volatility at lag of Month
            gsl_matrix_set (Moments, 20+ModuloMoments, t, MALVolatility (MarketBidAskTrue, t, 3*Month, 0+Modulo)); // Volatility at lag of 3*Month
            gsl_matrix_set (Moments, 22+ModuloMoments, t, MALVolatility (MarketBidAskTrue, t, 6*Month, 0+Modulo)); // Volatility at lag of 6*Month
            gsl_matrix_set (Moments, 24+ModuloMoments, t, MALVolatility (MarketBidAskTrue, t, Year, 0+Modulo)); // Volatility at lag of Year
            gsl_matrix_set (Moments, 26+ModuloMoments, t, gsl_matrix_get (MarketBidAskTrue, 10+Modulo, t)); // Volumes returns
            gsl_matrix_set (Moments, 27+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, Week, 10+Modulo)); // AC of volumes returns at lag of Week
            gsl_matrix_set (Moments, 28+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, 2*Week, 10+Modulo)); // AC of volumes returns at lag of 2*Week
            gsl_matrix_set (Moments, 29+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, Month, 10+Modulo)); // AC of volumes returns at lag of Month
            gsl_matrix_set (Moments, 30+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, 3*Month, 10+Modulo)); // AC of volumes returns at lag of 3*Month
            gsl_matrix_set (Moments, 31+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, 6*Month, 10+Modulo)); // AC of volumes returns at lag of 6*Month
            gsl_matrix_set (Moments, 32+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, Year, 10+Modulo)); // AC of volumes returns at lag of Year
            
            
            
            // AC-logreturns of [t-∆, t] & [t-2∆, t-∆] for ∆={w, 2w, 3w, m}
            gsl_matrix_set(Moments, 33+ModuloMoments, t, MALAutoCorrelation (Moments, t, 3*Week, 0+ModuloMoments)); // AC of log-returns at lag of Week
            // AC-logreturns of [t-∆, t] & [t-∆-∂,t-∂] for ∂={1, 2, 3, 4, w} and ∆={w, 2w, 3w, m}
            gsl_matrix_set(Moments, 34+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, Week, 1, 0+ModuloMoments)); // AC of log-returns at lag of Week shifted by 1
            gsl_matrix_set(Moments, 35+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, Week, 2, 0+ModuloMoments)); // AC of log-returns at lag of Week shifted by 2
            gsl_matrix_set(Moments, 36+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, Week, 3, 0+ModuloMoments)); // AC of log-returns at lag of Week shifted by 3
            gsl_matrix_set(Moments, 37+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, Week, 4, 0+ModuloMoments)); // AC of log-returns at lag of Week shifted by 4
            gsl_matrix_set(Moments, 38+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, Week, Week, 0+ModuloMoments)); // AC of log-returns at lag of Week shifted by Week
            
            gsl_matrix_set(Moments, 39+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*1, 0+ModuloMoments)); // AC of log-returns at lag of 2Week shifted by 2
            gsl_matrix_set(Moments, 40+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*2, 0+ModuloMoments)); // AC of log-returns at lag of 2Week shifted by 4
            gsl_matrix_set(Moments, 41+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*3, 0+ModuloMoments)); // AC of log-returns at lag of 2Week shifted by 6
            gsl_matrix_set(Moments, 42+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*4, 0+ModuloMoments)); // AC of log-returns at lag of 2Week shifted by 8
            gsl_matrix_set(Moments, 43+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*Week, 0+ModuloMoments)); // AC of log-returns at lag of 2Week shifted by 2Week
            
            gsl_matrix_set(Moments, 44+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*1, 0+ModuloMoments)); // AC of log-returns at lag of 3Week shifted by 3
            gsl_matrix_set(Moments, 45+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*2, 0+ModuloMoments)); // AC of log-returns at lag of 3Week shifted by 6
            gsl_matrix_set(Moments, 46+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*3, 0+ModuloMoments)); // AC of log-returns at lag of 3Week shifted by 9
            gsl_matrix_set(Moments, 47+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*4, 0+ModuloMoments)); // AC of log-returns at lag of 3Week shifted by 12
            gsl_matrix_set(Moments, 48+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*Week, 0+ModuloMoments)); // AC of log-returns at lag of 3Week shifted by 3Week
            
            gsl_matrix_set(Moments, 49+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, Month, 4*1, 0+ModuloMoments)); // AC of log-returns at lag of Month shifted by 4
            gsl_matrix_set(Moments, 50+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, Month, 4*2, 0+ModuloMoments)); // AC of log-returns at lag of Month shifted by 8
            gsl_matrix_set(Moments, 51+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, Month, 4*3, 0+ModuloMoments)); // AC of log-returns at lag of Month shifted by 12
            gsl_matrix_set(Moments, 52+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, Month, 4*4, 0+ModuloMoments)); // AC of log-returns at lag of Month shifted by 16
            gsl_matrix_set(Moments, 53+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, Month, Month, 0+ModuloMoments)); // AC of log-returns at lag of Month shifted by Month
            
            
        }; // closes t loop
    }; // closes j loop
    for (int j=0; j<NumberOfStocks; j++) {
        int Modulo=12*j; // Size1 of MarketBidAskTrue above
        int ModuloMoments=54*j; // Size1 of Moments above
        for (int t=0; t<Time; t++) {
            gsl_matrix_set (Moments, 15+ModuloMoments, t, MALAutoCorrelation (Moments, t, Week, 14+Modulo)); // AC of volatility at lag of Week for lag of Week
            gsl_matrix_set (Moments, 17+ModuloMoments, t, MALAutoCorrelation (Moments, t, 2*Week, 16+Modulo)); // AC of volatility at lag of 2*Week for lag of 2*Week
            gsl_matrix_set (Moments, 19+ModuloMoments, t, MALAutoCorrelation (Moments, t, Month, 18+Modulo)); // AC of volatility at lag of Month for lag of Month
            gsl_matrix_set (Moments, 21+ModuloMoments, t, MALAutoCorrelation (Moments, t, 3*Month, 20+Modulo)); // AC of volatility at lag of 3*Month for lag of 2*Month
            gsl_matrix_set (Moments, 23+ModuloMoments, t, MALAutoCorrelation (Moments, t, 6*Month, 22+Modulo)); // AC of volatility at lag of 6*Month for lag of 6*Month
            gsl_matrix_set (Moments, 25+ModuloMoments, t, MALAutoCorrelation (Moments, t, Year, 24+Modulo)); // AC of volatility at lag of Year for lag of Year
        }; // closes t loop
    }; // closes j loop
    
    // Populating Moments2, which is like Moments but without the first LearningPhase time steps
    for (int j=0; j<54*NumberOfStocks; j++) {
        for (int t=0; t<Time-LearningPhase; t++) {
            gsl_matrix_set(Moments2, j, t, gsl_matrix_get(Moments, j, t+LearningPhase));
        };
    };
    
    /*
    // Meta order impact
    // We want to see at a given time t where a large order is sent (say once per year for each simulation) the metaorder impact on prices i1(t)=sigma_prices[t->t+tau]/sigma_prices[t-tau->t] (and likewise i2(t) on volumes) as a function of the metaorder size j(t)=Q_order(t)/Q_total(t), sigma being the standard deviation (of prices or volumes) over an interval of size tau, and Q being the stock quantity. For each of the S=20 simulation runs, we record and output the last 5 such tuples (j, i1(tau=w), i1(tau=2w), i1(tau=m), i1(tau=3m), i1(tau=6m), i2(tau=w), i2(tau=2w), i2(tau=m), i2(tau=3m), i2(tau=6m)), disregarding the first 5 because of being not yet learned by the agents. We do 4 such S=20 simulation runs for values of j=5%, 10%, 25%, 50%. The metaorder is sent during each run once every year by a random (yet non-bankrupt) agent whose stock holding is suddenly increased to these values (5%, 10%, 25%, 50%) of the total share outstanding at time t and then revert back to its value at time t+1.
    gsl_matrix * MetaorderImpactResults = gsl_matrix_calloc (16, int(MetaorderInjectionTimes.size()));
    int Tau=Week;
    for (int k=0; k<int(MetaorderInjectionTimes.size()); k++) {
        // Quantities
        gsl_matrix_set (MetaorderImpactResults, 0, k, 100*gsl_matrix_get (MarketBidAskTrue2, 10, MetaorderInjectionTimes[k])/SharesOutstanding); // j(t)=Q_order(t)/Q_total(t)
        //Prices
        Tau=Week; gsl_matrix_set (MetaorderImpactResults, 1, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 0)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 0)); // i1(t)=sigma_prices[t->t+tau]/sigma_prices[t-tau->t] for tau=we
        Tau=2*Week; gsl_matrix_set (MetaorderImpactResults, 2, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 0)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 0)); // i1(t)=sigma_prices[t->t+tau]/sigma_prices[t-tau->t] for tau=2we
        Tau=3*Week; gsl_matrix_set (MetaorderImpactResults, 3, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 0)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 0)); // i1(t)=sigma_prices[t->t+tau]/sigma_prices[t-tau->t] for tau=mo
        Tau=4*Week; gsl_matrix_set (MetaorderImpactResults, 4, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 0)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 0)); // i1(t)=sigma_prices[t->t+tau]/sigma_prices[t-tau->t] for tau=3mo
        Tau=5*Week; gsl_matrix_set (MetaorderImpactResults, 5, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 0)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 0)); // i1(t)=sigma_prices[t->t+tau]/sigma_prices[t-tau->t] for tau=6mo
        Tau=6*Week; gsl_matrix_set (MetaorderImpactResults, 6, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 0)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 0)); // i1(t)=sigma_prices[t->t+tau]/sigma_prices[t-tau->t] for tau=6mo
        Tau=7*Week; gsl_matrix_set (MetaorderImpactResults, 7, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 0)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 0)); // i1(t)=sigma_prices[t->t+tau]/sigma_prices[t-tau->t] for tau=6mo
        //Volumes
        Tau=Week; gsl_matrix_set (MetaorderImpactResults, 8, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 10)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 10)); // i1(t)=sigma_volumes[t->t+tau]/sigma_volumes[t-tau->t] for tau=we
        Tau=2*Week; gsl_matrix_set (MetaorderImpactResults, 9, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 10)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 10)); // i1(t)=sigma_volumes[t->t+tau]/sigma_volumes[t-tau->t] for tau=2we
        Tau=3*Week; gsl_matrix_set (MetaorderImpactResults, 10, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 10)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 10)); // i1(t)=sigma_volumes[t->t+tau]/sigma_volumes[t-tau->t] for tau=mo
        Tau=4*Week; gsl_matrix_set (MetaorderImpactResults, 11, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 10)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 10)); // i1(t)=sigma_volumes[t->t+tau]/sigma_volumes[t-tau->t] for tau=3mo
        Tau=5*Week; gsl_matrix_set (MetaorderImpactResults, 12, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 10)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 10)); // i1(t)=sigma_volumes[t->t+tau]/sigma_volumes[t-tau->t] for tau=6mo
        Tau=6*Week; gsl_matrix_set (MetaorderImpactResults, 13, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 10)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 10)); // i1(t)=sigma_volumes[t->t+tau]/sigma_volumes[t-tau->t] for tau=6mo
        Tau=7*Week; gsl_matrix_set (MetaorderImpactResults, 14, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 10)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 10)); // i1(t)=sigma_volumes[t->t+tau]/sigma_volumes[t-tau->t] for tau=6mo
        // Injection times
        gsl_matrix_set (MetaorderImpactResults, 15, k, MetaorderInjectionTimes[k]); // Times of metaorder injection
    };
    */
    
    
    
    
    // Bull market: -20% then +20% then -20% =>
    // Crash: sudden -20% => return Pt-20% and then no Pt+10% after a month
    // Recession: sudden -20% with no recovery =>
    gsl_matrix * SystemicRisk = gsl_matrix_calloc (15*NumberOfStocks, Time-LearningPhase); // Without LearningPhase
    gsl_matrix * Prices = gsl_matrix_calloc (1, Time-LearningPhase); for (int t=0; t<Time-LearningPhase; t++) {gsl_matrix_set (Prices, 0, t, gsl_matrix_get (ReflexiveValues, 0, t+LearningPhase));};
    for (int j=0; j<NumberOfStocks; j++) {
        int Modulo=15*j; int Modulo2=54*j; double Crash=0; double TrendCounter=0;
        for (int t=0; t<Time-LearningPhase; t++) {
            gsl_matrix_set (SystemicRisk, 0+Modulo, t, gsl_matrix_get (ReflexiveValues, j, t+LearningPhase)); // Price ($)
            gsl_matrix_set (SystemicRisk, 1+Modulo, t, gsl_matrix_get (BiasedDelta, j, t+LearningPhase)); // Percentage of average of differences between biased and market prices (%)
            gsl_matrix_set (SystemicRisk, 2+Modulo, t, gsl_matrix_get (TrueDelta, j, t+LearningPhase)); // Percentage of average of differences between biased and market prices (%)
            double FormalSpread=gsl_matrix_get(HighestBid, j, t+LearningPhase) - gsl_matrix_get(LowestAsk, j, t+LearningPhase); // Formal spread definition ($)
            gsl_matrix_set (SystemicRisk, 3+Modulo, t, 100*FormalSpread/gsl_matrix_get (ReflexiveValues, j, t+LearningPhase)); // Formal spread as a percentage of market price (%)
            if (t==0) {gsl_matrix_set (SystemicRisk, 4+Modulo, t, 0);}
            else {gsl_matrix_set (SystemicRisk, 4+Modulo, t, -100+100*gsl_matrix_get (ReflexiveValues, j, t+LearningPhase)/gsl_matrix_get (ReflexiveValues, j, t-1+LearningPhase));}; // Percentage of returns -100+100*P(t)/P(t-1) in
            gsl_matrix_set (SystemicRisk, 5+Modulo, t, gsl_matrix_get (TotalStockQuantityTraded, j, t+LearningPhase)); // Quantity of stock traded at time t in basis points of total quantity of stock RPR
            gsl_matrix_set (SystemicRisk, 6+Modulo, t, gsl_matrix_get (Bankruptcies, 0, t+LearningPhase)); // Bankruptcy pct
            gsl_matrix_set (SystemicRisk, 7+Modulo, t, gsl_matrix_get (Bankruptcies, 1, t+LearningPhase)); // Market performance in % since Jan 1st
            gsl_matrix_set (SystemicRisk, 8+Modulo, t, gsl_matrix_get (Moments, 14+Modulo2, t+LearningPhase)); // Week normalized volatility
            gsl_matrix_set (SystemicRisk, 9+Modulo, t, gsl_matrix_get (Moments, 18+Modulo2, t+LearningPhase)); // Month normalized volatility
            gsl_matrix_set (SystemicRisk, 10+Modulo, t, gsl_matrix_get (Moments, 22+Modulo2, t+LearningPhase)); // 6 month normalized volatility
            double Mean1=0; for (int k=t-2*Week; k<t-Week; k++) {Mean1+=gsl_matrix_get (ReflexiveValues, j, k+LearningPhase)/Week;};
            double Mean2=0; for (int k=t-Week; k<t; k++) {Mean2+=gsl_matrix_get (ReflexiveValues, j, k+LearningPhase)/Week;};
            gsl_matrix_set (SystemicRisk, 11+Modulo, t, -100+100*Mean2/Mean1); // -100+100*<P(t-w,t)/P(t-2w,t-w)> Allows to see crashes (<-20% in the distribution)
            if (gsl_matrix_get (SystemicRisk, 11+Modulo, t)<=-20) {Crash=1;} else {Crash=0;}; gsl_matrix_set (SystemicRisk, 12+Modulo, t, Crash); // Crash
            if (t>0) {
                if ((gsl_matrix_get (SystemicRisk, 11+Modulo, t-1)>0) && (gsl_matrix_get (SystemicRisk, 11+Modulo, t)<0)) {TrendCounter=0;}; // Trend reversion
                if ((gsl_matrix_get (SystemicRisk, 11+Modulo, t-1)<0) && (gsl_matrix_get (SystemicRisk, 11+Modulo, t)>0)) {TrendCounter=0;}; // Trend reversion
                if ((gsl_matrix_get (SystemicRisk, 11+Modulo, t-1)>0) && (gsl_matrix_get (SystemicRisk, 11+Modulo, t)>0)) {TrendCounter+=1;}; // Trend reversion
                if ((gsl_matrix_get (SystemicRisk, 11+Modulo, t-1)<0) && (gsl_matrix_get (SystemicRisk, 11+Modulo, t)<0)) {TrendCounter-=1;}; // Trend reversion
            };
            gsl_matrix_set (SystemicRisk, 13+Modulo, t, TrendCounter);
        }; // closes t loop
        
        
        PlotGSLMatrix(Prices, "Prices.xls", 1);
        //gsl_matrix * Hawkes = HistoricalIntensity(Prices, 0);
        //for (int t=0; t<Time-LearningPhase; t++) {gsl_matrix_set(SystemicRisk, 14+Modulo, t, gsl_matrix_get (Hawkes, 0, t));};
        
        /*
        vector<double> V, GoF; for (int t=0; t<Time-LearningPhase; t++) {V.push_back(gsl_matrix_get (Prices, 0, t)); GoF.push_back(1);};
        vector<double> * pV= & V;
        //const double Rescale=10; // Rescale factor // POSTDOC2022
        //const double Threshold=0.1; // Threshold (0.1) // POSTDOC2022
        //const size_t NumIteration=10000000; // Number of iterations of the minimizer (stops before once a minimum is found) // POSTDOC2022
        double R; double * pR = & R; // Hawkes process parameters (does not need to be initialized)
        //vector<double> & pGoF = GoF; // Goodness of fit // POSTDOC2022
        //vector<double> * pRes = Intensity (pV, Rescale, Threshold, NumIteration, pR, pGoF); // POSTDOC2022
        vector<double> * pRes = 0;
        for (int t=0; t<Time-LearningPhase; t++) {gsl_matrix_set (SystemicRisk, 14, t, (*pRes)[t]);}; // Hawkes process intensity in column 9
        
        for (int t=0; t<Time-LearningPhase; t++) { // Hawkes process intensity in column 9
            if ((*pRes)[t]>10) { // Maximum value of Hawkes at 10
                gsl_matrix_set (SystemicRisk, 14, t, 10);
            }
            else if ((*pRes)[t]<0) { // Minimum value of Hawkes at 0
                gsl_matrix_set (SystemicRisk, 14, t, 0);
            };
        };
        
        pV=0; pR=0; pRes=0;
        //delete pV; delete pR; delete pRes;
        */
        
        
        
        
    }; // closes j loop
    
    
    
    
    // Volumes of j stocks, VolumesDistributions of j stocks, AutoCorrLagPVolumes of j stocks
    ReturnsDistributions = GSLDistribution(MarketBidAskTrue, 20);
    VolumesDistributions = GSLDistribution(TotalStockQuantityTraded, 1000);
    AutoCorrLagPVolumes = AutocorrelationLagP (TotalStockQuantityTraded);
    
    // Most successful and less successful agents
    vector<double> SummedCapitalEvolution;
    for (int i=0; i<NumberOfAgents; i++) {
        SummedCapitalEvolution.push_back(0);
        //for (int t=Time/2; t<Time; t++) {SummedCapitalEvolution[i]+=gsl_matrix_get(CapitalEvolution, i, t);};
        for (int t=0; t<Time; t++) {SummedCapitalEvolution[i]+=gsl_matrix_get(CapitalEvolution, i, t);};
    }; // closes i loop
    vector<int> DescendingRank;
    for (int k=0; k<NumberOfAgents; k++) {
        double Bar=0; int BarX=0;
        for (int i=0; i<NumberOfAgents; i++) {
            if ((SummedCapitalEvolution[i]>=Bar) && (SummedCapitalEvolution[i]>=0)) {Bar=SummedCapitalEvolution[i]; BarX=i;};
        }; // closes i loop
        SummedCapitalEvolution[BarX]=-1;
        DescendingRank.push_back(BarX);
    }; // closes k loop
    int AgentNumberPercentile=Percentile*NumberOfAgents/100;
    if (AgentNumberPercentile<=5) {AgentNumberPercentile=5;};
    gsl_matrix * MostSuccessful = gsl_matrix_alloc (AgentNumberPercentile, Time);
    gsl_matrix * LessSuccessful = gsl_matrix_alloc (AgentNumberPercentile, Time);
    for (int t=0; t<Time; t++) {
        for (int i=0; i<AgentNumberPercentile; i++) {
            gsl_matrix_set (MostSuccessful, i, t, gsl_matrix_get(CapitalEvolution, DescendingRank[i], t));};
        for (int i=0; i<AgentNumberPercentile; i++) {
            gsl_matrix_set (LessSuccessful, i, t, gsl_matrix_get(CapitalEvolution, DescendingRank[i+NumberOfAgents-AgentNumberPercentile], t));};
    }; // closes t loop
    
    
    //Plot the distribution of each parameter of the AgentNumberPercentile% most/less successful agents (cf. Versatility, Gesture, and all NEBs: NEB1, NEBLossAversion, NEB3, NEB4, NEBPositivity, NEB6, NEB1p, NEB3p, NEB4p, NEBPositivity, NEBLearningRate, NEB8, NEB9).
    gsl_matrix * MostSuccessfulParameters = gsl_matrix_alloc (7, AgentNumberPercentile);
    gsl_matrix * LessSuccessfulParameters = gsl_matrix_alloc (7, AgentNumberPercentile);
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (MostSuccessfulParameters, 0, i, Market[DescendingRank[i]].Future);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (LessSuccessfulParameters, 0, i, Market[DescendingRank[i+NumberOfAgents-AgentNumberPercentile]].Future);};
    gsl_matrix_set (MostSuccessfulParameters, 1, 0, 0); // To make sure the distributions will be between 0 and 100!! (X1<=X<X2)
    gsl_matrix_set (MostSuccessfulParameters, 1, 1, 100); // To make sure the distributions will be between 0 and 100!!
    gsl_matrix_set (LessSuccessfulParameters, 1, 0, 0); // To make sure the distributions will be between 0 and 100!!
    gsl_matrix_set (LessSuccessfulParameters, 1, 1, 100); // To make sure the distributions will be between 0 and 100!!
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (MostSuccessfulParameters, 1, i, Market[DescendingRank[i]].Reflexivity*100);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (LessSuccessfulParameters, 1, i, Market[DescendingRank[i+NumberOfAgents-AgentNumberPercentile]].Reflexivity*100);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (MostSuccessfulParameters, 2, i, Market[DescendingRank[i]].TradingWindow);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (LessSuccessfulParameters, 2, i, Market[DescendingRank[i+NumberOfAgents-AgentNumberPercentile]].TradingWindow);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (MostSuccessfulParameters, 3, i, Market[DescendingRank[i]].Gesture);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (LessSuccessfulParameters, 3, i, Market[DescendingRank[i+NumberOfAgents-AgentNumberPercentile]].Gesture);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (MostSuccessfulParameters, 4, i, Market[DescendingRank[i]].History);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (LessSuccessfulParameters, 4, i, Market[DescendingRank[i+NumberOfAgents-AgentNumberPercentile]].History);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (MostSuccessfulParameters, 5, i, Market[DescendingRank[i]].NEBLearningRate);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (LessSuccessfulParameters, 5, i, Market[DescendingRank[i+NumberOfAgents-AgentNumberPercentile]].NEBLearningRate);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (MostSuccessfulParameters, 6, i, Market[DescendingRank[i]].NEB);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (LessSuccessfulParameters, 6, i, Market[DescendingRank[i+NumberOfAgents-AgentNumberPercentile]].NEB);};
    
    
    
    // Ploting the distributions of these parameters
    gsl_matrix * MostSuccessfulParametersDistributions = gsl_matrix_alloc (AgentNumberPercentile, 2*7);
    gsl_matrix * LessSuccessfulParametersDistributions = gsl_matrix_alloc (AgentNumberPercentile, 2*7);
    MostSuccessfulParametersDistributions = GSLDistribution(MostSuccessfulParameters, 10);
    LessSuccessfulParametersDistributions = GSLDistribution(LessSuccessfulParameters, 10);
    //int pBots, int pNEBLossAversion24, int pNEBDelayDiscounting56, int pNEB89, int pNEBLossAversion23456, int pNEBLossAversion2489, int pNEBDelayDiscounting5689
    // Zipf's Law and Pareto's Principle
    gsl_matrix * SortedNAV = gsl_matrix_alloc (1, NumberOfAgents);
    gsl_matrix * ZipfDistribution = gsl_matrix_alloc (2, 50);
    gsl_matrix * Pareto = gsl_matrix_alloc (2, NumberOfAgents); // Percentage of agents who own a percentage of total market capital (80% should own around 20% according to Pareto's Principle)
    double TotalCapital=0;
    for (int i=0; i<NumberOfAgents; i++) {
        TotalCapital+=Market[DescendingRank[i]].Capital();
    };
    double TempTotalCapital=0;
    double ParetoPercentage=0;
    for (int i=0; i<NumberOfAgents; i++) {
        gsl_matrix_set (SortedNAV, 0, i, Market[DescendingRank[i]].Capital());
        TempTotalCapital+=100*(Market[DescendingRank[i]].Capital())/TotalCapital;
        ParetoPercentage=100*i/NumberOfAgents;
        gsl_matrix_set (Pareto, 0, i, ParetoPercentage);
        gsl_matrix_set (Pareto, 1, i, TempTotalCapital);
    };
    ZipfDistribution = GSLDistribution(SortedNAV, 50); // Distribution
    gsl_matrix * MAINDistribution = GSLDistribution(MarketBidAskTrue2, 50); // Distribution of MAIN
    gsl_matrix * MostSuccessfulDistribution = GSLDistribution(MostSuccessfulParameters, 50); // Distribution of MostSuccessfulParameters
    gsl_matrix * LessSuccessfulDistribution = GSLDistribution(LessSuccessfulParameters, 50); // Distribution of LessSuccessfulParameters
    gsl_matrix * MomentsDistribution = GSLDistribution(Moments2, 50); // Distribution of Moments
    
    
    // Output // XYZ
    ofstream outputMain(Machine+"MAIN.xls", ofstream::app);
    outputMain << "Market($)" << '\t' << "True($)" << '\t' << "Bid($)" << '\t' << "Ask($)" << '\t' << "Return" << '\t' << "Spread(%)" << '\t' << "Volume(bsp)" << '\t' << "Banruptcy(%)" << '\t' << "LogReturns" << '\t' << "AbsLogReturns" << '\t' << "Volumes(bsp)" << '\t' << "MarketPerfYTD(%)" << endl;
    //PlotGSLMatrix(MarketBidAskTrue2, "MAIN.xls", 1); // Contains info for all j stocks: market, average bids, average asks, top bids, top asks, True, daily returns pct, spread pct, volume bsp
    PlotGSLMatrix(MAINDistribution, "MAINDistribution.xls", 1);
    outputMain.close();
    //PlotGSLMatrix(MarketBidAskTrue, "MarketBidAskTrueReturnspctSpreadpctVolumebspBankrupcypct.xls", 1); // Contains info for all j stocks: market, average bids, average asks, top bids, top asks, True, daily returns pct, spread pct, volume bsp
    //PlotGSLMatrix(Spread, "MSpread.xls", 1); // Contains info for all j stocks
    //PlotGSLMatrix(ReturnsDistributions, "ReturnsDistributions.xls", 1);
    //PlotGSLMatrix(TotalStockQuantityTraded, "MVolumes.xls", 1);
    //PlotGSLMatrix(VolumesDistributions, "VolumesDistributions.xls", 1);
    //PlotGSLMatrix(AutoCorrLagPVolumes, "AutoCorrLagPVolumes.xls", 1);
    //PlotGSLMatrix(MostSuccessful, "MostSuccessful.xls", 1);
    //PlotGSLMatrix(LessSuccessful, "LessSuccessful.xls", 1);
    
    ofstream outputSystemic(Machine+"SystemicRisk.xls", ofstream::app);
    outputSystemic << "Market($)" << '\t' << "Biased/Pt(%)" << '\t' << "True/Pt(%)" << '\t' << "FormalSpread(%)" << '\t' << "Return(%)" << '\t' << "Volume(bsp)" << '\t' << "Banruptcy(%)" << '\t' << "MarketPerfYTD(%)" << '\t' << "w-Volatility" << '\t' << "m-Volatility" << '\t' << "6m-Volatility" << '\t' << "P[t-w,t]/P[t-2w,t-w](%)" << '\t' << "Crash" << '\t' << "TrendCounter" << '\t' << "Hawkes" << endl;
    //PlotGSLMatrix(SystemicRisk, "SystemicRisk.xls", 1); // Contains info for all j stocks: market, average bids, average asks, top bids, top asks, True, daily returns pct, spread pct, volume bsp
    outputSystemic.close();
    
    ofstream outputSystemicDistribution(Machine+"SystemicRiskDistribution.xls", ofstream::app);
    outputSystemicDistribution<< "Market($)" << '\t' << " " << '\t' << "Biased/Pt(%)" << '\t' << " " << '\t' << "True/Pt(%)" << '\t' << " " << '\t' << "FormalSpread(%)" << '\t' << " " << '\t' << "Return(%)" << '\t' << " " << '\t' << "Volume(bsp)" << '\t' << " " << '\t' << "Banruptcy(%)" << '\t' << " " << '\t' << "MarketPerfYTD(%)" << '\t' << " " << '\t' << "w-Volatility" << '\t' << " " << '\t' << "m-Volatility" << '\t' << " " << '\t' << "6m-Volatility" << '\t' << " " << '\t' << "P[t-w,t]/P[t-2w,t-w](%)" << '\t' << " " << '\t' << "Crash" << '\t' << " " << '\t' << "TrendCounter" << '\t' << " " << '\t' << "Hawkes" << endl;
    PlotGSLMatrix(GSLDistribution(SystemicRisk, 50), "SystemicRiskDistribution.xls", 1);
    outputSystemicDistribution.close();
    
    
    /*
    //MetaorderImpactResults
    //We want to see at a given time t where a large order is sent (say once per year for each simulation) the metaorder impact on prices i1(t)=sigma_prices[t->t+tau]/sigma_prices[t-tau->t] (and likewise i2(t) on volumes) as a function of the metaorder size j(t)=Q_order(t)/Q_total(t), sigma being the standard deviation (of prices or volumes) over an interval of size tau, and Q being the stock quantity. For each of the S=20 simulation runs, we record and output the last 5 such tuples (j, i1(tau=w), i1(tau=2w), i1(tau=m), i1(tau=3m), i1(tau=6m), i2(tau=w), i2(tau=2w), i2(tau=m), i2(tau=3m), i2(tau=6m)), disregarding the first 5 because of being not yet learned by the agents. We do 4 such S=20 simulation runs for values of j=5%, 10%, 25%, 50%. The metaorder is sent during each run once every year by a random (yet non-bankrupt) agent whose stock holding is suddenly increased to these values (5%, 10%, 25%, 50%) of the total share outstanding at time t and then revert back to its value at time t+1.
    ofstream outputMetaorderImpactResults(Machine+"MetaorderImpact_size" + IntToString(MetaorderImpact) + "pct.xls", ofstream::app);
    outputMetaorderImpactResults<< "Q/Qtot(%)" << '\t' << "Psig(t,t+w)/Psig(t-w,t)(%)" << '\t' << "Psig(t,t+2w)/Psig(t-2w,t)(%)" << '\t' << "Psig(t,t+3w)/Psig(t-3w,t)(%)" << '\t' << "Psig(t,t+4w)/Psig(t-4w,t)(%)" << '\t' << "Psig(t,t+5w)/Psig(t-5w,t)(%)" << '\t' << "Psig(t,t+6w)/Psig(t-6w,t)(%)" << '\t' << "Psig(t,t+7w)/Psig(t-7w,t)(%)" << '\t' << "Vsig(t,t+w)/Vsig(t-w,t)(%)" << '\t' << "Vsig(t,t+2w)/Vsig(t-2w,t)(%)" << '\t' << "Vsig(t,t+3w)/Vsig(t-3w,t)(%)" << '\t' << "Vsig(t,t+4w)/Vsig(t-4w,t)(%)" << '\t' << "Vsig(t,t+5w)/Vsig(t-5w,t)(%)" << '\t' << "Vsig(t,t+6w)/Vsig(t-6w,t)(%)" << '\t' << "Vsig(t,t+7w)/Vsig(t-7w,t)(%)" << '\t' << "t" << endl;
    PlotGSLMatrix(MetaorderImpactResults, "MetaorderImpact_size" + IntToString(MetaorderImpact) + "pct.xls", 1);
    outputMetaorderImpactResults.close();
    */
    
    ofstream outputSuccessful(Machine+"MostLessSuccessfulParameters.xls", ofstream::app);
    outputSuccessful<< "Future" << '\t' << "Reflexivity(%)" << '\t' << "TradingWindow" << '\t' << "Gesture" << '\t' << "History" << '\t' << "NEBLearningRate" << '\t' << "NEB" << endl;
    PlotGSLMatrix(MostSuccessfulParameters, "MostLessSuccessfulParameters.xls", 1);
    PlotGSLMatrix(LessSuccessfulParameters, "MostLessSuccessfulParameters.xls", 1);
    PlotGSLMatrix(MostSuccessfulDistribution, "MostLessSuccessfulParametersDistribution.xls", 1);
    PlotGSLMatrix(LessSuccessfulDistribution, "MostLessSuccessfulParametersDistribution.xls", 1);
    
    //PlotGSLMatrix(CapitalEvolution, "CapitalEvolution.xls", 1);
    
    outputSuccessful.close();
    //PlotGSLMatrix(MostSuccessfulParametersDistributions, "MostSuccessfulParametersDistributions.xls", 1);
    //PlotGSLMatrix(LessSuccessfulParametersDistributions, "LessSuccessfulParametersDistributions.xls", 1);
    PlotGSLMatrix(SortedNAV, "Zipf&Distribution.xls", 1);
    PlotGSLMatrix(ZipfDistribution, "Zipf&Distribution.xls", 1);
    ofstream outputMoments(Machine+"Moments.xls", ofstream::app);
    outputMoments << "Log-return" << '\t' << "AC-1w Log-return" << '\t' << "AC-2w Log-return" << '\t' << "AC-m Log-return" << '\t' << "AC-3m Log-return" << '\t' << "AC-6m Log-return" << '\t' << "AC-y Log-return" << '\t' << "Abs-log-return" << '\t' << "AC-1w Abs-log-return" << '\t' << "AC-2w Abs-log-return" << '\t' << "AC-m Abs-log-return" << '\t' << "AC-3m Abs-log-return" << '\t' << "AC-6m Abs-log-return" << '\t' << "AC-y Abs-log-return" << '\t' << "w-Volatility" << '\t' << "AC-w w-Volatility" << '\t' << "2w-Volatility" << '\t' << "AC-2w 2w-Volatility" << '\t' << "m-Volatility" << '\t' << "AC-m m-Volatility" << '\t' << "3m-Volatility" << '\t' << "AC-3m 3m-Volatility" << '\t' << "6m-Volatility" << '\t' << "AC-6m 6m-Volatility" << '\t' << "y-Volatility" << '\t' << "AC-y y-Volatility" << '\t' << "Volumes(bsp)" << '\t' << "AC-1w Volume" << '\t' << "AC-2w Volume" << '\t' << "AC-m Volume" << '\t' << "AC-3m Volume" << '\t' << "AC-6m Volume" << '\t' << "AC-y Volume" << '\t' << "AC-3w Log-return" << '\t' << "b1AC-1w Log-return" << '\t' << "b2AC-1w Log-return" << '\t' << "b3AC-1w Log-return" << '\t' << "b4AC-1w Log-return" << '\t' << "bwAC-1w Log-return" << '\t' << "b2AC-2w Log-return" << '\t' << "b4AC-2w Log-return" << '\t' << "b6AC-2w Log-return" << '\t' << "b8AC-2w Log-return" << '\t' << "b2wAC-2w Log-return" << '\t' << "b3AC-3w Log-return" << '\t' << "b6AC-3w Log-return" << '\t' << "b9AC-3w Log-return" << '\t' << "b12AC-3w Log-return" << '\t' << "b3wAC-3w Log-return" << '\t' << "b4AC-m Log-return" << '\t' << "b8AC-m Log-return" << '\t' << "b12AC-m Log-return" << '\t' << "b16AC-m Log-return" << '\t' << "bmAC-m Log-return" << endl;
    
    //PlotGSLMatrix(Moments2, "Moments.xls", 1);
    ofstream outputMomentsDistribution(Machine+"MomentsDistribution.xls", ofstream::app);
    outputMomentsDistribution << "Log-return" << '\t' << " " << '\t' << "AC-1w Log-return" << '\t' << " " << '\t' << "AC-2w Log-return" << '\t' << " " << '\t' << "AC-m Log-return" << '\t' << " " << '\t' << "AC-3m Log-return" << '\t' << " " << '\t' << "AC-6m Log-return" << '\t' << " " << '\t' << "AC-y Log-return" << '\t' << " " << '\t' << "Abs-log-return" << '\t' << " " << '\t' << "AC-1w Abs-log-return" << '\t' << " " << '\t' << "AC-2w Abs-log-return" << '\t' << " " << '\t' << "AC-m Abs-log-return" << '\t' << " " << '\t' << "AC-3m Abs-log-return" << '\t' << " " << '\t' << "AC-6m Abs-log-return" << '\t' << " " << '\t' << "AC-y Abs-log-return" << '\t' << " " << '\t' << "w-Volatility" << '\t' << " " << '\t' << "AC-w w-Volatility" << '\t' << " " << '\t' << "2w-Volatility" << '\t' << " " << '\t' << "AC-2w 2w-Volatility" << '\t' << " " << '\t' << "m-Volatility" << '\t' << " " << '\t' << "AC-m m-Volatility" << '\t' << " " << '\t' << "3m-Volatility" << '\t' << " " << '\t' << "AC-3m 3m-Volatility" << '\t' << " " << '\t' << "6m-Volatility" << '\t' << " " << '\t' << "AC-6m 6m-Volatility" << '\t' << " " << '\t' << "y-Volatility" << '\t' << " " << '\t' << "AC-y y-Volatility" << '\t' << " " << '\t' << "Volumes(bsp)" << '\t' << " " << '\t' << "AC-1w Volume" << '\t' << " " << '\t' << "AC-2w Volume" << '\t' << " " << '\t' << "AC-m Volume" << '\t' << " " << '\t' << "AC-3m Volume" << '\t' << " " << '\t' << "AC-6m Volume" << '\t' << " " << '\t' << "AC-y Volume" << '\t' << " " << '\t' << "AC-3w Log-return" << '\t' << " " << '\t' << "b1AC-1w Log-return" << '\t' << " " << '\t' << "b2AC-1w Log-return" << '\t' << " " << '\t' << "b3AC-1w Log-return" << '\t' << " " << '\t' << "b4AC-1w Log-return" << '\t' << " " << '\t' << "bwAC-1w Log-return" << '\t' << " " << '\t' << "b2AC-2w Log-return" << '\t' << " " << '\t' << "b4AC-2w Log-return" << '\t' << " " << '\t' << "b6AC-2w Log-return" << '\t' << " " << '\t' << "b8AC-2w Log-return" << '\t' << " " << '\t' << "b2wAC-2w Log-return" << '\t' << " " << '\t' << "b3AC-3w Log-return" << '\t' << " " << '\t' << "b6AC-3w Log-return" << '\t' << " " << '\t' << "b9AC-3w Log-return" << '\t' << " " << '\t' << "b12AC-3w Log-return" << '\t' << " " << '\t' << "b3wAC-3w Log-return" << '\t' << " " << '\t' << "b4AC-m Log-return" << '\t' << " " << '\t' << "b8AC-m Log-return" << '\t' << " " << '\t' << "b12AC-m Log-return" << '\t' << " " << '\t' << "b16AC-m Log-return" << '\t' << " " << '\t' << "bmAC-m Log-return" << endl;
    PlotGSLMatrix(MomentsDistribution, "MomentsDistribution.xls", 1);
    //PlotGSLMatrix(Pareto, "Pareto.xls", 1);
    //cout << "Produced outputs..." << endl;
    
    
    
    // Memory freeing // JJJ10
    gsl_matrix_free(Moments);
    gsl_matrix_free(MarketBidAskTrue);
    gsl_matrix_free(MarketBidAskTrue2);
    gsl_matrix_free(Prices);
    gsl_matrix_free(ReflexiveValues);
    gsl_matrix_free(AverageBids);
    gsl_matrix_free(AverageAsks);
    gsl_matrix_free(AverageTopBids);
    gsl_matrix_free(AverageTopAsks);
    gsl_matrix_free(MedianTopBids);
    gsl_matrix_free(MedianTopAsks);
    gsl_matrix_free(HighestBid);
    gsl_matrix_free(LowestAsk);
    gsl_matrix_free(GSLTrueValues);
    gsl_matrix_free(Spread);
    gsl_matrix_free(TotalStockQuantityTraded);
    gsl_matrix_free(VolumesDistributions);
    gsl_matrix_free(ReturnsDistributions);
    gsl_matrix_free(AutoCorrLagPVolumes);
    gsl_matrix_free(CapitalEvolution);
    gsl_matrix_free(Bankruptcies);
    gsl_matrix_free(BiasedDelta);
    gsl_matrix_free(TrueDelta);
    gsl_matrix_free(MomentsDistribution);
    gsl_matrix_free(ZipfDistribution);
    gsl_matrix_free(LessSuccessfulDistribution);
    gsl_matrix_free(MostSuccessfulDistribution);
    gsl_matrix_free(MAINDistribution);
    //gsl_matrix_free(MetaorderImpactResults);
    
    for (int i=0; i<NumberOfAgents; i++) {gsl_matrix_free(Market[i].BiasedValues);};
    
    
    
    outputLog.close();
    cout << "Outputted .xls files..." << endl;
    
    // BACKUP OF MATRIX RESULTS
    vector<gsl_matrix *> MatrixResults;
    MatrixResults.push_back(Moments2);
    MatrixResults.push_back(MostSuccessfulParameters);
    MatrixResults.push_back(LessSuccessfulParameters);
    MatrixResults.push_back(SortedNAV);
    MatrixResults.push_back(SystemicRisk);
    for (int k=0; k<int(PDs.size()); k++) {MatrixResults.push_back(PDs[k]);}; // PDistances
    
    
    Market.erase (Market.begin(), Market.end()); // JJJ10
    
    cout << "DONE" << endl;
    
    return MatrixResults;
};
















// Each run of MarketSimulator() outputs vector<gsl_matrix *> MatrixResults. We then compute S of these runs and save them in a vector of length S called vector<vector<gsl_matrix *> > MultiSim. Now function Bollinger() gets that MultiSim and accesses the m matrix of all its MatrixResults simulations (m=0 corresponds to MarketBidAskTrue, m=1 to MostSuccessfulParameters, m=2 to LessSuccessfulParameters, m=3 to SortedNAV). It then accesses all the cxr elements of all S matrix m and outputs their mean in a matrix of same dimension than matrix m.
vector<gsl_matrix *> Bollinger (vector<vector<gsl_matrix *> > MultiSim, int m) {
    vector<gsl_matrix * > Res;
    gsl_matrix * ResMean = gsl_matrix_calloc (MultiSim[0][m]->size1, MultiSim[0][m]->size2);
    gsl_matrix * ResMeanPlusSigma = gsl_matrix_calloc (MultiSim[0][m]->size1, MultiSim[0][m]->size2);
    gsl_matrix * ResMeanMinusSigma = gsl_matrix_calloc (MultiSim[0][m]->size1, MultiSim[0][m]->size2);
    gsl_matrix * ResMeanPlus2Sigma = gsl_matrix_calloc (MultiSim[0][m]->size1, MultiSim[0][m]->size2);
    gsl_matrix * ResMeanMinus2Sigma = gsl_matrix_calloc (MultiSim[0][m]->size1, MultiSim[0][m]->size2);
    int MSsize=int(MultiSim.size());
    for (int c=0; c<int(MultiSim[0][m]->size1); c++) {
        for (int r=0; r<int(MultiSim[0][m]->size2); r++) {
            double Mean=0; for (int s=0; s<MSsize; s++) {Mean+=gsl_matrix_get (MultiSim[s][m], c, r)/MSsize;};
            gsl_matrix_set (ResMean, c, r, Mean);
            double Variance=0; for (int s=0; s<MSsize; s++) {Variance+=(Mean - gsl_matrix_get (MultiSim[s][m], c, r)) * (Mean - gsl_matrix_get (MultiSim[s][m], c, r))/MSsize;};
            gsl_matrix_set (ResMeanPlusSigma, c, r, Mean+sqrt(Variance));
            gsl_matrix_set (ResMeanMinusSigma, c, r, Mean-sqrt(Variance));
            gsl_matrix_set (ResMeanPlusSigma, c, r, Mean+2*sqrt(Variance));
            gsl_matrix_set (ResMeanMinusSigma, c, r, Mean-2*sqrt(Variance));
        }; // closes r loop
    }; // closes c loop
    string A = "BollingerMean_m"; string B = IntToString(m); string C = ".xls"; A+=B+C; PlotGSLMatrix(ResMean, A.c_str(), 1);
    A = "BollingerMeanPlusSigma_m"; B = IntToString(m); C = ".xls"; A+=B+C; PlotGSLMatrix(ResMeanPlusSigma, A.c_str(), 1);
    A = "BollingerMeanMinusSigma_m"; B = IntToString(m); C = ".xls"; A+=B+C; PlotGSLMatrix(ResMeanMinusSigma, A.c_str(), 1);
    A = "BollingerMeanPlus2Sigma_m"; B = IntToString(m); C = ".xls"; A+=B+C; PlotGSLMatrix(ResMeanPlus2Sigma, A.c_str(), 1);
    A = "BollingerMeanMinus2Sigma_m"; B = IntToString(m); C = ".xls"; A+=B+C; PlotGSLMatrix(ResMeanMinus2Sigma, A.c_str(), 1);
    Res.push_back(ResMean); Res.push_back(ResMeanPlusSigma); Res.push_back(ResMeanMinusSigma); Res.push_back(ResMeanPlus2Sigma); Res.push_back(ResMeanMinus2Sigma);
    return Res;
}; // closes Bollinger()



// Computing mean, variance, skewness, kurtosis, median, q5, q25, q75, q95 of every column of a matrix M
gsl_matrix * Momenta (gsl_matrix * M) {
    int Size1=int(M->size1);
    int Size2=int(M->size2);
    gsl_matrix * Res = gsl_matrix_alloc (2*Size1, 9);
    for (int i=0; i<Size1; i++) {
        double Mean=0; double Variance=0; double Stdev=0; double Skewness=0; double Kurtosis=0;
        for (int j=0; j<Size2; j++) {Mean+=gsl_matrix_get (M, i, j)/Size2;};
        for (int j=0; j<Size2; j++) {Variance+=(gsl_matrix_get (M, i, j)-Mean)*(gsl_matrix_get (M, i, j)-Mean)/Size2;}; Stdev=sqrt(Variance);
        for (int j=0; j<Size2; j++) {
            Skewness+=pow(((gsl_matrix_get (M, i, j)-Mean)/Stdev), 3);
            Kurtosis+=pow(((gsl_matrix_get (M, i, j)-Mean)/Stdev), 4);
        }; // closes j loop
        vector<double> V; for (int j=0; j<Size2; j++) {V.push_back(gsl_matrix_get (M, i, j));};
        sort (V.begin(), V.begin()+Size2); // Sorting V ascendingly
        double Median=V[50*Size2/100];
        double q5=V[5*Size2/100];
        double q25=V[25*Size2/100];
        double q75=V[75*Size2/100];
        double q95=V[95*Size2/100];
        // Filling the matrix result Res
        gsl_matrix_set (Res, 1+2*i, 0, Mean);
        gsl_matrix_set (Res, 1+2*i, 1, Variance);
        gsl_matrix_set (Res, 1+2*i, 2, Skewness);
        gsl_matrix_set (Res, 1+2*i, 3, Kurtosis);
        gsl_matrix_set (Res, 1+2*i, 4, Median);
        gsl_matrix_set (Res, 1+2*i, 5, q5);
        gsl_matrix_set (Res, 1+2*i, 6, q25);
        gsl_matrix_set (Res, 1+2*i, 7, q75);
        gsl_matrix_set (Res, 1+2*i, 8, q95);
        V.erase(V.begin(), V.end());
    }; // closes i loop
    return Res;
}; // closes Momenta()



// Each run of MarketSimulator() outputs vector<gsl_matrix *> MatrixResults. We then compute S of these runs and save them in a vector of length S called vector<vector<gsl_matrix *> > MultiSim. Now function JointDistributions() gets that MultiSim and accesses the m matrix of all its MatrixResults simulations. It then accesses the [c,r]-element of all S matrices m and outputs their joint distributions over all S matrices and output a matrix similar to m but with twice its number of columns (each column has a range for 10 bins between Xmin and Xmax, and another for the number of counts per bin)
void JointDistributions (vector<vector<gsl_matrix *> > MultiSim, int m, int Precision) {
    int MSsize = int(MultiSim.size());
    int Size1=int(MultiSim[0][m]->size1);
    int Size2=int(MultiSim[0][m]->size2);
    gsl_matrix * Res = gsl_matrix_alloc (Size1, (MSsize*Size2));
    for (int s=0; s<MSsize; s++) {
        for (int c=0; c<Size1; c++) {
            for (int r=0; r<Size2; r++) {
                gsl_matrix_set (Res, c, s*Size2+r, gsl_matrix_get (MultiSim[s][m], c, r));
            }; // closes r loop
        }; // closes c loop
    }; // closes s loop
    //PlotGSLMatrix(Res, "StackedRes.xls", 1);
    gsl_matrix * JointDistributions = GSLDistribution(Res, Precision);
    //string A = "JointDistributions_m"; string B = IntToString(m); string C = ".xls"; A+=B+C;
    string A = "JointDistributions.xls";
    ofstream outputJointDistribution(Machine+"JointDistributions.xls", ofstream::app);
    if (m==0) {
        outputJointDistribution<< "Moments" << endl;
        outputJointDistribution << "Log-return" << '\t' << " " << '\t' << "AC-1w Log-return" << '\t' << " " << '\t' << "AC-2w Log-return" << '\t' << " " << '\t' << "AC-m Log-return" << '\t' << " " << '\t' << "AC-3m Log-return" << '\t' << " " << '\t' << "AC-6m Log-return" << '\t' << " " << '\t' << "AC-y Log-return" << '\t' << " " << '\t' << "Abs-log-return" << '\t' << " " << '\t' << "AC-1w Abs-log-return" << '\t' << " " << '\t' << "AC-2w Abs-log-return" << '\t' << " " << '\t' << "AC-m Abs-log-return" << '\t' << " " << '\t' << "AC-3m Abs-log-return" << '\t' << " " << '\t' << "AC-6m Abs-log-return" << '\t' << " " << '\t' << "AC-y Abs-log-return" << '\t' << " " << '\t' << "w-Volatility" << '\t' << " " << '\t' << "AC-w w-Volatility" << '\t' << " " << '\t' << "2w-Volatility" << '\t' << " " << '\t' << "AC-2w 2w-Volatility" << '\t' << " " << '\t' << "m-Volatility" << '\t' << " " << '\t' << "AC-m m-Volatility" << '\t' << " " << '\t' << "3m-Volatility" << '\t' << " " << '\t' << "AC-3m 3m-Volatility" << '\t' << " " << '\t' << "6m-Volatility" << '\t' << " " << '\t' << "AC-6m 6m-Volatility" << '\t' << " " << '\t' << "y-Volatility" << '\t' << " " << '\t' << "AC-y y-Volatility" << '\t' << " " << '\t' << "Volumes(bsp)" << '\t' << " " << '\t' << "AC-1w Volume" << '\t' << " " << '\t' << "AC-2w Volume" << '\t' << " " << '\t' << "AC-m Volume" << '\t' << " " << '\t' << "AC-3m Volume" << '\t' << " " << '\t' << "AC-6m Volume" << '\t' << " " << '\t' << "AC-y Volume" << '\t' << " " << '\t' << "AC-3w Log-return" << '\t' << " " << '\t' << "b1AC-1w Log-return" << '\t' << " " << '\t' << "b2AC-1w Log-return" << '\t' << " " << '\t' << "b3AC-1w Log-return" << '\t' << " " << '\t' << "b4AC-1w Log-return" << '\t' << " " << '\t' << "bwAC-1w Log-return" << '\t' << " " << '\t' << "b2AC-2w Log-return" << '\t' << " " << '\t' << "b4AC-2w Log-return" << '\t' << " " << '\t' << "b6AC-2w Log-return" << '\t' << " " << '\t' << "b8AC-2w Log-return" << '\t' << " " << '\t' << "b2wAC-2w Log-return" << '\t' << " " << '\t' << "b3AC-3w Log-return" << '\t' << " " << '\t' << "b6AC-3w Log-return" << '\t' << " " << '\t' << "b9AC-3w Log-return" << '\t' << " " << '\t' << "b12AC-3w Log-return" << '\t' << " " << '\t' << "b3wAC-3w Log-return" << '\t' << " " << '\t' << "b4AC-m Log-return" << '\t' << " " << '\t' << "b8AC-m Log-return" << '\t' << " " << '\t' << "b12AC-m Log-return" << '\t' << " " << '\t' << "b16AC-m Log-return" << '\t' << " " << '\t' << "bmAC-m Log-return" << endl;
    }
    else if (m==1) {
        outputJointDistribution<< "MostSuccessfulParameters" << endl;
        outputJointDistribution<< "Future" << '\t' << " " << '\t' << "Reflexivity(%)" << '\t' << " " << '\t' << "TradingWindow" << '\t' << " " << '\t' << "Gesture" << '\t' << " " << '\t' << "History" << '\t' << " " << '\t' << "Epsilon" << '\t' << " " << '\t' << "NEB" << endl;
    }
    else if (m==2) {
        outputJointDistribution<< "LessSuccessfulParameters" << endl;
        outputJointDistribution<< "Future" << '\t' << " " << '\t' << "Reflexivity(%)" << '\t' << " " << '\t' << "TradingWindow" << '\t' << " " << '\t' << "Gesture" << '\t' << " " << '\t' << "History" << '\t' << " " << '\t' << "Epsilon" << '\t' << " " << '\t' << "NEB" << endl;
    }
    else if (m==3) {
        outputJointDistribution<< "Pareto" << endl;
        outputJointDistribution<< "NAV" << '\t' << " " << '\t' << endl;
    }
    else if (m==4) {
        outputJointDistribution << "Systemics" << endl;
        outputJointDistribution<< "Market($)" << '\t' << " " << '\t' << "Biased/Pt(%)" << '\t' << " " << '\t' << "True/Pt(%)" << '\t' << " " << '\t' << "FormalSpread(%)" << '\t' << " " << '\t' << "Return(%)" << '\t' << " " << '\t' << "Volume(bsp)" << '\t' << " " << '\t' << "Banruptcy(%)" << '\t' << " " << '\t' << "MarketPerfYTD(%)" << '\t' << " " << '\t' << "w-Volatility" << '\t' << " " << '\t' << "m-Volatility" << '\t' << " " << '\t' << "6m-Volatility" << '\t' << " " << '\t' << "P[t-w,t]/P[t-2w,t-w](%)" << '\t' << " " << '\t' << "Crash" << '\t' << " " << '\t' << "TrendCounter" << '\t' << " " << '\t' << "Hawkes" << endl;
    }
    else {
        outputJointDistribution << "PDistances_y" << m-5 << endl;
        outputJointDistribution<< "BestToBest" << '\t' << " " << '\t' << "BestToAll" << '\t' << " " << '\t' << "BestToWorst" << '\t' << " " << '\t' << "WorstToAll" << '\t' << " " << '\t' << "WorstToWorst" << endl;
    };
    
    PlotGSLMatrix(JointDistributions, A.c_str(), 1);
    PlotGSLMatrix(Momenta(Res), A.c_str(), 1);
    
}; // closes JointDistributions()






// Function pauses code for a time T in s
void Pause (int T) {
    int Start=int(time(0));
    while (int(time(0))-Start<T) {continue;};
}







// Takes a string name of the csv file as input (N2 can be something like "AAPL.csv") and returns a vector<double> for a matrix of dimension MxN at batch m (in a file of stacked matrices). For Input() to work we must 1- pre-specify the dimensions of the input, 2- the input must be converted from xls to csv (by changing the name file, not via LibreOffice)
gsl_matrix * InputOLD (int M, int N, string S1, int m) {
    gsl_matrix * Res = gsl_matrix_alloc (M, N);
    string S2=Machine; S2+=S1; ifstream FileInput (S2);
    string Line;
    int j = -1 + m*(N+3);
    while (getline (FileInput, Line)) {
        j++;
        istringstream LineStream(Line);
        string Item;
        int i = -1;
        while (getline (LineStream, Item, '\t')) {
            //while (getline (LineStream, Item, ',')) {
            i++;
            if (j<N) {gsl_matrix_set (Res, i, j-m*(N+3), atof(Item.c_str()));};
        }
    }
    string S3="M"; string S3a=IntToString(m); string S3b=".txt"; S3+=S3a+S3b;
    PlotGSLMatrix(Res, S3.c_str(), 1);
    FileInput.close();
    return Res;
};




// Takes a string name of the csv file as input (N2 can be something like "AAPL.csv") and returns a vector<double> for a matrix of dimension MxN at batch m (in a file of stacked matrices). For Input() to work we must 1- pre-specify the dimensions of the input, 2- the input must be converted from xls to csv (by changing the name file, not via LibreOffice)
gsl_matrix * Input (int M, int N, string S1, int m, int Letter, int Opt, int Period) {
    gsl_matrix * Res = gsl_matrix_alloc (M, N);
    string S2=Machine; S2+=S1; ifstream FileInput (S2);
    string Line;
    //for (int j=m*(N+2); j<m*(N+2)+Period; j++) {
    //for (int j=m*(N+2)+N-Period; j<m*(N+2)+N; j++) {
    for (int j=0; j<m*(N+2)+N; j++) {
        getline (FileInput, Line);
        istringstream LineStream(Line);
        string Item;
        int i = -1;
        if (Opt==0) {
            while (getline (LineStream, Item, '\t')) {
                i++;
                //gsl_matrix_set (Res, i, j-(m*(N+2)+N-Period), atof(Item.c_str()));
                if (j>=m*(N+2)) {
                    gsl_matrix_set (Res, i, j-(m*(N+2)), atof(Item.c_str()));
                };
            }
        }
        else {
            while (getline (LineStream, Item, ',')) {
                i++;
                gsl_matrix_set (Res, i, j-(m*(N+2)+N-Period), atof(Item.c_str()));
            }
        }
    }
    gsl_matrix * Res2 = gsl_matrix_alloc (M, Period);
    for (int i=0; i<M; i++) {
        for (int j=N-Period; j<N; j++) {
            gsl_matrix_set (Res2, i, j-N+Period, gsl_matrix_get (Res, i, j));
        };
    };
    string S3="/OUTPUT/M"; string S3a=IntToString(m+Letter); string S3b=".txt"; S3+=S3a+S3b;
    PlotGSLMatrix(Res2, S3.c_str(), 1);
    FileInput.close();
    return Res2;
};




// Takes a string name of the csv file as input (N2 can be something like "AAPL.csv") and returns a vector<double> for a matrix of dimension MxN at batch m (in a file of stacked matrices). For Input() to work we must 1- pre-specify the dimensions of the input, 2- the input must be converted from xls to csv (by changing the name file, not via LibreOffice)
gsl_matrix * CSVInput (int M, int N, string S1, int Opt) {
    gsl_matrix * Res = gsl_matrix_alloc (M, N);
    string S2=Machine; S2+=S1; ifstream FileInput (S2);
    string Line;
    for (int j=0; j<=N; j++) {
        getline (FileInput, Line);
        istringstream LineStream(Line);
        string Item;
        int i = -1;
        if (Opt==0) {
            while (getline (LineStream, Item, '\t')) {
                i++;
                gsl_matrix_set (Res, i, j, atof(Item.c_str()));
            }
        }
        else {
            while (getline (LineStream, Item, ',')) {
                i++;
                if ((j==0) || (i==0)) {continue;};
                gsl_matrix_set (Res, i-1, j-1, atof(Item.c_str()));
            }
        }
    }
    // Inverting the matrix column elements
    gsl_matrix * Res2 = gsl_matrix_alloc (M, N);
    for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++) {
            gsl_matrix_set (Res2, i, j, gsl_matrix_get (Res, i, N-j-1));
        };
    };
    string S3="/OUTPUT/M_"; string S3b=".txt"; S3+=S1+S3b; PlotGSLMatrix(Res, S3.c_str(), 1);
    FileInput.close();
    return Res2;
};


// Computes the factorial of an integer
double Factorial (int n) {
    double Res=1;
    for (int i=1; i<=n; i++) {Res*=i;}
    return Res;
}
//In a set of n elements, the number of possible subsets of k elements (when non-ordered and without repetition) is given by the number of k-combination C(n,k)=n!/(k!(n-k)!)
double Combination (int n, int k) {
    double Res=Factorial(n)/(Factorial(k) * Factorial(n-k));
    return Res;
}
//Number of all poissble k-combinations of a set of n elements
double TotalCombination (int n) {
    double Res=n;
    for (int k=2; k<=n; k++) {Res+=Combination(n,k);}
    return Res;
}
void TotalCombinationMatrix (int N) {
    gsl_matrix * Res = gsl_matrix_calloc (1, N+1);
    for (int k=1; k<=N; k++) {
        gsl_matrix_set (Res, 0, k, TotalCombination (k));
        cout << "Combination at N=" << k << "/" << N << endl;
    }
    PlotGSLMatrix(Res, "TotalCombinationMatrix.csv", 1);
};
//TotalCombinationMatrix(100); exit(0);





void HPMarketSimulatorAll (int HPI, int HPGesture, double HPTrueMu, int HPAccuracy, int HPTime, int HPLiquidationFloor, string TypeNEB, string Leader, int ClusterLimit, int pNEB, double Rate, int S) {
    vector<vector<gsl_matrix *> > MultiSim;
    // int NEB0=50; int NEB1=13; int NEBLossAversion=13; int NEB3=12; int NEB4=12; // resp. bots, conservatism bias, loss aversion, positivity bias
    // HPAccuracy=9, 10, 11, 13 => probability of 0.2, 0.75, 1, 2 per year of bubble burst
    int J=1; double r=0.01; // int T=LearningPhase+2875; int LearningPhase=1000;
    for (int s=0; s<S; s++) {
        //MarketSimulator (int NumberOfAgents, int NumberOfStocks, int Time, double Rate, string Plot, string PDCondition, string TypeNEB, int HPGesture, double HPTrueMu, int HPAccuracy, int LiquidationFloor, string LeaderType, int ClusterLimit, int s)
        vector<gsl_matrix *> MSim = MarketSimulator (HPI, J, HPTime, r, "Off", "PDOff", TypeNEB, HPGesture, HPTrueMu, HPAccuracy, HPLiquidationFloor, Leader, ClusterLimit, pNEB, Rate, s, -1, 0);
        MultiSim.push_back(MSim); //delete MSim[0]; MSim[0]=NULL; MSim.erase(MSim.begin(), MSim.end());
    };
    for (int k=0; k<int(MultiSim[0].size()); k++) {JointDistributions (MultiSim, k, 100);};
    // Memory freeing
    for (int k=0; k<int(MultiSim[0].size()); k++) {
        for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][k]);};
    };
    MultiSim.erase(MultiSim.begin(), MultiSim.end()); // Simulated moments distribution
    /*
     JointDistributions (MultiSim, 0, 100); // Moments2
     JointDistributions (MultiSim, 1, 100); // MostSuccessfullParameters
     JointDistributions (MultiSim, 2, 100); // LessSuccessfullParameters
     JointDistributions (MultiSim, 3, 100); // SortedNAV
     JointDistributions (MultiSim, 4, 100); // SystemicRisk
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][0]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][1]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][2]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][3]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][4]);}; // JJJ10
     MultiSim.erase(MultiSim.begin(), MultiSim.end()); // Simulated moments distribution
     */
};
void HPMarketSimulatorAll2 (int HPI, int HPGesture, double HPTrueMu, int HPAccuracy, int HPTime, int HPLiquidationFloor, string TypeNEB, string Leader, int ClusterLimit, int pNEB, string PD, double Rate, int S) {
    vector<vector<gsl_matrix *> > MultiSim;
    // int NEB0=50; int NEB1=13; int NEBLossAversion=13; int NEB3=12; int NEB4=12; // resp. bots, conservatism bias, loss aversion, positivity bias
    // HPAccuracy=9, 10, 11, 13 => probability of 0.2, 0.75, 1, 2 per year of bubble burst
    int J=1; double r=0.01; // int T=LearningPhase+2875; int LearningPhase=1000;
    for (int s=0; s<S; s++) {
        //MarketSimulator (int NumberOfAgents, int NumberOfStocks, int Time, double Rate, string Plot, string PDCondition, string TypeNEB, int HPGesture, double HPTrueMu, int HPAccuracy, int LiquidationFloor, string LeaderType, int ClusterLimit, int s)
        vector<gsl_matrix *> MSim = MarketSimulator (HPI, J, HPTime, r, "Off", PD, TypeNEB, HPGesture, HPTrueMu, HPAccuracy, HPLiquidationFloor, Leader, ClusterLimit, pNEB, Rate, s, -1, 0);
        MultiSim.push_back(MSim); //delete MSim[0]; MSim[0]=NULL; MSim.erase(MSim.begin(), MSim.end());
    };
    for (int k=0; k<int(MultiSim[0].size()); k++) {JointDistributions (MultiSim, k, 100);};
    // Memory freeing
    for (int k=0; k<int(MultiSim[0].size()); k++) {
        for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][k]);};
    };
    MultiSim.erase(MultiSim.begin(), MultiSim.end()); // Simulated moments distribution
    /*
     JointDistributions (MultiSim, 0, 100); // Moments2
     JointDistributions (MultiSim, 1, 100); // MostSuccessfullParameters
     JointDistributions (MultiSim, 2, 100); // LessSuccessfullParameters
     JointDistributions (MultiSim, 3, 100); // SortedNAV
     JointDistributions (MultiSim, 4, 100); // SystemicRisk
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][0]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][1]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][2]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][3]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][4]);}; // JJJ10
     MultiSim.erase(MultiSim.begin(), MultiSim.end()); // Simulated moments distribution
     */
};
void MarketSimulatorTrunk (int HPI, int HPGesture, double HPTrueMu, int HPAccuracy, int HPTime, int HPLiquidationFloor, string TypeNEB, string Leader, int ClusterLimit, int pNEB, double Rate, int S, int Trunk) {
    vector<vector<gsl_matrix *> > MultiSim;
    // int NEB0=50; int NEB1=13; int NEBLossAversion=13; int NEB3=12; int NEB4=12; // resp. bots, conservatism bias, loss aversion, positivity bias
    // HPAccuracy=9, 10, 11, 13 => probability of 0.2, 0.75, 1, 2 per year of bubble burst
    int J=1; double r=0.01; // int T=LearningPhase+2875; int LearningPhase=1000;
    for (int s=0; s<S; s++) {
        //MarketSimulator (int NumberOfAgents, int NumberOfStocks, int Time, double Rate, string Plot, string PDCondition, string TypeNEB, int HPGesture, double HPTrueMu, int HPAccuracy, int LiquidationFloor, string LeaderType, int ClusterLimit, int s)
        vector<gsl_matrix *> MSim = MarketSimulator (HPI, J, HPTime, r, "On", "PDOff", TypeNEB, HPGesture, HPTrueMu, HPAccuracy, HPLiquidationFloor, Leader, ClusterLimit, pNEB, Rate, s, Trunk, 0);
        MultiSim.push_back(MSim); //delete MSim[0]; MSim[0]=NULL; MSim.erase(MSim.begin(), MSim.end());
    };
    for (int k=0; k<int(MultiSim[0].size()); k++) {JointDistributions (MultiSim, k, 100);};
    // Memory freeing
    for (int k=0; k<int(MultiSim[0].size()); k++) {
        for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][k]);};
    };
    MultiSim.erase(MultiSim.begin(), MultiSim.end()); // Simulated moments distribution
    /*
     JointDistributions (MultiSim, 0, 100); // Moments2
     JointDistributions (MultiSim, 1, 100); // MostSuccessfullParameters
     JointDistributions (MultiSim, 2, 100); // LessSuccessfullParameters
     JointDistributions (MultiSim, 3, 100); // SortedNAV
     JointDistributions (MultiSim, 4, 100); // SystemicRisk
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][0]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][1]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][2]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][3]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][4]);}; // JJJ10
     MultiSim.erase(MultiSim.begin(), MultiSim.end()); // Simulated moments distribution
     */
};

void HPMarketSimulatorAllPD (int HPI, int HPGesture, double HPTrueMu, int HPAccuracy, int HPTime, int HPLiquidationFloor, string TypeNEB, string Leader, int ClusterLimit, int pNEB, double Rate, int S) {
    vector<vector<gsl_matrix *> > MultiSim;
    // int NEB0=50; int NEB1=13; int NEBLossAversion=13; int NEB3=12; int NEB4=12; // resp. bots, conservatism bias, loss aversion, positivity bias
    // HPAccuracy=9, 10, 11, 13 => probability of 0.2, 0.75, 1, 2 per year of bubble burst
    int J=1; double r=0.01; // int T=LearningPhase+2875; int LearningPhase=1000;
    for (int s=0; s<S; s++) {
        //MarketSimulator (int NumberOfAgents, int NumberOfStocks, int Time, double Rate, string Plot, string PDCondition, string TypeNEB, int HPGesture, double HPTrueMu, int HPAccuracy, int LiquidationFloor, string LeaderType, int ClusterLimit, int s)
        vector<gsl_matrix *> MSim = MarketSimulator (HPI, J, HPTime, r, "Off", "PDOn", TypeNEB, HPGesture, HPTrueMu, HPAccuracy, HPLiquidationFloor, Leader, ClusterLimit, pNEB, Rate, s, -1, 0);
        MultiSim.push_back(MSim); //delete MSim[0]; MSim[0]=NULL; MSim.erase(MSim.begin(), MSim.end());
    };
    for (int k=0; k<int(MultiSim[0].size()); k++) {JointDistributions (MultiSim, k, 100);};
    // Memory freeing
    for (int k=0; k<int(MultiSim[0].size()); k++) {
        for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][k]);};
    };
    MultiSim.erase(MultiSim.begin(), MultiSim.end()); // Simulated moments distribution
    /*
     JointDistributions (MultiSim, 0, 100); // Moments2
     JointDistributions (MultiSim, 1, 100); // MostSuccessfullParameters
     JointDistributions (MultiSim, 2, 100); // LessSuccessfullParameters
     JointDistributions (MultiSim, 3, 100); // SortedNAV
     JointDistributions (MultiSim, 4, 100); // SystemicRisk
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][0]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][1]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][2]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][3]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][4]);}; // JJJ10
     MultiSim.erase(MultiSim.begin(), MultiSim.end()); // Simulated moments distribution
     */
};

void HPMarketSimulator (int T, int S) {HPMarketSimulatorAll (500, 0, 0.10, 10, T, 50, "Classic", "None", 1, 0, 1, S);};

void HPMarketSimulator2 (string TypeNEB, string Leader, int ClusterLimit, int S) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "Classic", Leader, ClusterLimit, 0, 1, S);};

void HPMarketSimulatorNEB (string TypeNEB, int S) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, TypeNEB, "None", 1, 0, 1, S);};

void HPMarketSimulatorDummy () {HPMarketSimulatorAll (500, 0, 0.10, 10, 1000+1000, 50, "Classic", "None", 1, 0, 1, 3);};





void HPMarketSimulatorAllOB (int HPI, int HPGesture, double HPTrueMu, int HPAccuracy, int HPTime, int HPLiquidationFloor, string TypeNEB, string Leader, int ClusterLimit, int pNEB, double Rate, int S, int MetaOrderImpact) {
    vector<vector<gsl_matrix *> > MultiSim;
    int J=1; double r=0.01; // int T=LearningPhase+2875; int LearningPhase=1000;
    for (int s=0; s<S; s++) {
        //MarketSimulator (int NumberOfAgents, int NumberOfStocks, int Time, double Rate, string Plot, string PDCondition, string TypeNEB, int HPGesture, double HPTrueMu, int HPAccuracy, int LiquidationFloor, string LeaderType, int ClusterLimit, int s)
        vector<gsl_matrix *> MSim = MarketSimulator (HPI, J, HPTime, r, "Off", "PDOff", TypeNEB, HPGesture, HPTrueMu, HPAccuracy, HPLiquidationFloor, Leader, ClusterLimit, pNEB, Rate, s, -1, MetaOrderImpact);
        MultiSim.push_back(MSim); //delete MSim[0]; MSim[0]=NULL; MSim.erase(MSim.begin(), MSim.end());
    };
    for (int k=0; k<int(MultiSim[0].size()); k++) {JointDistributions (MultiSim, k, 100);};
    // Memory freeing
    for (int k=0; k<int(MultiSim[0].size()); k++) {
        for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][k]);};
    };
    MultiSim.erase(MultiSim.begin(), MultiSim.end()); // Simulated moments distribution
};



/*
 // Memory freeing check
 gsl_matrix * Coucou1 (gsl_matrix * M) {
 gsl_matrix * Res = gsl_matrix_calloc (100, 100);
 for (int i=0; i<int(M->size1); i++) {
 for (int j=0; j<int(M->size2); j++) {gsl_matrix_set (Res, i, j, gsl_matrix_get (M, i, j));};
 };
 return Res;
 };
 gsl_matrix * Coucou2 (gsl_matrix * M) {
 gsl_matrix * Res = gsl_matrix_calloc (100, 100);
 for (int i=0; i<int(M->size1); i++) {
 for (int j=0; j<int(M->size2); j++) {gsl_matrix_set (Res, i, j, gsl_matrix_get (M, i, j));};
 };
 return Res;
 };
 gsl_matrix * MemoryFreeingTest () {
 gsl_matrix * A = gsl_matrix_calloc (100, 100);
 for (int i=0; i<int(A->size1); i++) {
 for (int j=0; j<int(A->size2); j++) {gsl_matrix_set (A, i, j, i*j);};
 };
 gsl_matrix * Temp = Coucou1(A); gsl_matrix * Res = Coucou2(Temp);
 gsl_matrix_free(A); gsl_matrix_free(Temp);
 return Res;
 };
 */




/*
 vector<double> CheckTrueGeneration (double HPTrueMu, double HPLeash) {
 vector<double> Res;
 // INITIALIZING AGENTS
 int NumberOfAgents=10; int NumberOfStocks=1; int Time=1*Year; int HPAccuracy=10;
 gsl_rng * r = make_rng(); gsl_rng_set(r, static_cast<unsigned long int>(time(0)));
 gsl_matrix * GSLTrueValues = gsl_matrix_calloc (NumberOfStocks, Time);
 //gsl_matrix * Spikes = gsl_matrix_calloc (1, NumberOfStocks);
 gsl_matrix * TrueBiased = gsl_matrix_calloc (NumberOfAgents+1, Time); // For j=0 only
 vector<Agent> Market;
 for (int i=0; i<NumberOfAgents; i++) {
 Agent TempAgent; TempAgent.BiasedValues = gsl_matrix_calloc (NumberOfStocks, Time);
 for (int j=0; j<NumberOfStocks; j++) {
 TempAgent.Leash.push_back(HPLeash*0.001*gsl_rng_uniform (r));
 TempAgent.LeashVol.push_back(0.05*gsl_rng_uniform (r));
 TempAgent.Accuracy.push_back(HPAccuracy+gsl_rng_uniform (r));
 };
 Market.push_back(TempAgent);
 };
 // GENERATING TRUE VALUES
 vector<double> Gen = STLRandom (NumberOfStocks, "Uniform", "NoPlot");
 double AnnualCounter=0; double AnnualAmplitude=0;
 vector<double> VAnnualCounter, VAnnualAmplitude;
 for (int j=0; j<NumberOfStocks; j++) {
 vector<double> S = PoissonRandomWalk (100, Week+HPTrueMu*3*Month*Gen[j], int(0.1*Time/250+0.9*Gen[j]*Time/250), Time, 1000*Gen[j]);
 for (int t=0; t<Time; t++) {
 gsl_matrix_set (GSLTrueValues, j, t, S[t]);
 gsl_matrix_set (TrueBiased, 0, t, S[t]);
 };
 int Counter=0; double Amplitude=0;
 for (int t=1; t<Time; t++) {
 if (S[t]-S[t-1]>0.001) {Counter+=1; Amplitude+=100*(S[t]-S[t-1])/S[t];};
 };
 AnnualCounter+=Counter/(Time/Year);
 AnnualAmplitude+=Amplitude/Counter;
 VAnnualCounter.push_back(Counter/(Time/Year));
 VAnnualAmplitude.push_back(Amplitude/Counter);
 //cout << "Counter=" << Counter/(Time/Year) << ", Amplitude=" << Amplitude/Counter << endl;
 }; // closes j loop
 AnnualCounter/=NumberOfStocks;
 AnnualAmplitude/=NumberOfStocks;
 cout << "HPTrueMu=" << HPTrueMu << ", HPLeash=" << HPLeash << ": ";
 cout << "AnnualCounter=" << AnnualCounter << ", AnnualAmplitude=" << AnnualAmplitude << "%";
 //PlotGSLMatrix(Spikes, "Spikes.csv", 1);
 vector<vector<double> > TrueValues = GSLMatrixToSTLMatrix(GSLTrueValues);
 // GENERATING BIASED VALUES
 double BiasedDistance=0; vector<double> VBiasedDistance;
 for (int i=0; i<NumberOfAgents; i++) {
 for (int j=0; j<NumberOfStocks; j++) {
 vector<double> C = CointegratedWalk(TrueValues[j], Market[i].Leash[j], Market[i].LeashVol[j], Market[i].Accuracy[j]); // Master, Leash, LeashVolatility, Accuracy
 for (int t=0; t<Time; t++) {
 gsl_matrix_set(Market[i].BiasedValues, j, t, C[t]);
 gsl_matrix_set(TrueBiased, i+1, t, C[t]);
 BiasedDistance+=100*((gsl_matrix_get (GSLTrueValues, j, t)-C[t])/gsl_matrix_get (GSLTrueValues, j, t))/(Time*NumberOfStocks*NumberOfAgents);
 VBiasedDistance.push_back(BiasedDistance);
 };
 C.clear();
 }; // closes j loop
 }; // closes i loop
 cout << ", BiasedDistance=" << BiasedDistance << "%" << endl;
 PlotGSLMatrix(TrueBiased, "TrueBiased.csv", 1);
 sort(VAnnualCounter.begin(), VAnnualCounter.end()); // sort() is by ascending order
 sort(VAnnualAmplitude.begin(), VAnnualAmplitude.end()); // sort() is by ascending order
 sort(VBiasedDistance.begin(), VBiasedDistance.end()); // sort() is by ascending order
 double VarAnnualCounter=0; for (int k=0; k<int(VAnnualCounter.size()); k++) {VarAnnualCounter+=(VAnnualCounter[k]-AnnualCounter)*(VAnnualCounter[k]-AnnualCounter)/int(VAnnualCounter.size());};
 double VarAnnualAmplitude=0; for (int k=0; k<int(VAnnualAmplitude.size()); k++) {VarAnnualAmplitude+=(VAnnualAmplitude[k]-AnnualAmplitude)*(VAnnualAmplitude[k]-AnnualAmplitude)/int(VAnnualAmplitude.size());};
 double VarBiasedDistance=0; for (int k=0; k<int(VBiasedDistance.size()); k++) {VarBiasedDistance+=(VBiasedDistance[k]-BiasedDistance)*(VBiasedDistance[k]-BiasedDistance)/int(VBiasedDistance.size());};
 Res.push_back(AnnualCounter); // Mean
 Res.push_back(AnnualAmplitude); // Mean
 Res.push_back(BiasedDistance); // Mean
 Res.push_back(VAnnualCounter[int(VAnnualCounter.size()/2)]); // Median
 Res.push_back(VAnnualAmplitude[int(VAnnualAmplitude.size()/2)]); // Median
 Res.push_back(VBiasedDistance[int(VBiasedDistance.size()/2)]); // Median
 Res.push_back(sqrt(VarAnnualCounter)); // Stdev
 Res.push_back(sqrt(VarAnnualAmplitude)); // Stdev
 Res.push_back(sqrt(VarBiasedDistance)); // Stdev
 return Res;
 }; // closes CheckTrueGeneration()
 void CheckTrueGenerationAll () {
 gsl_matrix * Res = gsl_matrix_calloc (9, 1); // For j=0 only
 //double Mu[3] = {0.10, 0.10, 0.10}; double Leash[3] = {1, 1, 1};
 double Mu[1] = {0.10}; double Leash[1] = {1};
 for (int k1=0; k1<1; k1++) {
 for (int k2=0; k2<1; k2++) {
 vector<double> V = CheckTrueGeneration (Mu[k1], Leash[k2]);
 gsl_matrix_set(Res, 0, k1*3+k2, V[0]);
 gsl_matrix_set(Res, 1, k1*3+k2, V[1]);
 gsl_matrix_set(Res, 2, k1*3+k2, V[2]);
 gsl_matrix_set(Res, 3, k1*3+k2, V[3]);
 gsl_matrix_set(Res, 4, k1*3+k2, V[4]);
 gsl_matrix_set(Res, 5, k1*3+k2, V[5]);
 gsl_matrix_set(Res, 6, k1*3+k2, V[6]);
 gsl_matrix_set(Res, 7, k1*3+k2, V[7]);
 gsl_matrix_set(Res, 8, k1*3+k2, V[8]);
 }; // closes k2 loop
 }; // closes k1 loop
 PlotGSLMatrix(Res, "Res.csv", 1);
 };
 */









// MAIN PROGRAM - MAIN PROGRAM - MAIN PROGRAM - MAIN PROGRAM - MAIN PROGRAM
// SPECIFY "string Machine" JUSTE BELOW "using namespace std" ACCORDING TO LOCAL PATHWAY
// DOWNLOAD GNU SCIENTIFIC LIBRARY GSL-v2.5 at https://www.gnu.org/software/gsl/

int main () {
    
    

    // CALIBRATION
    // Parameters
    //int HPI=500; // Nb of agents
    //int HPGesture=0; // Agent gesture amplification
    //double HPTrueMu=0.10; // Parameter in the Poisson random walk of the fundamental values
    int HPAccuracy=10; // Probability in bubble burst in the fundamental values
    int HPTime=1453+1000; // Nb of time steps (from 27/09/2018 to 27/09/2022: 1453 steps on crypto (BTC.csv) for 355 steps per trading year, 1010 steps on LSE (UB45.L) for 252 steps per trading year
    //int HPLiquidationFloor=50; // Agent drawdown threshold
    string TypeNEB="Algorithmic"; // 100% bots, 0% "human" agents
    string Leader="NoCluster"; // No herding input
    int ClusterLimit=1; // Threshold of the number of herding agents (1 for "off")
    int pNEB=0;
    double Rate=1; // Interest rate
    int S=20; // Nb of simulation runs



    HPMarketSimulatorAll (500, 1.5, 0.1, 11, 1453+1000, 50, "Algorithmic", "NoCluster", 1, 0, 1, S); // CALIB TO BINANCE 
    


    // Out of 530 coins for Binance: 153 BTC and 20 USDT (stitched), and xxx 57 and 14 USDT (unstitched)
    vector<Coin> PF_coin = readStockData_coin("/Users/johann/Desktop/Work/Postdoc2022/cryptotick/BINANCE_OUTPUT_STITCHED/");
    vector<vector<Coin> > PFTwin = PFTrainingTesting_coin (PF_coin, int(PF_coin.size()/2)); // LOADS FULL LONDON STOCK EXCHANGE RAW DATA IN TRAINING AND TESTING SETS
    JointDistributions (PortfolioMultiOutput_coin (PFTwin[0], 0), 0, 100);
    JointDistributions (PortfolioMultiOutput_coin (PFTwin[1], 0), 0, 100);
    
    
    // FULL HYPERPARAMETER SPACE OPTIMIZATION: 1200 combinations
    /*
    for (int LiquidationFloor=10; LiquidationFloor<=90; LiquidationFloor+=20) {
        for (double TrueMu=0.1; TrueMu<=1.5; TrueMu+=0.2) {
            for (double Gesture=1.0; Gesture<=3.0; Gesture+=0.5) {
                for (int I=500; I<=5500; I+=1000) {
                    // This outputs to "/Users/johann/Documents/SYMBA/JointDistributions.xls"
                    HPMarketSimulatorAll (I, Gesture, TrueMu, HPAccuracy, HPTime, LiquidationFloor, TypeNEB, Leader, ClusterLimit, pNEB, Rate, S);
                };
            };
        };
    };
    
    
    
    return 0;
}

