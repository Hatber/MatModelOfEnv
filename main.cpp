#include <iostream>
#include <vector>

#include <cmath>
#include <algorithm>
#include <numeric>

#include <boost/lexical_cast.hpp>

#include "lib/gnuplot_i.hpp"

#include "include/TaskLoader.h"
#include "lib/fftw/fftw3.h"

using namespace std;

const int pointCount = 1024;
const string taskFilePath = "../resource/Task.txt";

const size_t rangesCount = 100;
const size_t polynomialMaxDegree = 38; // GnuPlot constraint 38

double calcAverageValue(const vector<double>& heights);
double calcRootMeanSquareValue(const vector<double>& heights);

vector< double > getAproximateCoefficients(const size_t polynomialDegree, const vector<double>& x, const vector<double>& y);
void GaussSolve(vector< vector< double > >& a, vector< double >& b, vector< double >& x, const size_t n);
double leastSquaresMethod(const vector<double>& polynom, const vector<double>& x, const vector<double>& y);
double calcPolynomValue(const vector<double>& polynom, double xValue);
std::string makePolynom(const vector<double>& polynom);

int main()
{
    // ---- LOAD TASK ---- //
    ITaskLoader *loader = TaskLoaderFactory::makeLoader(TXT, taskFilePath.c_str());
    if(loader->load(pointCount)) {
        cerr << "Task Load error!" << endl;
        return -1;
    }

    vector<double>& years = loader->getYears();
    vector<double>& heights = loader->getHeights();


    // **** PART I **** \\
    // ---- Сalculation the moments of series ---- //
    double averageValue = calcAverageValue(heights);
    cout << "Average Value: " << averageValue << endl;

    double rootMeanSquareValue = calcRootMeanSquareValue(heights);
    cout << "Root Mean Square Value: " << rootMeanSquareValue << endl;


    // ---- Construction of the distribution function ---- //
    double minHeight, maxHeight;
    minHeight = *std::min_element(heights.begin(), heights.end());
    maxHeight = *std::max_element(heights.begin(), heights.end());
    cout << "Min Value of water height: " << minHeight << endl;
    cout << "Max Value of water height: " << maxHeight << endl;

    vector<double> normalizedHeight(pointCount, 0);
    normalizedHeight.resize(pointCount);
    for(size_t i = 0; i < pointCount; i++) {
        normalizedHeight[i] = (heights[i] - minHeight) / (maxHeight - minHeight);
    }

    vector<double> rangeEntering(rangesCount);
    for(size_t heightIt = 0; heightIt < pointCount; heightIt++) {
        for(size_t range = 0; range < rangesCount; range++) {
            if(normalizedHeight[heightIt] < ((double)1/rangesCount)*(range+1)) {
                rangeEntering[range] += 1./pointCount;
            }
        }
    }

    vector<double> ranges(rangesCount, 0);
    for(size_t range = 0; range < rangesCount-1; range++) {
        ranges[range+1] = ranges[range] + (double)1/rangesCount;
    }

    Gnuplot distributionFunctionPlot;
    distributionFunctionPlot.set_title("Distribution function");
    distributionFunctionPlot.set_xlabel("Heights");
    distributionFunctionPlot.set_style("lines");
    distributionFunctionPlot.plot_xy(ranges, rangeEntering);

    // ---- Construction of the density distribution ---- //
    vector<double> densityDistribution(rangesCount);
    for(size_t range = 0; range < rangesCount-1; range++) {
        densityDistribution[range] = fabs(rangeEntering[range] - rangeEntering[range+1]);
    }

    double minDensity, maxDensity;
    minDensity = *std::min_element(densityDistribution.begin(), densityDistribution.end());
    maxDensity = *std::max_element(densityDistribution.begin(), densityDistribution.end());
    for(size_t range = 0; range < rangesCount; range++) {
        densityDistribution[range] = (densityDistribution[range] - minDensity) / (maxDensity - minDensity);
    }

    Gnuplot densityDistributionPlot;
    densityDistributionPlot.set_title("Normalized Density distribution");
    densityDistributionPlot.set_style("impulses");
    densityDistributionPlot.set_ylabel("f(x)");
    densityDistributionPlot.set_xlabel("Heights");
    densityDistributionPlot.plot_xy(ranges, densityDistribution);

    // ---- Approximation of a theoretical dependences ---- //
    vector<double> polynom, deviations;
    for(size_t degree = 1; degree < polynomialMaxDegree; degree++) {
        polynom = getAproximateCoefficients(degree, ranges, densityDistribution);
        deviations.push_back(leastSquaresMethod(polynom, ranges, densityDistribution));
    }

    size_t bestDegreeDeviation = 0;
    double bestDeviation = deviations.front();
    for(size_t i = 0; i < deviations.size(); i++) {
        if(bestDeviation > deviations[i]) {
            bestDeviation = deviations[i];
            bestDegreeDeviation = i;
        }
    }

    polynom = getAproximateCoefficients(bestDegreeDeviation, ranges, densityDistribution);
    cout << "Best Deviation: " << bestDegreeDeviation << " polynom degree";

    Gnuplot aproximateDensityDistributionPlot;
    aproximateDensityDistributionPlot.set_title("Aproximate Density Distribution");
    aproximateDensityDistributionPlot.set_style("lines");
    aproximateDensityDistributionPlot.set_ylabel("f(x)");
    aproximateDensityDistributionPlot.set_xlabel("Heights");
    aproximateDensityDistributionPlot.set_xrange(0., 0.98);

    aproximateDensityDistributionPlot.reset_plot();
    aproximateDensityDistributionPlot.plot_xy(ranges, densityDistribution);
    aproximateDensityDistributionPlot.set_pointsize(10.0);
    aproximateDensityDistributionPlot.plot_equation(makePolynom(polynom));


    //**** PART II ****\\
    // ---- Stationarity of the series ---- //
    fftw_complex *in, *out;
    fftw_plan p_forward, p_backward;

    in = ( fftw_complex* ) fftw_malloc(sizeof (fftw_complex ) * pointCount);
    out = ( fftw_complex* ) fftw_malloc(sizeof (fftw_complex ) * pointCount);

    for(size_t i = 0; i < pointCount; i++) {
        in[i][0] = normalizedHeight[i];
        in[i][1] = 0;
    }

    p_forward = fftw_plan_dft_1d(pointCount, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p_forward);

    for(size_t i = 0; i < pointCount; i++) {
        out[i][0] = pow(out[i][0], 2);
        out[i][1] = pow(out[i][1], 2);
    }

    p_backward = fftw_plan_dft_1d (pointCount, out, in, FFTW_BACKWARD, FFTW_ESTIMATE );
    fftw_execute(p_backward);

    vector<double> ACF;
    ACF.resize(pointCount);
    for(size_t i = 0; i < pointCount; i++) {
        ACF[i] = in[i][0]/in[0][0];
    }

    fftw_destroy_plan(p_forward);
    fftw_destroy_plan(p_backward);
    fftw_free(in);
    fftw_free(out);

    Gnuplot ACFPlot;
    ACFPlot.set_title("ACF");
    ACFPlot.plot_x(ACF);


    

    cin.get();

    return 0;
}


double calcAverageValue(const vector<double>& heights) {
    return  std::accumulate(heights.begin(), heights.end(), 0.) / heights.size();
}

double calcRootMeanSquareValue(const vector<double>& heights) {
    double averageValue = calcAverageValue(heights);

    double summ = 0;
    size_t heightCount = heights.size();

    for(size_t heightIt = 0; heightIt < heightCount; heightIt++) {
        summ += pow(heights[heightIt] - averageValue, 2.);
    }

    return summ/heightCount;
}

vector< double > getAproximateCoefficients(const size_t polynomialDegree, const vector<double>& x, const vector<double>& y) {
    vector< vector< double > > equationCoefficients;
    vector<double> freeMembers;

    freeMembers.resize(polynomialDegree);
    equationCoefficients.resize(polynomialDegree);
    for(size_t i = 0; i< polynomialDegree; i++) {
        equationCoefficients[i].resize(polynomialDegree);
    }

    size_t pointCount = x.size();
    for(size_t i = 0; i < polynomialDegree; i++) {
        for(size_t j = 0; j < polynomialDegree; j++) {
            for(size_t pointIt = 0; pointIt < pointCount; pointIt++) {
                equationCoefficients[i][j] += pow(x[pointIt], j + i);
            }
        }
        for(size_t pointIt = 0; pointIt < pointCount; pointIt++) {
            freeMembers[i] += pow(x[pointIt], i)*y[pointIt];
        }
    }

    vector< double > result;
    result.resize(polynomialDegree);

    GaussSolve(equationCoefficients, freeMembers, result, polynomialDegree);

    return result;
}

void GaussSolve(vector< vector< double > >& a, vector< double >& b, vector< double >& x, const size_t n) {
    double tmpValue;
    int i,j,k,z;
    double dblLeadElement;

    for(i=0; i<n; i++) {
        dblLeadElement=a[i][i];

        double tmpMax = dblLeadElement;
        int tmpMaxNumber = i;

        for (z=i; z<n; z++) {
            if (a[z][i]>tmpMax) {
                tmpMax = a[z][i];
                tmpMaxNumber=z;
            }
        }

        for (z=i; z<n; z++) {
            tmpValue = a[i][z];
            a[i][z] = a[tmpMaxNumber][z];
            a[tmpMaxNumber][z] = tmpValue;
        }
        tmpValue = b[i];
        b[i] = b[tmpMaxNumber];
        b[tmpMaxNumber] = tmpValue;

        dblLeadElement = tmpMax;

        for(j=i; j<n; j++) {
            a[i][j]/=dblLeadElement;
        }
        b[i]/=dblLeadElement;

        for(k=i+1; k<n; k++) {
            double dblToDivide=a[k][i]/a[i][i];
            for(z=i; z<n; z++) {
                a[k][z]-=a[i][z]*dblToDivide;
            }
            b[k]-=b[i]*dblToDivide;
        }
    }

    x[n-1]=b[n-1];

    for(k=n-2; k>=0; k--) {
        double sum=b[k];

        for(j=k+1; j<n; j++) {
            sum-=a[k][j]*x[j];
        }
        x[k]=sum;
    }
}

double leastSquaresMethod(const vector<double>& polynom, const vector<double>& x, const vector<double>& y) {
    double deviation = 0;
    size_t pointCount = x.size();

    for(size_t i = 0; i < pointCount; i++) {
        deviation += powl(y[i] - calcPolynomValue(polynom, x[i]), 2);
    }

    return deviation;
}

double calcPolynomValue(const vector<double>& polynom, double xValue) {
    double value = 0;
    size_t polynomSize = polynom.size();

    for(size_t i = 0; i < polynomSize; i++) {
        value += powl(xValue, i) * polynom[i];
    }

    return value;
}

std::string makePolynom(const vector<double>& polynom) {
    std::string result = "";
    size_t polynomSize = polynom.size();

    for(size_t i = 0; i < polynomSize; i++) {
        result += boost::lexical_cast<string>(polynom[i]) + "*(x**" + boost::lexical_cast<string>(i) + ")";
        if(i!=polynomSize-1) { result+="+"; }
    }

    return result;
}
