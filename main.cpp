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

const double alpha = 2;
const double betta = 2*alpha;

const int pointCount = 1024;
const string taskFilePath = "../resource/Barinova.txt";

const size_t rangesCount = 100;
const size_t polynomialMaxDegree = 9; // GnuPlot constraint 38

const string pngPrefix = "Result/";
const string imageFormat = "png";

const double a = 0.95;

double calcAverageValue(const vector<double>& heights);
double calcRootMeanSquareValue(const vector<double>& heights);

vector< double > getAproximateCoefficients(const size_t polynomialDegree, const vector<double>& x, const vector<double>& y);
void GaussSolve(vector< vector< double > >& a, vector< double >& b, vector< double >& x, const size_t n);
double leastSquaresMethod(const vector<double>& polynom, const vector<double>& x, const vector<double>& y);
double calcPolynomValue(const vector<double>& polynom, double xValue);
std::string makePolynom(const vector<double>& polynom);

size_t factorial(size_t n) {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

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

    Gnuplot primaryFilePlot;
    primaryFilePlot.set_title("Primary File");
    primaryFilePlot.set_xlabel("Years");
    primaryFilePlot.set_ylabel("Heights");
    primaryFilePlot.set_style("lines");
    primaryFilePlot.cmd("set terminal " + imageFormat + " size 1024,768");
    primaryFilePlot.cmd("set output \"" + pngPrefix + "Primary File." + imageFormat + "\"");
    primaryFilePlot.plot_xy(years, heights);


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
    distributionFunctionPlot.cmd("set terminal " + imageFormat + " size 1024,768");
    distributionFunctionPlot.cmd("set output \"" + pngPrefix + "Distribution function." + imageFormat + "\"");
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
    densityDistributionPlot.cmd("set terminal " + imageFormat + " size 1024,768");
    densityDistributionPlot.cmd("set output \"" + pngPrefix + "Density Distribution." + imageFormat + "\"");
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

    aproximateDensityDistributionPlot.cmd("set terminal " + imageFormat + " size 1024,768");
    aproximateDensityDistributionPlot.cmd("set output \"" + pngPrefix + "Aproximate Density Distribution." + imageFormat + "\"");
    aproximateDensityDistributionPlot.plot_xy(ranges, densityDistribution);

    aproximateDensityDistributionPlot.cmd("set terminal " + imageFormat + " size 1024,768");
    aproximateDensityDistributionPlot.cmd("set output \"" + pngPrefix + "Aproximate Density Distribution." + imageFormat + "\"");
    aproximateDensityDistributionPlot.set_pointsize(10.0);
    aproximateDensityDistributionPlot.plot_equation(makePolynom(polynom));


    //**** PART II ****\\
    // ---- Stationarity of the series ---- //
    fftw_complex *in, *out;
    fftw_plan p_forward, p_backward;

    in  = ( fftw_complex* ) fftw_malloc(sizeof (fftw_complex ) * pointCount);
    out = ( fftw_complex* ) fftw_malloc(sizeof (fftw_complex ) * pointCount);

    for(size_t i = 0; i < pointCount; i++) {
        in[i][0] = heights[i];
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

    ACFPlot.set_xrange(0, pointCount);
    ACFPlot.cmd("set terminal " + imageFormat + " size 1024,768");
    ACFPlot.cmd("set output \"" + pngPrefix + "ACF." + imageFormat + "\"");

    ACFPlot.plot_x(ACF);

    //**** PART III ****\\
    // ---- Calculate the cumulative frequency ---- //
    vector<double> cumulativeFrequency;
    for(size_t i = 0; i < rangesCount; i++) {
        cumulativeFrequency.push_back(rangeEntering[rangesCount-i-1]);
    }

    Gnuplot cumulativeFrequencyPlot;
    cumulativeFrequencyPlot.set_title("Cumulative Frequency");

    cumulativeFrequencyPlot.set_xrange(0, 1.);
    cumulativeFrequencyPlot.set_style("lines");
    cumulativeFrequencyPlot.cmd("set terminal " + imageFormat + " size 1024,768");
    cumulativeFrequencyPlot.cmd("set output \"" + pngPrefix + "Cumulative Frequency." + imageFormat + "\"");
    cumulativeFrequencyPlot.plot_xy(ranges, cumulativeFrequency);

    // ---- Stationarity of the cumulative frequency ----
    in  = ( fftw_complex* ) fftw_malloc(sizeof (fftw_complex ) * rangesCount);
    out = ( fftw_complex* ) fftw_malloc(sizeof (fftw_complex ) * rangesCount);

    for(size_t i = 0; i < rangesCount; i++) {
        in[i][0] = heights[i];
        in[i][1] = 0;
    }

    p_forward = fftw_plan_dft_1d(rangesCount, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p_forward);

    for(size_t i = 0; i < rangesCount; i++) {
        out[i][0] = pow(out[i][0], 2);
        out[i][1] = pow(out[i][1], 2);
    }

    p_backward = fftw_plan_dft_1d (rangesCount, out, in, FFTW_BACKWARD, FFTW_ESTIMATE );
    fftw_execute(p_backward);

    vector<double> ACF2;
    ACF2.resize(rangesCount);
    for(size_t i = 0; i < rangesCount; i++) {
        ACF2[i] = in[i][0]/in[0][0];
    }

    fftw_destroy_plan(p_forward);
    fftw_destroy_plan(p_backward);
    fftw_free(in);
    fftw_free(out);

    Gnuplot ACF2Plot;
    ACF2Plot.set_title("ACF2");

    ACF2Plot.set_xrange(0, rangesCount);
    ACF2Plot.set_style("lines");
    ACF2Plot.cmd("set terminal " + imageFormat + " size 1024,768");
    ACF2Plot.cmd("set output \"" + pngPrefix + "ACF2." + imageFormat + "\"");

    ACF2Plot.plot_x(ACF2);

    //**** PART IV ****\\
    // ---- Long term prognosis excess ---- //
    Gnuplot excessPlot;
    excessPlot.set_title("Excess");


    double v = cumulativeFrequency[a*rangesCount];
    const double T = 400;
    double e = 2.7;
    for(size_t m = 1; m < 13; m++) {
        string vStr = boost::lexical_cast<string>(v);
        string mStr = boost::lexical_cast<string>(m);
        string eStr = boost::lexical_cast<string>(e);
        string equation = "((" + vStr + "*x)" +
                "**" + mStr + ")/"  + boost::lexical_cast<string>(factorial(m)) +
                "*" + eStr + "**(-" + vStr + "*x)";

        excessPlot.set_style("lines");
        excessPlot.set_xrange(0, T);
        excessPlot.cmd("set terminal " + imageFormat + " size 1024,768");
        excessPlot.cmd("set output \"" + pngPrefix + "Excess." + imageFormat + "\"");
        excessPlot.plot_equation(equation);
    }

    //**** PART V ****\\
    // ---- Short-term forecast of the water level ---- //
    //У меня получилось стационарны 146 точек
    const size_t stationaryPointsCount = 71;
    vector< vector< double > > x(stationaryPointsCount, vector< double >(stationaryPointsCount, 0));
    vector< double > xn(stationaryPointsCount*2, 0);

    std::copy(normalizedHeight.begin() + stationaryPointsCount-1,
              normalizedHeight.begin() + stationaryPointsCount*2,
              xn.begin());
    for(size_t i = 0; i < stationaryPointsCount; i++) {
        for(size_t j = 0; j< stationaryPointsCount; j++) {
            x[i][j] = normalizedHeight[j+i];
        }
    }

    vector< double > result(stationaryPointsCount, 0);
    GaussSolve(x, xn, result, stationaryPointsCount);

    for(size_t j = 0; j< stationaryPointsCount; j++) {
        //cout << result[j] << endl;
    }

    std::copy(normalizedHeight.begin() + stationaryPointsCount,
              normalizedHeight.begin() + stationaryPointsCount*2,
              xn.begin());

    for(size_t i = 0; i < stationaryPointsCount; i++) {
        for(size_t j = 0; j< stationaryPointsCount; j++) {
            xn[i + stationaryPointsCount] += result[j] * normalizedHeight[i + j + stationaryPointsCount*2];
        }
    }


    vector< double > forecast(stationaryPointsCount);
    vector< double > realData(stationaryPointsCount);

    std::copy(xn.begin() + stationaryPointsCount, xn.end(), forecast.begin());
    std::copy(
                normalizedHeight.begin() + stationaryPointsCount*2+2,
                normalizedHeight.begin() + stationaryPointsCount*3+2,
                realData.begin()
             );

    Gnuplot forecastPlot;
    forecastPlot.set_title("Forecast");

    forecastPlot.set_style("lines");
    forecastPlot.cmd("set terminal " + imageFormat + " size 1024,768");
    forecastPlot.cmd("set output \"" + pngPrefix + "Forecast." + imageFormat + "\"");
    forecastPlot.plot_x(forecast);

    forecastPlot.cmd("set terminal " + imageFormat + " size 1024,768");
    forecastPlot.cmd("set output \"" + pngPrefix + "Forecast." + imageFormat + "\"");
    forecastPlot.plot_x(realData);

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
