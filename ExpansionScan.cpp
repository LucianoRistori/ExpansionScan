/*
==============================================================================
File:        ExpansionScan.cpp
Version:     1.3.0
Author:      Luciano Ristori

Description:
------------
ExpansionScan analyzes thermal expansion effects by comparing two sets of
measured 3D points (e.g. warm vs cold scans) acquired with a CMM.

For each matched point, the program:
  • Computes raw displacements (warm → cold)
  • Fits a global expansion model
  • Computes residual displacements after the fit
  • Produces ROOT plots (histograms and XY scatter plots)
  • Optionally writes a detailed CSV summary

The program is designed to be run from a dataset-specific *results directory*.

Input and output path handling:
-------------------------------
• Input files:
    - Absolute paths are used as given
    - Relative paths are interpreted relative to $CMM_DATA if set,
      otherwise relative to the current working directory

• Output files (ROOT, PNG, CSV):
    - Absolute paths are used as given
    - Relative paths are interpreted relative to $CMM_RESULTS/ExpansionScan
      if CMM_RESULTS is set
    - Otherwise, relative paths refer to the current working directory

Command-line usage:
-------------------
ExpansionScan <warm.csv> <cold.csv> [output.root] [options]

Positional arguments:
  warm.csv        CSV file with warm measurements
  cold.csv        CSV file with cold measurements
  output.root     Optional ROOT output file
                  (default: ExpansionScan.root)

Options:
  --labels            Draw point labels on XY scatter plots
  --dT <value>        Specify temperature difference (ΔT)
  --csv <file.csv>    Write a CSV summary of per-point results

CSV output:
-----------
When --csv is specified, a CSV file is written with one row per point,
containing:

  label,
  warm_x_mm, warm_y_mm, warm_z_mm,
  cold_x_mm, cold_y_mm, cold_z_mm,
  dx_before_um, dy_before_um, dr_before_um,
  dx_after_um,  dy_after_um,  dr_after_um

Formatting conventions:
  • Coordinates in millimeters are written with 3 decimal digits
  • Displacements in micrometers are written as integers

Notes:
------
• Plotting behavior and --labels semantics are unchanged from previous versions
• CSV generation does not alter any analysis or plotting logic
• Executables and build artifacts are not tracked in git

==============================================================================
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <tuple>
#include <cstdio>

#include "Points.h"

// ROOT includes
#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h" 
#include "TMarker.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TColor.h"
#include "TROOT.h"
#include "TPad.h"
#include "TPaveStats.h"


using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;

#include <cstdlib>
#include <filesystem>

///////////////////////////////////////////////
// Helper functions                          //
///////////////////////////////////////////////
//
// Resolve input path:
//  - absolute → unchanged
//  - relative → prepend CMM_DATA if set, else unchanged
//
static std::string resolveInputPath(const std::string& path)
{
    if (path.empty()) return path;

    std::filesystem::path p(path);
    if (p.is_absolute()) return path;

    const char* base = std::getenv("CMM_DATA");
    if (base && *base) {
        return (std::filesystem::path(base) / p).string();
    }

    return path;
}

// Resolve output path:
static std::string resolveOutputPath(const std::string& path)
{
    if (path.empty()) return path;

    std::filesystem::path p(path);
    if (p.is_absolute()) return path;

    const char* base = std::getenv("CMM_RESULTS");
    if (base && *base) {
        return (std::filesystem::path(base) / "ExpansionScan" / p).string();
    }

    // No CMM_RESULTS: write to current working directory
    return path;
}


//------------------------------------------------------------------------------
// Compile-time constant: maximum arrow length (in mm) on the plots
//   (Option B: relatively large arrows → 500 mm max normalized length.)
//------------------------------------------------------------------------------
constexpr double MAX_ARROW_LEN = 100.0;

//------------------------------------------------------------------------------
// Simple struct to hold displacement info
//------------------------------------------------------------------------------
struct Disp {
    double dx;
    double dy;
    double r;
};

//------------------------------------------------------------------------------
// Simple statistics container
//------------------------------------------------------------------------------
struct Stats {
    double mean  = 0.0;
    double sigma = 0.0;  // standard deviation
    size_t n     = 0;
};

//------------------------------------------------------------------------------
// Compute mean and sigma for a vector of values
//------------------------------------------------------------------------------
Stats computeStats(const vector<double>& v)
{
    Stats s;
    s.n = v.size();
    if (s.n == 0) return s;

    double sum = 0.0;
    double sum2 = 0.0;
    for (double x : v) {
        sum  += x;
        sum2 += x * x;
    }
    s.mean = sum / s.n;
    double var = sum2 / s.n - s.mean * s.mean;
    if (var < 0.0) var = 0.0;
    s.sigma = std::sqrt(var);
    return s;
}

//------------------------------------------------------------------------------
// Compute similarity transform (translation + uniform scale + rotation)
// mapping (xw,yw) → (xc,yc).
//
// The model is:
//   [xc_i]   [Cx]   + s * R(theta) * ([xw_i] - [cw])
//   [yc_i] = [Cy]         where cw=(cx,cy), R is 2×2 rotation matrix,
//   and cx,cy (resp. Cx,Cy) are the centroids of warm (resp. cold) points.
//
// On success, returns true and fills (tx,ty,s,theta):
//   (tx,ty) is the translation to apply AFTER rotation and scaling, i.e.
//     [xw’]   [tx]   + s * R(theta) * [xw]
//     [yw’] = [ty]                   [yw]
//
// so that (xw’,yw’) best matches (xc,yc) in least-squares sense.
//
// We implement the classic 2D Procrustes/Umeyama-style solution.
//------------------------------------------------------------------------------
bool fitSimilarity2D(const vector<double>& xw,
                     const vector<double>& yw,
                     const vector<double>& xc,
                     const vector<double>& yc,
                     double& tx, double& ty,
                     double& s,  double& theta)
{
    size_t n = xw.size();
    if (n == 0 || xc.size() != n || yc.size() != n || yw.size() != n) {
        return false;
    }

    // centroids
    double cxw = 0.0, cyw = 0.0;
    double cxc = 0.0, cyc = 0.0;
    for (size_t i = 0; i < n; ++i) {
        cxw += xw[i];
        cyw += yw[i];
        cxc += xc[i];
        cyc += yc[i];
    }
    cxw /= static_cast<double>(n);
    cyw /= static_cast<double>(n);
    cxc /= static_cast<double>(n);
    cyc /= static_cast<double>(n);

    // centered coordinates and scatter terms
    double Sxx  = 0.0;
    double Sxy  = 0.0;
    double Sx_y = 0.0;

    for (size_t i = 0; i < n; ++i) {
        double Xw = xw[i] - cxw;
        double Yw = yw[i] - cyw;
        double Xc = xc[i] - cxc;
        double Yc = yc[i] - cyc;

        Sxx  += Xw * Xw + Yw * Yw;
        Sxy  += Xw * Xc + Yw * Yc;
        Sx_y += Xw * Yc - Yw * Xc;  // cross term for rotation
    }

    if (Sxx <= 0.0) {
        return false;
    }

    // scale
    double denom = Sxx;
    double norm  = std::sqrt(Sxy * Sxy + Sx_y * Sx_y);
    if (norm <= 0.0) {
        // purely translational? set scale=1 and zero rotation
        s     = 1.0;
        theta = 0.0;
    } else {
        s = norm / denom;
        double cos_th = Sxy  / (s * denom);
        double sin_th = Sx_y / (s * denom);
        // guard numerically
        if (cos_th >  1.0) cos_th =  1.0;
        if (cos_th < -1.0) cos_th = -1.0;
        theta = std::atan2(sin_th, cos_th);
    }

    double cos_t = std::cos(theta);
    double sin_t = std::sin(theta);

    // translation from centroid relation
    double tx_cent = cxc - s * (cos_t * cxw - sin_t * cyw);
    double ty_cent = cyc - s * (sin_t * cxw + cos_t * cyw);

    tx = tx_cent;
    ty = ty_cent;

    return true;
}

//------------------------------------------------------------------------------
// Compute displacements between two sets of points
//------------------------------------------------------------------------------
void computeDisplacements(const vector<double>& x1,
                          const vector<double>& y1,
                          const vector<double>& x2,
                          const vector<double>& y2,
                          vector<Disp>& disp)
{
    size_t n = std::min({x1.size(), y1.size(), x2.size(), y2.size()});
    disp.resize(n);
    for (size_t i = 0; i < n; ++i) {
        double dx = x2[i] - x1[i];
        double dy = y2[i] - y1[i];
        disp[i].dx = dx;
        disp[i].dy = dy;
        disp[i].r  = std::sqrt(dx * dx + dy * dy);
    }
}

//------------------------------------------------------------------------------
// Normalize arrow lengths so that max |(dx,dy)| → MAX_ARROW_LEN
//------------------------------------------------------------------------------
void normalizeArrows(vector<Disp>& disp)
{
    double maxR = 0.0;
    for (const auto& d : disp) {
        if (d.r > maxR) maxR = d.r;
    }
    if (maxR <= 0.0) return;

    double factor = MAX_ARROW_LEN / maxR;
    for (auto& d : disp) {
        d.dx *= factor;
        d.dy *= factor;
        d.r  *= factor; // not strictly needed for drawing, but consistent
    }
}

//------------------------------------------------------------------------------
// Adjust StatBox size and font for better readability
//------------------------------------------------------------------------------
void tuneStatBox(TH1* h,
                 double x1 = 0.60, double y1 = 0.60,
                 double x2 = 0.90, double y2 = 0.90,
                 double textSize = 0.035)
{
    if (!h) return;

    TPaveStats* st =
        dynamic_cast<TPaveStats*>(h->FindObject("stats"));
    if (!st) return;

    st->SetX1NDC(x1);
    st->SetY1NDC(y1);
    st->SetX2NDC(x2);
    st->SetY2NDC(y2);
    st->SetTextSize(textSize);
}



//------------------------------------------------------------------------------
// Main
//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    bool drawLabels = false;
    bool haveDeltaT = false;
    double deltaT   = 0.0;

   string warmFile, coldFile, outRoot("ExpansionScan.root"), csvFile;

    if (argc < 3) {
        cerr << "Usage:\n"
             << "  " << argv[0]
             << " [--labels] [--dT value] warm.txt cold.txt [out.root]\n";
        return 1;
    }

	// Parse arguments
	vector<string> positional;
	for (int i = 1; i < argc; ++i) {
		string arg = argv[i];
	
		if (arg == "--labels") {
			drawLabels = true;
		}
		else if (arg == "--dT") {
			if (i + 1 >= argc) {
				cerr << "Error: --dT requires a numeric value.\n";
				return 1;
			}
			try {
				deltaT = std::stod(argv[++i]);
				haveDeltaT = true;
			} catch (...) {
				cerr << "Error: invalid value for --dT.\n";
				return 1;
			}
		}
		
		
		else if (arg == "--csv") {
			if (i + 1 >= argc) {
				cerr << "Error: --csv requires a filename.\n";
				return 1;
			}
			csvFile = argv[++i];
		}
	
		else if (!arg.empty() && arg[0] == '-') {
			cerr << "Error: unknown option '" << arg << "'.\n";
			return 1;
		}
		else {
			positional.push_back(arg);
		}
	}


    if (positional.size() < 2) {
        cerr << "Error: need warm and cold files.\n";
        return 1;
    }

    warmFile = positional[0];
    coldFile = positional[1];

    if (positional.size() >= 3) {
        outRoot = positional[2];
    }
      
    warmFile = resolveInputPath(warmFile);
	coldFile = resolveInputPath(coldFile);

	outRoot = resolveOutputPath(outRoot);

	if (!csvFile.empty()) {
    csvFile = resolveOutputPath(csvFile);
}
  

    cout << "\n====================================\n";
    cout << " ExpansionScan v0.4\n";
    cout << " Built: " << __DATE__ << " " << __TIME__ << "\n";
    cout << "====================================\n";
    cout << "Warm file : " << warmFile << "\n";
    cout << "Cold file : " << coldFile << "\n";
    cout << "Out ROOT  : " << outRoot << "\n";
    cout << "Labels    : " << (drawLabels ? "ON" : "OFF") << "\n";
    if (haveDeltaT) {
        cout << "dT        : " << deltaT << "\n";
    }

    // Read points (label + X,Y,Z; Z ignored)
    // We ask for 3 coordinates, but only the first two are used.
    vector<Point> warmPts = readPoints(warmFile, 3);
    vector<Point> coldPts = readPoints(coldFile, 3);

    if (warmPts.empty() || coldPts.empty()) {
        cerr << "Error: one of the input files yielded no points.\n";
        return 1;
    }

    size_t n = std::min(warmPts.size(), coldPts.size());
    if (warmPts.size() != coldPts.size()) {
        cerr << "Warning: different point counts (warm=" << warmPts.size()
             << ", cold=" << coldPts.size() << "); using first " << n
             << " points.\n";
    }

    vector<double> xw(n), yw(n), xc(n), yc(n);
    for (size_t i = 0; i < n; ++i) {
        xw[i] = warmPts[i].coords[0];
        yw[i] = warmPts[i].coords[1];
        xc[i] = coldPts[i].coords[0];
        yc[i] = coldPts[i].coords[1];
    }

    cout << "Using " << n << " common points.\n";

    // BEFORE-fit displacements
    vector<Disp> dispBefore;
    computeDisplacements(xw, yw, xc, yc, dispBefore);

    // Fit similarity transform
    double tx = 0.0, ty = 0.0, s = 1.0, theta = 0.0;
    bool okFit = fitSimilarity2D(xw, yw, xc, yc, tx, ty, s, theta);

    if (!okFit) {
        cerr << "Error: similarity fit failed.\n";
        return 1;
    }

    double theta_deg = theta * 180.0 / M_PI;

    cout << std::fixed << std::setprecision(6);
    cout << "\nBest-fit similarity transform (warm → cold):\n";
    cout << "  scale s      = " << s << "\n";
    cout << "  rotation θ   = " << theta_deg << " deg\n";
    cout << "  translation  = (" << tx << ", " << ty << ") mm\n";

    if (haveDeltaT && std::fabs(deltaT) > 0.0) {
        double alpha = (s - 1.0) / deltaT;
        cout << "  alpha ≈ (s - 1)/dT = " << alpha << " per degree\n";
    }

    // Apply transform to warm points
    vector<double> xwFit(n), ywFit(n);
    double cos_t = std::cos(theta);
    double sin_t = std::sin(theta);

    for (size_t i = 0; i < n; ++i) {
        double x = xw[i];
        double y = yw[i];
        double xr = cos_t * x - sin_t * y;
        double yr = sin_t * x + cos_t * y;
        xwFit[i] = tx + s * xr;
        ywFit[i] = ty + s * yr;
    }

    // AFTER-fit displacements
    vector<Disp> dispAfter;
    computeDisplacements(xwFit, ywFit, xc, yc, dispAfter);

    //------------------------------------------------------------------------
    // Textual summary of displacements
    //------------------------------------------------------------------------
    vector<double> dxB(n), dyB(n), rB(n);
    vector<double> dxA(n), dyA(n), rA(n);
    for (size_t i = 0; i < n; ++i) {
        dxB[i] = dispBefore[i].dx;
        dyB[i] = dispBefore[i].dy;
        rB[i]  = dispBefore[i].r;

        dxA[i] = dispAfter[i].dx;
        dyA[i] = dispAfter[i].dy;
        rA[i]  = dispAfter[i].r;
    }

    Stats sDxB = computeStats(dxB);
    Stats sDyB = computeStats(dyB);
    Stats sRB  = computeStats(rB);

    Stats sDxA = computeStats(dxA);
    Stats sDyA = computeStats(dyA);
    Stats sRA  = computeStats(rA);

    auto computeRMS = [](const vector<double>& v) {
        if (v.empty()) return 0.0;
        double sum2 = 0.0;
        for (double x : v) sum2 += x * x;
        return std::sqrt(sum2 / static_cast<double>(v.size()));
    };

    double rmsBefore = computeRMS(rB);
    double rmsAfter  = computeRMS(rA);

    cout << "\nDisplacement summary (mm):\n";
    cout << "  BEFORE fit (N = " << n << "):\n";
    cout << "    dX: mean = " << sDxB.mean << "  sigma = " << sDxB.sigma << "\n";
    cout << "    dY: mean = " << sDyB.mean << "  sigma = " << sDyB.sigma << "\n";
    cout << "     R: mean = " << sRB.mean  << "  sigma = " << sRB.sigma
         << "  RMS = " << rmsBefore << "\n";

    cout << "  AFTER fit (N = " << n << "):\n";
    cout << "    dX: mean = " << sDxA.mean << "  sigma = " << sDxA.sigma << "\n";
    cout << "    dY: mean = " << sDyA.mean << "  sigma = " << sDyA.sigma << "\n";
    cout << "     R: mean = " << sRA.mean  << "  sigma = " << sRA.sigma
         << "  RMS = " << rmsAfter << "\n";
			 
	// --- CSV output (additive, no side effects) ---
	//
	// Assumptions (already true in your code):
	//   - vectors warmPts and coldPts exist and are matched 1:1
	//   - dispBefore[i] and dispAfter[i] are already fully computed
	//   - Point has fields: label (string), coords[0]=x, coords[1]=y, coords[2]=z
	//   - Units: coords in mm, displacements in mm
	//
	// This block must be placed AFTER dispBefore / dispAfter are computed,
	// and BEFORE any TApplication::Run() or ROOT event loop.
	
	if (!csvFile.empty()) {
		std::ofstream csv(csvFile);
		if (!csv) {
			std::cerr << "Error: cannot open CSV file '" << csvFile << "'\n";
			return 1;
		}
		csv << std::fixed << std::setprecision(3);
	
		// Header
		csv << "label,"
			<< "x_warm_mm,y_warm_mm,z_warm_mm,"
			<< "x_cold_mm,y_cold_mm,z_cold_mm,"
			<< "dx_before_um,dy_before_um,r_before_um,"
			<< "dx_after_um,dy_after_um,r_after_um\n";
	
		// Rows
		for (size_t i = 0; i < warmPts.size(); ++i) {
			const auto& w = warmPts[i];
			const auto& c = coldPts[i];
			csv << w.label << ","
			<< w.coords[0] << "," << w.coords[1] << "," << w.coords[2] << ","
			<< c.coords[0] << "," << c.coords[1] << "," << c.coords[2] << ","
			<< static_cast<long>(std::lround(dispBefore[i].dx * 1000.0)) << ","
			<< static_cast<long>(std::lround(dispBefore[i].dy * 1000.0)) << ","
			<< static_cast<long>(std::lround(dispBefore[i].r  * 1000.0)) << ","
			<< static_cast<long>(std::lround(dispAfter[i].dx  * 1000.0)) << ","
			<< static_cast<long>(std::lround(dispAfter[i].dy  * 1000.0)) << ","
			<< static_cast<long>(std::lround(dispAfter[i].r   * 1000.0)) << "\n";
		}
	
		csv.close();
	}


    // Prepare ROOT
    int fakeArgc = 0;
    char* fakeArgv[] = { (char*)"app", nullptr };
    TApplication app("app", &fakeArgc, fakeArgv);

    gStyle->SetOptStat(1110);   // show Entries, Mean, StdDev
    gStyle->SetStatFontSize(0.05);   // default ~0.025
	gStyle->SetStatFormat("6.3f");

    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(64);

    TFile outF(outRoot.c_str(), "RECREATE");

    //======================================================================
    // Histograms: BEFORE and AFTER
    //======================================================================
    auto [minDxB_it, maxDxB_it] = std::minmax_element(
        dispBefore.begin(), dispBefore.end(),
        [](const Disp& a, const Disp& b){ return a.dx < b.dx; });
    auto [minDyB_it, maxDyB_it] = std::minmax_element(
        dispBefore.begin(), dispBefore.end(),
        [](const Disp& a, const Disp& b){ return a.dy < b.dy; });
    auto [minRB_it,  maxRB_it ] = std::minmax_element(
        dispBefore.begin(), dispBefore.end(),
        [](const Disp& a, const Disp& b){ return a.r  < b.r;  });

    auto [minDxA_it, maxDxA_it] = std::minmax_element(
        dispAfter.begin(), dispAfter.end(),
        [](const Disp& a, const Disp& b){ return a.dx < b.dx; });
    auto [minDyA_it, maxDyA_it] = std::minmax_element(
        dispAfter.begin(), dispAfter.end(),
        [](const Disp& a, const Disp& b){ return a.dy < b.dy; });
    auto [minRA_it,  maxRA_it ] = std::minmax_element(
        dispAfter.begin(), dispAfter.end(),
        [](const Disp& a, const Disp& b){ return a.r  < b.r;  });

	// Fixed histogram bin width: 1 µm = 0.001 mm
	constexpr double BIN_WIDTH_MM = 0.001;

	// Helper: nice range + fixed physical bin width
	auto niceRangeFixedBin = [](double mn, double mx, double binw)
	{
		if (mx <= mn) mx = mn + binw;

		// Add small margins
		double span = mx - mn;
		mn -= 0.05 * span;
		mx += 0.05 * span;

		// Snap to bin grid so bin CENTERS are exact multiples of binw
		double firstCenter = std::floor(mn / binw) * binw;
		double lastCenter  = std::ceil (mx / binw) * binw;

		double xmin = firstCenter - 0.5 * binw;
		double xmax = lastCenter  + 0.5 * binw;

		int nBins = static_cast<int>(
			std::round((xmax - xmin) / binw)
		);

		return std::make_tuple(xmin, xmax, nBins);
	};

	// Compute ranges and bin counts (BEFORE)
	double dxB_min, dxB_max; int nBinsDxB;
	double dyB_min, dyB_max; int nBinsDyB;
	double rB_min,  rB_max;  int nBinsRB;

	std::tie(dxB_min, dxB_max, nBinsDxB) =
		niceRangeFixedBin(minDxB_it->dx, maxDxB_it->dx, BIN_WIDTH_MM);
	std::tie(dyB_min, dyB_max, nBinsDyB) =
		niceRangeFixedBin(minDyB_it->dy, maxDyB_it->dy, BIN_WIDTH_MM);
	std::tie(rB_min, rB_max, nBinsRB) =
		niceRangeFixedBin(minRB_it->r, maxRB_it->r, BIN_WIDTH_MM);

	// Compute ranges and bin counts (AFTER)
	double dxA_min, dxA_max; int nBinsDxA;
	double dyA_min, dyA_max; int nBinsDyA;
	double rA_min,  rA_max;  int nBinsRA;

	std::tie(dxA_min, dxA_max, nBinsDxA) =
		niceRangeFixedBin(minDxA_it->dx, maxDxA_it->dx, BIN_WIDTH_MM);
	std::tie(dyA_min, dyA_max, nBinsDyA) =
		niceRangeFixedBin(minDyA_it->dy, maxDyA_it->dy, BIN_WIDTH_MM);
	std::tie(rA_min, rA_max, nBinsRA) =
		niceRangeFixedBin(minRA_it->r, maxRA_it->r, BIN_WIDTH_MM);

	// Histogram declarations (Fill loop stays exactly as in v0.5)
	TH1D* hDxBefore = new TH1D("hDxBefore","dX before fit;dX [mm];Counts",
							   nBinsDxB, dxB_min, dxB_max);
	TH1D* hDyBefore = new TH1D("hDyBefore","dY before fit;dY [mm];Counts",
							   nBinsDyB, dyB_min, dyB_max);
	TH1D* hRBefore  = new TH1D("hRBefore","R before fit;R [mm];Counts",
							   nBinsRB, rB_min, rB_max);

	TH1D* hDxAfter  = new TH1D("hDxAfter","dX after fit;dX [mm];Counts",
							   nBinsDxA, dxA_min, dxA_max);
	TH1D* hDyAfter  = new TH1D("hDyAfter","dY after fit;dY [mm];Counts",
							   nBinsDyA, dyA_min, dyA_max);
	TH1D* hRAfter   = new TH1D("hRAfter","R after fit;R [mm];Counts",
							   nBinsRA, rA_min, rA_max);
	
	double xfontSize = 0.030;						   
	hDxBefore->GetXaxis()->SetLabelSize(xfontSize);
	hDyBefore->GetXaxis()->SetLabelSize(xfontSize);
	hRBefore ->GetXaxis()->SetLabelSize(xfontSize);

	hDxAfter ->GetXaxis()->SetLabelSize(xfontSize);
	hDyAfter ->GetXaxis()->SetLabelSize(xfontSize);
	hRAfter  ->GetXaxis()->SetLabelSize(xfontSize);

    for (size_t i = 0; i < n; ++i) {
        hDxBefore->Fill(dispBefore[i].dx);
        hDyBefore->Fill(dispBefore[i].dy);
        hRBefore->Fill (dispBefore[i].r);

        hDxAfter->Fill (dispAfter[i].dx);
        hDyAfter->Fill (dispAfter[i].dy);
        hRAfter->Fill  (dispAfter[i].r);
    }
   

    // Canvas for histograms: 2 rows (before, after) × 3 columns (dX, dY, R)
    TCanvas* cH = new TCanvas("cHists", "ExpansionScan: displacements", 1200, 800);
    cH->Divide(3, 2, 0.01, 0.01); // small gaps between pads

	cH->cd(1);
	hDxBefore->SetLineColor(kBlue+2);
	hDxBefore->Draw("HIST");
	gPad->Update();
	tuneStatBox(hDxBefore, 0.60, 0.75, 0.90, 0.90);
	gPad->Modified();
	gPad->Update();
	
	cH->cd(2);
	hDyBefore->SetLineColor(kBlue+2);
	hDyBefore->Draw("HIST");
	gPad->Update();
	tuneStatBox(hDyBefore, 0.60, 0.75, 0.90, 0.90);
	gPad->Modified();
	gPad->Update();

	cH->cd(3);
	hRBefore->SetLineColor(kBlue+2);
	hRBefore->Draw("HIST");
	gPad->Update();
	tuneStatBox(hRBefore, 0.60, 0.75, 0.90, 0.90);
	gPad->Modified();
	gPad->Update();
	
	cH->cd(4);
	hDxAfter->SetLineColor(kBlue+2);
	hDxAfter->Draw("HIST");
	gPad->Update();
	tuneStatBox(hDxAfter, 0.60, 0.75, 0.90, 0.90);
	gPad->Modified();
	gPad->Update();
	
	cH->cd(5);
	hDyAfter->SetLineColor(kBlue+2);
	hDyAfter->Draw("HIST");
	gPad->Update();
	tuneStatBox(hDyAfter, 0.60, 0.75, 0.90, 0.90);
	gPad->Modified();
	gPad->Update();
	
	cH->cd(6);
	hRAfter->SetLineColor(kBlue+2);
	hRAfter->Draw("HIST");
	gPad->Update();
	tuneStatBox(hRAfter, 0.60, 0.75, 0.90, 0.90);
	gPad->Modified();
	gPad->Update();

    cH->Write("ExpansionScanHists");

    //======================================================================
    // Scatter plots with arrows
    //======================================================================
    // Normalize arrows separately for BEFORE and AFTER
    vector<Disp> arrowsBefore = dispBefore;
    vector<Disp> arrowsAfter  = dispAfter;
    normalizeArrows(arrowsBefore);
    normalizeArrows(arrowsAfter);

    // Determine XY ranges from warm positions
    double xmin = xw[0], xmax = xw[0];
    double ymin = yw[0], ymax = yw[0];
    for (size_t i = 1; i < n; ++i) {
        if (xw[i] < xmin) xmin = xw[i];
        if (xw[i] > xmax) xmax = xw[i];
        if (yw[i] < ymin) ymin = yw[i];
        if (yw[i] > ymax) ymax = yw[i];
    }
    double dxRange = xmax - xmin;
    double dyRange = ymax - ymin;
    if (dxRange <= 0) dxRange = 1.0;
    if (dyRange <= 0) dyRange = 1.0;
    xmin -= 0.05 * dxRange; xmax += 0.05 * dxRange;
    ymin -= 0.05 * dyRange; ymax += 0.05 * dyRange;

    // COMMON color scale for BEFORE and AFTER (so they are directly comparable)
    double Rb_min = rB_min, Rb_max = rB_max;
    double Ra_min = rA_min, Ra_max = rA_max;
    double Rmin_global = std::min(Rb_min, Ra_min);
    double Rmax_global = std::max(Rb_max, Ra_max);
    if (Rmax_global <= Rmin_global) {
        Rmax_global = Rmin_global + 1e-6;
    }

    int nCol = gStyle->GetNumberOfColors();

    TCanvas* cS = new TCanvas("cScat", "ExpansionScan: scatter + arrows", 1400, 650);
    // Add some spacing between pads so they don't visually overlap
    cS->Divide(2, 1, 0.03, 0.03);

    // BEFORE
    cS->cd(1);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.06);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.12);

    TH2D* frameB = new TH2D("frameB", "Before fit;X [mm];Y [mm]",
                            100, xmin, xmax, 100, ymin, ymax);
    frameB->SetStats(false);
    frameB->Draw("AXIS");

    {
        TLatex t;
        t.SetNDC(true);
        t.SetTextSize(0.045);
        t.SetTextAlign(13); // left-top
        t.DrawLatex(0.16, 0.92, "Before fit");
    }

    for (size_t i = 0; i < n; ++i) {
        double xx = xw[i];
        double yy = yw[i];
        double R  = dispBefore[i].r;

        double norm = (R - Rmin_global) / (Rmax_global - Rmin_global);
        if (norm < 0.0) norm = 0.0;
        if (norm > 1.0) norm = 1.0;
        int ci = gStyle->GetColorPalette(int(norm * (nCol - 1)));

        TMarker* m = new TMarker(xx, yy, 20);
        m->SetMarkerColor(ci);
        m->SetMarkerSize(1.2);
        m->Draw("SAME");

        // arrow from warm to warm+normalized displacement
        double x2 = xx + arrowsBefore[i].dx;
        double y2 = yy + arrowsBefore[i].dy;
		
		TArrow* arr = new TArrow(xx, yy, x2, y2, 0.0);  // 0.0 → no head
		arr->SetLineColor(ci);
		arr->SetLineWidth(2);
		arr->Draw("SAME");

        if (drawLabels) {
            double Rmm = dispBefore[i].r;
            char buf[64];
            std::snprintf(buf, sizeof(buf), "%.3f", Rmm);
            TLatex* t = new TLatex(xx, yy + 0.02 * (ymax - ymin), buf);
            t->SetTextAlign(21);
            t->SetTextSize(0.025);
            t->SetTextColor(kBlack);
            t->Draw("SAME");
        }
    }

    // AFTER
    cS->cd(2);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.06);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.12);

    TH2D* frameA = new TH2D("frameA", "After fit;X [mm];Y [mm]",
                            100, xmin, xmax, 100, ymin, ymax);
    frameA->SetStats(false);
    frameA->Draw("AXIS");

    {
        TLatex t;
        t.SetNDC(true);
        t.SetTextSize(0.045);
        t.SetTextAlign(13);
        t.DrawLatex(0.16, 0.92, "After fit");
    }

    for (size_t i = 0; i < n; ++i) {
        // Still plotting at original warm coordinates, as per design
        double xx = xw[i];
        double yy = yw[i];
        double R  = dispAfter[i].r;

        double norm = (R - Rmin_global) / (Rmax_global - Rmin_global);
        if (norm < 0.0) norm = 0.0;
        if (norm > 1.0) norm = 1.0;
        int ci = gStyle->GetColorPalette(int(norm * (nCol - 1)));

        TMarker* m = new TMarker(xx, yy, 20);
        m->SetMarkerColor(ci);
        m->SetMarkerSize(1.2);
        m->Draw("SAME");

        double x2 = xx + arrowsAfter[i].dx;
        double y2 = yy + arrowsAfter[i].dy;
        TArrow* arr = new TArrow(xx, yy, x2, y2, 0.03, ">");
        arr->SetLineColor(ci);
        arr->SetLineWidth(2);
        arr->Draw("SAME");

        if (drawLabels) {
            double Rmm = dispAfter[i].r;
            char buf[64];
            std::snprintf(buf, sizeof(buf), "%.3f", Rmm);
            TLatex* t = new TLatex(xx, yy + 0.02 * (ymax - ymin), buf);
            t->SetTextAlign(21);
            t->SetTextSize(0.025);
            t->SetTextColor(kBlack);
            t->Draw("SAME");
        }
    }

    cS->Write("ExpansionScanScatter");

    // Save PNG snapshots (before interactive session)
    string base = outRoot.substr(0, outRoot.find_last_of('.'));
    string pngH = base + "_hists.png";
    string pngS = base + "_scatter.png";

    cH->SaveAs(pngH.c_str());
    cS->SaveAs(pngS.c_str());

    cout << "\nRunning interactive ROOT...\n";
    app.Run();

    // After the interactive session, write and close the ROOT file
    outF.Write();
    outF.Close();

    cout << "\nSaved ROOT file : " << outRoot << "\n";
    cout << "Saved PNG (hists): " << pngH << "\n";
    cout << "Saved PNG (scat) : " << pngS << "\n";

    return 0;
}
