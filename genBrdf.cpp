/*

Copyright 2018 Lingqi Yan

This file is part of WaveOpticsBrdf.

WaveOpticsBrdf is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

WaveOpticsBrdf is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with WaveOpticsBrdf.  If not, see <https://www.gnu.org/licenses/>.

*/

#include <iostream>
#include <unistd.h>
#include <chrono>
#include <string>
#include "waveBrdf.h"
#include "spectrum.h"

using namespace std;
using namespace Eigen;

double deg2rad(double deg) {
    return deg / 180.0 * M_PI;
}

int main(int argc, char **argv) {
    srand(time(NULL));

    char *heightfieldFilename = NULL;
    double texelWidth = 1.0;
    double vertScale = 1.0;

    double mu_x = 0.0;
    double mu_y = 0.0;
    double sigma_p = 10.0;

    string method = "Wave";
    int sampleNum = 10000000;
    string diffModel = "OHS";
    double lambda = 0.0;

    // Resolution of the BRDF incident direction in degrees.
    double theta_res_deg = 5.0;
    double phi_res_deg = 30.0;

    char *outputFilename = NULL;
    int outputResolution = 256;

    double theta_i = 0.0;
    double phi_i = 0.0;
    bool rgb = false;  // Output BRDF in RGB or not.
    bool full = false;  // Output full BRDF or not.

    int opt;
    while ((opt = getopt(argc, argv, ":i:o:l:r:e:w:v:p:m:g:d:n:u:f:s:t:qz")) != -1) {
        switch (opt) {
            case 'i': heightfieldFilename = optarg; break;      // heightfield filename.
            case 'w': texelWidth = atof(optarg); break;         // The width of a texel in microns on the heightfield.
            case 'v': vertScale = atof(optarg); break;          // The vertical scaling factor of the heightfield.

            case 'p': sigma_p = atof(optarg); break;            // Size (1 sigma) of the Gaussian footprint.

            case 'm': method = optarg; break;                   // Method. Choose between "Geom" and "Wave".
            case 'n': sampleNum = atoi(optarg); break;          // Number of binning samples. Only valid for geometric optics.
            case 'd': diffModel = optarg; break;                // Diffraction model. Choose between "OHS", "GHS", "ROHS", "RGHS" and "Kirchhoff". And expect only subtle difference.
            case 'l': lambda = atof(optarg); break;             // Wavelength in microns. Once set, single wavelength mode is ON.

            case 'u': theta_i = atof(optarg); break;            // Incoming light's coordinate in spherical coordinates in degrees.
            case 'f': phi_i = atof(optarg); break;              // Incoming light's coordinate in spherical coordinates in degrees.

            case 's': theta_res_deg = atof(optarg); break;      // Resolution of the BRDF incident direction in degrees.
            case 't': phi_res_deg = atof(optarg); break;        // Resolution of the BRDF incident direction in degrees.

            case 'q': rgb = true; break;                // Output BRDF in RGB or not.
            case 'z': full = true; break;               // Output full BRDF or not.

            case 'o': outputFilename = optarg; break;           // output filename.
            case 'r': outputResolution = atoi(optarg); break;   // output resolution.
        }
    }

    if (heightfieldFilename == NULL) {
        cout << "A heightfield file must be specified." << endl;
        return -1;
    }

    if (outputFilename == NULL) {
        cout << "An output file must be specified." << endl;
        return -1;
    }

    printf("RGB mode: %d\n", rgb);

    SpectrumInit();

    EXRImage heightfieldImage(heightfieldFilename);
    Heightfield heightfield(&heightfieldImage, texelWidth, vertScale);

    Query query;

    // Set the center of the Gaussian footprint to the center of the heightfield.
    mu_x = heightfieldImage.width * texelWidth / 2.0;
    mu_y = heightfieldImage.height * texelWidth / 2.0;
    query.mu_p = Vector2f(mu_x, mu_y);
    printf("Center of the Gaussian footprint: (%f, %f)\n", mu_x, mu_y);

    query.sigma_p = sigma_p;
    query.lambda = lambda;

    if (method == "Geom") {
        GeometricBrdf geometricBrdf(&heightfield, sampleNum);

        // Convert spherical coordinates to Cartesian coordinates.
        theta_i = deg2rad(theta_i);
        phi_i = deg2rad(phi_i);
        query.omega_i = Vector2f(sin(theta_i) * cos(phi_i), sin(theta_i) * sin(phi_i));
        printf("Incoming light direction from sph (%f, %f): (x,y) = (%f, %f)\n", theta_i, phi_i, query.omega_i(0), query.omega_i(1));

        Float *brdfImage = geometricBrdf.genBrdfImage(query, outputResolution, rgb);
        EXRImage::writeImage(brdfImage, outputFilename, outputResolution, outputResolution);
        delete[] brdfImage;
    } else if (method == "Wave") {
        WaveBrdfAccel waveBrdfAccel(&heightfield, diffModel);

        if (!full) {

            Float *brdfImage = waveBrdfAccel.genBrdfImage(query, outputResolution, rgb);
            if (rgb) {
                EXRImage::writeImage(brdfImage, outputFilename, outputResolution, outputResolution);
            } else {
                EXRImage::writeRawSingleLayer(brdfImage, outputFilename, outputResolution, outputResolution, theta_i,
                                              phi_i);
            }
        } else {
            // Measure the BRDF with different incident angles and write to a single file.
            int num_theta = int(90.0 / theta_res_deg) + 1;
            int num_phi = int(360.0 / phi_res_deg);
            printf("Number of θ: %d, per %.2f°, Number of φ: %d, per %.2f°\n", num_theta, num_phi, theta_res_deg, phi_res_deg);

            float theta[num_theta];
            float phi[num_phi];

            for (int i = 0; i < num_phi; i++) {
                phi[i] = i * 2.0 * M_PI / num_phi;
            }

            for (int i = 0; i < num_theta; i++) {
                theta[i] = i * M_PI * 0.5 / num_theta;
            }

            Float **brdf_data = new Float*[num_theta * num_phi];
            for (int i = 0; i < num_phi; i++) {
                for (int j = 0; j < num_theta; j++) {
                    auto t = theta[j];
                    auto p = phi[i];
                    printf("Generating BRDF: θ%.2f, φ%.2f...\n", t, p);
                    query.omega_i = Vector2f(sin(t) * cos(p), sin(t) * sin(p));
                    brdf_data[i * num_theta + j] = waveBrdfAccel.genBrdfImage(query, outputResolution, rgb);
                }
            }

            EXRImage::writeRawMultiLayers(brdf_data, outputFilename, outputResolution, outputResolution, theta, phi, num_phi, num_theta);

            for (int i = 0; i < num_theta * num_phi; i++) {
                delete[] brdf_data[i];
            }
            delete[] brdf_data;
        }
    } else if (method == "GeomNdf") {
        GeometricBrdf geometricBrdf(&heightfield, sampleNum);
        Float *ndfImage = geometricBrdf.genNdfImage(query, outputResolution);
        EXRImage::writeImage(ndfImage, outputFilename, outputResolution, outputResolution);
        delete[] ndfImage;
    }

    return 0;
}
