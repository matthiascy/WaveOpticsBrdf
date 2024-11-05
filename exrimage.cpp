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

#include "exrimage.h"
#include "spectrum.h"
#include <OpenEXR/ImfMultiPartOutputFile.h>
#include <OpenEXR/ImfOutputPart.h>
#include <OpenEXR/ImfHeader.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfOutputFile.h>
#include <format>
#include <regex>
#include <string>

Float A_inv[16][16] = {{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                       {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                       {-3, 3, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                       {2, -2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                       {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
                       {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
                       {0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0},
                       {0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0},
                       {-3, 0, 3, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0},
                       {0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, -2, 0, -1, 0},
                       {9, -9, -9, 9, 6, 3, -6, -3, 6, -6, 3, -3, 4, 2, 2, 1},
                       {-6, 6, 6, -6, -3, -3, 3, 3, -4, 4, -2, 2, -2, -2, -1, -1},
                       {2, 0, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
                       {0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0},
                       {-6, 6, 6, -6, -4, -2, 4, 2, -3, 3, -3, 3, -2, -1, -2, -1},
                       {4, -4, -4, 4, 2, 2, -2, -2, 2, -2, 2, -2, 1, 1, 1, 1}};

double rad2deg(double rad) {
    return rad / M_PI * 180.0;
}

EXRImage::EXRImage(const char *filename) {
    readImage(filename);
}

void EXRImage::readImage(const char *filename) {
    RgbaInputFile file(filename);
    Imath::Box2i dw = file.dataWindow();
    width = dw.max.x - dw.min.x + 1;
    height = dw.max.y - dw.min.y + 1;
    image.resizeErase(height, width);
    file.setFrameBuffer(&image[0][0] - dw.min.x - dw.min.y * width, 1, width);
    file.readPixels(dw.min.y, dw.max.y);
}

void EXRImage::computeCoeff(Float *alpha, const Float *x) {
    memset(alpha, 0, sizeof(Float) * 16);
    for (int i = 0; i < 16; i++)
        for (int j = 0; j < 16; j++)
            alpha[i] += A_inv[i][j] * x[j];
}

// Bicubic interpolation
Float EXRImage::getValue(Float x, Float y) {
    return getValueUV(x / height, y / width);
}

// Bicubic interpolation
Float EXRImage::getValueUV(Float u, Float v) {
    Float x = u * height;
    Float y = v * width;
    int x1 = (int) floor(x);
    int y1 = (int) floor(y);
    int x2 = x1 + 1;
    int y2 = y1 + 1;

    Float a[16];
    Float xp[16] = {hp(x1, y1), hp(x2, y1), hp(x1, y2), hp(x2, y2),
                    hpx(x1, y1), hpx(x2, y1), hpx(x1, y2), hpx(x2, y2),
                    hpy(x1, y1), hpy(x2, y1), hpy(x1, y2), hpy(x2, y2),
                    hpxy(x1, y1), hpxy(x2, y1), hpxy(x1, y2), hpxy(x2, y2)};

    computeCoeff(a, xp);

    Float coeffA[4][4] = {{a[0], a[4], a[8], a[12]},
                          {a[1], a[5], a[9], a[13]},
                          {a[2], a[6], a[10], a[14]},
                          {a[3], a[7], a[11], a[15]}};
    
    Float h = 0.0f;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            h += coeffA[i][j] * pow(x - x1, (Float)(i)) * pow(y - y1, (Float)(j));
        }

    return h;
}

void EXRImage::writeImage(const Float *image, const char *filename, int outputHeight, int outputWidth) {
    // Clear image
    Rgba *pixels = new Rgba[outputHeight * outputWidth];

    // Write to image
    for (int i = 0; i < outputHeight; i++)
        for (int j = 0; j < outputWidth; j++) {
            pixels[i * outputWidth + j].r = image[(i * outputWidth + j) * 3 + 0];
            pixels[i * outputWidth + j].g = image[(i * outputWidth + j) * 3 + 1];
            pixels[i * outputWidth + j].b = image[(i * outputWidth + j) * 3 + 2];
            pixels[i * outputWidth + j].a = 1.0f;
        }

    RgbaOutputFile file(filename, outputHeight, outputWidth, WRITE_RGBA);
    file.setFrameBuffer(pixels, 1, outputWidth);
    file.writePixels(outputHeight);

    delete[] pixels;
}

string construct_layer_name(float theta, float phi) {
    string theta_str = std::format("{:.2f}", rad2deg(theta));
    auto idx = theta_str.find('.');
    theta_str.replace(idx, 1, "_");
    string phi_str = std::format("{:.2f}", rad2deg(phi));
    idx = phi_str.find('.');
    phi_str.replace(idx, 1, "_");

    return std::format("θ{}.φ{}", theta_str, phi_str);
}

void reorganise_data(const Float *data, Float *data_copy, int outputHeight, int outputWidth) {
    for (int k = 0; k < SPECTRUM_SAMPLES; k++)
        for (int i = 0; i < outputHeight; i++)
            for (int j = 0; j < outputWidth; j++)
                data_copy[(k * outputHeight + i) * outputWidth + j] = data[(i * outputWidth + j) * SPECTRUM_SAMPLES + k];
}

void EXRImage::writeRawSingleLayer(const Float *data, const char *filename, int outputHeight, int outputWidth, float theta, float phi) {
    Float *data_copy = new Float[SPECTRUM_SAMPLES * outputHeight * outputWidth];
    reorganise_data(data, data_copy, outputHeight, outputWidth);

    Header header(outputWidth, outputHeight);
    string theta_str = std::format("{:.2f}", theta);
    header.setName(construct_layer_name(theta, phi));
    header.setType("scanlineimage");
    FrameBuffer frameBuffer;

    // Each channel is a certain wavelength.
    for (int k = 0; k < SPECTRUM_SAMPLES; k++) {
        // in nanometers
        double lambda = (k + 0.5) / SPECTRUM_SAMPLES * (0.83 - 0.36) + 0.36;
        char channelName[256];
        sprintf(channelName, "%.0f nm", lambda * 1000.0);
        header.channels().insert(channelName, Channel(FLOAT));
        frameBuffer.insert(channelName, Slice(FLOAT, (char *) &data_copy[k * outputHeight * outputWidth], sizeof(Float), sizeof(Float) * outputWidth));
    }

    MultiPartOutputFile file(filename, &header, 1);
    OutputPart file_part(file, 0);
    file_part.setFrameBuffer(frameBuffer);
    file_part.writePixels(outputHeight);

    delete[] data_copy;
}

void EXRImage::writeRawMultiLayers(const Float* const*data, const char *filename, int outputHeight, int outputWidth, float* theta, float* phi, int n_phi, int n_theta) {
    auto n = n_phi * n_theta;
    Header headers[n];
    FrameBuffer framebuffers[n];
    Float *data_copy = new Float[n * SPECTRUM_SAMPLES * outputHeight * outputWidth];
    memset(data_copy, 0, sizeof(Float) * n * SPECTRUM_SAMPLES * outputHeight * outputWidth);

    for (int i = 0; i < n_phi; i++) {
        for (int j = 0; j < n_theta; j++) {
            auto idx = i * n_theta + j;
            headers[idx] = Header(outputWidth, outputHeight);
            headers[idx].setName(construct_layer_name(theta[j], phi[i]));
            // printf("Layer name: %s\n", headers[idx].name().c_str());
            headers[idx].setType("scanlineimage");

            reorganise_data(data[idx], &data_copy[idx * SPECTRUM_SAMPLES * outputHeight * outputWidth], outputHeight, outputWidth);

            // EXRImage::writeRawSingleLayer(data[idx], headers[idx].name().c_str(), outputHeight, outputWidth, theta[j], phi[i]);

            // Each channel is a certain wavelength.
            for (int k = 0; k < SPECTRUM_SAMPLES; k++) {
                // in nanometers
                double lambda = (k + 0.5) / SPECTRUM_SAMPLES * (0.83 - 0.36) + 0.36;
                char channelName[256];
                sprintf(channelName, "%.0f nm", lambda * 1000.0);
                headers[idx].channels().insert(channelName, Channel(FLOAT));
                framebuffers[idx].insert(channelName, Slice(FLOAT, (char *) &data_copy[(idx * SPECTRUM_SAMPLES + k) * outputHeight * outputWidth], sizeof(Float), sizeof(Float) * outputWidth));
            }
        }
    }


    MultiPartOutputFile file(filename, headers, n);

    for (int i = 0; i < n; i++) {
        OutputPart file_part(file, i);
        file_part.setFrameBuffer(framebuffers[i]);
        file_part.writePixels(outputHeight);
    }

    delete[] data_copy;
}