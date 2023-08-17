/************************************************/
/* Encoder.cc, (c) Rene Puchinger               */
/************************************************/

#include <cstdlib>
#include <cmath>
#include "Exception.h"
#include "FileIOStream.h"
#include "BitIOStream.h"
#include "Encoder.h"
#include <iterator>
#include <algorithm>
#include <random>
#include <iostream>
#include <iomanip>

Encoder::Encoder()
{
    im1 = NULL;
    im2 = NULL;
}

Encoder::~Encoder()
{
    if (im1 != NULL)
        delete im1;
    if (im2 != NULL)
        delete im2;
}

int Encoder::sub_dc(int component)
{
    float dc = 0;
    for (int y = 0; y < im1->get_size_y() + im1->get_extra_y(); y++)
        for (int x = 0; x < im1->get_size_x() + im1->get_extra_x(); x++)
            dc += im1->at(component, x, y);
    dc /= (im1->get_size_x() + im1->get_extra_x()) * (im1->get_size_y() + im1->get_extra_y());
    dc = floor(dc);
    for (int y = 0; y < im1->get_size_y() + im1->get_extra_y(); y++)
        for (int x = 0; x < im1->get_size_x() + im1->get_extra_x(); x++)
            im1->at(component, x, y) -= dc;
    return (int)dc;
}

void Encoder::dwt_row(int component, int row, int num_items)
{
    const float alpha = -1.586134342;
    const float beta = -0.05298011854;
    const float gamma = 0.8829110762;
    const float delta = 0.44355068522;
    const float xi = 1.149604398;

    for (int x = 0; x < num_items / 2 - 1; x++)
        im2->at(component, num_items / 2 + x, row) = im1->at(component, 2 * x + 1, row) + alpha * (im1->at(component, 2 * x, row) + im1->at(component, 2 * x + 2, row));
    im2->at(component, num_items - 1, row) = im1->at(component, num_items - 1, row) + 2 * alpha * im1->at(component, num_items - 2, row);

    im2->at(component, 0, row) = im1->at(component, 0, row) + beta * (im2->at(component, num_items / 2, row) + im2->at(component, num_items / 2 + 1, row));
    for (int x = 1; x < num_items / 2; x++)
        im2->at(component, x, row) = im1->at(component, 2 * x, row) + beta * (im2->at(component, num_items / 2 + x, row) + im2->at(component, num_items / 2 + x - 1, row));

    for (int x = 0; x < num_items / 2 - 1; x++)
        im2->at(component, num_items / 2 + x, row) += gamma * (im2->at(component, x, row) + im2->at(component, x + 1, row));
    im2->at(component, num_items - 1, row) += gamma * (im2->at(component, num_items / 2 - 1, row) + im2->at(component, num_items / 2 - 2, row));

    im2->at(component, 0, row) += delta * (im2->at(component, num_items / 2, row) + im2->at(component, num_items / 2 + 1, row));
    for (int x = 1; x < num_items / 2; x++)
        im2->at(component, x, row) += delta * (im2->at(component, num_items / 2 + x, row) + im2->at(component, num_items / 2 + x - 1, row));

    for (int x = 0; x < num_items / 2; x++)
    {
        im2->at(component, x, row) *= xi;
        im2->at(component, num_items / 2 + x, row) /= xi;
    }
}

void Encoder::dwt_col(int component, int col, int num_items)
{
    const float alpha = -1.586134342;
    const float beta = -0.05298011854;
    const float gamma = 0.8829110762;
    const float delta = 0.44355068522;
    const float xi = 1.149604398;

    for (int y = 0; y < num_items / 2 - 1; y++)
        im1->at(component, col, num_items / 2 + y) = im2->at(component, col, 2 * y + 1) + alpha * (im2->at(component, col, 2 * y) + im2->at(component, col, 2 * y + 2));
    im1->at(component, col, num_items - 1) = im2->at(component, col, num_items - 1) + 2 * alpha * im2->at(component, col, num_items - 2);

    im1->at(component, col, 0) = im2->at(component, col, 0) + beta * (im1->at(component, col, num_items / 2) + im1->at(component, col, num_items / 2 + 1));
    for (int y = 1; y < num_items / 2; y++)
        im1->at(component, col, y) = im2->at(component, col, 2 * y) + beta * (im1->at(component, col, num_items / 2 + y) + im1->at(component, col, num_items / 2 + y - 1));

    for (int y = 0; y < num_items / 2 - 1; y++)
        im1->at(component, col, num_items / 2 + y) += gamma * (im1->at(component, col, y) + im1->at(component, col, y + 1));
    im1->at(component, col, num_items - 1) += gamma * (im1->at(component, col, num_items / 2 - 1) + im1->at(component, col, num_items / 2 - 2));

    im1->at(component, col, 0) += delta * (im1->at(component, col, num_items / 2) + im1->at(component, col, num_items / 2 + 1));
    for (int y = 1; y < num_items / 2; y++)
        im1->at(component, col, y) += delta * (im1->at(component, col, num_items / 2 + y) + im1->at(component, col, num_items / 2 + y - 1));

    for (int y = 0; y < num_items / 2; y++)
    {
        im1->at(component, col, y) *= xi;
        im1->at(component, col, num_items / 2 + y) /= xi;
    }
}

void Encoder::dwt2(int component, int size_x, int size_y)
{
    for (int y = 0; y < size_y; y++)
        dwt_row(component, y, size_x);
    for (int x = 0; x < size_x; x++)
        dwt_col(component, x, size_y);
}

void Encoder::dwt2full(int component)
{
    int size_x = im1->get_size_x() + im1->get_extra_x();
    int size_y = im1->get_size_y() + im1->get_extra_y();
    for (int i = 0; i < num_stages; i++)
    {
        dwt2(component, size_x, size_y);
        size_x /= 2;
        size_y /= 2;
    }
}

void Encoder::normalize(int component)
{
    for (int y = 0; y < im1->get_size_y() + im1->get_extra_y(); y++)
        for (int x = 0; x < im1->get_size_x() + im1->get_extra_x(); x++)
        {
            if (im1->at(component, x, y) >= 0)
                im1->at(component, x, y) = floor(im1->at(component, x, y));
            else
                im1->at(component, x, y) = -floor(fabs(im1->at(component, x, y)));
        }
}

void Encoder::encode(char *fn_in, char *fn_out, float bit_rate)
{
    im1 = new Image();
    im1->load(fn_in, &num_stages);
    // std::vector<float> data;
    // data.resize(im1->get_size_x() * im1->get_size_y());
    // std::mt19937 gen;
    // gen.seed(0);
    // std::uniform_real_distribution<float> dis(0, 1);
    // std::generate(data.begin(), data.end(), [&dis, &gen]()
    //               { return dis(gen); });
    // std::transform(data.begin(), data.end(), data.begin(), [](float value)
    //                { return value * 255.f; });
    // im1->Y = data.data();
    im1->Cb = nullptr;
    im1->Cr = nullptr;

    for (int y = 0; y < im1->get_size_y(); y++)
    {
        for (int x = 0; x < im1->get_size_x(); x++)
        {
            int t = (int)trunc(im1->at(0, x, y));
            // std::cout << "0x" << std::hex << t << std::dec << ", ";
            std::cout << t << ", ";
        }
        std::cout << std::endl;
    }

    float *Y = new float[]{
        4765, -871.5, -5.187, 493.8, 45.35, -156.5, -92.17, -49.51, -28.18, -49.14, -40.73, 5.897, 9.047, -5.99, -111.7, 14.25, -27.69, -6.088, -15.36, -8.036, -9.983, -12.59, -2.983, -0.3529, 0.9839, 2.323, 0.5499, 11.91, -16.17, 5.088, -10.17, 0.4448,
        -476.4, -1033, 57.24, -441.4, -25.16, 31.19, -206.7, -242.5, -2.782, 1.768, -0.08953, 24.75, -5.571, -18.05, 51.49, 2.032, 1.508, 0.2646, -0.444, -0.4754, -0.5755, 0.9333, -3.997, 3.492, -0.04992, 0.2774, 0.3492, -0.08407, -0.3482, -1.822, 2.682, 8.415,
        -675.4, -309.8, -17.68, -552, -56.87, -32.41, 202.7, -91.86, -5.58, -4.789, -5.649, -27.06, -20.74, -7.355, 0.1976, 24.2, -1.018, -0.8362, 1.045, -0.2505, -0.4954, 0.2663, 1.772, -3.65, -11.56, 9.109, 0.7137, -2.421, 1.096, 10.05, 3.26, -2.64,
        348.7, -1428, -18.85, -347.3, 83.41, 44.21, -89.37, 132.5, -4.818, -6.66, -2.53, 22.39, -27.68, 90.76, 53.3, 18.48, 2.824, -0.5578, 0.2577, 1.03, 1.519, 0.4675, 4.514, 17.52, 5.995, -3.523, 1.827, -1.156, -22.7, -4.468, -0.5315, 3.672,
        -328.9, 20.86, -27.64, 103.3, -42.81, -34.44, 131.8, 288.7, 7.806, 51.78, 73.43, 34.69, 165.9, 133.8, 8.011, -29.5, -0.9664, -0.6209, -0.0101, -0.5325, 1.288, -3.872, -7.598, -7.027, -0.7516, -3.347, -7.777, 25.43, 1.503, -2.759, 7.835, 1.017,
        409.3, 1.101, 130.2, -705.7, 6.807, -2.443, 175.6, 139.6, -6.299, 19.16, 26.25, -25.61, -35.58, -11.49, -6.045, -7.902, 2.555, -0.3656, -1.319, -0.6755, 2.298, -5.921, -1.28, -9.375, 13.42, -7.446, -6.415, -0.5507, -25.49, -5.024, -3.065, 0.528,
        527.4, 139.6, 145.1, 70.5, 62.28, -41.39, -331.7, -48.42, 74.85, -3.485, -18.64, 12.65, 19.17, 151, 95.31, -59.19, 0.5885, 1.111, 0.02895, 4.731, 4.151, 1.312, -13.24, -1.923, 11.37, -1.409, -16.19, 15.6, 5.972, 3.737, 0.2859, 0.9307,
        205.5, -5.407, -67.76, -279.3, 60.27, 73.88, -116.5, 55.62, -2.931, 3.43, 9.885, 13.11, 4.668, -0.9589, -161.2, -57.1, -1.647, 0.5475, 0.4203, 0.6646, -0.09409, -4.987, 4.377, 5.114, 6.067, 1.784, -2.237, 0.8905, 0.17, -1.967, 4.479, -3.396,
        -91.07, 39.15, -25.9, 11.18, 0.07828, 52.12, 6.977, 12.2, 41.08, 14.85, 4.322, -2.413, -2.8, 9.312, -40.15, 34.67, 2.747, 0.8912, -0.2081, -13.3, -12.69, -8.435, -10.97, -1.671, 29.58, -12.66, 6.803, -9.322, -0.7772, -1.008, -0.8723, -3.855,
        -142, 55.88, -0.3284, -30.45, 8.477, 66.6, -61.57, 25.14, 54.57, -1.47, -12.12, 5.853, -26.94, -32.61, -46.45, -40.65, -5.477, -0.1006, 0.9391, 9.982, -0.4072, -11.88, 6.249, 8.627, -2.136, 1.591, 0.2882, -1.297, 0.289, -1.484, -1.454, -8.544,
        -9.564, 81.95, -7.636, 13.65, 31.72, -84.98, 173, 101, -1.686, -3.108, -21.65, 17.33, -0.7871, 41.56, -30.37, -40.11, -0.2894, -0.9004, -0.6926, -12.35, 0.8284, 22.21, -0.1467, 1.644, 1.23, 15.91, -2.41, 3.646, 1.111, -0.7268, -0.5127, -0.7971,
        -6.887, 63.83, -185.2, 50.72, 79.86, -219, 4.434, 30.23, -1.491, 3.08, -33.49, 24.38, -87.29, -3.836, -93.35, -4.012, 0.1302, 0.4912, -9.634, 3.366, 7.19, 5.176, 2.056, -0.3612, -3.988, -30.44, -1.496, -1.787, -0.8129, 0.9477, -0.5242, -0.8197,
        -26.85, 55.49, 26.52, 179.5, 100.3, 44.95, 182.4, -29.98, 5.602, -13.04, -22.18, 112.5, -5.983, -46.76, 7.541, -21.41, -6.225, -4.334, -24.83, 1.671, 17.75, 22.86, 2.568, -4.615, 3.503, -13.35, 11.6, 6.801, 10.62, -14.06, -7.723, -27.74,
        47.81, 61.16, -80.43, 72.61, 5.358, -93.7, -46.68, 53.36, 7.02, 3.546, -29.83, -65.28, -6.898, 17.6, 2.609, 6.16, 0.8536, 0.7215, 6.209, -0.5273, -1.835, -16.76, 2.273, -2.918, -1.208, 1.275, 3.333, 11.03, 3.015, -5.026, -0.9705, -1.835,
        -177.4, -11.37, 9.161, 34, 0.8473, -89.77, -25.13, -38.47, -23.05, 13.52, 5.811, 1.831, -11.39, -42.96, -4.242, -25.47, -0.5935, 0.9248, 4.737, 5.929, 5.215, 10.1, -3.632, -2.943, 0.7763, 0.8212, 0.4111, 6.639, -3.102, -1.055, 7.186, 2.338,
        14.99, 28.15, 42.3, 77.07, -4.236, -50.77, -81.03, 96.48, 12.71, -9.808, 39.59, -28.71, 3.196, 5.877, -48, 80.73, -8.722, -1.819, -5.619, -0.9512, -6.87, 2.055, -24.83, -10.88, -0.1872, -0.8055, 0.4968, 0.4856, -4.152, -26.64, 12.11, 7.947,
        -14.24, 1.285, 9.709, 0.4253, -2.127, 3.428, 0.904, -2.697, 0.3477, 1.847, 9.899, -0.2747, -6.775, -39.61, 4.277, 9.758, -4, 3, -0.1875, 0.9062, 0.5, -1.281, 0.2188, 0.2188, -0.09375, 0.4375, 0.625, 0.1875, 3.531, 1.094, -0.375, 4.219,
        -24.75, 0.2829, 11.59, 1.669, -0.2171, 0.2459, -2.012, -0.121, -2.21, -1.282, 1.054, -1.495, -0.07897, 23.91, 31.08, 13.28, -0.6875, 0.2812, -0.2812, -0.0625, 0.2188, 0.625, 0.03125, 0.4688, 0.6875, 0.1563, -0.1875, -0.1562, -0.125, 0.4687, 8.375, -3.875,
        -29.29, -5.14, 13.84, 0.9172, -1.686, 6.301, 0.3258, 1.236, 1.973, -18.4, 1.384, -7.075, 4.823, -0.4406, 16.62, 0.03228, 0.5937, 0.25, 0.0625, 0.4063, -0.2188, -0.8125, 1.89E-15, 0.7188, -2.625, -12.09, 4.35E-15, -1.469, 0.1875, -2.812, -7.781, -0.0625,
        21.47, -2.508, 12.46, 0.2267, -16.7, 0.02615, -1.836, 0.8834, -0.4849, 1.569, 31.35, 1.799, 7.299, -35.22, 5.49, -0.8627, -2.188, -0.03125, -0.6875, 1.031, -2.531, -0.5313, 1.75, -3.469, 1.375, 1.969, 2.812, 0.4375, 2.562, 6.75, -1.125, -13.28,
        13.3, -1.836, 12.44, 3.457, 9.681, -2.251, -2.852, 3.13, -1.587, -1.898, 19.33, -13.08, -4.981, 15.32, -2.107, -1.4, -0.09375, -0.0625, -0.09375, -0.2188, -0.1562, -0.8125, 2.719, 2.5, -0.9375, 2.188, -8.312, 5.625, 2.437, 1.969, 10.28, 1.125,
        12.97, -3.656, 14.6, -4.667, 9.278, -3.683, 2.019, 1.335, -1.714, 1.334, 3.607, 2.664, 9.507, -0.4945, -9.55, -1.528, -0.0312, 0.2187, 0.4063, -0.1562, 2.312, 0.625, -2.125, 0.6875, -2.219, -0.4688, -1.75, -1.188, -8.406, 1.125, 2.813, -0.09375,
        9.033, -4.068, 15.93, 2.323, -14.89, -2.329, 2.458, 32.28, -6.175, 3.835, 14.41, 5.372, -4.434, 7.368, -1.011, 0.2752, -0.3437, -0.09375, 0.1563, 0.125, -2.656, -1.969, 1.031, -3.906, -5.281, 0.3438, -1.781, -5.281, -3.813, 4.125, -0.03125, 0.2187,
        10.08, -3.07, 14.11, 6.132, 2.341, -5.367, -3.832, -5.471, 13.81, 0.6134, 11.71, 71.58, -4.013, -12.53, 0.6195, 0.3269, 0.3437, 0.2188, 0.375, -0.2188, -1.563, 8.75, 2.187, 12.69, 4.219, -2.438, -0.4062, 2.031, -0.6875, 1.625, -0.8438, 0.3125,
        9.752, -4.828, 14.92, 19.58, 4.388, -35.48, -18.66, -3.349, -11.62, -3.342, 18.99, 11.71, -4.099, -6.946, 0.03792, -1.46, -0.6563, -0.5313, -0.03125, -3, 15.63, 0.5625, 4.625, -0.1563, -8.406, 4.062, 3.25, 1.344, -2.125, 0.2812, 0.2812, 0.3125,
        25.57, -7.237, 12.51, -0.09905, -13.04, 46.9, 46.5, -2.672, -3.063, -12.75, -0.9568, -15.69, 17.53, -0.6827, -4.561, -3.26, 0.6562, 0.4062, -0.1875, 3.625, 5.281, 2.125, 7.594, 1.906, -0.125, -0.125, 1, 0.5625, -0.625, -0.875, -1.125, 1.625,
        27.13, -6.735, 9.946, 0.7146, -43.43, -12.81, 3.969, -1.738, 4.637, 1.733, 22.25, -23.77, 3.307, -0.02164, 24.6, 0.3495, 0.75, -0.2187, 0.0625, -3.437, -0.9375, 18.16, 5.563, -0.9062, -0.25, 1.219, 1.625, -0.7813, -0.5, 0.2812, -0.3438, 0.09375,
        33.4, -4.071, 24.27, 8.741, -7.234, -16.38, 2.868, 12.19, -2.747, -1.733, 13.22, -4.425, 1.449, 1.312, 20.44, -2.345, -0.6875, -0.2187, 9.531, -2.75, 12.31, -21.22, 0.75, 1.344, 0.25, -0.1563, -0.9687, 6.66E-16, -0.9063, 0.4375, 0.3437, -0.8125,
        -15.25, -7.402, 18.19, -9.1, -10.44, -28.31, 3.402, 3.685, -0.07084, 3.845, -8.35, -13.03, 7.536, 0.152, -18.02, -19.64, 1.625, -0.375, 10.53, 3.344, -3.5, -11.97, -0.4063, -1.438, 0.6563, 0.8437, -9.781, -0.03125, -1.906, 1.406, -0.09375, -7.906,
        59.24, -8.916, -15.15, 9.099, 8.446, -30.24, -6.677, -37.76, -0.4617, -0.8695, -1.714, -37.48, -0.1264, -11.87, 18.29, -9.267, -0.9062, -0.0625, -1.938, -2.937, 6.031, 1.094, -0.375, 2.031, -0.4688, -0.5625, 0.4688, -14.66, 0.0625, -3.406, 2.094, -3.656,
        39.06, -1.141, -11.2, 7.817, 4.287, 22.64, -13.78, -12.79, -0.9268, 0.04607, 1.238, 22.13, 2.12, -5.886, -17.03, -2.213, 0.09375, 0.09375, -2.25, -0.3438, 1.813, 13.66, 0.5625, 2.938, -0.9375, -0.0625, -0.0625, 1.625, 0.8437, -1, -0.25, -5.219,
        -2.83, 12.1, 2.453, 7.136, -2.65, -5.009, 0.5901, -10.9, 1.162, -0.5426, -0.6655, -2.171, -0.9408, 20.57, 9.893, 3.079, 1.344, -0.9062, 4.219, -1.062, -0.8437, 0.4063, -1.5, 3, -0.03125, 0.5, 0.6563, 1.063, 1.719, -7.437, 7.187, -1.125};
    num_stages = 5;

    im1->Y = Y;
    im1->size_x = 32;
    im1->size_y = 32;

    im2 = new Image(im1->is_color(), im1->get_size_x(), im1->get_size_y(), im1->get_extra_x(), im1->get_extra_y());
    FileOutputStream *fout = new FileOutputStream(fn_out);
    BitOutputStream *bout = new BitOutputStream(fout);
    /* write IMS header */
    bout->put_bits('I', 8);
    bout->put_bits('M', 8);
    bout->put_bits('S', 8);
    bout->put_bits(num_stages, 6);
    bout->put_bits(im1->get_size_x(), 12);
    bout->put_bits(im1->get_size_y(), 12);
    bout->put_bits(im1->get_extra_x(), 10);
    bout->put_bits(im1->get_extra_y(), 10);
    bout->put_bit((int)im1->is_color());
    SPIHT_Encoder *spiht_enc = new SPIHT_Encoder(im1, num_stages);
    if (im1->is_color())
    {
        int bits0 = (int)ceil(bit_rate * im1->get_size() * 0.6);
        int bits1 = (int)ceil(bit_rate * im1->get_size() * 0.2);
        int bits2 = (int)ceil(bit_rate * im1->get_size() * 0.2);
        bout->put_bits(bits0, 29);
        bout->put_bits(bits1, 29);
        bout->put_bits(bits2, 29);
        int dc0 = sub_dc(0);
        int dc1 = sub_dc(1);
        int dc2 = sub_dc(2);
        bout->put_bits(dc0, 8);
        bout->put_bits(dc1, 8);
        bout->put_bits(dc2, 8);
        dwt2full(0);
        dwt2full(1);
        dwt2full(2);
        normalize(0);
        normalize(1);
        normalize(2);
        spiht_enc->encode(0, bits0, bout);
        spiht_enc->encode(1, bits1, bout);
        spiht_enc->encode(2, bits2, bout);
    }
    else
    {
        int bits0 = (int)ceil(bit_rate * im1->get_size());
        // bout->put_bits(bits0, 29);
        // int dc0 = sub_dc(0);
        // bout->put_bits(dc0, 8);
        // dwt2full(0);
        // normalize(0);
        spiht_enc->encode(0, bits0, bout);
    }
    bout->flush_bits();
    bout->flush();
    delete spiht_enc;
    delete bout;
    delete fout;
}
