// --------------------------------------------------------------------------------------------------
// EE569 Homework_3 algorithm functions
// Date : 3/27/2018
// Name : Shanlins Sun
// USCID: 8376995433
// email: shanlins@usc.edu
// --------------------------------------------------------------------------------------------------
// implementations of algorithms
// --------------------------------------------------------------------------------------------------

#include "HW3_Shanlins.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <queue>
#include <vector>
#include <iterator>
#include <valarray>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

unsigned char*** SpaceForRGB(unsigned char*** Img_Rgb, int Img_Height, int Img_Width)
{
    /* Allocate a 3D array for RGB image -- array[Height][Width][3] */
    Img_Rgb = new unsigned char **[Img_Height]();
    for (int i = 0; i < Img_Height; i++)
    {
        Img_Rgb[i] = new unsigned char*[Img_Width]();
        for (int j = 0; j < Img_Width; j++)
        {
            Img_Rgb[i][j] = new unsigned char[3]();
            for (int k = 0; k < 3; k++)
            {
                Img_Rgb[i][j][k] = 0;
            }
        }
    }
    return Img_Rgb;
}

unsigned char*** Convert1dTo3d(unsigned char* Img_1d, int Img_Height, int Img_Width)
{
    /* Convert 1D image data to 3D */
    unsigned char*** Img_3d = new unsigned char **[Img_Height]();
    for (int i = 0; i < Img_Height; i++)
    {
        Img_3d[i] = new unsigned char *[Img_Width]();
        for (int j = 0; j < Img_Width; j++)
        {
            Img_3d[i][j] = new unsigned char[3]();
            for (int k = 0; k < 3; k++)
            {
                Img_3d[i][j][k] = 0;
            }
        }
    }
    int m = 0;
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            Img_3d[i][j][0] = Img_1d[m];
            Img_3d[i][j][1] = Img_1d[m + 1];
            Img_3d[i][j][2] = Img_1d[m + 2];
            m += 3;
        }
    }
    return Img_3d;
}

unsigned char* Convert3dTo1d(unsigned char*** Img_3d, int Img_Height, int Img_Width)
{
    unsigned char* Img_1d = new unsigned char[Img_Height*Img_Width * 3]();
    int k = 0;
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            Img_1d[k] = Img_3d[i][j][0];
            Img_1d[k + 1] = Img_3d[i][j][1];
            Img_1d[k + 2] = Img_3d[i][j][2];
            k += 3;
        }
    }
    return Img_1d;
}

unsigned char** SpaceFor2D(unsigned char** Img_2d, int Img_Height, int Img_Width)
{
    /* Allocate a 2D array for gray image -- array[Height][Width] */
    Img_2d = new unsigned char *[Img_Height]();
    for (int i = 0; i < Img_Height; i++)
    {
        Img_2d[i] = new unsigned char[Img_Width]();
        for (int j = 0; j < Img_Width; j++)
        {
            Img_2d[i][j] = 0;
        }
    }
    return Img_2d;
}

int** SpaceFor2DInt(int** Img_2d, int Img_Height, int Img_Width)
{
    /* Allocate a 2D array for gray image -- array[Height][Width] */
    Img_2d = new int *[Img_Height]();
    for (int i = 0; i < Img_Height; i++)
    {
        Img_2d[i] = new int[Img_Width]();
        for (int j = 0; j < Img_Width; j++)
        {
            Img_2d[i][j] = 0;
        }
    }
    return Img_2d;
}

double** SpaceForNor2D(double** Img_2d, int Img_Height, int Img_Width)
{
    /* Allocate a normalized 2D array for gray image -- array[Height][Width] */
    Img_2d = new double *[Img_Height]();
    for (int i = 0; i < Img_Height; i++)
    {
        Img_2d[i] = new double[Img_Width]();
        for (int j = 0; j < Img_Width; j++)
        {
            Img_2d[i][j] = 0.0;
        }
    }
    return Img_2d;
}

unsigned char** Convert1dTo2d(unsigned char* Img_1d, int Img_Height, int Img_Width)
{
    unsigned char** Img_2d = new unsigned char *[Img_Height]();
    for (int i = 0; i < Img_Height; i++)
    {
        Img_2d[i] = new unsigned char[Img_Width]();
        for (int j = 0; j < Img_Width; j++)
        {
            Img_2d[i][j] = 0;
        }
    }

    int m = 0;
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            Img_2d[i][j] = Img_1d[m];
            m += 1;
        }
    }
    return Img_2d;
}

unsigned char* Convert2dTo1d(unsigned char** Img_2d, int Img_Height, int Img_Width)
{
    unsigned char* Img_1d = new unsigned char[Img_Height*Img_Width]();
    int k = 0;
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            Img_1d[k] = Img_2d[i][j];
            k += 1;
        }
    }
    return Img_1d;
}

double* Convert2dTo1dD(double** Img_2d_D, int Img_Height, int Img_Width)
{
    double* Img_1d_D = new double[Img_Height*Img_Width]();
    int k = 0;
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            Img_1d_D[k] = Img_2d_D[i][j];
            k += 1;
        }
    }
    return Img_1d_D;
}

void Raw_Read(char* path, unsigned char* Img_1d, long int data_size)
{
    FILE* file = fopen(path, "rb");

    if (file == NULL)
    {
        cout << "Fail to open the file " << path << endl;
        exit(1);
    }
    fread(Img_1d, sizeof(unsigned char), data_size, file);
    fclose(file);
    cout << "Loading " << path << "... Successfully" << endl;
}

void Raw_Write(char* path, unsigned char* Img_1d, long int data_size)
{
    FILE* file = fopen(path, "wb");

    if (file == NULL)
    {
        cout << "Fail to write the image into the file " << path << endl;
        exit(1);
    }
    fwrite(Img_1d, sizeof(unsigned char), data_size, file);
    fclose(file);
    cout << "Saving " << path << "... Successfully" << endl;
}

void Hist_Write(char* path, unsigned long Hist[256])
{
    ofstream file(path);
    if (!file)
    {
        cout << "Fail to write into " << path << endl;
        exit(1);
    }

    for (int i = 0; i < 256; i++)
    {
        file << Hist[i] << endl;
    }
    file.close();
    cout << "Writing " << path << "... Successfully" << endl;
}

void PSNR_Write(char* path, double Psnr_Red, double Psnr_Green, double Psnr_Blue, double Psnr)
{
    ofstream file(path);
    if (!file)
    {
        cout << "Fail to write into " << path << endl;
        exit(1);
    }

    file << Psnr_Red << endl;
    file << Psnr_Green << endl;
    file << Psnr_Blue << endl;
    file << Psnr << endl;

    file.close();
    cout << "Writing " << path << "... Successfully" << endl;
}

void Gray_Lightness(unsigned char* Img_1d_input, unsigned char* Img_1d_output, int image_size)
{
    int lightest_pixel = 0;
    int darkest_pixel = 0;
    int a = 0;
    int b = 0;
    int c = 0;
    for (int i = 0; i < image_size; i++)
    {
        a = (int)(Img_1d_input[3 * i]);
        b = (int)(Img_1d_input[3 * i + 1]);
        c = (int)(Img_1d_input[3 * i + 2]);
        lightest_pixel = max(max(a, b), c);
        darkest_pixel = min(min(a, b), c);
        Img_1d_output[i] = (unsigned char)(((double)(lightest_pixel)+(double)(darkest_pixel)+1) / 2);
    }
}

void Gray_Average(unsigned char* Img_1d_input, unsigned char* Img_1d_output, int image_size)
{
    for (int i = 0; i < image_size; i++)
    {
        Img_1d_output[i] = (unsigned char)(((double)(Img_1d_input[3 * i]) + (double)(Img_1d_input[3 * i + 1]) + (double)(Img_1d_input[3 * i + 2]) + 1.5) / 3);
    }
}

void Gray_Luminosity(unsigned char* Img_1d_input, unsigned char* Img_1d_output, int image_size)
{
    for (int i = 0; i < image_size; i++)
    {
        Img_1d_output[i] = (unsigned char)((double)(Img_1d_input[3 * i]) * 0.21 + (double)(Img_1d_input[3 * i + 1]) * 0.72 + (double)(Img_1d_input[3 * i + 2]) * 0.07 + 0.5);
    }
}

void RGB2CMY(unsigned char* Img_1d_input, unsigned char* Img_1d_output, int image_size)
{
    for (int i = 0; i < image_size; i++)
    {
        Img_1d_output[i] = 255 - Img_1d_input[i];
    }
}

void ChannelSplit(unsigned char* Img_1d_input, unsigned char* Img_1d_output, int image_size, int channel)
{
    for (int i = 0; i < image_size; i++)
    {
        Img_1d_output[i] = Img_1d_input[3 * i + channel];
    }
}

void ChannelCombine(unsigned char* Img_1d_red, unsigned char* Img_1d_green, unsigned char* Img_1d_blue, unsigned char* Img_1d_color, int image_size)
{
    for (int i = 0; i < image_size; i++)
    {
        Img_1d_color[3 * i] = Img_1d_red[i];
        Img_1d_color[3 * i + 1] = Img_1d_green[i];
        Img_1d_color[3 * i + 2] = Img_1d_blue[i];
    }
}

void Resize3d(unsigned char*** Img_3d_input, unsigned char*** Img_3d_output, int img_input_height, int img_input_width, int img_output_height, int img_output_width)
{
    double x = 0;
    double y = 0;
    int x_input_1 = 0;
    int x_input_2 = 0;
    int y_input_1 = 0;
    int y_input_2 = 0;
    int x_output = 0;
    int y_output = 0;
    int a = 0;
    int b = 0;
    int c = 0;
    int d = 0;
    double x_del = 0.00;
    double y_del = 0.00;
    double pixel_val = 0.00;
    int pixel = 0;

    /* Allocate 3D dynamic space for augmented input image */
    unsigned char*** Img_3d = NULL;
    Img_3d = SpaceForRGB(Img_3d, img_input_height + 1, img_input_width + 1);
    for (int i = 0; i < img_input_height; i++)
    {
        for (int j = 0; j < img_input_width; j++)
        {
            Img_3d[i][j][0] = Img_3d_input[i][j][0];
            Img_3d[i][j][1] = Img_3d_input[i][j][1];
            Img_3d[i][j][2] = Img_3d_input[i][j][2];
        }
    }

    for (x_output = 0; x_output < img_output_height; x_output++)
    {
        for (y_output = 0; y_output < img_output_width; y_output++)
        {
            x = (double)((img_input_height - 1)*x_output) / (double)(img_output_height - 1);
            x_input_1 = (int)(x);
            x_input_2 = x_input_1 + 1;
            y = (double)((img_input_width - 1)*y_output) / (double)(img_output_width - 1);
            y_input_1 = (int)(y);
            y_input_2 = y_input_1 + 1;
            x_del = x - x_input_1;
            y_del = y - y_input_1;

            a = Img_3d[x_input_1][y_input_1][0];
            b = Img_3d[x_input_1][y_input_1 + 1][0];
            c = Img_3d[x_input_1 + 1][y_input_1][0];
            d = Img_3d[x_input_1 + 1][y_input_1 + 1][0];
            pixel_val = (a*(1 - y_del) + b * y_del)*(1 - x_del) + (c*(1 - y_del) + d * y_del)*x_del;
            pixel = (int)(pixel_val + 0.5);
            Img_3d_output[x_output][y_output][0] = pixel;

            a = Img_3d[x_input_1][y_input_1][1];
            b = Img_3d[x_input_1][y_input_1 + 1][1];
            c = Img_3d[x_input_1 + 1][y_input_1][1];
            d = Img_3d[x_input_1 + 1][y_input_1 + 1][1];
            pixel_val = (a*(1 - y_del) + b * y_del)*(1 - x_del) + (c*(1 - y_del) + d * y_del)*x_del;
            pixel = (int)(pixel_val + 0.5);
            Img_3d_output[x_output][y_output][1] = pixel;

            a = Img_3d[x_input_1][y_input_1][2];
            b = Img_3d[x_input_1][y_input_1 + 1][2];
            c = Img_3d[x_input_1 + 1][y_input_1][2];
            d = Img_3d[x_input_1 + 1][y_input_1 + 1][2];
            pixel_val = (a*(1 - y_del) + b * y_del)*(1 - x_del) + (c*(1 - y_del) + d * y_del)*x_del;
            pixel = (int)(pixel_val + 0.5);
            Img_3d_output[x_output][y_output][2] = pixel;
        }
    }
}

void Hist1d(unsigned long hist[256], unsigned char *Img_gray, int Img_Height, int Img_Width)
{
    for (int i = 0; i < 256; i++)
    {
        hist[i] = 0;
    }

    for (int j = 0; j < Img_Height*Img_Width; j++)
    {
        hist[(int)(Img_gray[j])] += 1;
    }
}

void Hist2d(unsigned long hist[256], unsigned char **Img_gray, int Img_Height, int Img_Width)
{
    for (int i = 0; i < 256; i++)
    {
        hist[i] = 0;
    }

    for (int j = 0; j < Img_Height; j++)
    {
        for (int k = 0; k < Img_Width; k++)
        {
            hist[(int)(Img_gray[j][k])] += 1;
        }
    }
}

void CumHist1d(unsigned long hist[256], unsigned long cum_hist[256])
{
    for (int i = 0; i < 256; i++)
    {
        cum_hist[i] = 0;
    }
    cum_hist[0] = hist[0];
    for (int j = 1; j < 256; j++)
    {
        cum_hist[j] += cum_hist[j - 1] + hist[j];
    }
}

void Equ_CMF(unsigned long hist_input[256], unsigned long gray_value_transfer[256], unsigned char** Img_gray, unsigned char** Img_gray_equ, int Img_Height, int Img_Width)
{
    int cum_pixel = 0;
    int gray_value = 0;
    double cum_prob = 0;
    for (int i = 0; i < 256; i++)
    {
        gray_value_transfer[i] = 0;
    }
    for (int i = 0; i < 256; i++)
    {
        cum_pixel += hist_input[i];
        cum_prob = double(cum_pixel) / (Img_Height*Img_Width);
        gray_value = (int)(cum_prob * 255);
        gray_value_transfer[i] = gray_value;
    }
    for (int j = 0; j < Img_Height; j++)
    {
        for (int k = 0; k < Img_Width; k++)
        {
            Img_gray_equ[j][k] = (unsigned char)gray_value_transfer[(int)(Img_gray[j][k])];
        }
    }
}

void Equ_BF(unsigned char* Img_gray, unsigned char* Img_gray_equ, int Img_Height, int Img_Width)
{
    /* Initialize Img_gray_equ */
    for (int i = 0; i < Img_Height*Img_Width; i++)
    {
        Img_gray_equ[i] = 0;
    }
    /* Sort the original image pixel by the grayscales */
    unsigned int *sort_index = new unsigned int[Img_Height*Img_Width]();
    long temp_index = 0;
    int temp_value = 0;
    for (int i = 0; i < Img_Height*Img_Width; i++)
    {
        sort_index[i] = i;
    }
    for (int j = 0; j < Img_Height*Img_Width; j++)
    {
        for (int k = j + 1; k < Img_Height*Img_Width; k++)
        {
            if (Img_gray[k] < Img_gray[j])
            {
                temp_value = Img_gray[j];
                Img_gray[j] = Img_gray[k];
                Img_gray[k] = temp_value;
                temp_index = sort_index[j];
                sort_index[j] = sort_index[k];
                sort_index[k] = temp_index;
            }
        }
    }

    /* Generate the rand number */
    srand((unsigned)time(NULL));
    int vol_rand = (Img_Height * Img_Width) % 256;
    unsigned int *rand_index = new unsigned int[vol_rand]();
    for (int p = 0; p < vol_rand; p++)
    {
        rand_index[p] = (rand() % (Img_Height*Img_Width)) + 0;
        Img_gray_equ[rand_index[p]] = (unsigned char)(255);
    }

    /* Distribute gray values */
    unsigned char* std_bf = new unsigned char[Img_Height*Img_Width - vol_rand]();
    for (int r = 0; r < 256; r++)
    {
        for (int q = 0; q < Img_Height*Img_Width / 256; q++)
        {
            std_bf[(Img_Height*Img_Width / 256) * r + q] = r;
        }
    }
    int t = 0;
    for (int s = 0; s < Img_Height*Img_Width; s++)
    {
        if (Img_gray_equ[sort_index[s]] == 0) {
            Img_gray_equ[sort_index[s]] = std_bf[t];
            t += 1;
        }
        else {
            Img_gray_equ[sort_index[s]] = std_bf[t];
        }
    }
}

void Color_Reduction(unsigned long cum_hist[256], unsigned long hist[256], unsigned char* Img_input, unsigned char* Img_reduced, int color_reduced, int Img_Height, int Img_Width)
{
    unsigned int seg_cum = Img_Height * Img_Width / color_reduced;
    unsigned long* seg_value = new unsigned long[color_reduced + 1];
    seg_value[0] = 0;
    seg_value[color_reduced] = 255;
    int j = 1;
    for (int i = 1; j <= color_reduced - 1; i++)
    {
        int a = cum_hist[i - 1];
        if (cum_hist[i - 1] >= seg_cum * j)
        {
            seg_value[j] = i - 1;
            j += 1;
        }
    }
    unsigned long color_reduced_trans_func[256];
    for (int i = 0; i < color_reduced; i++)
    {
        int sum = 0;
        for (int m = seg_value[i] + 1; m <= (int)(seg_value[i + 1]); m++)
        {
            sum += (m * hist[m]);
        }
        int mean = sum / (cum_hist[seg_value[i + 1]] - cum_hist[seg_value[i]]);
        for (int m = seg_value[i] + 1; m <= (int)(seg_value[i + 1]); m++)
        {
            color_reduced_trans_func[m] = mean;
        }
    }
    color_reduced_trans_func[0] = color_reduced_trans_func[1];
    for (int i = 0; i < Img_Height*Img_Width; i++)
    {
        Img_reduced[i] = (unsigned char)(color_reduced_trans_func[(int)(Img_input[i])]);
    }
}

void ExtendBoundary(unsigned char** Img_2d, unsigned char** Img_2d_Extend, int Extend_Scale, int Img_Height, int Img_Width)
{
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            Img_2d_Extend[i + Extend_Scale / 2][j + Extend_Scale / 2] = Img_2d[i][j];
        }
    }
    for (int i = 0; i < (Extend_Scale / 2); i++)
    {
        for (int j = (Extend_Scale / 2); j < (Extend_Scale / 2 + Img_Width); j++)
        {
            Img_2d_Extend[i][j] = Img_2d_Extend[Extend_Scale - 1 - i][j];
        }
    }
    for (int i = Extend_Scale / 2 + Img_Height; i < (Extend_Scale + Img_Height - 1); i++)
    {
        for (int j = (Extend_Scale / 2); j < (Extend_Scale / 2 + Img_Width); j++)
        {
            Img_2d_Extend[i][j] = Img_2d_Extend[Extend_Scale + 2 * Img_Height - 3 - i][j];
        }
    }
    for (int i = 0; i < (Img_Height + Extend_Scale - 1); i++)
    {
        for (int j = 0; j < (Extend_Scale / 2); j++)
        {
            Img_2d_Extend[i][j] = Img_2d_Extend[i][Extend_Scale - 1 - j];
        }
    }
    for (int i = 0; i < (Img_Height + Extend_Scale - 1); i++)
    {
        for (int j = (Extend_Scale / 2 + Img_Width); j < (Extend_Scale + Img_Width - 1); j++)
        {
            Img_2d_Extend[i][j] = Img_2d_Extend[i][Extend_Scale + 2 * Img_Width - 3 - j];
        }
    }
}

void RemoveExtend(unsigned char** Img_2d, unsigned char** Img_2d_Extend, int Extend_Scale, int Img_Height, int Img_Width)
{
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            Img_2d[i][j] = Img_2d_Extend[i + Extend_Scale / 2][j + Extend_Scale / 2];
        }
    }
}

int Mask_Most_Freq(unsigned char** Img_Input, int Mask_Size, unsigned long count[256])
{
    for (int i = 0; i < 256; i++)
    {
        count[i] = 0;
    }
    for (int i = 0; i < Mask_Size; i++)
    {
        for (int j = 0; j < Mask_Size; j++)
        {
            count[(int)(Img_Input[i][j])]++;
        }
    }
    int maxCount = (int)(count[0]);
    int gray_value = 0;
    for (int i = 1; i < 256; i++)
    {
        if ((int)(count[i]) > maxCount)
            maxCount = (int)(count[i]);
    }
    for (int i = 0; i < 256; i++)
    {
        if ((int)(count[i]) == maxCount)
            gray_value = i;
    }
    return gray_value;
}

void OilPaint(unsigned char** Img_input, unsigned char** Img_oil, int Img_Height, int Img_Width, int Mask_Size)
{
    int size_half = Mask_Size / 2;
    unsigned char** Img_mask_2d = NULL;
    Img_mask_2d = SpaceFor2D(Img_mask_2d, Mask_Size, Mask_Size);
    unsigned long hist[256];
    for (int i = 0; i < 256; i++)
    {
        hist[i] = 0;
    }
    for (int i = size_half; i < size_half + Img_Height; i++)
    {
        for (int j = size_half; j < size_half + Img_Width; j++)
        {
            for (int p = -size_half; p < size_half + 1; p++)
            {
                for (int q = -size_half; q < size_half + 1; q++)
                {
                    Img_mask_2d[p + size_half][q + size_half] = Img_input[i + p][j + q];
                }
            }
            int gray_value = Mask_Most_Freq(Img_mask_2d, Mask_Size, hist);
            Img_oil[i][j] = (unsigned char)(gray_value);
        }
    }
}

void SqueezeHist(unsigned char* Img_input, unsigned char* Img_output, int Img_size, int lower_bound, int upper_bound)
{
    double c = (upper_bound - lower_bound + 1) / 256.0;
    for (int i = 0; i < Img_size; i++)
    {
        int value = (int)(lower_bound + c * (1 + Img_input[i]) + (rand() % 101) / 100.0);
        if (value > 255)
        {
            value = 255;
        }
        Img_output[i] = (unsigned char)(value);
    }
}

void Mirror(unsigned char** Img_input, unsigned char** Img_output, int Img_Height, int Img_Width)
{
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            Img_output[i][j] = Img_input[i][Img_Width - 1 - j];
        }
    }
}

void MedianFilter(unsigned char** Img_input, unsigned char** Img_output, int Mask_Size, int Img_Height, int Img_Width)
{
    int size_half = Mask_Size / 2;
    int temp = 0;
    unsigned char* Img_mask_1d = new unsigned char[2 * Mask_Size - 1]();
    for (int i = size_half; i < size_half + Img_Height; i++)
    {
        for (int j = size_half; j < size_half + Img_Width; j++)
        {
            int k = 0;
            for (int p = -size_half; p < size_half + 1; p++)
            {
                Img_mask_1d[k] = Img_input[i + p][j];
                k += 1;
            }
            for (int q = -size_half; q < size_half + 1; q++)
            {
                Img_mask_1d[k] = Img_input[i][j + q];
                k += 1;
            }
            Img_mask_1d[k] = Img_input[i][j];
            for (int s = 0; s < 2 * Mask_Size - 1; s++)
            {
                for (int r = s + 1; r < 2 * Mask_Size - 1; r++)
                {
                    if ((int)(Img_mask_1d[s]) >(int)(Img_mask_1d[r]))
                    {
                        temp = (int)(Img_mask_1d[s]);
                        Img_mask_1d[s] = Img_mask_1d[r];
                        Img_mask_1d[r] = (unsigned char)(temp);
                    }
                }
            }
            Img_output[i - size_half][j - size_half] = Img_mask_1d[Mask_Size];
        }
    }
}

double CalPSNR(unsigned char* Img_denoise, unsigned char* Img_std, int Img_size)
{
    double MSE = 0.0;
    for (int i = 0; i < Img_size; i++)
    {
        MSE += pow((double)(Img_denoise[i] - Img_std[i]), 2.0);
    }
    MSE = MSE / Img_size;
    double PSNR = 10.0 * log10(255 * 255 / MSE);

    return PSNR;
}

void BilateralFilter(unsigned char** Img_input, unsigned char** Img_output, int Img_Height, int Img_Width, double sigma_color, double sigma_dist)
{
    int mask_size = (int)(6 * sigma_dist + 1);
    int half_size = (int)(3 * sigma_dist);
    double exp_pow_dist = -0.5 / (sigma_dist*sigma_dist);
    double exp_pow_color = -0.5 / (sigma_color*sigma_color);
    double **mask_dist = NULL;
    mask_dist = new double *[mask_size]();
    for (int i = 0; i < mask_size; i++)
    {
        mask_dist[i] = new double[mask_size]();
        for (int j = 0; j < mask_size; j++)
        {
            mask_dist[i][j] = exp(exp_pow_dist * (pow(half_size - i, 2) + pow(half_size - j, 2)));
        }
    }

    double** color_dist = NULL;
    color_dist = new double *[256]();
    for (int i = 0; i < 256; i++)
    {
        color_dist[i] = new double[256]();
        for (int j = 0; j < 256; j++)
        {
            color_dist[i][j] = exp(exp_pow_color * pow(abs(i - j), 2));
        }
    }

    double **mask_bilateral_2d = NULL;
    mask_bilateral_2d = new double *[mask_size]();
    for (int i = 0; i < mask_size; i++)
    {
        mask_bilateral_2d[i] = new double[mask_size]();
        for (int j = 0; j < mask_size; j++)
        {
            mask_bilateral_2d[i][j] = 0;
        }
    }
    unsigned char** Img_mask_2d = NULL;
    Img_mask_2d = SpaceFor2D(Img_mask_2d, mask_size, mask_size);
    for (int i = half_size; i < half_size + Img_Height; i++) {
        for (int j = half_size; j < half_size + Img_Width; j++) {
            double sum = 0;
            double weight_sum = 0;
            for (int p = -half_size; p < half_size + 1; p++) {
                for (int q = -half_size; q < half_size + 1; q++) {
                    Img_mask_2d[p + half_size][q + half_size] = Img_input[i + p][j + q];
                    mask_bilateral_2d[p + half_size][q + half_size] = mask_dist[p + half_size][q + half_size] *
                                                                      color_dist[Img_input[i + p][j +
                                                                                                  q]][Img_input[i][j]];
                }
            }
            for (int s = 0; s < mask_size; s++) {
                for (int r = 0; r < mask_size; r++) {
                    sum += (double) (Img_mask_2d[s][r]) * mask_bilateral_2d[s][r];
                    weight_sum += mask_bilateral_2d[s][r];
                }
            }
            Img_output[i - half_size][j - half_size] = (unsigned char) (sum / weight_sum + 0.5);
        }
    }
}

int Interpolation_1(unsigned char** Img_input_2d, double ii, double jj)
{
    int i_1 = (int)(ii);
    int i_2 = i_1 + 1;
    int j_1 = (int)(jj);
    int j_2 = j_1 + 1;
    double i_del = ii - i_1;
    double j_del = jj - j_1;

    int a = (int)Img_input_2d[i_1][j_1];
    int b = (int)Img_input_2d[i_1][j_2];
    int c = (int)Img_input_2d[i_2][j_1];
    int d = (int)Img_input_2d[i_2][j_2];
    double pixel_val = (a*(1 - j_del) + b * j_del)*(1 - i_del) + (c*(1 - j_del) + d * j_del)*i_del;
    int pixel = (int)(pixel_val + 0.5);
    return pixel;
}

void SquMapDisc(unsigned char** Img_input_2d, unsigned char** Img_output_2d, int Img_Height, int Img_Width)
{
    double u = 0.0;
    double v = 0.0;
    double x = 0.0;
    double y = 0.0;
    double ii = 0;
    double jj = 0;
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            u = i - ((double)(Img_Height-1))/2;
            v = j - ((double)(Img_Width-1))/2;
            if (u * u + v * v > (double)(Img_Height-1)*(Img_Width-1)/4) {
                continue;
            } else {
                if (abs(u) > abs(v)) {
                    x = sqrt(u * u + v * v)*(abs(u)/u);
                    y = abs(x * v / u)*(abs(v)/v);
                } else {
                    y = sqrt(u * u + v * v)*(abs(v)/v);
                    x = abs(y * u / v)*(abs(u)/u);
                }
            }
            ii = x + ((double)(Img_Height-1))/2;
            jj = y + ((double)(Img_Width-1))/2;
            Img_output_2d[i][j] = (unsigned char)Interpolation_1(Img_input_2d, ii, jj);
        }
    }
}

void SquMapDisc_1(unsigned char** Img_input_2d, unsigned char** Img_output_2d, int Img_Height, int Img_Width)
{
    double u = 0.0;
    double v = 0.0;
    double x = 0.0;
    double y = 0.0;
    double ii = 0;
    double jj = 0;
    double ee = pow(2.0,-126);
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            u = i - ((double)(Img_Height-1))/2;
            v = j - ((double)(Img_Width-1))/2;
            if (u * u + v * v > (double)(Img_Height-1)*(Img_Width-1)/4) {
                continue;
            } else {
                if (abs(u) > abs(v)) {
                    x = sqrt(u * u + v * v)*(abs(u)/u);
                    y = sqrt(u * u + v * v)*atan(v/abs(u))*M_2_PI*2;
                } else {
                    y = sqrt(u * u + v * v)*(abs(v)/v);
                    x = sqrt(u * u + v * v)*atan(u/abs(v)+ee)*M_2_PI*2;
                }
            }
            ii = x + ((double)(Img_Height-1))/2;
            jj = y + ((double)(Img_Width-1))/2;
            Img_output_2d[i][j] = (unsigned char)Interpolation_1(Img_input_2d, ii, jj);
        }
    }
}

void SquMapDisc_2(unsigned char** Img_input_2d, unsigned char** Img_output_2d, int Img_Height, int Img_Width)
{
    double u = 0.0;
    double v = 0.0;
    double x = 0.0;
    double y = 0.0;
    double ii = 0.0;
    double jj = 0.0;
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            u = 2*i/(double)(Img_Height-1) - 1;
            v = 2*j/(double)(Img_Width-1) - 1;
            if (u * u + v * v > 1) {
                continue;
            } else {
                x = (0.5*sqrt(2+u*u-v*v+2*sqrt(2)*u)-0.5*sqrt(2+u*u-v*v-2*sqrt(2)*u))*(double)(Img_Height-1)/2;
                y = (0.5*sqrt(2-u*u+v*v+2*sqrt(2)*v)-0.5*sqrt(2-u*u+v*v-2*sqrt(2)*v))*(double)(Img_Width-1)/2;
            }
            ii = x + ((double)(Img_Height-1))/2;
            jj = y + ((double)(Img_Width-1))/2;
            Img_output_2d[i][j] = (unsigned char)Interpolation_1(Img_input_2d, ii, jj);
        }
    }
}

void SquMapDisc_3(unsigned char** Img_input_2d, unsigned char** Img_output_2d, int Img_Height, int Img_Width)
{
    double u = 0.0;
    double v = 0.0;
    double x = 0.0;
    double y = 0.0;
    double ii = 0.0;
    double jj = 0.0;
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            u = i - ((double)(Img_Height-1))/2;
            v = j - ((double)(Img_Width-1))/2;
            if (u * u + v * v > (double)(Img_Height-1)*(Img_Width-1)/4) {
                continue;
            } else {
                x = u;
                y = (Img_Height-1)*v/(2*sqrt((double)(Img_Height-1)*(Img_Width-1)/4 - u*u));
            }
            ii = x + ((double)(Img_Height-1))/2;
            jj = y + ((double)(Img_Width-1))/2;
            Img_output_2d[i][j] = (unsigned char)Interpolation_1(Img_input_2d, ii, jj);
        }
    }
}

void DiscMapSqu_2(unsigned char** Img_input_2d, unsigned char** Img_output_2d, int Img_Height, int Img_Width)
{
    double u = 0.0;
    double v = 0.0;
    double x = 0.0;
    double y = 0.0;
    double ii = 0.0;
    double jj = 0.0;
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            x = 2*i/(double)(Img_Height-1) - 1;
            y = 2*j/(double)(Img_Width-1) - 1;
            u = x*sqrt(1-(y*y)/2)*(double)(Img_Height-1)/2;
            v = y*sqrt(1-(x*x)/2)*(double)(Img_Width-1)/2;
            ii = u + ((double)(Img_Height-1))/2;
            jj = v + ((double)(Img_Width-1))/2;
            Img_output_2d[i][j] = (unsigned char)Interpolation_1(Img_input_2d, ii, jj);
        }
    }
}

void DiscMapSqu_3(unsigned char** Img_input_2d, unsigned char** Img_output_2d, int Img_Height, int Img_Width)
{
    double u = 0.0;
    double v = 0.0;
    double x = 0.0;
    double y = 0.0;
    double ii = 0.0;
    double jj = 0.0;
    for (int i = 0; i < Img_Height; i++) {
        for (int j = 0; j < Img_Width; j++) {
            x = i - ((double)(Img_Height-1))/2;
            y = j - ((double)(Img_Width-1))/2;
            u = x;
            v = y*sqrt((double)(Img_Height-1)*(Img_Width-1)/4-x*x)/sqrt((double)(Img_Height-1)*(Img_Width-1)/4);
            ii = u + ((double) (Img_Height - 1)) / 2;
            jj = v + ((double) (Img_Width - 1)) / 2;

            if (ii >= Img_Height-1)
            {
                ii -= 1;
            }
            if (jj >= Img_Width-1)
            {
                jj -= 1;
            }
            Img_output_2d[i][j] = (unsigned char) Interpolation_1(Img_input_2d, ii, jj);
        }
    }
}

void FillBackground(unsigned char** Img_Left, unsigned char** Img_Back, unsigned char** Img_Flag, int Height_mov, int Width_mov, int Img_Height, int Img_Width)
{
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            Img_Back[i+Height_mov][j+Width_mov] = Img_Left[i][j];
            Img_Flag[i+Height_mov][j+Width_mov] = 1;
        }
    }
}

void Stitch_Left(unsigned char** Img_Input_2d, unsigned char** Img_Output_2d_Pan1, unsigned char** Img_Flag, int Img_Height, int Img_Width, int Pan_Height, int Pan_Width)
{
    float y1 = 564;
    float x1 = 213;

    float y2 = 656;
    float x2 = 388;

    float y3 = 681;
    float x3 = 650;

    float y4 = 560;
    float x4 = 744;

    float Y1 = 318;
    float X1 = 40;

    float Y2 = 403;
    float X2 = 208;

    float Y3 = 430;
    float X3 = 485;

    float Y4 = 313;
    float X4 = 557;

    float aa[8][8] = {{x1, y1, 1, 0,  0,  0, -x1 * X1, -y1 * X1},
                      {x2, y2, 1, 0,  0,  0, -x2 * X2, -y2 * X2},
                      {x3, y3, 1, 0,  0,  0, -x3 * X3, -y3 * X3},
                      {x4, y4, 1, 0,  0,  0, -x4 * X4, -y4 * X4},
                      {0,  0,  0, x1, y1, 1, -x1 * Y1, -y1 * Y1},
                      {0,  0,  0, x2, y2, 1, -x2 * Y2, -y2 * Y2},
                      {0,  0,  0, x3, y3, 1, -x3 * Y3, -y3 * Y3},
                      {0,  0,  0, x4, y4, 1, -x4 * Y4, -y4 * Y4}};

    float bb[8][1] = {{X1},
                      {X2},
                      {X3},
                      {X4},
                      {Y1},
                      {Y2},
                      {Y3},
                      {Y4}};

    Mat A = Mat(8, 8, CV_32FC1, aa);
    Mat B = Mat(8, 1, CV_32FC1, bb);

    Mat p(8, 1, CV_32FC1);

    solve(A, B, p);

    float X, Y, x, y;
    int pixel_value;
    float a = p.at<float>(0,0);
    float b = p.at<float>(1,0);
    float c = p.at<float>(2,0);
    float d = p.at<float>(3,0);
    float e = p.at<float>(4,0);
    float f = p.at<float>(5,0);
    float g = p.at<float>(6,0);
    float h = p.at<float>(7,0);

    for (int i = 0; i < Pan_Height; i++)
    {
        for (int j = 0; j < Pan_Width; j++)
        {
            x = (float)i;
            y = (float)j;
            X = (a * x + b * y + c) / (g * x + h * y + 1);
            Y = (d * x + e * y + f) / (g * x + h * y + 1);
            if ((X < Img_Height-1) && (Y < Img_Width-1) && (X > 0) && (Y > 0))
            {
                pixel_value = Interpolation_1(Img_Input_2d, (double)X, (double)Y);
                if (Img_Flag[i][j] == 1) {
                    Img_Output_2d_Pan1[i][j] = 0.5*(pixel_value+Img_Output_2d_Pan1[i][j]);
                } else {
                    Img_Output_2d_Pan1[i][j] = pixel_value;
                    Img_Flag[i][j] = 1;
                }
            }
        }
    }
}

void Stitch_Right(unsigned char** Img_Input_2d, unsigned char** Img_Output_2d_Pan2, unsigned char** Img_Flag, int Img_Height, int Img_Width, int Pan_Height, int Pan_Width)
{
    float y1 = 762;
    float x1 = 422;

    float y2 = 732;
    float x2 = 549;

    float y3 = 762;
    float x3 = 719;

    float y4 = 938;
    float x4 = 750;

    float Y1 = 51;
    float X1 = 239;

    float Y2 = 18;
    float X2 = 374;

    float Y3 = 55;
    float X3 = 546;

    float Y4 = 229;
    float X4 = 554;

    float aa[8][8] = {{x1, y1, 1, 0,  0,  0, -x1 * X1, -y1 * X1},
                      {x2, y2, 1, 0,  0,  0, -x2 * X2, -y2 * X2},
                      {x3, y3, 1, 0,  0,  0, -x3 * X3, -y3 * X3},
                      {x4, y4, 1, 0,  0,  0, -x4 * X4, -y4 * X4},
                      {0,  0,  0, x1, y1, 1, -x1 * Y1, -y1 * Y1},
                      {0,  0,  0, x2, y2, 1, -x2 * Y2, -y2 * Y2},
                      {0,  0,  0, x3, y3, 1, -x3 * Y3, -y3 * Y3},
                      {0,  0,  0, x4, y4, 1, -x4 * Y4, -y4 * Y4}};

    float bb[8][1] = {{X1},
                      {X2},
                      {X3},
                      {X4},
                      {Y1},
                      {Y2},
                      {Y3},
                      {Y4}};

    Mat A = Mat(8, 8, CV_32FC1, aa);
    Mat B = Mat(8, 1, CV_32FC1, bb);

    Mat p(8, 1, CV_32FC1);

    solve(A, B, p);

    float X, Y, x, y;
    int pixel_value;
    float a = p.at<float>(0,0);
    float b = p.at<float>(1,0);
    float c = p.at<float>(2,0);
    float d = p.at<float>(3,0);
    float e = p.at<float>(4,0);
    float f = p.at<float>(5,0);
    float g = p.at<float>(6,0);
    float h = p.at<float>(7,0);

    for (int i = 0; i < Pan_Height; i++)
    {
        for (int j = 0; j < Pan_Width; j++)
        {
            x = (float)i;
            y = (float)j;
            X = (a * x + b * y + c) / (g * x + h * y + 1);
            Y = (d * x + e * y + f) / (g * x + h * y + 1);
            if ((X < Img_Height-1) && (Y < Img_Width-1) && (X > 0) && (Y > 0))
            {
                pixel_value = Interpolation_1(Img_Input_2d, (double)X, (double)Y);
                if (Img_Flag[i][j] == 1) {
                    Img_Output_2d_Pan2[i][j] = 0.5*(pixel_value+Img_Output_2d_Pan2[i][j]);
                } else {
                    Img_Output_2d_Pan2[i][j] = pixel_value;
                    Img_Flag[i][j] = 1;
                }
            }
        }
    }
}

Mat SloveAxB(float x1, float x2, float x3, float x4, float X1, float X2, float X3, float X4, float y1, float y2, float y3, float y4, float Y1, float Y2, float Y3, float Y4)
{
    float aa[8][8] = {{x1, y1, 1, 0,  0,  0, -x1 * X1, -y1 * X1},
                  {x2, y2, 1, 0,  0,  0, -x2 * X2, -y2 * X2},
                  {x3, y3, 1, 0,  0,  0, -x3 * X3, -y3 * X3},
                  {x4, y4, 1, 0,  0,  0, -x4 * X4, -y4 * X4},
                  {0,  0,  0, x1, y1, 1, -x1 * Y1, -y1 * Y1},
                  {0,  0,  0, x2, y2, 1, -x2 * Y2, -y2 * Y2},
                  {0,  0,  0, x3, y3, 1, -x3 * Y3, -y3 * Y3},
                  {0,  0,  0, x4, y4, 1, -x4 * Y4, -y4 * Y4}};

    float bb[8][1] = {{X1},
                  {X2},
                  {X3},
                  {X4},
                  {Y1},
                  {Y2},
                  {Y3},
                  {Y4}};

    Mat A = Mat(8, 8, CV_32FC1, aa);
    Mat B = Mat(8, 1, CV_32FC1, bb);

    Mat p(8, 1, CV_32FC1);

    solve(A, B, p);

    return p;

}

void Stitch_Left_5Points(unsigned char** Img_Input_2d, unsigned char** Img_Output_2d_Pan2, unsigned char** Img_Flag, int Img_Height, int Img_Width, int Pan_Height, int Pan_Width)
{
    float y1 = 564;
    float x1 = 213;

    float y2 = 656;
    float x2 = 388;

    float y3 = 681;
    float x3 = 650;

    float y4 = 560;
    float x4 = 744;

    float y5 = 619;
    float x5 = 648;

    float Y1 = 318;
    float X1 = 40;

    float Y2 = 403;
    float X2 = 208;

    float Y3 = 430;
    float X3 = 485;

    float Y4 = 313;
    float X4 = 557;

    float Y5 = 363;
    float X5 = 472;

    Mat p1 = SloveAxB(x1, x2, x3, x4, X1, X2, X3, X4, y1, y2, y3, y4, Y1, Y2, Y3, Y4);
    Mat p2 = SloveAxB(x5, x2, x3, x4, X5, X2, X3, X4, y5, y2, y3, y4, Y5, Y2, Y3, Y4);
    Mat p3 = SloveAxB(x1, x5, x3, x4, X1, X5, X3, X4, y1, y5, y3, y4, Y1, Y5, Y3, Y4);
    Mat p4 = SloveAxB(x1, x2, x5, x4, X1, X2, X5, X4, y1, y2, y5, y4, Y1, Y2, Y5, Y4);
    Mat p5 = SloveAxB(x1, x2, x3, x4, X5, X2, X3, X5, y1, y2, y3, y5, Y1, Y2, Y3, Y5);

    Mat p = p1+p2+p3+p4+p5;

    float X, Y, x, y;
    int pixel_value;
    float a = p.at<float>(0,0)/5;
    float b = p.at<float>(1,0)/5;
    float c = p.at<float>(2,0)/5;
    float d = p.at<float>(3,0)/5;
    float e = p.at<float>(4,0)/5;
    float f = p.at<float>(5,0)/5;
    float g = p.at<float>(6,0)/5;
    float h = p.at<float>(7,0)/5;

    for (int i = 0; i < Pan_Height; i++)
    {
        for (int j = 0; j < Pan_Width; j++)
        {
            x = (float)i;
            y = (float)j;
            X = (a * x + b * y + c) / (g * x + h * y + 1);
            Y = (d * x + e * y + f) / (g * x + h * y + 1);
            if ((X < Img_Height-1) && (Y < Img_Width-1) && (X > 0) && (Y > 0))
            {
                pixel_value = Interpolation_1(Img_Input_2d, (double)X, (double)Y);
                if (Img_Flag[i][j] == 1) {
                    Img_Output_2d_Pan2[i][j] = 0.5*(pixel_value+Img_Output_2d_Pan2[i][j]);
                } else {
                    Img_Output_2d_Pan2[i][j] = pixel_value;
                    Img_Flag[i][j] = 1;
                }
            }
        }
    }
}

void StitchMultiPoints(unsigned char** Img_Input_2d_left, unsigned char** Img_Input_2d_right, unsigned char** Img_Output_2d_Pan, unsigned char** Img_Flag, int Img_Height, int Img_Width, int Pan_Height, int Pan_Width)
{
    float a = 0.6713;
    float b = -0.1867;
    float c = -11.1897;
    float d = 0.0250;
    float e = 0.4054;
    float f = -29.0901;
    float g = 0.000071919;
    float h = -0.00063355;

    float X, Y, x, y;
    int pixel_value;

    for (int i = 0; i < Pan_Height; i++)
    {
        for (int j = 0; j < Pan_Width; j++)
        {
            x = (float)i;
            y = (float)j;
            X = (a * x + b * y + c) / (g * x + h * y + 1);
            Y = (d * x + e * y + f) / (g * x + h * y + 1);
            if ((X < Img_Height-1) && (Y < Img_Width-1) && (X > 0) && (Y > 0))
            {
                pixel_value = Interpolation_1(Img_Input_2d_left, (double)X, (double)Y);
                if (Img_Flag[i][j] == 1) {
                    Img_Output_2d_Pan[i][j] = 0.5*(pixel_value+Img_Output_2d_Pan[i][j]);
                } else {
                    Img_Output_2d_Pan[i][j] = pixel_value;
                    Img_Flag[i][j] = 1;
                }
            }
        }
    }

    a = 1.6996;
    b = 0.1400;
    c = -458.6130;
    d = 0.0242;
    e = 1.7202;
    f = -1241.4;
    g = 0.00016441;
    h = 0.00061866;

    for (int i = 0; i < Pan_Height; i++)
    {
        for (int j = 0; j < Pan_Width; j++)
        {
            x = (float)i;
            y = (float)j;
            X = (a * x + b * y + c) / (g * x + h * y + 1);
            Y = (d * x + e * y + f) / (g * x + h * y + 1);
            if ((X < Img_Height-1) && (Y < Img_Width-1) && (X > 0) && (Y > 0))
            {
                pixel_value = Interpolation_1(Img_Input_2d_right, (double)X, (double)Y);
                if (Img_Flag[i][j] == 1) {
                    Img_Output_2d_Pan[i][j] = 0.5*(pixel_value+Img_Output_2d_Pan[i][j]);
                } else {
                    Img_Output_2d_Pan[i][j] = pixel_value;
                    Img_Flag[i][j] = 1;
                }
            }
        }
    }
}

void Dither_F_Thre(unsigned char** Img_Input_2d, unsigned char** Img_Output_2d, int Img_Height, int Img_Width)
{
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            if ((int)Img_Input_2d[i][j] < 127.5) {
                Img_Output_2d[i][j] = 0;
            } else {
                Img_Output_2d[i][j] = 255;
            }
        }
    }
}

void Dither_F_Thre_4(unsigned char** Img_Input_2d, unsigned char** Img_Output_2d, int Img_Height, int Img_Width)
{
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            if ((int)Img_Input_2d[i][j] < 45) {
                Img_Output_2d[i][j] = 0;
            } else if ((int)Img_Input_2d[i][j] < 85) {
                Img_Output_2d[i][j] = 85;
            } else if ((int)Img_Input_2d[i][j] < 120) {
                Img_Output_2d[i][j] = 170;
            } else {
                Img_Output_2d[i][j] = 255;
            }
        }
    }
}

void Dither_R_Thre(unsigned char** Img_Input_2d, unsigned char** Img_Output_2d, int Img_Height, int Img_Width)
{
    int irand;
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            irand = rand() % 256;
            if ((int)Img_Input_2d[i][j] < irand) {
                Img_Output_2d[i][j] = 0;
            } else {
                Img_Output_2d[i][j] = 255;
            }
        }
    }
}

void Dither_R_Thre_4(unsigned char** Img_Input_2d, unsigned char** Img_Output_2d, int Img_Height, int Img_Width)
{
    int irand_1;
    int irand_2;
    int irand_3;
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            irand_1 = rand() % 60;
            irand_2 = rand() % 40 + 60;
            irand_3 = rand() % 156 + 100;

            if ((int)Img_Input_2d[i][j] < irand_1) {
                Img_Output_2d[i][j] = 0;
            } else if ((int)Img_Input_2d[i][j] < irand_2) {
                Img_Output_2d[i][j] = 85;
            } else if ((int)Img_Input_2d[i][j] < irand_3) {
                Img_Output_2d[i][j] = 170;
            } else {
                Img_Output_2d[i][j] = 255;
            }
        }
    }
}

void Normalize(unsigned char** Img_Input_2d, double** Img_Input_nor_2d, int Img_Height, int Img_Width)
{
    int c;
    double m;
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            c = (int)(Img_Input_2d[i][j]);
            m = (double)(c)/255.0;
            Img_Input_nor_2d[i][j] = m;
        }
    }
}

void Dither_D_Mat(double** Img_Input_nor_2d, unsigned char** Img_Output_2d, int Img_Height, int Img_Width, int size)
{
    if (size == 2)
    {
        double M[2][2] = {{0, 2},
                          {3, 1}};
        for (int i = 0 ; i < Img_Height; i++)
        {
            for (int j = 0; j < Img_Width; j++)
            {
                if (Img_Input_nor_2d[i][j] > (M[i%2][j%2]+0.5)/4)
                {
                    Img_Output_2d[i][j] = 255;
                } else {
                    Img_Output_2d[i][j] = 0;
                }
            }
        }
    }
    if (size == 4)
    {
        double M[4][4] = {{0, 8, 2, 10},
                          {12, 4, 14, 6},
                          {3, 11, 1, 9},
                          {15, 7, 13, 5}};
        for (int i = 0 ; i < Img_Height; i++)
        {
            for (int j = 0; j < Img_Width; j++)
            {
                if (Img_Input_nor_2d[i][j] > (M[i%4][j%4]+0.5)/16)
                {
                    Img_Output_2d[i][j] = 255;
                } else {
                    Img_Output_2d[i][j] = 0;
                }
            }
        }
    }
    if (size == 8)
    {
        double M[8][8] = {{0, 48, 12, 60, 3, 51, 15, 63},
                          {32, 16, 44, 28, 35, 19, 47, 31},
                          {8, 56, 4, 52, 11, 59, 7, 55},
                          {40, 24, 36, 20, 43, 27, 39, 23},
                          {2, 50, 14, 62, 1, 49, 13, 61},
                          {34, 18, 46, 30, 33, 17, 45, 29},
                          {10, 58, 6, 54, 9, 57, 5, 53},
                          {42, 26, 38, 22, 41, 25, 37, 21}};
        for (int i = 0 ; i < Img_Height; i++)
        {
            for (int j = 0; j < Img_Width; j++)
            {
                if (Img_Input_nor_2d[i][j] > (M[i%8][j%8]+0.5)/64)
                {
                    Img_Output_2d[i][j] = 255;
                } else {
                    Img_Output_2d[i][j] = 0;
                }
            }
        }
    }
}

void Dither_D_Mat_4(unsigned char** Img_Input_2d, unsigned char** Img_Output_2d, int Img_Height, int Img_Width, int size)
{
    double pixel;

    if (size == 2)
    {
        double M[2][2] = {{0, 2},
                          {3, 1}};
        for (int i = 0 ; i < Img_Height; i++)
        {
            for (int j = 0; j < Img_Width; j++)
            {
                if (Img_Input_2d[i][j] < 60) {
                    pixel = (double)((int)Img_Input_2d[i][j])/60;
                    if (pixel > (M[i % 2][j % 2] + 0.5) / 4) {
                        Img_Output_2d[i][j] = 85;
                    } else {
                        Img_Output_2d[i][j] = 0;
                    }
                } else if (Img_Input_2d[i][j] < 100) {
                    pixel = (double)((int)Img_Input_2d[i][j]-60)/40;
                    if (pixel > (M[i % 2][j % 2] + 0.5) / 4) {
                        Img_Output_2d[i][j] = 170;
                    } else {
                        Img_Output_2d[i][j] = 85;
                    }
                } else {
                    pixel = (double)((int)Img_Input_2d[i][j]-100)/155;
                    if (pixel > (M[i % 2][j % 2] + 0.5) / 4) {
                        Img_Output_2d[i][j] = 255;
                    } else {
                        Img_Output_2d[i][j] = 170;
                    }
                }
            }
        }
    }
    if (size == 4)
    {
        double M[4][4] = {{0, 8, 2, 10},
                          {12, 4, 14, 6},
                          {3, 11, 1, 9},
                          {15, 7, 13, 5}};
        for (int i = 0 ; i < Img_Height; i++)
        {
            for (int j = 0; j < Img_Width; j++)
            {
                if (Img_Input_2d[i][j] < 60) {
                    pixel = (double)((int)Img_Input_2d[i][j])/60;
                    if (pixel > (M[i % 4][j % 4] + 0.5) / 16) {
                        Img_Output_2d[i][j] = 85;
                    } else {
                        Img_Output_2d[i][j] = 0;
                    }
                } else if (Img_Input_2d[i][j] < 100) {
                    pixel = (double)((int)Img_Input_2d[i][j]-60)/40;
                    if (pixel > (M[i % 4][j % 4] + 0.5) / 16) {
                        Img_Output_2d[i][j] = 170;
                    } else {
                        Img_Output_2d[i][j] = 85;
                    }
                } else {
                    pixel = (double)((int)Img_Input_2d[i][j]-100)/155;
                    if (pixel > (M[i % 4][j % 4] + 0.5) / 16) {
                        Img_Output_2d[i][j] = 255;
                    } else {
                        Img_Output_2d[i][j] = 170;
                    }
                }
            }
        }
    }
    if (size == 8)
    {
        double M[8][8] = {{0, 48, 12, 60, 3, 51, 15, 63},
                          {32, 16, 44, 28, 35, 19, 47, 31},
                          {8, 56, 4, 52, 11, 59, 7, 55},
                          {40, 24, 36, 20, 43, 27, 39, 23},
                          {2, 50, 14, 62, 1, 49, 13, 61},
                          {34, 18, 46, 30, 33, 17, 45, 29},
                          {10, 58, 6, 54, 9, 57, 5, 53},
                          {42, 26, 38, 22, 41, 25, 37, 21}};
        for (int i = 0 ; i < Img_Height; i++)
        {
            for (int j = 0; j < Img_Width; j++)
            {
                if (Img_Input_2d[i][j] < 60) {
                    pixel = (double)((int)Img_Input_2d[i][j])/60;
                    if (pixel > (M[i % 8][j % 8] + 0.5) / 64) {
                        Img_Output_2d[i][j] = 85;
                    } else {
                        Img_Output_2d[i][j] = 0;
                    }
                } else if (Img_Input_2d[i][j] < 100) {
                    pixel = (double)((int)Img_Input_2d[i][j]-60)/40;
                    if (pixel > (M[i % 8][j % 8] + 0.5) / 64) {
                        Img_Output_2d[i][j] = 170;
                    } else {
                        Img_Output_2d[i][j] = 85;
                    }
                } else {
                    pixel = (double)((int)Img_Input_2d[i][j]-100)/155;
                    if (pixel > (M[i % 8][j % 8] + 0.5) / 64) {
                        Img_Output_2d[i][j] = 255;
                    } else {
                        Img_Output_2d[i][j] = 170;
                    }
                }
            }
        }
    }
}

void CovertCharToDouble(unsigned char** Img_Input_2d, double** Img_Output_2d, int Img_Height, int Img_Width)
{
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            Img_Output_2d[i][j] = (double)Img_Input_2d[i][j];
        }
    }
}

void FsErrDiffusion(double** Img_Input_2d, unsigned char** Img_Output_2d, int Img_Height, int Img_Width)
{
    float fs[3][3] = {{0, 0, 0},
                      {0, 0, 7},
                      {3, 5, 1}};
    Mat FS = Mat(3, 3, CV_32FC1, fs);

    Mat ERR = Mat::ones(Size(3,3), CV_32FC1);

    float err = 0;

    for (int i = 1; i < Img_Height+1; i++)
    {
        for (int j = 1; j < Img_Width+1; j++)
        {
            Mat ERR = Mat::ones(Size(3,3), CV_32FC1);
            double a = Img_Input_2d[i][j];
            if (Img_Input_2d[i][j] < 127) {
                Img_Output_2d[i][j] = 0;
            } else {
                Img_Output_2d[i][j] = 255;
            }
            err = Img_Input_2d[i][j] - Img_Output_2d[i][j];
            ERR = ERR * err;
            ERR = ERR .mul (FS) / 16;
            for (int p = -1; p <= 1; p++)
            {
                for (int q = -1; q <= 1; q++)
                {
                    Img_Input_2d[i+p][j+q] = (double)ERR.at<float>(p+1,q+1) + Img_Input_2d[i+p][j+q];
                }
            }
        }
    }
}

void MyErrDiffusion(double** Img_Input_2d, unsigned char** Img_Output_2d, int Img_Height, int Img_Width)
{
    float fs[3][3] = {{0, 0, 0},
                      {0, 0, 2},
                      {1, 2, 1}};
    Mat FS = Mat(3, 3, CV_32FC1, fs);

    Mat ERR = Mat::ones(Size(3,3), CV_32FC1);

    float err = 0;

    for (int i = 1; i < Img_Height+1; i++)
    {
        for (int j = 1; j < Img_Width+1; j++)
        {
            Mat ERR = Mat::ones(Size(3,3), CV_32FC1);
            double a = Img_Input_2d[i][j];
            if (Img_Input_2d[i][j] < 127) {
                Img_Output_2d[i][j] = 0;
            } else {
                Img_Output_2d[i][j] = 255;
            }
            err = Img_Input_2d[i][j] - Img_Output_2d[i][j];
            ERR = ERR * err;
            ERR = ERR .mul (FS) / 6;
            for (int p = -1; p <= 1; p++)
            {
                for (int q = -1; q <= 1; q++)
                {
                    Img_Input_2d[i+p][j+q] = (double)ERR.at<float>(p+1,q+1) + Img_Input_2d[i+p][j+q];
                }
            }
        }
    }
}

int Rendering(double** Img_Red, double** Img_Green, double** Img_Blue, int i, int j)
{
    double R = Img_Red[i][j];
    double G = Img_Green[i][j];
    double B = Img_Blue[i][j];

    if (R+G > 255) {
        if (G+B > 255) {
            if (R+G+B > 510) {
                return 1;
            } else {
                return 2;
            }
        } else {
            return 3;
        }
    } else if (G+B <= 255) {
        if (R+G+B <= 255) {
            return 4;
        } else {
            return 5;
        }
    } else {
        return 6;
    }
}

double Dist_2(double a, double b, double c, double R, double G, double B)
{
    double dist = pow(a-R, 2) + pow(b-G, 2) + pow(c-B, 2);
    return dist;
}

int MinPos(double a, double b, double c, double d)
{
    double array[] = { a, b, c, d};
    int option = distance(array, min_element(array, array + sizeof(array)/sizeof(array[0])));

    return option;
}

void CMYW(double** Img_Red, double** Img_Green, double** Img_Blue,
          unsigned char** Img_Output_Red, unsigned char** Img_Output_Green, unsigned char** Img_Output_Blue,
          int i, int j)
{
    double R = Img_Red[i][j];
    double G = Img_Green[i][j];
    double B = Img_Blue[i][j];

    double dist_c = Dist_2(0.0, 255.0, 255.0, R, G, B);
    double dist_m = Dist_2(255.0, 0.0, 255.0, R, G, B);
    double dist_y = Dist_2(255.0, 255.0, 0.0, R, G, B);
    double dist_w = Dist_2(255.0, 255.0, 255.0, R, G, B);

    int option = MinPos(dist_c, dist_m, dist_y, dist_w);

    if (option == 0)
    {
        Img_Output_Red[i][j] = 0;
        Img_Output_Green[i][j] = 255;
        Img_Output_Blue[i][j] = 255;
    }

    if (option == 1)
    {
        Img_Output_Red[i][j] = 255;
        Img_Output_Green[i][j] = 0;
        Img_Output_Blue[i][j] = 255;
    }

    if (option == 2)
    {
        Img_Output_Red[i][j] = 255;
        Img_Output_Green[i][j] = 255;
        Img_Output_Blue[i][j] = 0;
    }

    if (option == 3)
    {
        Img_Output_Red[i][j] = 255;
        Img_Output_Green[i][j] = 255;
        Img_Output_Blue[i][j] = 255;
    }
}

void MYGC(double** Img_Red, double** Img_Green, double** Img_Blue,
          unsigned char** Img_Output_Red, unsigned char** Img_Output_Green, unsigned char** Img_Output_Blue,
          int i, int j)
{
    double R = Img_Red[i][j];
    double G = Img_Green[i][j];
    double B = Img_Blue[i][j];

    double dist_m = Dist_2(255.0, 0.0, 255.0, R, G, B);
    double dist_y = Dist_2(255.0, 255.0, 0.0, R, G, B);
    double dist_g = Dist_2(0.0, 255.0, 0.0, R, G, B);
    double dist_c = Dist_2(0.0, 255.0, 255.0, R, G, B);

    int option = MinPos(dist_m, dist_y, dist_g, dist_c);

    if (option == 0)
    {
        Img_Output_Red[i][j] = 255;
        Img_Output_Green[i][j] = 0;
        Img_Output_Blue[i][j] = 255;
    }

    if (option == 1)
    {
        Img_Output_Red[i][j] = 255;
        Img_Output_Green[i][j] = 255;
        Img_Output_Blue[i][j] = 0;
    }

    if (option == 2)
    {
        Img_Output_Red[i][j] = 0;
        Img_Output_Green[i][j] = 255;
        Img_Output_Blue[i][j] = 0;
    }

    if (option == 3)
    {
        Img_Output_Red[i][j] = 0;
        Img_Output_Green[i][j] = 255;
        Img_Output_Blue[i][j] = 255;
    }
}

void RGMY(double** Img_Red, double** Img_Green, double** Img_Blue,
          unsigned char** Img_Output_Red, unsigned char** Img_Output_Green, unsigned char** Img_Output_Blue,
          int i, int j)
{
    double R = Img_Red[i][j];
    double G = Img_Green[i][j];
    double B = Img_Blue[i][j];

    double dist_r = Dist_2(255, 0.0, 0.0, R, G, B);
    double dist_g = Dist_2(0.0, 255.0, 0.0, R, G, B);
    double dist_m = Dist_2(255.0, 0.0, 255.0, R, G, B);
    double dist_y = Dist_2(255.0, 255.0, 0.0, R, G, B);

    int option = MinPos(dist_r, dist_g, dist_m, dist_y);

    if (option == 0)
    {
        Img_Output_Red[i][j] = 255;
        Img_Output_Green[i][j] = 0;
        Img_Output_Blue[i][j] = 0;
    }

    if (option == 1)
    {
        Img_Output_Red[i][j] = 0;
        Img_Output_Green[i][j] = 255;
        Img_Output_Blue[i][j] = 0;
    }

    if (option == 2)
    {
        Img_Output_Red[i][j] = 255;
        Img_Output_Green[i][j] = 0;
        Img_Output_Blue[i][j] = 255;
    }

    if (option == 3)
    {
        Img_Output_Red[i][j] = 255;
        Img_Output_Green[i][j] = 255;
        Img_Output_Blue[i][j] = 0;
    }
}

void KRGB(double** Img_Red, double** Img_Green, double** Img_Blue,
          unsigned char** Img_Output_Red, unsigned char** Img_Output_Green, unsigned char** Img_Output_Blue,
          int i, int j)
{
    double R = Img_Red[i][j];
    double G = Img_Green[i][j];
    double B = Img_Blue[i][j];

    double dist_k = Dist_2(0.0, 0.0, 0.0, R, G, B);
    double dist_r = Dist_2(255, 0.0, 0.0, R, G, B);
    double dist_g = Dist_2(0.0, 255.0, 0.0, R, G, B);
    double dist_b = Dist_2(0.0, 0.0, 255.0, R, G, B);

    int option = MinPos(dist_k, dist_r, dist_g, dist_b);

    if (option == 0)
    {
        Img_Output_Red[i][j] = 0;
        Img_Output_Green[i][j] = 0;
        Img_Output_Blue[i][j] = 0;
    }

    if (option == 1)
    {
        Img_Output_Red[i][j] = 255;
        Img_Output_Green[i][j] = 0;
        Img_Output_Blue[i][j] = 0;
    }

    if (option == 2)
    {
        Img_Output_Red[i][j] = 0;
        Img_Output_Green[i][j] = 255;
        Img_Output_Blue[i][j] = 0;
    }

    if (option == 3)
    {
        Img_Output_Red[i][j] = 0;
        Img_Output_Green[i][j] = 0;
        Img_Output_Blue[i][j] = 255;
    }
}

void RGBM(double** Img_Red, double** Img_Green, double** Img_Blue,
          unsigned char** Img_Output_Red, unsigned char** Img_Output_Green, unsigned char** Img_Output_Blue,
          int i, int j)
{
    double R = Img_Red[i][j];
    double G = Img_Green[i][j];
    double B = Img_Blue[i][j];

    double dist_r = Dist_2(255, 0.0, 0.0, R, G, B);
    double dist_g = Dist_2(0.0, 255.0, 0.0, R, G, B);
    double dist_b = Dist_2(0.0, 0.0, 255.0, R, G, B);
    double dist_m = Dist_2(255.0, 0.0, 255.0, R, G, B);

    int option = MinPos(dist_r, dist_g, dist_b, dist_m);

    if (option == 0)
    {
        Img_Output_Red[i][j] = 255;
        Img_Output_Green[i][j] = 0;
        Img_Output_Blue[i][j] = 0;
    }

    if (option == 1)
    {
        Img_Output_Red[i][j] = 0;
        Img_Output_Green[i][j] = 255;
        Img_Output_Blue[i][j] = 0;
    }

    if (option == 2)
    {
        Img_Output_Red[i][j] = 0;
        Img_Output_Green[i][j] = 0;
        Img_Output_Blue[i][j] = 255;
    }

    if (option == 3)
    {
        Img_Output_Red[i][j] = 255;
        Img_Output_Green[i][j] = 0;
        Img_Output_Blue[i][j] = 255;
    }
}

void CMGB(double** Img_Red, double** Img_Green, double** Img_Blue,
          unsigned char** Img_Output_Red, unsigned char** Img_Output_Green, unsigned char** Img_Output_Blue,
          int i, int j)
{
    double R = Img_Red[i][j];
    double G = Img_Green[i][j];
    double B = Img_Blue[i][j];

    double dist_c = Dist_2(0.0, 255.0, 255.0, R, G, B);
    double dist_m = Dist_2(255.0, 0.0, 255.0, R, G, B);
    double dist_g = Dist_2(0.0, 255.0, 0.0, R, G, B);
    double dist_b = Dist_2(0.0, 0.0, 255.0, R, G, B);

    int option = MinPos(dist_c, dist_m, dist_g, dist_b);

    if (option == 0)
    {
        Img_Output_Red[i][j] = 0;
        Img_Output_Green[i][j] = 255;
        Img_Output_Blue[i][j] = 0;
    }

    if (option == 1)
    {
        Img_Output_Red[i][j] = 255;
        Img_Output_Green[i][j] = 0;
        Img_Output_Blue[i][j] = 255;
    }

    if (option == 2)
    {
        Img_Output_Red[i][j] = 0;
        Img_Output_Green[i][j] = 255;
        Img_Output_Blue[i][j] = 0;
    }

    if (option == 3)
    {
        Img_Output_Red[i][j] = 0;
        Img_Output_Green[i][j] = 0;
        Img_Output_Blue[i][j] = 255;
    }
}

void MBVQ(double** Img_Red, double** Img_Green, double** Img_Blue, unsigned char** Img_Output_Red,
          unsigned char** Img_Output_Green, unsigned char** Img_Output_Blue, int Img_Height, int Img_Width)
{
    float fs[3][3] = {{0, 0, 0},
                      {0, 0, 7},
                      {3, 5, 1}};
    Mat FS = Mat(3, 3, CV_32FC1, fs);

    int option = 0;

    float err_r = 0;
    float err_g = 0;
    float err_b = 0;

    for (int i = 1; i < Img_Height+1; i++)
    {
        for (int j = 1; j < Img_Width+1; j++)
        {
            Mat ERR_R = Mat::ones(Size(3,3), CV_32FC1);
            Mat ERR_G = Mat::ones(Size(3,3), CV_32FC1);
            Mat ERR_B = Mat::ones(Size(3,3), CV_32FC1);
            option = Rendering(Img_Red, Img_Green, Img_Blue, i, j);
            switch (option) {
                case 1: CMYW(Img_Red, Img_Green, Img_Blue, Img_Output_Red, Img_Output_Green, Img_Output_Blue, i, j); break;
                case 2: MYGC(Img_Red, Img_Green, Img_Blue, Img_Output_Red, Img_Output_Green, Img_Output_Blue, i, j); break;
                case 3: RGMY(Img_Red, Img_Green, Img_Blue, Img_Output_Red, Img_Output_Green, Img_Output_Blue, i, j); break;
                case 4: KRGB(Img_Red, Img_Green, Img_Blue, Img_Output_Red, Img_Output_Green, Img_Output_Blue, i, j); break;
                case 5: RGBM(Img_Red, Img_Green, Img_Blue, Img_Output_Red, Img_Output_Green, Img_Output_Blue, i, j); break;
                case 6: CMGB(Img_Red, Img_Green, Img_Blue, Img_Output_Red, Img_Output_Green, Img_Output_Blue, i, j); break;
            }

            err_r = Img_Red[i][j] - Img_Output_Red[i][j];
            ERR_R = ERR_R * err_r;
            ERR_R = ERR_R .mul (FS) / 16;
            for (int p = -1; p <= 1; p++)
            {
                for (int q = -1; q <= 1; q++)
                {
                    Img_Red[i+p][j+q] = (double)ERR_R.at<float>(p+1,q+1) + Img_Red[i+p][j+q];
                }
            }

            err_g = Img_Green[i][j] - Img_Output_Green[i][j];
            ERR_G = ERR_G * err_g;
            ERR_G = ERR_G .mul (FS) / 16;
            for (int p = -1; p <= 1; p++)
            {
                for (int q = -1; q <= 1; q++)
                {
                    Img_Green[i+p][j+q] = (double)ERR_G.at<float>(p+1,q+1) + Img_Green[i+p][j+q];
                }
            }

            err_b = Img_Blue[i][j] - Img_Output_Blue[i][j];
            ERR_B = ERR_B * err_b;
            ERR_B = ERR_B .mul (FS) / 16;
            for (int p = -1; p <= 1; p++)
            {
                for (int q = -1; q <= 1; q++)
                {
                    Img_Blue[i+p][j+q] = (double)ERR_B.at<float>(p+1,q+1) + Img_Blue[i+p][j+q];
                }
            }
        }
    }
}

void JJNErrDiffusion(double** Img_Input_2d, unsigned char** Img_Output_2d, int Img_Height, int Img_Width)
{
    float jjn[5][5] = {{0, 0, 0, 0, 0},
                       {0, 0, 0, 0, 0},
                       {0, 0, 0, 7, 5},
                       {3, 5, 7, 5, 3},
                       {1, 3, 5, 3, 1}};
    Mat JJN = Mat(5, 5, CV_32FC1, jjn);

    Mat ERR = Mat::ones(Size(5,5), CV_32FC1);

    float err = 0;

    for (int i = 2; i < Img_Height+2; i++)
    {
        for (int j = 2; j < Img_Width+2; j++)
        {
            Mat ERR = Mat::ones(Size(5,5), CV_32FC1);
            if (Img_Input_2d[i][j] < 127) {
                Img_Output_2d[i][j] = 0;
            } else {
                Img_Output_2d[i][j] = 255;
            }
            err = Img_Input_2d[i][j] - Img_Output_2d[i][j];
            ERR = ERR * err;
            ERR = ERR .mul (JJN) / 48;
            for (int p = -2; p <= 2; p++)
            {
                for (int q = -2; q <= 2; q++)
                {
                    Img_Input_2d[i+p][j+q] = (double)ERR.at<float>(p+2,q+2) + Img_Input_2d[i+p][j+q];
                }
            }
        }
    }
}

void StErrDiffusion(double** Img_Input_2d, unsigned char** Img_Output_2d, int Img_Height, int Img_Width)
{
    float jjn[5][5] = {{0, 0, 0, 0, 0},
                       {0, 0, 0, 0, 0},
                       {0, 0, 0, 8, 4},
                       {2, 4, 8, 4, 2},
                       {1, 2, 4, 2, 1}};
    Mat JJN = Mat(5, 5, CV_32FC1, jjn);

    Mat ERR = Mat::ones(Size(5,5), CV_32FC1);

    float err = 0;

    for (int i = 2; i < Img_Height+2; i++)
    {
        for (int j = 2; j < Img_Width+2; j++)
        {
            Mat ERR = Mat::ones(Size(5,5), CV_32FC1);
            if (Img_Input_2d[i][j] < 127) {
                Img_Output_2d[i][j] = 0;
            } else {
                Img_Output_2d[i][j] = 255;
            }
            err = Img_Input_2d[i][j] - Img_Output_2d[i][j];
            ERR = ERR * err;
            ERR = ERR .mul (JJN) / 48;
            for (int p = -2; p <= 2; p++)
            {
                for (int q = -2; q <= 2; q++)
                {
                    Img_Input_2d[i+p][j+q] = (double)ERR.at<float>(p+2,q+2) + Img_Input_2d[i+p][j+q];
                }
            }
        }
    }
}

void Binary2Bool(unsigned char** Img_Input, int** Img_Output, int Img_Height, int Img_Width)
{
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            Img_Output[i][j] = ((int)Img_Input[i][j])/255;
        }
    }
}

void Bool2Binary(int** Img_Input, unsigned char** Img_Output, int Img_Height, int Img_Width)
{
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            Img_Output[i][j] = (unsigned char)(Img_Input[i][j]*255);
        }
    }
}

void Shrinking(unsigned char** Img_Input, unsigned char** Img_Output, int Img_Height, int Img_Width)
{
    /* Conditional mark pattern */
    int mark_s_1[54] = {1000000,10000,100,1,
                        10000000,100000,1000,10,
                        11000000,1100000,110000,11000,1100,110,11,10000001,
                        11000001,1110000,11100,111,
                        10110000,10100001,1101000,11000010,
                        11100000,111000,1110,10000011,
                        10110001,1101100,
                        11110000,11100001,1111000,111100,11110,1111,10000111,11000011,
                        11110001,1111100,11111,11000111,
                        11100011,11111000,111110,10001111,
                        11110011,11100111,11111100,11111001,111110,11111,10011111,11001111};
    /* Unconditional mark pattern*/
    /*int mark_s_2[102] = {1000000,10000,
                         10,10000000,
                         11000000,1100000,110000,11000,1100,110,11,10000001,
                         1101000,10110000,10100001,11000010,
                         1100100,11000100,11100100,11001,110001,111001,1001100,1000110,1001110,10011,10010001,10010011,
                         111000,11111111,
                         10101000,10111100,11101001,10001010,11001011,10011110,101010,1111010,101111,10100010,10100111,11110010,
                         11101001,10001010,101010,10100010,
                         1010001,1010010,1010011,1010100,1010101,1010110,1010111,11111001,11111010,11111011,11111100,11111101,11111110,11111111,
                         1010100,10010100,11010100,10101,1010101,10010101,11010101,1111110,10111110,11111110,111111,1111111,10111111,11111111,
                         10101,100101,110101,1000101,1010101,1100101,1110101,10011111,10101111,10111111,11001111,11011111,11101111,11111111,
                         1000101,1001001,1001101,1010001,1010101,1011001,1011101,11100111,11101011,11101111,11110011,11110111,11111011,11111111};*/

    int mark_s_2[108] = {1000000,10000,
                         10,10000000,
                         11000000,1100000,110000,11000,1100,110,11,10000001,
                         1101000,10110000,10100001,11000010,
                         1100100,11000100,11100100,11001,110001,111001,1001100,1000110,1001110,10011,10010001,10010011,
                         111000,11111111,
                         10101000,10111100,11101001,10001010,11001011,10011110,101010,1111010,101111,10100010,10100111,11110010,
                         11101001,10001010,101010,10100010,
                         1010001,1010010,1010011,1010101,1010110,1010111,11111010,11111011,11111100,11111101,11111110,11111111,
                         10010100,11010100,10101,1010101,10010101,11010101,10111110,11111110,111111,1111111,10111111,11111111,
                         100101,110101,1000101,1010101,1100101,1110101,10101111,10111111,11001111,11011111,11101111,11111111,
                         1001001,1001101,1010001,1010101,1011001,1011101,11101011,11101111,11110011,11110111,11111011,11111111,
                         10110101,10100100,1101101,101001,1011011,1001010,10010010,11010110,11110000,10010111,101110,1111,1011,10110100};

    int B_mask[Img_Height][Img_Width];
    int M_mask[Img_Height][Img_Width];

    int b_pixel = 0;
    int m_pixel = 0;

    /* Transfer unsigned char array into bool array */
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            B_mask[i][j] = ((int)Img_Input[i][j])/255;
        }
    }

    for (int t = 0; t < 40; t++)
    {
        memset(M_mask, 0, sizeof(M_mask));    // Initialize M_mask

        /* 1_stage match */
        for (int i = 1; i < Img_Height-1; i++)
        {
            for (int j = 1; j < Img_Width-1; j++)
            {
                if (B_mask[i][j] == 1)
                {
                    b_pixel = 10000000* B_mask[i][j+1]+
                              1000000* B_mask[i-1][j+1]+
                              100000* B_mask[i-1][j]+
                              10000* B_mask[i-1][j-1]+
                              1000* B_mask[i][j-1]+
                              100 * B_mask[i+1][j-1]+
                              10* B_mask[i+1][j]+
                              B_mask[i+1][j+1];
                    for (int k = 0; k < 54; k++)
                    {
                        if (b_pixel == mark_s_1[k])
                            M_mask[i][j] = 1;
                    }
                }
            }
        }

        for (int p = 1; p < Img_Height-1; p++)
        {
            for (int q = 1; q < Img_Width-1; q++)
            {
                if (M_mask[p][q] == 1)
                {
                    m_pixel = 10000000* M_mask[p][q+1]+
                              1000000* M_mask[p-1][q+1]+
                              100000* M_mask[p-1][q]+
                              10000* M_mask[p-1][q-1]+
                              1000* M_mask[p][q-1]+
                              100 * M_mask[p+1][q-1]+
                              10* M_mask[p+1][q]+
                              M_mask[p+1][q+1];
                    for (int k = 0; k < 108; k++)
                    {
                        if (m_pixel == mark_s_2[k])
                        {
                            B_mask[p][q] = 1;
                            break;
                        } else {
                            B_mask[p][q] = 0;
                        }
                    }
                }
            }
        }
    }

    /* Transfer bool array into unsigned char array */
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            Img_Output[i][j] = (unsigned char)(B_mask[i][j]*255);
        }
    }

}

int CheckEightConnected(int** Img, int i, int j)
{
    int flag = 0;
    for (int p = -1; p < 2; p++)
    {
        for (int q = -1; q < 2; q++)
        {
            if ((Img[i+p][j+q] == 1) && (p != 0 || p != 0))
            {
                flag = 1;
                return flag;
            }
        }
    }
    return flag;
}

void InitializeInt(int** Mask, int Height, int Width)
{
    for (int i = 0; i < Height; i++)
    {
        for (int j = 0; j < Width; j++)
        {
            Mask[i][j] = 0;
        }
    }
}

unsigned char** ShrinkHist(unsigned char** Img_S, unsigned char** Img_B, int Img_Height, int Img_Width, int num)
{
    int pixel[num*2];
    int index = 0;

    unsigned char** Img_C = NULL;
    Img_C = SpaceFor2D(Img_C, Img_Height, Img_Width);

    for (int i = 0; i < Img_Height; i++) {
        for (int j = 0; j < Img_Width; j++) {
            if (Img_S[i][j] == 255)
            {
                pixel[index] = i;
                pixel[index + 1] = j;
                index += 2;
            }
        }
    }

    /* Conditional mark pattern */
    int mark_s_1[54] = {1000000,10000,100,1,
                        10000000,100000,1000,10,
                        11000000,1100000,110000,11000,1100,110,11,10000001,
                        11000001,1110000,11100,111,
                        10110000,10100001,1101000,11000010,
                        11100000,111000,1110,10000011,
                        10110001,1101100,
                        11110000,11100001,1111000,111100,11110,1111,10000111,11000011,
                        11110001,1111100,11111,11000111,
                        11100011,11111000,111110,10001111,
                        11110011,11100111,11111100,11111001,111110,11111,10011111,11001111};
    /* Unconditional mark pattern*/
    int mark_s_2[108] = {1000000,10000,
                         10,10000000,
                         11000000,1100000,110000,11000,1100,110,11,10000001,
                         1101000,10110000,10100001,11000010,
                         1100100,11000100,11100100,11001,110001,111001,1001100,1000110,1001110,10011,10010001,10010011,
                         111000,11111111,
                         10101000,10111100,11101001,10001010,11001011,10011110,101010,1111010,101111,10100010,10100111,11110010,
                         11101001,10001010,101010,10100010,
                         1010001,1010010,1010011,1010101,1010110,1010111,11111010,11111011,11111100,11111101,11111110,11111111,
                         10010100,11010100,10101,1010101,10010101,11010101,10111110,11111110,111111,1111111,10111111,11111111,
                         100101,110101,1000101,1010101,1100101,1110101,10101111,10111111,11001111,11011111,11101111,11111111,
                         1001001,1001101,1010001,1010101,1011001,1011101,11101011,11101111,11110011,11110111,11111011,11111111,
                         10110101,10100100,1101101,101001,1011011,1001010,10010010,11010110,11110000,10010111,101110,1111,1011,10110100};

    int** B_mask = NULL;
    B_mask = SpaceFor2DInt(B_mask, Img_Height, Img_Width);
    int** M_mask = NULL;
    M_mask = SpaceFor2DInt(M_mask, Img_Height, Img_Width);

    int b_pixel = 0;
    int m_pixel = 0;

    /* Transfer unsigned char array into bool array */
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            B_mask[i][j] = ((int)Img_B[i][j])/255;
        }
    }

    for (int t = 0; t < 20; t++)
    {
        InitializeInt(M_mask, Img_Height, Img_Width);    // Initialize M_mask

        int x = 0;
        while(x < num*2)
        {
            if (CheckEightConnected(B_mask, pixel[x], pixel[x+1]))
            {
                Img_C[pixel[x]][pixel[x+1]] += 1;
            }
            x += 2;
        }

        /* 1_stage match */
        for (int i = 1; i < Img_Height-1; i++)
        {
            for (int j = 1; j < Img_Width-1; j++)
            {
                if (B_mask[i][j] == 1)
                {
                    b_pixel = 10000000* B_mask[i][j+1]+
                              1000000* B_mask[i-1][j+1]+
                              100000* B_mask[i-1][j]+
                              10000* B_mask[i-1][j-1]+
                              1000* B_mask[i][j-1]+
                              100 * B_mask[i+1][j-1]+
                              10* B_mask[i+1][j]+
                              B_mask[i+1][j+1];
                    for (int k = 0; k < 54; k++)
                    {
                        if (b_pixel == mark_s_1[k])
                            M_mask[i][j] = 1;
                    }
                }
            }
        }

        for (int p = 1; p < Img_Height-2; p++)
        {
            for (int q = 1; q < Img_Width-2; q++)
            {
                if (M_mask[p][q] == 1)
                {
                    m_pixel = 10000000* M_mask[p][q+1]+
                              1000000* M_mask[p-1][q+1]+
                              100000* M_mask[p-1][q]+
                              10000* M_mask[p-1][q-1]+
                              1000* M_mask[p][q-1]+
                              100 * M_mask[p+1][q-1]+
                              10* M_mask[p+1][q]+
                              M_mask[p+1][q+1];
                    for (int k = 0; k < 108; k++)
                    {
                        if (m_pixel == mark_s_2[k])
                        {
                            B_mask[p][q] = 1;
                            break;
                        } else {
                            B_mask[p][q] = 0;
                        }
                    }
                }
            }
        }

    }

    return Img_C;

}

int CountNum(unsigned char** Img, int Img_Height, int Img_Width)
{
    int num = 0;
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            if (Img[i][j] == 255)
                num += 1;
        }
    }
    cout << "The total number of stars is " << num << endl;
    return num;
}

void SizeHist2d(unsigned int hist[20], unsigned char **Img_C, int Img_Height, int Img_Width)
{
    for (int i = 0; i < 20; i++)
    {
        hist[i] = 0;
    }

    for (int j = 0; j < Img_Height; j++)
    {
        for (int k = 0; k < Img_Width; k++)
        {
            if (Img_C[j][k] > 0)
                hist[(int)(Img_C[j][k])] += 1;
        }
    }
}

void SizeHist_Write(char* path, unsigned int Hist[20])
{
    ofstream file(path);
    if (!file)
    {
        cout << "Fail to write into " << path << endl;
        exit(1);
    }

    for (int i = 0; i < 20; i++)
    {
        file << Hist[i] << endl;
    }
    file.close();
    cout << "Writing " << path << "... Successfully" << endl;
}

void Thining(unsigned char** Img_Input, unsigned char** Img_Output, int Img_Height, int Img_Width)
{
    /* Conditional mark pattern */
    int mark_s_1[46] = {10100000,101000,1010,10000010,
                        11000001,1110000,11100,111,
                        10110000,10100001,1101000,11000010,
                        11100000,111000,1110,10000011,
                        10110001,1101100,
                        11110000,11100001,1111000,111100,11110,1111,10000111,11000011,
                        11110001,1111100,11111,11000111,
                        11100011,11111000,111110,10001111,
                        11110011,11100111,11111100,11111001,111110,11111,10011111,11001111,
                        11110111,11111101,111111,11011111};
    /* Unconditional mark pattern*/
    /*int mark_s_2[102] = {1000000,10000,
                     10,10000000,
                     11000000,1100000,110000,11000,1100,110,11,10000001,
                     1101000,10110000,10100001,11000010,
                     1100100,11000100,11100100,11001,110001,111001,1001100,1000110,1001110,10011,10010001,10010011,
                     111000,11111111,
                     10101000,10111100,11101001,10001010,11001011,10011110,101010,1111010,101111,10100010,10100111,11110010,
                     11101001,10001010,101010,10100010,
                     1010001,1010010,1010011,1010100,1010101,1010110,1010111,11111001,11111010,11111011,11111100,11111101,11111110,11111111,
                     1010100,10010100,11010100,10101,1010101,10010101,11010101,1111110,10111110,11111110,111111,1111111,10111111,11111111,
                     10101,100101,110101,1000101,1010101,1100101,1110101,10011111,10101111,10111111,11001111,11011111,11101111,11111111,
                     1000101,1001001,1001101,1010001,1010101,1011001,1011101,11100111,11101011,11101111,11110011,11110111,11111011,11111111};*/
    int mark_s_2[116] = {1000000,10000,
                         10,10000000,
                         11000000,1100000,110000,11000,1100,110,11,10000001,
                         1101000,10110000,10100001,11000010,
                         1100100,11000100,11100100,11001,110001,111001,1001100,1000110,1001110,10011,10010001,10010011,
                         111000,11111111,
                         10101000,10111100,11101001,10001010,11001011,10011110,101010,1111010,101111,10100010,10100111,11110010,
                         11101001,10001010,101010,10100010,
                         1010001,1010010,1010011,1010101,1010110,1010111,11111010,11111011,11111100,11111101,11111110,11111111,
                         10010100,11010100,10101,1010101,10010101,11010101,10111110,11111110,111111,1111111,10111111,11111111,
                         100101,110101,1000101,1010101,1100101,1110101,10101111,10111111,11001111,11011111,11101111,11111111,
                         1001001,1001101,1010001,1010101,1011001,1011101,11101011,11101111,11110011,11110111,11111011,11111111,
                         10110101,10100100,1101101,101001,1011011,1001010,10010010,11010110,11111000,11110000,11111001,1111001,
                         11010011,11100001,11010010,11110100,1111000,11100101,11100111,11100011,11000011,10010111};

    int B_mask[Img_Height][Img_Width];
    int M_mask[Img_Height][Img_Width];

    int b_pixel = 0;
    int m_pixel = 0;

    /* Transfer unsigned char array into bool array */
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            B_mask[i][j] = 1-(((int)Img_Input[i][j])/255);
        }
    }

    for (int t = 0; t < 25; t++)
    {
        memset(M_mask, 0, sizeof(M_mask));    // Initialize M_mask

        /* 1_stage match */
        for (int i = 1; i < Img_Height-2; i++)
        {
            for (int j = 1; j < Img_Width-2; j++)
            {
                if (B_mask[i][j] == 1)
                {
                    b_pixel = 10000000* B_mask[i][j+1]+
                              1000000* B_mask[i-1][j+1]+
                              100000* B_mask[i-1][j]+
                              10000* B_mask[i-1][j-1]+
                              1000* B_mask[i][j-1]+
                              100 * B_mask[i+1][j-1]+
                              10* B_mask[i+1][j]+
                              B_mask[i+1][j+1];
                    for (int k = 0; k < 46; k++)
                    {
                        if (b_pixel == mark_s_1[k])
                            M_mask[i][j] = 1;
                    }
                }
            }
        }

        for (int p = 1; p < Img_Height-2; p++)
        {
            for (int q = 1; q < Img_Width-2; q++)
            {
                if (M_mask[p][q] == 1)
                {
                    m_pixel = 10000000* M_mask[p][q+1]+
                              1000000* M_mask[p-1][q+1]+
                              100000* M_mask[p-1][q]+
                              10000* M_mask[p-1][q-1]+
                              1000* M_mask[p][q-1]+
                              100 * M_mask[p+1][q-1]+
                              10* M_mask[p+1][q]+
                              M_mask[p+1][q+1];
                    for (int k = 0; k < 116; k++)
                    {
                        if (m_pixel == mark_s_2[k])
                        {
                            B_mask[p][q] = 1;
                            break;
                        } else {
                            B_mask[p][q] = 0;
                        }
                    }
                }
            }
        }
    }

    /* Transfer bool array into unsigned char array */
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            Img_Output[i][j] = (unsigned char)(B_mask[i][j]*255);
        }
    }

}

void Skeleton(unsigned char** Img_Input, unsigned char** Img_Output, int Img_Height, int Img_Width)
{
    /* Conditional mark pattern */
    int mark_s_1[40] = {10100000,101000,1010,10000010,
                        11000001,1110000,11100,111,
                        11110000,11100001,1111000,111100,11110,1111,10000111,11000011,
                        11110001,1111100,11111,11000111,
                        11100011,11111000,111110,10001111,
                        11110011,11100111,11111100,11111001,111110,11111,10011111,11001111,
                        11110111,11111101,111111,11011111,
                        11111011,11111110,10111111,11101111};
    /* Unconditional mark pattern*/
    int mark_s_2[88] = {1,100,1000000,10000,
                         10,10000000,1000,100000,
                         10100000,101000,10000010,1010,
                         111000,11111111,11111111,10000011,
                         10101000,11111111,101010,11111111,10001010,11111111,10100010,11111111,
                         1010001,1010010,1010011,1010100,1010101,1010110,1010111,11111001,11111010,11111011,11111100,11111101,11111110,11111111,
                         1010100,10010100,11010100,10101,1010101,10010101,11010101,1111110,10111110,11111110,111111,1111111,10111111,11111111,
                         10101,100101,110101,1000101,1010101,1100101,1110101,10011111,10101111,10111111,11001111,11011111,11101111,11111111,
                         1000101,1001001,1001101,1010001,1010101,1011001,1011101,11100111,11101011,11101111,11110011,11110111,11111011,11111111,
                         10100100,10110101,101001,1101101,1001010,1011011,10010010,11010110};

    int B_mask[Img_Height][Img_Width];
    int M_mask[Img_Height][Img_Width];

    int b_pixel = 0;
    int m_pixel = 0;

    /* Transfer unsigned char array into bool array */
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            B_mask[i][j] = 1-(((int)Img_Input[i][j])/255);
        }
    }

    for (int t = 0; t < 40; t++)
    {
        memset(M_mask, 0, sizeof(M_mask));    // Initialize M_mask

        /* 1_stage match */
        for (int i = 1; i < Img_Height-2; i++)
        {
            for (int j = 1; j < Img_Width-2; j++)
            {
                if (B_mask[i][j] == 1)
                {
                    b_pixel = 10000000* B_mask[i][j+1]+
                              1000000* B_mask[i-1][j+1]+
                              100000* B_mask[i-1][j]+
                              10000* B_mask[i-1][j-1]+
                              1000* B_mask[i][j-1]+
                              100 * B_mask[i+1][j-1]+
                              10* B_mask[i+1][j]+
                              B_mask[i+1][j+1];
                    for (int k = 0; k < 40; k++)
                    {
                        if (b_pixel == mark_s_1[k])
                            M_mask[i][j] = 1;
                    }
                }
            }
        }

        for (int p = 1; p < Img_Height-2; p++)
        {
            for (int q = 1; q < Img_Width-2; q++)
            {
                if (M_mask[p][q] == 1)
                {
                    m_pixel = 10000000* M_mask[p][q+1]+
                              1000000* M_mask[p-1][q+1]+
                              100000* M_mask[p-1][q]+
                              10000* M_mask[p-1][q-1]+
                              1000* M_mask[p][q-1]+
                              100 * M_mask[p+1][q-1]+
                              10* M_mask[p+1][q]+
                              M_mask[p+1][q+1];
                    for (int k = 0; k < 88; k++)
                    {
                        if (m_pixel == mark_s_2[k])
                        {
                            B_mask[p][q] = 1;
                            break;
                        } else {
                            B_mask[p][q] = 0;
                        }
                    }
                }
            }
        }
    }

    /* Transfer bool array into unsigned char array */
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            Img_Output[i][j] = (unsigned char)(B_mask[i][j]*255);
        }
    }

}

void Inverse(unsigned char** Img_Input, unsigned char** Img_Output, int Img_Height, int Img_Width)
{
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            Img_Output[i][j] = 255 - Img_Input[i][j];
        }
    }
}

unsigned char** Erosion(unsigned char** Img, int Img_Height, int Img_Width)
{
    int flag;

    unsigned char **temp = NULL;
    temp = SpaceFor2D(temp, Img_Height, Img_Width);
    unsigned char **Img_Output = NULL;
    Img_Output = SpaceFor2D(Img_Output, Img_Height, Img_Width);
    for (int i =0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            temp[i][j] = Img[i][j];
            Img_Output[i][j] = Img[i][j];
        }
    }

    for (int i = 1; i < Img_Height-1; i++)
    {
        for (int j = 1; j < Img_Width-1; j++)
        {
            flag = 1;
            for (int p = -1; p < 2; p++)
            {
                for (int q = -1; q < 2; q++)
                {
                    if (temp[i+p][j+q] == 0)
                    {
                        flag = 0;
                        break;
                    }
                }
                if (flag == 0)
                    break;
            }
            if (flag == 0) {
                Img_Output[i][j] = 0;
            } else {
                Img_Output[i][j] = 255;
            }
        }
    }
    free(temp);
    return Img_Output;
}

unsigned char** Dilation(unsigned char** Img, int Img_Height, int Img_Width)
{
    int flag;

    unsigned char **temp = NULL;
    temp = SpaceFor2D(temp, Img_Height, Img_Width);
    unsigned char **Img_Output = NULL;
    Img_Output = SpaceFor2D(Img_Output, Img_Height, Img_Width);
    for (int i =0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            temp[i][j] = Img[i][j];
            Img_Output[i][j] = Img[i][j];
        }
    }

    for (int i = 1; i < Img_Height-1; i++)
    {
        for (int j = 1; j < Img_Width-1; j++)
        {
            flag = 1;
            for (int p = -1; p < 2; p++)
            {
                for (int q = -1; q < 2; q++)
                {
                    if (temp[i+p][j+q] == 255)
                    {
                        flag = 0;
                        break;
                    }
                }
                if (flag == 0)
                    break;
            }
            if (flag == 0) {
                Img_Output[i][j] = 255;
            } else {
                Img_Output[i][j] = 0;
            }
        }
    }
    free(temp);
    return Img_Output;
}

int CheckNull(unsigned char** Img, int Img_Height, int Img_Width)
{
    int flag = 1;

    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            if (Img[i][j] == 255)
            {
                flag = 0;
                return flag;
            }
        }
    }
    return flag;
}

unsigned char** Union(unsigned char** Img_a, unsigned char** Img_b, int Img_Height, int Img_Width)
{
    unsigned char **Img_U = NULL;
    Img_U = SpaceFor2D(Img_U, Img_Height, Img_Width);

    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            Img_U[i][j] = (unsigned char)(max((int)(Img_a[i][j]),(int)(Img_b[i][j])));
        }
    }

    return Img_U;
}

unsigned char** Deduct(unsigned char** Img_a, unsigned char** Img_b, int Img_Height, int Img_Width)
{
    int val;

    unsigned char **Img_D = NULL;
    Img_D = SpaceFor2D(Img_D, Img_Height, Img_Width);

    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            val = (int)Img_a[i][j] - (int)Img_b[i][j];
            Img_D[i][j] = (unsigned char)(max(0, val));
        }
    }

    return Img_D;
}

unsigned char** Skeleton_2(unsigned char** Img, int Img_Height, int Img_Width)
{
    unsigned char** Img_2 = NULL;
    Img_2 = SpaceFor2D(Img_2, Img_Height, Img_Width);

    unsigned char** Img_i = NULL;
    Img_i = SpaceFor2D(Img_i, Img_Height, Img_Width);

    unsigned char** Img_new = NULL;
    Img_new = SpaceFor2D(Img_new, Img_Height, Img_Width);

    Img = Erosion(Img, Img_Height, Img_Width);

    Img_2 = Erosion(Img, Img_Height, Img_Width);
    Img_2 = Dilation(Img_2, Img_Height, Img_Width);

    Img_i = Deduct(Img, Img_2, Img_Height, Img_Width);

    while(CheckNull(Img, Img_Height, Img_Width) == 0)
    {

        Img_new = Union(Img_new, Img_i, Img_Height, Img_Width);

        Img = Erosion(Img, Img_Height, Img_Width);

        Img_2 = Erosion(Img, Img_Height, Img_Width);
        Img_2 = Dilation(Img_2, Img_Height, Img_Width);

        Img_i = Deduct(Img, Img_2, Img_Height, Img_Width);
    }

    return Img_new;
}

void TraceLabel(int x, int y, int* parent)
{
    int i = x;
    int j = y;

    while(0 != parent[i])
        i = parent[i];
    while(0 != parent[j])
        j = parent[j];
    if(i != j)
        parent[i] = j;
}

void TwoPass(unsigned char** Img, unsigned char** Label, int Img_Height, int Img_Width)
{
    int max_size = Img_Height * Img_Width;
    int *parent = new int[max_size]();
    memset(parent, 0, max_size * sizeof(int));
    int label_small;
    int label_big;
    int label_temp;
    int label_last = 0;
    int label_1;
    int label_2;
    int label_3;
    int label_4;

    /* First Pass */
    for (int i = 1; i < Img_Height-1; i++)
    {
        for (int j = 1; j < Img_Width-1; j++)
        {
            if (Img[i][j] == 255)
            {
                label_small = 0;
                label_big = 0;
                if (Img[i-1][j+1] != 0 || Img[i-1][j] != 0 || Img[i-1][j-1] != 0)
                {
                    if (Img[i-1][j+1] == 255) {
                        label_small = Label[i-1][j+1];
                    } else if (Img[i-1][j] == 255) {
                        label_small = Label[i-1][j];
                    } else if (Img[i-1][j-1] == 255) {
                        label_small = Label[i-1][j-1];
                    }
                    if (Img[i][j-1] == 255) {
                        label_big = Label[i][j-1];
                        if (label_big < label_small)
                        {
                            label_temp = label_small;
                            label_small = label_big;
                            label_big = label_temp;
                        }
                        Label[i][j] = label_small;
                        TraceLabel(label_big, label_small, parent);
                    } else {
                        Label[i][j] = label_small;
                    }
                } else if (Img[i][j-1] == 255) {
                    Label[i][j] = Label[i][j-1];
                } else {
                    label_last += 1;
                    Label[i][j] = label_last;
                }
            }
        }
    }

    int label;

    /* Second Pass */
    for (int i = 1; i < Img_Height-1; i++)
    {
        for (int j = 1; j < Img_Width - 1; j++)
        {
            if (Label[i][j] != 0)
            {
                label = Label[i][j];
                while(0 != parent[label])
                    label = parent[label];
                Label[i][j] = label;
            }
        }
    }
}

void PosConnectedRegion(unsigned char** Img_Input, int Img_Height, int Img_Width, int pos[4])
{
    int flag = 0;
    int val;
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            if (Img_Input[i][j] != 0)
            {
                flag = 1;
                val = Img_Input[i][j];
                break;
            }
        }
        if (flag == 1)
        {
            break;
        }
    }

    int left = Img_Width;
    int right = 0;
    int up = Img_Height;
    int bottom = 0;

    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            if (Img_Input[i][j] == val)
            {
                if (i < up)
                {
                    up = i;
                }
                if (i > bottom)
                {
                    bottom = i;
                }
                if (j < left)
                {
                    left = j;
                }
                if (j > right)
                {
                    right = j;
                }
            }
        }
    }

    pos[0] = up;
    pos[1] = bottom;
    pos[2] = left;
    pos[3] = right;

    for(int i = up; i <= bottom; i++)
    {
        for (int j = left; j <= right; j++)
        {
            Img_Input[i][j] = 0;
        }
    }
}

void Crop(unsigned char** Img_Crop, unsigned char** Img_Input_Connected_2d, int pos[4])
{
    for (int i = pos[0]; i <= pos[1]; i++)
    {
        for (int j = pos[2]; j <= pos[3]; j++)
        {
            if (Img_Input_Connected_2d[i][j] != 0)
            Img_Crop[i-pos[0]][j-pos[2]] = 255;
        }
    }
}

int JigsawFeature(unsigned char** Img_Crop, int Img_Height, int Img_Width)
{
    /* Right */
    int num_right = 0;
    int flag_right = 0;
    for (int i = 0; i < Img_Height; i++) {
        if (Img_Crop[i][Img_Width - 1] != 0) {
            num_right += 1;
        }
    }
    if (num_right > 40) {
        flag_right = 2;
    } else if (num_right > 20) {
        flag_right = 1;
    }

    /* Top */
    int num_top = 0;
    int flag_top = 0;
    for (int j = 0; j < Img_Width; j++) {
        if (Img_Crop[0][j] != 0) {
            num_top += 1;
        }
    }
    if (num_top > 40) {
        flag_top = 2;
    } else if (num_top > 20) {
        flag_top = 1;
    }

    /* Left */
    int num_left = 0;
    int flag_left = 0;
    for (int i = 0; i < Img_Height; i++) {
        if (Img_Crop[i][0] != 0) {
            num_left += 1;
        }
    }
    if (num_left > 40) {
        flag_left = 2;
    } else if (num_left > 20) {
        flag_left = 1;
    }

    /* Bottom */
    int num_bottom = 0;
    int flag_bottom = 0;
    for (int j = 0; j < Img_Width; j++) {
        if (Img_Crop[Img_Height - 1][j] != 0) {
            num_bottom += 1;
        }
    }
    if (num_bottom > 40) {
        flag_bottom = 2;
    } else if (num_bottom > 20) {
        flag_bottom = 1;
    }

    int feature_1 = 1000 * flag_right + 100 * flag_top + 10 * flag_left + flag_bottom;
    int feature_2 = 1000 * flag_top + 100 * flag_left + 10 * flag_bottom + flag_right;
    int feature_3 = 1000 * flag_left + 100 * flag_bottom + 10 * flag_right + flag_top;
    int feature_4 = 1000 * flag_bottom + 100 * flag_right + 10 * flag_top + flag_left;

    int feature = max(max(max(feature_1, feature_2), feature_3), feature_4);
    return feature;
}

void MatchJigsaw(int* feature, int num)
{
    for (int i = 1; i < num; i++)
    {
        queue<int>match_queue;
        match_queue.push(i+1);
        for (int j = i + 1; j < num; j++)
        {
            if (feature[i] == feature[j] && feature[i] != 3)
            {
                match_queue.push(j+1);
                feature[j] = 3;
            }
        }
        if (match_queue.size() != 1)
        {
            cout << endl << "Group of matched jigsaws: ";
            int queue_size = match_queue.size();
            for (int q = 0; q < queue_size; q++)
            {
                cout << match_queue.front() << ", ";
                match_queue.pop();
            }
        }
    }
    cout << endl;
}

int* SaveCrop(unsigned char** Img_Label_2d, unsigned char** Img_Input_Connected_2d, int Img_Height, int Img_Width, int Img_Channel, int num)
{
    int pos[4];
    int Crop_Height;
    int Crop_Width;
    int *feature = new int[num]();

    for (int i = 1; i <= num; i++)
    {
        PosConnectedRegion(Img_Label_2d, Img_Height, Img_Width, pos);

        Crop_Height = pos[1] - pos[0] + 1;
        Crop_Width = pos[3] - pos[2] + 1;
        unsigned char** Img_Crop = NULL;
        Img_Crop = SpaceFor2D(Img_Crop, Crop_Height, Crop_Width);
        Crop(Img_Crop, Img_Input_Connected_2d, pos);

        unsigned char *Img_Crop_1d = new unsigned char[Img_Channel * Crop_Height * Crop_Width]();
        Img_Crop_1d = Convert2dTo1d(Img_Crop, Crop_Height, Crop_Width);

        string path = "/Users/shanlinsun/Desktop/result//p3d/crop_";
        path = path + to_string(i) + ".raw";

        char *Crop_Write_Path = const_cast<char *>(path.c_str());

        Raw_Write(Crop_Write_Path, Img_Crop_1d, Crop_Height * Crop_Width * Img_Channel);
        cout << "Saving cropped image height: ";
        cout << pos[1] - pos[0] + 1 << endl;
        cout << "Saving cropped image width: ";
        cout << pos[3] - pos[2] + 1 << endl;

        feature[i-1] = JigsawFeature(Img_Crop, Crop_Height, Crop_Width);
    }
    return feature;
}

Mat FilterGenerate_5(int option)
{
    float e5[1][5] = {{-1.0/6, -2.0/6, 0, 2.0/6, 1.0/6}};
    float s5[1][5] = {{-1.0/4, 0, 2.0/4, 0, -1.0/4}};
    float w5[1][5] = {{-1.0/6, 2.0/6, 0, -2.0/6, 1.0/6}};
    Mat E5 = Mat(1, 5, CV_32FC1, e5);
    Mat S5 = Mat(1, 5, CV_32FC1, s5);
    Mat W5 = Mat(1, 5, CV_32FC1, w5);
    Mat M(5, 5, CV_32FC1, Scalar(1));

    switch (option) {
        case 1:
            M = E5.t() * E5;
            break;
        case 2:
            M = E5.t() * S5;
            break;
        case 3:
            M = E5.t() * W5;
            break;
        case 4:
            M = S5.t() * E5;
            break;
        case 5:
            M = S5.t() * S5;
            break;
        case 6:
            M = S5.t() * W5;
            break;
        case 7:
            M = W5.t() * E5;
            break;
        case 8:
            M = W5.t() * S5;
            break;
        case 9:
            M = W5.t() * W5;
            break;
    }

    return M;

}

void FeatureExtraction_5(unsigned char** Img_Input, int Height, int Width, int Extend, vector< vector<double> > &Feature, int img_num)
{
    for (int num = 1; num <= 9; num++)
    {
        Mat M = FilterGenerate_5(num).clone();
        double sum = 0;
        for (int i = Extend/2; i < Extend/2+Height; i++)
        {
            for (int j = Extend/2; j <Extend/2+Width; j++)
            {
                double center_val = 0;
                for (int p = -2; p < 3; p++)
                {
                    for (int q = -2; q < 3; q++) {
                        center_val += (double) M.at<float>(p + 2, q + 2) * (double) Img_Input[i + p][j + q];
                        sum += abs(center_val);
                    }
                }
            }
        }
        Feature.at(img_num-1).at(num-1) = sum / (Height*Width);
    }
}

Mat FilterGenerate_3(int option)
{
    float l3[1][3] = {{1.0/6, 2.0/6, 1.0/6}};
    float e3[1][3] = {{-1.0/2, 0, 1.0/2}};
    float s3[1][3] = {{-1.0/2, 2.0/2, -1.0/2}};
    Mat L3 = Mat(1, 3, CV_32FC1, l3);
    Mat E3 = Mat(1, 3, CV_32FC1, e3);
    Mat S3 = Mat(1, 3, CV_32FC1, s3);
    Mat M(3, 3, CV_32FC1, Scalar(1));

    switch (option) {
        case 1:
            M = L3.t() * L3;
            break;
        case 2:
            M = L3.t() * E3;
            break;
        case 3:
            M = L3.t() * S3;
            break;
        case 4:
            M = E3.t() * L3;
            break;
        case 5:
            M = E3.t() * E3;
            break;
        case 6:
            M = E3.t() * S3;
            break;
        case 7:
            M = S3.t() * L3;
            break;
        case 8:
            M = S3.t() * E3;
            break;
        case 9:
            M = S3.t() * S3;
            break;
    }

    return M;

}

void FeatureExtraction_3(unsigned char** Img_Input, unsigned char** Img_Output, int Height, int Width, int num, vector< vector<double> > &Feature, int Extend, double mean)
{
    Mat M = FilterGenerate_3(num).clone();

    for (int i = Extend/2; i < Extend/2+Height; i++)
    {
        for (int j = Extend/2; j < Extend/2+Width; j++)
        {
            double center_val = 0;
            for (int p = -1; p < 2; p++)
            {
                for (int q = -1; q < 2; q++) {
                    center_val += abs((double) M.at<float>(p + 1, q + 1) * (double) Img_Input[i + p][j + q]);
                }
            }
            Img_Output[i][j] = (unsigned char)(center_val-mean);
        }
    }
}

vector<double> Vector_Dist(vector<double> &x, vector< vector<double> > &mu)
{
    vector<double> dist(mu.size(), 0.0);
    for (int i = 0; i < mu.size(); i++)
    {
        double sum = 0.0;
        for (int j = 0; j < x.size(); j++)
        {
            sum += pow(mu.at(i).at(j)-x.at(j), 2);
        }
        dist.at(i) = sum;
    }
    return dist;
}

void K_Means(vector< vector<double> > &Feature, vector< vector<double> > &Mean_vec, vector<int> &label, int index[4])
{
    /* Initialize mean vector */
    for (int i = 0; i < Mean_vec.size(); i++)
    {
        for (int j = 0; j < Mean_vec.at(0).size(); j++)
            Mean_vec.at(i).at(j) = Feature.at(index[i]).at(j);
    }

    /* Build distance vector */
    vector<double> dist_vec(Feature.size(), 0);

    /* Save the old mean vector */
    vector<vector<double> >Mean_vec_Old(Mean_vec);

    for (int ita = 0; ita < 20 ; ita++)
    {
        for (int j = 0; j < Feature.size(); j++)
        {
            dist_vec = Vector_Dist(Feature.at(j), Mean_vec);
            vector<double>::iterator smallest = min_element(dist_vec.begin(), dist_vec.end());
            int pos = distance(dist_vec.begin(), smallest);
            label.at(j) = pos;
        }

        for (int i = 0; i < Mean_vec.size(); i++)
        {
            int cnt = 0;
            Mean_vec.at(i).assign(Mean_vec.at(i).size(),0);
            for (int p = 0; p < Feature.size(); p++)
            {
                if (label.at(p) == i)
                {
                    cnt += 1;
                    for (int q = 0; q < Feature.at(0).size(); q++)
                    {
                        Mean_vec.at(i).at(q) += Feature.at(p).at(q);
                    }
                }
            }
            for (int x = 0; x < Feature.at(0).size(); x++)
            {
                Mean_vec.at(i).at(x) = Mean_vec.at(i).at(x) / cnt;
            }
        }
        if (Mean_vec_Old == Mean_vec)
        {
            cout << ita << endl;
            break;
        }
    }
}

void AveEngery(unsigned char** Img_Input_Extend_2d, double** Img_Output_2d_D, int Img_Height, int Img_Width, int Extend, int mask_size)
{
    for (int i = Extend/2; i < Extend/2 + Img_Height; i++)
    {
        for (int j = Extend/2; j < Extend/2 + Img_Width; j++)
        {
            double sum = 0.0;
            for (int p = -mask_size/2; p <= mask_size/2; p++)
            {
                for (int q = -mask_size/2; q <= mask_size/2; q++)
                {
                    sum += (double)Img_Input_Extend_2d[i+p][j+q];
                }
            }
            Img_Output_2d_D[i-Extend/2][j-Extend/2] = sum / pow(mask_size,2);
        }
    }
}

double Mean(unsigned char** Img_Input, int Img_Height, int Img_Width)
{
    double sum = 0.0;
    for (int i = 0; i < Img_Height; i++)
    {
        for (int j = 0; j < Img_Width; j++)
        {
            sum += (double)(Img_Input[i][j]);
        }
    }
    return sum/(Img_Height*Img_Width);
}

void NorFeature(vector<vector<double> > &Feature)
{
    for (int i = 0; i < Feature.size(); i++)
    {
        double nor = Feature.at(i).at(0);
        for (int j = 0; j < Feature.at(0).size(); j++)
        {
            Feature.at(i).at(j) = Feature.at(i).at(j)/nor;
        }
    }
}

Mat FilterGenerate_25(int option)
{
    float l5[1][5] = {{1.0/16, 4.0/16, 6.0/16, 4.0/16, 1.0/16}};
    float e5[1][5] = {{-1.0/6, -2.0/6, 0, 2.0/6, 1.0/6}};
    float s5[1][5] = {{-1.0/4, 0, 2.0/4, 0, -1.0/4}};
    float w5[1][5] = {{-1.0/6, 2.0/6, 0, -2.0/6, 1.0/6}};
    float r5[1][5] = {{1.0/16, -4.0/16, 6.0/16, -4.0/16, 1.0/16}};
    Mat L5 = Mat(1, 5, CV_32FC1, l5);
    Mat E5 = Mat(1, 5, CV_32FC1, e5);
    Mat S5 = Mat(1, 5, CV_32FC1, s5);
    Mat W5 = Mat(1, 5, CV_32FC1, w5);
    Mat R5 = Mat(1, 5, CV_32FC1, r5);
    Mat M(5, 5, CV_32FC1, Scalar(1));

    switch (option) {
        case 1:
            M = L5.t() * L5; break;
        case 2:
            M = L5.t() * E5; break;
        case 3:
            M = L5.t() * S5; break;
        case 4:
            M = L5.t() * W5; break;
        case 5:
            M = L5.t() * R5; break;
        case 6:
            M = E5.t() * L5; break;
        case 7:
            M = E5.t() * E5; break;
        case 8:
            M = E5.t() * S5; break;
        case 9:
            M = E5.t() * W5; break;
        case 10:
            M = E5.t() * R5; break;
        case 11:
            M = S5.t() * L5; break;
        case 12:
            M = S5.t() * E5; break;
        case 13:
            M = S5.t() * S5; break;
        case 14:
            M = S5.t() * W5; break;
        case 15:
            M = S5.t() * R5; break;
        case 16:
            M = W5.t() * L5; break;
        case 17:
            M = W5.t() * E5; break;
        case 18:
            M = W5.t() * S5; break;
        case 19:
            M = W5.t() * W5; break;
        case 20:
            M = W5.t() * R5; break;
        case 21:
            M = R5.t() * L5; break;
        case 22:
            M = R5.t() * E5; break;
        case 23:
            M = R5.t() * S5; break;
        case 24:
            M = R5.t() * W5; break;
        case 25:
            M = R5.t() * R5; break;
    }

    return M;

}

void FeatureExtraction_25(unsigned char** Img_Input, unsigned char** Img_Output, int Height, int Width, int num, vector< vector<double> > &Feature, int Extend, double average)
{
    Mat M = FilterGenerate_25(num).clone();
    for (int i = Extend/2; i < Extend/2+Height; i++)
    {
        for (int j = Extend/2; j <Extend/2+Width; j++)
        {
            double center_val = 0;
            for (int p = -2; p < 3; p++)
            {
                for (int q = -2; q < 3; q++) {
                    center_val += abs((double) M.at<float>(p + 2, q + 2) * (double) Img_Input[i + p][j + q]);
                }
            }
            Img_Output[i][j] = (unsigned char)(center_val-average);
        }
    }
}

void Reduce_Dim(vector< vector<double> > &Feature, vector< vector<double> > &Feature_PCA)
{
    Mat frame(Feature.size(), Feature.at(0).size(), CV_64FC1);

    for (int i = 0; i < frame.rows; i++)
        for (int j = 0; j < frame.cols; j++)
            frame.at<double>(i,j) = Feature.at(i).at(j);

    int dim = Feature_PCA.at(0).size();

    PCA map(frame, Mat(), CV_PCA_DATA_AS_ROW, dim);

    Mat dst = map.project(frame);

    for (int i = 0; i < dst.rows; i++)
        for (int j = 0; j < dst.cols; j++)
            Feature_PCA.at(i).at(j) = dst.at<double>(i,j) ;
}

void SobelX(unsigned char** Img_Input, unsigned char** Img_Output, int Img_Height, int Img_Width)
{
    int sobel_x[3][3] = {{-1, 0, 1},{-2, 0, 2},{-1, 0, 1}};

    double** Img_temp = NULL;
    Img_temp = SpaceForNor2D(Img_temp, Img_Height+2, Img_Width+2);

    double min = 10000.0;
    double max = -10000.0;

    for (int i = 1; i < Img_Height + 1; i++)
    {
        for (int j = 1; j < Img_Width + 1; j++)
        {
            for (int p = -1; p < 2; p++)
            {
                for (int q = -1; q < 2; q++)
                {
                    Img_temp[i][j] += (double)(sobel_x[1 + p][1 + q] * Img_Input[i+p][j+q]);
                }
            }
            if (Img_temp[i][j] > max)
                max = Img_temp[i][j];
            if (Img_temp[i][j] < min)
                min = Img_temp[i][j];
        }
    }

    for (int i = 1; i < Img_Height + 1; i++)
    {
        for (int j = 1; j < Img_Width + 1; j++)
        {
            Img_Output[i][j] = (unsigned char)((Img_temp[i][j]-min)*255/(max-min));
        }
    }
}

void SobelY(unsigned char** Img_Input, unsigned char** Img_Output, int Img_Height, int Img_Width)
{
    int sobel_y[3][3] = {{1, 2, 1},{0, 0, 0},{-1, -2, -1}};

    double** Img_temp = NULL;
    Img_temp = SpaceForNor2D(Img_temp, Img_Height+2, Img_Width+2);

    double min = 10000;
    double max = -10000;

    for (int i = 1; i < Img_Height + 1; i++)
    {
        for (int j = 1; j < Img_Width + 1; j++)
        {
            for (int p = -1; p < 2; p++)
            {
                for (int q = -1; q < 2; q++)
                {
                    Img_temp[i][j] += (double)(sobel_y[1 + p][1 + q] * Img_Input[i+p][j+q]);
                }
            }
            if (Img_temp[i][j] > max)
                max = Img_temp[i][j];
            if (Img_temp[i][j] < min)
                min = Img_temp[i][j];
        }
    }

    for (int i = 1; i < Img_Height + 1; i++)
    {
        for (int j = 1; j < Img_Width + 1; j++)
        {
            Img_Output[i][j] = (unsigned char)((Img_temp[i][j]-min)*255/(max-min));
        }
    }
}

void SobelGradient(unsigned char** Img_Output, unsigned char** Img_X, unsigned char** Img_Y, int Height, int Width)
{
    double** Img_temp = NULL;
    Img_temp = SpaceForNor2D(Img_temp, Height, Width);

    double min = 10000;
    double max = -10000;

    for (int i = 0; i < Height; i++)
    {
        for (int j = 0; j < Width; j++)
        {
            Img_temp[i][j] = sqrt(pow((double)Img_X[i][j], 2) + pow((double)Img_Y[i][j], 2));
            if (Img_temp[i][j] > max)
                max = Img_temp[i][j];
            if (Img_temp[i][j] < min)
                min = Img_temp[i][j];
        }
    }

    for (int i = 0; i < Height; i++)
    {
        for (int j = 0; j < Width; j++)
        {
            Img_Output[i][j] = (unsigned char)((Img_temp[i][j]-min)*255/(max-min));
        }
    }
}

void GradientThrPer(unsigned char** Img, unsigned long cum_hist[256], int Height, int Width, double per)
{
    for (int i = 0; i < Height; i++)
    {
        for (int j = 0; j < Width; j++)
        {
            if ((double)cum_hist[Img[i][j]]/(double)cum_hist[255] >= per) {
                Img[i][j] = 255;
            } else {
                Img[i][j] = 0;
            }
        }
    }
}

void LoG(unsigned char** Img_Input, unsigned char** Img_Output, double a, int Height, int Width, int extend)
{
    double** Img_temp = NULL;
    Img_temp = SpaceForNor2D(Img_temp, Height+extend-1, Width+extend-1);

    double max= -10000;
    double min = 10000;

    for (int i = extend/2; i < Height + extend/2; i++)
    {
        for (int j = extend/2; j < Width + extend/2; j++)
        {
            for (int x = -extend/2; x < extend/2 + 1; x++)
            {
                for (int y = -extend/2; y < extend/2 + 1; y++)
                {
                    double d = -1/(M_PI*pow(a,4))*(1-(double)(x*x+y*y)/(2*a*a))*exp(-1*(double)(x*x+y*y)/(2*a*a));
                    Img_temp[i][j] += (double)(Img_Input[i+x][j+y])*d;
                }
            }
            if (Img_temp[i][j] > max)
                max = Img_temp[i][j];
            if (Img_temp[i][j] < min)
                min = Img_temp[i][j];
        }
    }

    for (int i = extend/2; i < Height + extend/2; i++)
    {
        for (int j = extend/2; j < Width + extend/2; j++)
        {
            if (Img_temp[i][j] <= 0) {
                Img_Output[i][j] = (unsigned char)((Img_temp[i][j]-min)*127/(0-min) + 0.5);
            } else {
                Img_Output[i][j] = (unsigned char)((Img_temp[i][j])*127/(max-0) + 128.5);
            }
        }
    }
}

void LoGEdge(unsigned char** Img, unsigned long cum_hist[256], int Height, int Width, double per)
{
    for (int i = 0; i < Height; i++)
    {
        for (int j = 0; j < Width; j++)
        {
            if (Img[i][j] > 130) {
                if ((cum_hist[Img[i][j]] - cum_hist[130]) <= per*(cum_hist[255]-cum_hist[130])) {
                    Img[i][j] = 192;
                } else {
                    Img[i][j] = 0;
                }
            } else if (Img[i][j] < 124) {
                if ((cum_hist[124] - cum_hist[Img[i][j]]) <= per*cum_hist[124]) {
                    Img[i][j] = 64;
                } else {
                    Img[i][j] = 0;
                }
            } else {
                Img[i][j] = 128;
            }
        }
    }
}

void LoGEdgeThr(unsigned char** Img, unsigned char** Img_Output, int Height, int Width)
{
    int extend = 3;
    unsigned char** Img_temp = NULL;
    Img_temp = SpaceFor2D(Img_temp, Height+extend-1, Width+extend-1);

    ExtendBoundary(Img, Img_temp, extend, Height, Width);

    for (int i = extend/2; i < Height+extend/2; i++)
    {
        for (int j = extend/2; j < Width+extend/2; j++)
        {
            int flag_high = 0;
            int flag_low = 0;
            for (int p = -extend/2; p < extend/2+1; p++)
            {
                for (int q = -extend/2; q < extend/2+1; q++)
                {
                    if (Img_temp[i+p][j+q] == 192)
                    {
                        flag_high = 1;
                    }
                    if (Img_temp[i+p][j+q] == 64)
                    {
                        flag_low = 1;
                    }
                    if (flag_low*flag_high == 1)
                    {
                        Img_Output[i-extend/2][j-extend/2] = 255;
                        break;
                    }
                }
                if (flag_low*flag_high == 1)
                    break;
            }
        }
    }
}

void GaussianFilter(unsigned char** Img_Input, unsigned char** Img_Output, int Height, int Width, int Extend)
{
    unsigned char** Img_Input_E = NULL;
    Img_Input_E = SpaceFor2D(Img_Input_E, Height+Extend-1, Width+Extend-1);

    ExtendBoundary(Img_Input, Img_Input_E, Extend, Height, Width);

    double a = (((double)Extend)/2 - 1) * 0.3 + 0.8;

    for (int i = Extend/2; i < Height + Extend/2; i++)
    {
        for (int j = Extend/2; j < Width + Extend/2; j++)
        {
            double val = 0;
            for (int p = -Extend/2; p < Extend/2 + 1; p++)
            {
                for (int q = -Extend/2; q < Extend/2 + 1; q++)
                {
                    val += (double)Img_Input_E[i+p][j+q]*exp((0-p*p-q*q)/(2*a*a))/(2*M_PI*a*a);
                }
            }
            Img_Output[i-Extend/2][j-Extend/2] = (unsigned char)val;
        }
    }
}

void LaplacianFilter(unsigned char** Img_Input, unsigned char** Img_Output, int Height, int Width)
{
    vector< vector<double> >Laplace {{0, -1, 0},{-1, 4, -1},{0, -1, 0}};

    unsigned char** Img_Input_E = NULL;
    Img_Input_E = SpaceFor2D(Img_Input_E, Height+2, Width+2);

    double** Img_temp = NULL;
    Img_temp = SpaceForNor2D(Img_temp, Height, Width);

    ExtendBoundary(Img_Input, Img_Input_E, 3, Height, Width);

    double max = -10000;
    double min = 10000;

    for (int i = 1; i < 1+Height; i++)
    {
        for (int j = 1; j < 1+Width; j++)
        {

            for (int p = -1; p < 2; p++)
            {
                for (int q = -1; q < 2; q++)
                {
                    Img_temp[i-1][j-1] += (double)Img_Input_E[i+p][j+q]*Laplace.at(p+1).at(q+1);
                }
            }
            if (Img_temp[i-1][j-1] > max)
                max = Img_temp[i-1][j-1];
            if (Img_temp[i-1][j-1] < min)
                min = Img_temp[i-1][j-1];
        }
    }

    for (int i = 0; i < Height; i++)
    {
        for (int j = 0; j < Width; j++)
        {
            if (Img_temp[i][j] <= 0) {
                Img_Output[i][j] = (unsigned char)((Img_temp[i][j]-min)*127/(0-min) + 0.5);
            } else {
                Img_Output[i][j] = (unsigned char)((Img_temp[i][j])*127/(max-0) + 128.5);
            }
        }
    }
}

