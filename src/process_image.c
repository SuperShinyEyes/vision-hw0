#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

/*Images are stored in CHW format*/
int intmax(int a, int b) 
{
    if (a > b) {
        return a;
    } else {
        return b;
    }
}

int intmin(int a, int b) 
{
    if (a < b) {
        return a;
    } else {
        return b;
    }
}


float get_pixel(image im, int x, int y, int c)
{
    // TODO Fill this in
    int x_clamped = intmin(intmax(x, 0), im.w - 1);
    int y_clamped = intmin(intmax(y, 0), im.h - 1);

    return im.data[(c * im.w * im.h) + im.w * y_clamped + x_clamped];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    int x_clamped = intmin(intmax(x, 0), im.w - 1);
    int y_clamped = intmin(intmax(y, 0), im.h - 1);

    im.data[(c * im.w * im.h) + im.w * y_clamped + x_clamped] = v;
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    // TODO Fill this in
    memcpy(copy.data, im.data, im.w * im.h * im.c * sizeof(float));
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    // TODO Fill this in
    const int image_area = im.w * im.h;
    for(int i = 0; i < image_area; i++) {
        gray.data[i] = 0.299 * im.data[i] + 0.587 * im.data[i + image_area] + 0.114 * im.data[i + 2 * image_area];
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    // TODO Fill this in CHW
    int start_i = c * im.h * im.w;
    int end_i = start_i + im.h * im.w;
    for(int i = start_i; i < end_i; i++) {
        im.data[i] += v;
    }
}

void clamp_image(image im)
{
    int end_i = 3 * im.h * im.w;
    for(int i = 0; i < end_i; i++) {
        im.data[i] = fmin(1, im.data[i]);
    }
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    int size = im.h * im.w;
    for(int i = 0; i < size; i++){
        float * ptr_R = im.data + i;
        float * ptr_G = im.data + size + i;
        float * ptr_B = im.data + size * 2 + i;
        float V = three_way_max(*ptr_R, *ptr_G, *ptr_B );
        float C = V - three_way_min(*ptr_R, *ptr_G, *ptr_B );
        float S = 0;
        if (V > 0) {
            S = C / V;
        }
        float H_prime = 0;
        if(C == 0) {
            H_prime = 0;
        }else if(V == *ptr_R) {
            H_prime = (*ptr_G - *ptr_B) / C;
        }else if(V == *ptr_G) {
            H_prime = (*ptr_B - *ptr_R) / C + 2;
        }else if(V == *ptr_B) {
            H_prime = (*ptr_R - *ptr_G) / C + 4;
        }

        float H = 0;
        if (H_prime < 0) {
            H = H_prime / 6 + 1;
        } else {
            H = H_prime / 6;
        }
        *ptr_R = H, *ptr_G = S, *ptr_B = V;
    }
}

void hsv_to_rgb(image im)
{
    int size = im.h * im.w;
    for(int i=0; i<size; ++i) {
        float* ptr_H = im.data + i;
        float* ptr_S = ptr_H + size;
        float* ptr_V = ptr_S + size;
        float H = *ptr_H, S = *ptr_S, V = *ptr_V;

        float H_prime = 6 * H, C = S * V, m = V - C;
        float R = 0, G = 0, B = 0;
        if (H_prime < 1) {
            R = V;
            G = H_prime * (V - m) + m;
            B = m;
        } else if (H_prime < 2) {
            R = m + (2 - H_prime) * (V - m);
            G = V;
            B = m;
        } else if (H_prime < 3) {
            R = m;
            G = V;
            B = (H_prime - 2) * (V - m) + m;
        } else if (H_prime < 4) {
            R = m;
            G = (4 - H_prime) * (V - m) + m ;
            B = V;
        } else if (H_prime < 5) {
            R = (H_prime - 4) * (V - m) + m ;
            G = m;
            B = V;
        } else {
            H_prime = (H - 1) * 6;
            R = V;
            G = m;
            B = - H_prime * (V - m) + m;
        }
        *ptr_H = R, *ptr_S = G, *ptr_V = B;
    }

    // TODO Fill this in
}
