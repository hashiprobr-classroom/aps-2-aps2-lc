#include <math.h>

#include "fourier.h"

void normalize(complex s[], int n) {
    for (int k = 0; k < n; k++) {
        s[k].a /= n;
        s[k].b /= n;
    }
}

void nft(complex s[], complex t[], int n, int sign) {
    for (int k = 0; k < n; k++) {
        t[k].a = 0;
        t[k].b = 0;

        for (int j = 0; j < n; j++) {
            double x = sign * 2 * PI * k * j / n;

            double cosx = cos(x);
            double sinx = sin(x);

            t[k].a += s[j].a * cosx - s[j].b * sinx;
            t[k].b += s[j].a * sinx + s[j].b * cosx;
        }
    }
}

void nft_forward(complex s[], complex t[], int n) {
    nft(s, t, n, -1);
}

void nft_inverse(complex t[], complex s[], int n) {
    nft(t, s, n, 1);
    normalize(s, n);
}

void separa_pares(complex s[], complex s_par[], int n){
    for (int i = 0; i < n/2; i++) { // i é o índice de s_par, i*2 é o índice de s
        s_par[i] =  s[i * 2];  //ex: s_par[1] = s[2], só os pares
    }
    return;
}

void separa_impares(complex s[], complex s_impar[], int n){
    for (int i = 0; i < n/2; i++){
        s_impar[i] = s[i * 2 + 1];
    }
    return;
}

void fft(complex s[], complex t[], int n, int sign) {
    if (n == 1){
        t[0] = s[0];
        return;
    }
    complex s_par[n/2];
    complex s_impar[n/2];
    complex E[n/2];
    complex O[n/2];
    
    separa_pares(s, s_par, n);
    separa_impares(s, s_impar, n);
    
    fft(s_par, E, n/2, sign);
    fft(s_impar, O, n/2, sign);

    for (int k  = 0; k < n/2; k ++){
        double x = sign * 2 * PI * k/n;
        complex w;
        w.a = cos(x);
        w.b = sin(x);

      
        double wO_a = w.a * O[k].a - w.b * O[k].b;
        double wO_b = w.a * O[k].b + w.b * O[k].a;

        t[k].a = E[k].a + wO_a;
        t[k].b = E[k].b + wO_b;

        t[k + n/2].a = E[k].a - wO_a;
        t[k + n/2].b = E[k].b - wO_b;
    }   
}

void fft_forward(complex s[], complex t[], int n) {
    fft(s, t, n, -1);
}

void fft_inverse(complex t[], complex s[], int n) {
    fft(t, s, n, 1);
    normalize(s, n);
}

void fft_forward_2d(complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
    complex in[MAX_SIZE];
    complex out[MAX_SIZE];

    for (int c = 0; c < width; c++) {
        for (int l = 0; l < height; l++) {
            in[l] = matrix[l][c];
        }

        fft_forward(in, out, height);

        for (int l = 0; l < height; l++) {
            matrix[l][c] = out[l];
        }
    }

    for (int l = 0; l < height; l++) {
        for (int c = 0; c < width; c++) {
            in[c] = matrix[l][c];
        }

        fft_forward(in, out, width);

        for (int c = 0; c < width; c++) {
            matrix[l][c] = out[c];
        }
    }
}


void fft_inverse_2d(complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
     complex in[MAX_SIZE];
    complex out[MAX_SIZE];

    for (int l = 0; l < height; l++) {
        for (int c = 0; c < width; c++) {
            in[c] = matrix[l][c];
        }

        fft_inverse(in, out, width);

        for (int c = 0; c < width; c++) {
            matrix[l][c] = out[c];
        }
    }

    for (int c = 0; c < width; c++) {
        for (int l = 0; l < height; l++) {
            in[l] = matrix[l][c];
        }

        fft_inverse(in, out, height);

        for (int l = 0; l < height; l++) {
            matrix[l][c] = out[l];
        }
    }
}
    

void filter(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height, int flip) {
    int center_x = width / 2;
    int center_y = height / 2;

    double variance = -2 * SIGMA * SIGMA;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int dx = center_x - (x + center_x) % width;
            int dy = center_y - (y + center_y) % height;

            double d = dx * dx + dy * dy;
            double g = exp(d / variance);

            if (flip) {
                g = 1 - g;
            }

            output[y][x].a = g * input[y][x].a;
            output[y][x].b = g * input[y][x].b;
        }
    }
}

void filter_lp(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 0);
}

void filter_hp(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 1);
}
