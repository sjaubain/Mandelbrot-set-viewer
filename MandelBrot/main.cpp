#include <exception>
#include <string>
#include <iostream>
#include <cmath>
#include <string>
#include <SDL2/SDL.h>
#include <immintrin.h>
using namespace std;

#define WIDTH 400
#define HEIGHT 400
#define N_COLORS 1024
#define N_THREADS 4

#define NUM_TURN 6
#define REFRESH_DELAY_MILLIS 20
#define ZOOM_IN_FACTOR 0.85
#define ZOOM_OUT_FACTOR 1 / ZOOM_IN_FACTOR
#define KERNEL_SIZE 3

int kernel[KERNEL_SIZE][KERNEL_SIZE] = {
    {1, 1, 1},
    {1, 4, 1},
    {1, 1, 1}
};

__m256i mask1 = _mm256_set1_epi32(0xFFFFFFFF);

__m256d _4 = _mm256_set1_pd(4.0);

/*
__m512i mask11 = _mm512_set1_epi32(0xFFFFFFFF);

__m512d _22 = _mm512_set1_pd(2);
__m512d _44 = _mm512_set1_pd(4);
*/

typedef struct args {
    int threadId;
    int startX;
    int endX;
} threadArgs;

int rgbToInt(int r, int g, int b);

void initColorMap();

void applyGaussianBlur();

void render();

void repaint();

int threadDraw(void* data);

int threadDrawSIMD(void* data);

//int threadDrawSIMD2(void* data);

int threadRender(void* data);

void displayThreadMessage(int threadId, string msg);

// shared between threads
double centerX = -0.63405500073356879653;
double centerY = -0.38243699932279362486;
double scaleFactor = 0.006;//0.0009765625;
int maxIter;
int turn;

/**
* the color is a 32 bits value composed
* of [0 r g b]. Each color component is an
* 8 bits sub-value. To get r, g or b, just
* shift right and mask with 0x000000FF
**/
int colorMap[N_COLORS];

int image[WIDTH][HEIGHT];

SDL_Window* screen;

SDL_Renderer* renderer;

SDL_Event event;

SDL_mutex* mutex;

SDL_cond* mustRestart;

SDL_sem* mustRender;

SDL_sem* nextTurn;

threadArgs argsTab[N_THREADS];

int restart = 0, quit = 0;

int main(int argc, char* argv[]) {

    SDL_Init(SDL_INIT_VIDEO);

    SDL_CreateWindowAndRenderer(WIDTH, HEIGHT, 0, &screen, &renderer);

    mutex = SDL_CreateMutex();

    mustRestart = SDL_CreateCond();

    mustRender = SDL_CreateSemaphore(0);

    nextTurn = SDL_CreateSemaphore(0);

    initColorMap();


    /*
    // a enlever plus tard
    SDL_Renderer* colorMapRenderer;
    SDL_Window* colorMapScreen;
    SDL_CreateWindowAndRenderer(1024, 100, 0, &colorMapScreen, &colorMapRenderer);
    for(int i = 0; i < 1024; ++i) {
        for(int j = 0; j < 40; ++j) {
            SDL_SetRenderDrawColor(colorMapRenderer, colorMap[i] >> 16, (colorMap[i] >> 8) & 0xFF, colorMap[i] & 0xFF, 255);
            SDL_RenderDrawPoint(colorMapRenderer, i, j);

            SDL_SetRenderDrawColor(colorMapRenderer, colorMap[i] >> 16, (colorMap[i] >> 8) & 0xFF, colorMap[i] & 0xFF, 255);
            SDL_RenderDrawPoint(colorMapRenderer, i, 60 + j);
        }
    }
    SDL_RenderPresent(colorMapRenderer);
    */

    // Init rendering thread
    SDL_Thread* renderThread = SDL_CreateThread(threadRender, "renderThread", NULL);

    // Init drawing threads
    for(int i = 0; i < N_THREADS; ++i) {

        argsTab[i] = {
            i,
            -1 * WIDTH / 2 + (WIDTH / N_THREADS) * i,
            -1 * WIDTH / 2 + (WIDTH / N_THREADS) * (i + 1)
        };

        SDL_Thread* thread = SDL_CreateThread(threadDrawSIMD, "drawingThread", (void*)(&argsTab[i]));
    }

    int time, oldX, oldY, newX, newY; // coordinate of mouse when clicked
    int buttonDown = 0;
    while(!quit) {

        time = SDL_GetTicks();

        while(SDL_PollEvent(&event)) {

            switch(event.type) {

            case SDL_MOUSEBUTTONDOWN:

                SDL_GetMouseState(&oldX, &oldY);
                buttonDown = 1;
                break;

            case SDL_MOUSEMOTION:

                if(buttonDown) {
                    SDL_GetMouseState(&newX, &newY);
                    centerX -= (newX - oldX) * scaleFactor;
                    centerY -= (newY - oldY) * scaleFactor;
                    oldX = newX;
                    oldY = newY;
                    repaint();
                }
                break;

            case SDL_MOUSEBUTTONUP:

                buttonDown = 0;
                break;

            case SDL_MOUSEWHEEL:

                if(event.wheel.y > 0) {
                    scaleFactor *= ZOOM_IN_FACTOR;
                } else {
                    scaleFactor *= ZOOM_OUT_FACTOR;
                }

                repaint();
                break;

            case SDL_QUIT:
                quit = 1;
                break;
            }
        }

        if(SDL_GetTicks() - time < REFRESH_DELAY_MILLIS) {
            SDL_Delay(REFRESH_DELAY_MILLIS - (SDL_GetTicks() - time));
        }
    }

    SDL_DestroyWindow(screen);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyMutex(mutex);
    SDL_DestroyCond(mustRestart);
    SDL_DestroySemaphore(mustRender);
    SDL_DestroySemaphore(nextTurn);
    SDL_Quit();

    return EXIT_SUCCESS;
}

void applyGaussianBlur() {

    int sumR, sumG, sumB;
    for(int i = KERNEL_SIZE / 2; i < WIDTH - KERNEL_SIZE / 2; ++i) {
        for(int j = 2; j < HEIGHT - 2; ++j) {
            sumR = sumG = sumB = 0;
            for(int k = i - KERNEL_SIZE / 2, m = 0; k <= i + KERNEL_SIZE / 2; ++k, ++m) {
                for(int l = j - KERNEL_SIZE / 2, n = 0; l <= j + KERNEL_SIZE / 2; ++l, ++n) {
                    sumR += (image[k][l] >> 16) * kernel[m][n] / 12.0;
                    sumG += ((image[k][l] >> 8) & 0xFF) * kernel[m][n] / 12.0;
                    sumB += (image[k][l] & 0xFF) * kernel[m][n] / 12.0;
                }
            }
            image[i][j] = rgbToInt((int)sumR, (int)sumG, (int)sumB);
        }
    }
}

void render() {

    //applyGaussianBlur();

    for(int i = 0; i < WIDTH; ++i) {
        for(int j = 0; j < HEIGHT; ++j) {
            SDL_SetRenderDrawColor(renderer, image[i][j] >> 16, (image[i][j] >> 8) & 0xFF, image[i][j] & 0xFF, 255);
            SDL_RenderDrawPoint(renderer, i, j);
        }
    }
    SDL_RenderPresent(renderer);
}

void repaint() {
    restart = 1;
    SDL_CondBroadcast(mustRestart);
}

void initColorMap() {

    for(int i = 0; i < 256; ++i) {
        if(i < 128) {
            colorMap[i] = rgbToInt(1.32 * i, 2 * i, 2 * i);
            colorMap[256 + i] = rgbToInt(42 + i / 2, 42 + i / 2, 255);
            colorMap[512 + i] = rgbToInt(170 + i / 3, 170 + i / 3, 255 - i);
            colorMap[768 + i] = rgbToInt(255, 255 - 2 * i, 0);
        } else {
            colorMap[i] = rgbToInt(170 - (i - 128), 255 - 1.66 * (i - 128), 255);
            colorMap[256 + i] = rgbToInt(42 + i / 2, 42 + i / 2, 255);
            colorMap[512 + i] = rgbToInt(170 + i / 3, 170 + i / 3, 255 - i);
            colorMap[768 + i] = rgbToInt(255 - 2 * (i - 128), 0, 0);
        }

        //colorMap[512 + i] = rgbToInt(170 + i / 3, 170 + i / 3, 255 - i);
    }
}

int rgbToInt(int r, int g, int b) {
    return ((r << 16) & 0x00FF0000 | (g << 8) & 0x0000FF00 | b & 0xFF);
}

int threadDraw(void* data) {

    int iter = 0;

    threadArgs* args = (threadArgs*)data;
    int id = (int)args -> threadId;

    while(1) {

        SDL_SemWait(nextTurn);

        // complex number a + bi
        double a, b, aTmp;
        for(int i = (int)args -> startX; i < (int)args -> endX; ++i) {

            if(restart)
                break;
            if(quit)
                return 0;

            double cX = centerX + i * scaleFactor;
            for(int j = -1 * HEIGHT / 2; j < HEIGHT / 2; ++j) {

                iter = a = b = 0;

                // zn = zn ^ 2 + c
                double cY = centerY + j * scaleFactor;
                do {
                    aTmp = a * a - b * b + cX;
                    b = a * b + a * b + cY;
                    a = aTmp;
                } while(a * a + b * b < 4 && ++iter < maxIter);

                // set black pixel if iter == maxIter (colorMap[0])
                image[i + WIDTH / 2][j + HEIGHT / 2] = colorMap[(iter != maxIter) * (iter % N_COLORS)];
            }
        }

        // emit rendering
        SDL_SemPost(mustRender);
    }

    return 0;
}

/*
int threadDrawSIMD2(void* data) {

    int iter = 0;

    threadArgs* args = (threadArgs*)data;
    int id = (int)args -> threadId;

    while(1) {

        SDL_SemWait(nextTurn);

        // complex number a + bi
        __m512i _iter;
        __m512d _atmp;
        __m512d _a;
        __m512d _b;
        __m512i _maxIter = _mm512_set1_epi32(maxIter);
        for(int i = (int)args -> startX; i < (int)args -> endX; ++i) {

            if(restart)
                break;
            if(quit)
                return 0;

            __m512d _cX = _mm512_set1_pd(centerX + i * scaleFactor);

            for(int j = -1 * HEIGHT / 2; j < HEIGHT / 2; j += 8) {

                _iter = _mm512_setzero_si512();
                _a = _mm512_setzero_pd();
                _b = _mm512_setzero_pd();

                __m512d _cY = _mm512_set4_pd(centerY + j * scaleFactor,
                                            centerY + (j + 1) * scaleFactor,
                                            centerY + (j + 2) * scaleFactor,
                                            centerY + (j + 3) * scaleFactor);

                // zn = zn ^ 2 + c
                __mmask8 cond1;
                __mmask8 cond2;
                do {

                    // aTmp = a * a - b * b + cX
                    _atmp = _mm512_sub_pd(_mm512_mul_pd(_a, _a), _mm512_mul_pd(_b, _b));
                    _atmp = _mm512_add_pd(_atmp, _cX);

                    // b = 2 * a * b + cY
                    _b = _mm512_mul_pd(_22, _mm512_mul_pd(_a, _b));
                    _b = _mm512_add_pd(_b, _cY);

                    // a = aTmp
                    _a = _atmp;

                    // cond1 : a^2 + b^2 < 4
                    cond1 = _mm512_cmp_pd_mask(_mm512_add_pd(_mm512_mul_pd(_a, _a), _mm512_mul_pd(_b, _b)), _44, _CMP_LT_OS);

                    // increments iter if cond1 is true for each packed double
                    _iter = _mm512_add_epi64(_iter, _mm512_setr_epi64(cond1 & 128,cond1 & 64,cond1 & 32,cond1 & 16,cond1 & 8,cond1 & 4,cond1 & 2,cond1 & 1));

                    // cond2 : iter < maxIter
                    cond2 = _mm512_cmp_epi64_mask(_iter, _maxIter, _CMP_LT_OS);

                } while(cond1 & cond2 != 0x0);

                long long iters[8];
                _mm512_store_epi32(iters, _iter);

                for(int k = 0; k < 8; ++k) {
                    // set black pixel if iter == maxIter (colorMap[0])
                    image[i + WIDTH / 2][j + k + HEIGHT / 2] = colorMap[(iters[8 - k - 1] != maxIter) * (iters[8 - k - 1] % N_COLORS)];
                }
            }
        }

        // emit rendering
        SDL_SemPost(mustRender);
    }

    return 0;
}
*/

int threadDrawSIMD(void* data) {

    int iter = 0;

    threadArgs* args = (threadArgs*)data;
    int id = (int)args -> threadId;

    while(1) {

        SDL_SemWait(nextTurn);

        // complex number a + bi
        __m256i _iter;
        __m256d _atmp;
        __m256d _a, _b, _ab, _aa, _bb;
        __m256i _maxIter = _mm256_set1_epi64x(maxIter);
        for(int i = (int)args -> startX; i < (int)args -> endX; ++i) {

            if(restart)
                break;
            if(quit)
                return 0;

            __m256d _cX = _mm256_set1_pd(centerX + i * scaleFactor);
            __m256d _cY;
            for(int j = -1 * HEIGHT / 2; j < HEIGHT / 2; j += 4) {

                _iter = _mm256_setzero_si256();
                _a = _mm256_setzero_pd();
                _b = _mm256_setzero_pd();

                _cY = _mm256_set_pd(centerY + j * scaleFactor,
                                            centerY + (j + 1) * scaleFactor,
                                            centerY + (j + 2) * scaleFactor,
                                            centerY + (j + 3) * scaleFactor);

                // zn = zn ^ 2 + c
                __m256d cond1;
                __m256i cond2;
                do {

                    // aTmp = a * a - b * b + cX
                    _atmp = _mm256_sub_pd(_aa = _mm256_mul_pd(_a, _a), _bb = _mm256_mul_pd(_b, _b));
                    _atmp = _mm256_add_pd(_atmp, _cX);

                    // b = 2 * a * b + cY
                    _b = _mm256_add_pd(_ab = _mm256_mul_pd(_a, _b), _ab);
                    _b = _mm256_add_pd(_b, _cY);

                    // a = aTmp
                    _a = _atmp;

                    // cond1 : a^2 + b^2 < 4
                    cond1 = _mm256_cmp_pd(_mm256_add_pd(_aa, _bb), _4, _CMP_LT_OS);

                    // increments iter if cond1 is true for each packed double
                    _iter = _mm256_add_epi64(_iter, _mm256_srli_epi64((__m256i)cond1, 63));

                    // cond2 : iter < maxIter
                    cond2 = _mm256_cmpgt_epi64(_maxIter, _iter);

                } while(!_mm256_testz_si256(_mm256_and_si256(__m256i(cond1), cond2), mask1));

                long long iters[4];
                _mm256_maskstore_epi64(iters, mask1, _iter);

                for(int k = 0; k < 4; ++k) {
                    // set black pixel if iter == maxIter (colorMap[0])
                    image[i + WIDTH / 2][j + k + HEIGHT / 2] = colorMap[(iters[4 - k - 1] != maxIter) * (iters[4 - k - 1] % N_COLORS)];
                }
            }
        }

        // emit rendering
        SDL_SemPost(mustRender);
    }

    return 0;
}

int threadRender(void* data) {

    while(1) {

        turn = 0;
        maxIter = (1 << (2 * turn + 6)) + 32;
        while(turn++ < NUM_TURN) {

            if(restart)
                break;
            if(quit)
                return 0;

            maxIter *= (1 << 2);
            int curTime = SDL_GetTicks();

            // wake all threads
            for(int t = 0; t < N_THREADS; ++t)
                SDL_SemPost(nextTurn);

            // wait for all threads to render image
            for(int t = 0; t < N_THREADS; ++t)
                SDL_SemWait(mustRender);

            //printf("center : ( %20.20f ; %20.20f )  scale factor : %20.20f\n", centerX, centerY, scaleFactor);
            cout << "time (ms) for turn " << turn << " : " << SDL_GetTicks() - curTime << endl;
            if(!restart)
                render();
        }

        SDL_LockMutex(mutex);
        if(!restart)
            SDL_CondWait(mustRestart, mutex);
        SDL_UnlockMutex(mutex);

        restart = 0;
    }

    return 0;
}

void displayThreadMessage(int threadId, string msg) {
    SDL_LockMutex(mutex);
    cout << "[" << threadId << "] : " << msg << endl;
    SDL_UnlockMutex(mutex);
}
