/*

File: example3.cpp

An example of the use of the Quasisampler class.
Reminder: This is a toy implementation, created
to aid in understanding how the system works.

This example generates a set of points for a
grayscale bitmap.

Usage: example3 image.pgm [magnitude_factor=200]

This is a toy (non-optimized) implementation of the importance sampling
technique proposed in the paper:
"Fast Hierarchical Importance Sampling with Blue Noise Properties",
by Victor Ostromoukhov, Charles Donohue and Pierre-Marc Jodoin,
to be presented at SIGGRAPH 2004.


Implementation by Charles Donohue,
Based on Mathematica code by Victor Ostromoukhov.
Universite de Montreal
05.08.04

*/
#include <numeric>
#include <algorithm>
#include <string.h>
#include <iostream>
#include <fstream>
#include "quasisampler_prototype.h"
#include <opencv2/opencv.hpp>
using namespace std;




class ImageQuasiSampler_ : public Quasisampler {
public:
    unsigned *data;
    unsigned w, h;

    ImageQuasiSampler_(char *filename, double mag) {
        // Load the grayscale image
        if (!loadPGM(filename, mag)) {
            data = 0;
            width = height = 0;
        }
    }

    // Simple PGM parser (Low fault tolerance)
    bool loadPGM(char *filename, double mag = 1.0) {
        if (!filename) {
            cerr << "Could not load PGM: no filename given." << endl;
            return false;
        }
        char buffer[80];
        ifstream infile(filename);
        if (!infile.is_open()) {
            cerr << "Could not open file: " << filename << endl;
            return false;
        }
        infile.getline(buffer, 80);
        if (strcmp(buffer, "P2")) {
//            cerr << "PGM file header not recognized (P2 type only)" << endl; return false;
        }

        do { infile.getline(buffer, 80); } while (buffer[0] == '#');
        w = atoi(buffer);
        char *tmp = strchr(buffer, ' ');
        tmp++; // skip whitespace
        h = atoi(tmp);

        do { infile.getline(buffer, 80); } while (buffer[0] == '#');
        unsigned maxval, _buffer;
        maxval = atoi(buffer);  // nb: not used.
        data = new unsigned [w * h];
        for (unsigned i = 0; i < w * h; i++) {
            infile >> data[i];
        }
        infile.close();

        if (mag != 1.0) for (unsigned i = 0; i < w * h; i++) data[i] = (unsigned) (data[i] * mag);

        width = w;
        height = h;
        return true;
    }


    unsigned getImportanceAt( Point2D pt ){
        // Nearest pixel sampling.
        return data[ w * ((unsigned)(pt.y)) + (unsigned)(pt.x) ];
    }
};


namespace HilbertCurve{


//rotate/flip a quadrant appropriately
void rot(int n, int *x, int *y, int rx, int ry) {
    if (ry == 0) {
        if (rx == 1) {
            *x = n-1 - *x;
            *y = n-1 - *y;
        }
        std::swap(*x, *y);
    }
}

int xy2d (int n, int x, int y) {
    int rx, ry, s, d=0;
    for (s=n/2; s>0; s >>= 1) {
        rx = (x & s) > 0;
        ry = (y & s) > 0;
        d += s * s * ((3 * rx) ^ ry);
        rot(n, &x, &y, rx, ry);
    }
    return d;
}

//convert d to (x,y)
void d2xy(int n, int d, int *x, int *y) {
    int rx, ry, s, t=d;
    *x = *y = 0;
    for (s=1; s<n; s <<= 1) {
        rx = 1 & (t/2);
        ry = 1 & (t ^ rx);
        rot(s, x, y, rx, ry);
        *x += s * rx;
        *y += s * ry;
        t /= 4;
    }
}

struct HilbertSequence{
    unsigned x = -1, y = 0, i = 0;
    void hilbert(int dir, int rot, int order) {
        if (order == 0) return;
        dir = dir + rot;
        hilbert(dir, -rot, order - 1);
        step(dir);
        dir = dir - rot;
        hilbert(dir, rot, order - 1);
        step(dir);
        hilbert(dir, rot, order - 1);
        dir = dir - rot;
        step(dir);
        hilbert(dir, -rot, order - 1);
    }

    virtual void func() = 0;

    void step(int dir) {
        static int dx[4] = {1, 0, -1, 0};
        static int dy[4] = {0, 1, 0, -1};
        x += dx[dir&3];
        y += dy[dir&3];
        func();
    }
};


class ImageQuasisampler : public ImageQuasiSampler_
{
public:
    ImageQuasisampler(char* filename, double mag): ImageQuasiSampler_(filename, mag){}

    unsigned getImportanceAt( Point2D pt ){
        // Nearest pixel sampling.
        return data[ w * ((unsigned)(pt.y)) + (unsigned)(pt.x) ];
    }

    std::vector<Point2D> hilbert_dithering(size_t k){
        unsigned _width = std::pow( 2, std::ceil(std::log2(double(w))) );
        unsigned _height = std::pow( 2, std::ceil(std::log2(double(h))) );
        _width = _height = std::max(_width, _height);

        struct CurveDither: public HilbertSequence{
            unsigned *data, w, h, error{0}, k;
            std::vector<Point2D> result;
            CurveDither(unsigned *data_, unsigned w_, unsigned h_, unsigned k_): data(data_), w(w_), h(h_), k(k_){
               unsigned sum = std::accumulate(data, data + w * h, 0);
               result.reserve(sum / k );
//               std::cerr << "reserve to " << sum / k  << std::endl;
            }
            void func() override {
                if(x >= w || y >= h) return ;
                auto& weight = data[w*y + x];
                if(weight+error > k) {
                    ((CurveDither*)(this))->result.emplace_back(x, y);
//                    std::cout << x << ' ' << y << "\n";
                }
                error = (weight+error)%k;
            }
        } dither_sampling(data, w, h, k);

        dither_sampling.hilbert(0, 1, std::round( std::log2(_width) ) );

        return dither_sampling.result;
    }


};

}




namespace MortenCurve {



// "Insert" a 0 bit after each of the 16 low bits of x
uint32_t Part1By1(uint32_t x) {
    x &= 0x0000ffff;                  // x = ---- ---- ---- ---- fedc ba98 7654 3210
    x = (x ^ (x <<  8)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x <<  4)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x <<  2)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x <<  1)) & 0x55555555; // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    return x;
}

// "Insert" two 0 bits after each of the 10 low bits of x
uint32_t Part1By2(uint32_t x) {
    x &= 0x000003ff;                  // x = ---- ---- ---- ---- ---- --98 7654 3210
    x = (x ^ (x << 16)) & 0xff0000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x ^ (x <<  8)) & 0x0300f00f; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x ^ (x <<  4)) & 0x030c30c3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x ^ (x <<  2)) & 0x09249249; // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
    return x;
}

uint32_t EncodeMorton2(uint32_t x, uint32_t y) { return (Part1By1(y) << 1) + Part1By1(x); }
// Inverse of Part1By1 - "delete" all odd-indexed bits
uint32_t Compact1By1(uint32_t x){
    x &= 0x55555555;                  // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    x = (x ^ (x >>  1)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x >>  2)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x >>  4)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x >>  8)) & 0x0000ffff; // x = ---- ---- ---- ---- fedc ba98 7654 3210
    return x;
}

// Inverse of Part1By2 - "delete" all bits not at positions divisible by 3
uint32_t Compact1By2(uint32_t x){
    x &= 0x09249249;                  // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
    x = (x ^ (x >>  2)) & 0x030c30c3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x ^ (x >>  4)) & 0x0300f00f; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x ^ (x >>  8)) & 0xff0000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x ^ (x >> 16)) & 0x000003ff; // x = ---- ---- ---- ---- ---- --98 7654 3210
    return x;
}

uint32_t DecodeMorton2X(uint32_t code) { return Compact1By1(code >> 0); }

uint32_t DecodeMorton2Y(uint32_t code) { return Compact1By1(code >> 1); }



class ImageQuasisampler : public ImageQuasiSampler_ {
public:
    ImageQuasisampler(char *filename, double mag) : ImageQuasiSampler_(filename, mag) {
        unsigned _width = std::pow( 2, std::ceil(std::log2(double(w))) );
        unsigned _height = std::pow( 2, std::ceil(std::log2(double(h))) );
        _width = _height = std::max(_width, _height);
        unsigned  *new_data = new unsigned[_width*_height];
        std::fill(new_data, new_data + _width*_height, 0);
        for(size_t x = 0; x < w; x++)
            for(size_t y = 0; y < h; y++)
                new_data[ EncodeMorton2(x, y) ] = data[y*w +x ];
        std::swap(new_data, data);
    }


    unsigned getImportanceAt( Point2D pt ){
        // Nearest pixel sampling.
        return data[EncodeMorton2(pt.x, pt.y)];
    }

};


}


typedef std::vector<Point2D> PointList;

int main(int argc, char* argv[])
{
    if (argc<2)
    {
        std::cout << "Usage: " << argv[0] << " image.pgm [magnitude_factor = 200]" << std::endl;
    }

    double mag_factor = 500.0;
    if (argc>2) mag_factor = atof(argv[2]);

    cv::Mat img = cv::imread(argv[1], 0);


    {
        // initialize sampler
        ImageQuasiSampler_ test( argv[1], mag_factor );

        // generate points
        uint64_t begin_t = clock();
        PointList points = test.getSamplingPoints();

        std::cout << "[baseline ] \t time cost: " <<  double(clock() - begin_t)/CLOCKS_PER_SEC * 1000  << "ms  "
                  << points.size()<< " pts/["<< test.w << 'x'  << test.h << "]" << std::endl; // 834x1193


        std::ofstream obj_file("baseline.obj");
        for ( PointList::iterator it=points.begin(); it!=points.end(); it++ )
            obj_file << "v " << it->x << ' ' << test.h-it->y << " 0\n";
        obj_file.close();

//        cv::cvtColor(img, img, cv::COLOR_GRAY2BGR );
//        for ( PointList::iterator it=points.begin(); it!=points.end(); it++ ) {
//            cv::circle(img, cv::Point2f(it->x,  it->y), 1, cv::Scalar(255, 0, 0));
//        }
//
//        cv::imshow("img", img);
//        cv::waitKey(0);
    }

    {
        // initialize sampler
        MortenCurve::ImageQuasisampler test( argv[1], mag_factor );

        // generate points
        uint64_t begin_t = clock();
        PointList points = test.getSamplingPoints();

        std::cout << "[morton indexed] \t time cost: " <<  double(clock() - begin_t)/CLOCKS_PER_SEC * 1000  << "ms  "
                  << points.size()<< " pts/["<< test.w << 'x'  << test.h << "]" << std::endl; // 834x1193


        std::ofstream obj_file("morton.obj");
        for ( PointList::iterator it=points.begin(); it!=points.end(); it++ )
            obj_file << "v " << it->x << ' ' << test.h-it->y << " 0\n";
        obj_file.close();

//        cv::cvtColor(img, img, cv::COLOR_GRAY2BGR );
//        for ( PointList::iterator it=points.begin(); it!=points.end(); it++ ) {
//            cv::circle(img, cv::Point2f(it->x,  it->y), 1, cv::Scalar(255, 0, 0));
//        }
//
//        cv::imshow("img", img);
//        cv::waitKey(0);
    }

    {
        // initialize sampler
        HilbertCurve::ImageQuasisampler test( argv[1], mag_factor );

        // generate points
        uint64_t begin_t = clock();
//        PointList points = test.getSamplingPoints();
//
//        std::cout << "[hilbert indexed]\t time cost: " <<  double(clock() - begin_t)/CLOCKS_PER_SEC * 1000  << "ms  "
//                  << points.size()<< " pts/["<< test.w << 'x'  << test.h << "]" << std::endl; // 834x1193
//
//        std::ofstream obj_file("2.obj");
//        for ( PointList::iterator it=points.begin(); it!=points.end(); it++ )
//            obj_file << "v " << it->x << ' ' << it->y << " 0\n";
//        obj_file.close();


        for(size_t s :{10, 100, 500, 1000, 2000, 5000}) {
            begin_t = clock();
            PointList points2 = test.hilbert_dithering(mag_factor*s);
            std::cout << "[hilbert dither]\t time cost: " <<  double(clock() - begin_t)/CLOCKS_PER_SEC * 1000  << "ms  "
                      << points2.size()<< " pts/["<< test.w << 'x'  << test.h << "]" << std::endl; // 834x1193

                    std::ofstream obj2_file(std::to_string(s)+".obj");
        for ( PointList::iterator it=points2.begin(); it!=points2.end(); it++ )
            obj2_file << "v " << it->x << ' ' << test.h-it->y << " 0\n";
        obj2_file.close();

        }

    }




//        std::cout << std::distance(points.begin(), it) << ": " << it->x << "," << it->y << std::endl;

    return 0;
}
