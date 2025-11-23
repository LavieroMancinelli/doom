#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include <stdexcept>
#include <cmath>
#include <windows.h>
using namespace std; 

const int IMAGE_HEIGHT = 40;
const int IMAGE_WIDTH = 80;
vector<vector<int>> image(IMAGE_HEIGHT, vector<int>(IMAGE_WIDTH, 0));

vector<vector<double>> cube = { {-1,-1,-1}, {1,-1,-1}, {1,-1,1}, {-1, -1, 1},
                                {-1,1,-1}, {1,1,-1}, {1,1,1}, {-1, 1, 1}};

vector<double> three_to_four(vector<double> &v) {
    if (v.size() != 3) throw invalid_argument("Vector must be size 3");
    return vector<double>{v[0], v[1], v[2], 1};
}


vector<vector<double>> matrix_mult(vector<vector<double>> &a, vector<vector<double>> &b) { // 2d b
    size_t ay = a.size(), ax = a[0].size(), by = b.size(), bx = b[0].size();
    vector<vector<double>> res(ay, vector<double>(bx, 0));
    // ax must == by, result is ay by bx
    for (size_t i = 0; i < ay; ++i) {
        for (size_t j = 0; j < bx; ++j) {
            for (size_t k = 0; k < ax; ++k) // columns in first ax
                res[i][j] += a[i][k] * b[k][j];
        }
    }
    return res;
}

vector<double> matrix_mult(const vector<vector<double>> &a, const vector<double> b) { // 1d b
    size_t ay = a.size(), ax = a[0].size(), by = b.size();
    vector<double> res(ay, 0);
    // ax must == by, result is ay by bx
    for (size_t i = 0; i < ay; ++i) {
        for (size_t k = 0; k < ax; ++k) // columns in first ax
            res[i] += a[i][k] * b[k];
    }
    return res;
}

vector<vector<double>> translate(double x, double y, double z) {
    return {{1, 0, 0, x},
            {0, 1, 0, y},
            {0, 0, 1, z},
            {0, 0, 0, 1}};
}

vector<vector<double>> projection_m = {{1, 0, 0, 0},
                                    {0, 1, 0, 0},
                                    {0, 0, 1, 1},
                                    {0, 0, 1, 0}};


vector<vector<int>> project(vector<vector<double>> points, vector<int> translation) { // returns vec of point pairs
    size_t size = points.size();
    vector<vector<double>> projected_points(size, vector<double>(4, 0));
    vector<vector<int>> out;
    for (size_t i = 0; i < size; ++i) {
        vector<double> ttf = three_to_four(points[i]);
        vector<double> translated = matrix_mult(translate(translation[0], translation[1], translation[2]), ttf);
        projected_points[i] = matrix_mult(projection_m, translated);
    
        double x = projected_points[i][0];
        double y = projected_points[i][1];
        double z = projected_points[i][2];
        double w = projected_points[i][3];
        double xr = double(x)/w, yr = double(y)/w; // x real
        int xs = floor((xr + 1.0) * IMAGE_WIDTH / 2); // x screen
        int ys = floor((-yr + 1.0) * IMAGE_HEIGHT / 2); // x screen
        out.push_back({xs, ys});
    }
        
    return out;

}

bool offscreen(int x, int y) {
    return x < 0 || x >= IMAGE_WIDTH || y < 0 || y >= IMAGE_HEIGHT;
}

void draw_line(vector<int> a, vector<int> b) { // bresenham's line algo
    int x0 = a[0], x1 = b[0], y0 = a[1], y1 = b[1];
    
    bool steep = abs(y1 - y0) > abs(x1 - x0);

    if (steep) {
        swap(x0, y0);
        swap(x1, y1);
    }
    if (x0 > x1) {
        swap(x0, x1);
        swap(y0, y1);
    }

    int dx = x1 - x0, dy = abs(y1 - y0);
    int error = dx / 2;
    int ystep = (y0 < y1) ? 1 : -1;
    double D = 2*dy - dx;
    double y = y0;

    for (double x = x0; x <= x1; ++x) {
        if (!steep) {
            if (!offscreen(x, y)) image[y][x] = 1;
        } else {
            if (!offscreen(y, x)) image[x][y] = 1;
        }
        
        error -= dy;
        if (error < 0) {
            y += ystep;
            error += dx;
        }
    }
}

void render(vector<vector<int>> image) {
    size_t y_size = image.size(), x_size = image[0].size();
    string output = "";
    for (size_t i = 0; i < y_size; ++i) {
        for (size_t j = 0; j < x_size; ++j) {
            if (image[i][j] == 1) {
                output += '#';
                
                /*
                if (i == 0 || i == y_size-1) output += '-';
                else if (j == 0 || j == x_size-1) output += '|';
                else if (image[i-1][j] == 1 || image[i+1][j] == 1) output += '|';
                else if (image[i-1][j+1] == 1 && image[i+1][j-1] == 1) output += '/';
                else if (image[i-1][j-1] == 1 && image[i+1][j+1] == 1) output += '\\\\';
                else  output += '-';
                */
                

            }
            else output += ' ';
        }
        output += '\n';
    }
    cout << "\x1b[H" << output;
}

int main() {
    cout << "\x1b[2J";
    cout << "\x1b[?25l";

    int s = 0;
    for (int i = 0; i < 1000; ++i) {
        int t = i % 6;
        if (i % 6 == 0) s = (s+1) % 4;
        for (auto &row : image) fill(row.begin(), row.end(), 0);
        vector<vector<int>> cube_points;
        if (s == 0) cube_points = project(cube, {-3+t, 3, 6});
        else if (s == 1) cube_points = project(cube, {3, 3-t, 6});
        else if (s == 2) cube_points = project(cube, {3-t, -3, 6});
        else if (s == 3) cube_points = project(cube, {-3, -3+t, 6});
        draw_line(cube_points[0], cube_points[1]);
        draw_line(cube_points[1], cube_points[2]);
        draw_line(cube_points[2], cube_points[3]);
        draw_line(cube_points[3], cube_points[0]);
        draw_line(cube_points[4], cube_points[5]);
        draw_line(cube_points[5], cube_points[6]);
        draw_line(cube_points[6], cube_points[7]);
        draw_line(cube_points[7], cube_points[4]);
        draw_line(cube_points[0], cube_points[4]);
        draw_line(cube_points[1], cube_points[5]);
        draw_line(cube_points[2], cube_points[6]);
        draw_line(cube_points[3], cube_points[7]);
        render(image);
        Sleep(100);
    }
    return 0;
}