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

vector<vector<double>> matrix_3x3to4x4(const vector<vector<double>> &v) {
    return {{v[0][0], v[0][1], v[0][2], 0},
            {v[1][0], v[1][1], v[1][2], 0},
            {v[2][0], v[2][1], v[2][2], 0},
            {0, 0, 0, 1},
    };
}


vector<vector<double>> matrix_mult(const vector<vector<double>> &a, const vector<vector<double>> &b) { // 2d b
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


vector<double> inverse(vector<int> v) {
    vector<double> out;
    for (auto& i : v)
        out.push_back(i * -1);
    return out;
}

// inverse = transpose / determinant
vector<vector<double>> inverse(vector<vector<double>> v) {
    size_t y = v.size(), x = v[0].size();
    double det =    (x == 2 && y == 2) ? (v[0][0]*v[1][1] - v[0][1]*v[1][0]) : 
                    (x == 3 && y == 3) ? (v[0][0]*v[1][1]*v[2][2] + v[0][1]*v[1][2]*v[2][0] + v[0][2]*v[1][0]*v[2][1]
                                         - v[0][2]*v[1][1]*v[2][0] - v[0][1]*v[1][0]*v[2][2] - v[0][0]*v[1][2]*v[2][1])
                                       : 1; // must be 2x2 or 3x3
    for (size_t i = 0; i < y; ++i) {
        for (size_t j = i; j < x; ++j) {
            swap(v[i][j],v[j][i]);
        }
    }
    for (size_t i = 0; i < y; ++i) {
        for (size_t j = 0; j < x; ++j) {
            v[i][j] /= det;
        }
    }
    return v;
}

void print_m (vector<vector<double>> v) {
    size_t y = v.size(), x = v[0].size();
    for (size_t i = 0; i < y; ++i) {
        cout << "{ ";
        for (size_t j = 0; j < x; ++j) {
            cout << v[i][j] << " ";
        }
        cout << "}\n";
    }
}

vector<double> camera_pos = {0, 0, 0};
vector<double> camera_rot = {0, 0, 0};
vector<vector<double>> camera_rot_m(vector<double> v) {
    double x = v[0], y= v[1], z=v[2];
    double cx = cos(x), cy = cos(y), cz = cos(z),
           sx = sin(x), sy = sin(y), sz = sin(z);
           
    return matrix_3x3to4x4({{ cy*cz, cy*sz, -sy },
            { sx*sy*cz - cx*sz, sx*sy*sz + cx*cz, sx*cy },
            { cx*sy*cz + sx*sz, cx*sy*sz - sx*cz, cx*cy }});
}

vector<vector<int>> project(vector<vector<double>> points, vector<double> translation, vector<double> c_pos, vector<double> c_rot) { // returns vec of point pairs
    size_t size = points.size();
    vector<vector<double>> projected_points(size, vector<double>(4, 0));
    vector<vector<int>> out;
    for (size_t i = 0; i < size; ++i) {
        vector<double> ttf = three_to_four(points[i]);
        vector<double> translated = matrix_mult(translate(translation[0], translation[1], translation[2]), ttf);
        // convert to camera space (mult with rot^-1 * trans^-1)
        vector<double> camerad = matrix_mult(matrix_mult(inverse(camera_rot_m(camera_rot)), inverse(translate(c_pos[0], c_pos[1], c_pos[2]))), translated);
        projected_points[i] = matrix_mult(projection_m, camerad);
    
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
    /*
    vector<vector<double>> m = {{1,2,3,4},
                                {5,6,7,8},
                                {1,2,3,4},
                                {5,6,7,8}};
    cout << "og:\n";
    print_m(m);
    cout << "inverse:\n";
    print_m(inverse(m));
    */

    
    cout << "\x1b[2J";
    cout << "\x1b[?25l";

    int s = 0;
    for (int i = 0; i < 1000; ++i) {
        int t = i % 6;
        if (i % 6 == 0) s = (s+1) % 4;
        for (auto &row : image) fill(row.begin(), row.end(), 0);
        vector<vector<int>> cube_points;
        vector<double> cube_trans = {0.0, 0.0, 5.0};
        camera_pos = {0, 0, 0};
        camera_rot = {0, 0, 0};
        /*
        if (s == 0) cube_trans = {-3.0+t, 3, 6.0};
        else if (s == 1) cube_trans = {3.0, 3.0-t, 6.0};
        else if (s == 2) cube_trans = {3.0-t, -3.0, 6.0};
        else if (s == 3) cube_trans = {-3.0, -3.0+t, 6.0};
        */
        if (s == 0 || s == 2) {camera_rot = {0, -0.6+t*0.2, 0.0}; camera_pos = {0, -1.2+t*0.4, 5.0};}
        else if (s == 1 || s == 3) {camera_rot = {0, 0.6-t*0.2, 0.0}; camera_pos = {0, -1.2+t*0.4, 5.0};}
        cube_points = project(cube, cube_trans, camera_pos, camera_rot);
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
        Sleep(250);
    }
        
    return 0;
}