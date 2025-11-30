#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <windows.h>
#include <unordered_map>
using namespace std; 

const int IMAGE_HEIGHT = 40;
const int IMAGE_WIDTH = 80;
int sprites_created = 0;
vector<vector<int>> image(IMAGE_HEIGHT, vector<int>(IMAGE_WIDTH, 0));
const int AMMO_MAX = 16;
int ammo = AMMO_MAX;
bool reloading = false;
const int RELOAD_DUR = 32;
int reload_counter = 0;

vector<vector<double>> cube = { {-1,-1,-1}, {1,-1,-1}, {1,-1,1}, {-1, -1, 1},
                                {-1,1,-1}, {1,1,-1}, {1,1,1}, {-1, 1, 1}};

                                
vector<vector<int>> enemy ={{0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,3,3,3,0,0,0,0},
                            {0,0,3,3,2,2,2,3,3,0,0},
                            {0,3,2,2,2,2,2,2,2,3,0},
                            {3,2,2,2,3,3,3,2,2,2,3},
                            {3,2,2,2,3,3,3,2,2,2,3},
                            {0,3,2,2,2,2,2,2,2,3,0},
                            {0,0,3,3,2,2,2,3,3,0,0},
                            {0,0,0,0,3,3,3,0,0,0,0},
                            {0,0,0,0,0,0,0,0,0,0,0}};
/*{{0,0,0,0,0,0,0,0,0,0,0,0},
                            {0,0,0,0,1,1,1,1,0,0,0,0},
                            {0,0,0,0,1,1,1,1,0,0,0,0},
                            {0,0,0,1,1,1,1,1,1,1,0,0},
                            {0,0,1,1,1,1,1,1,1,1,1,0},
                            {0,1,1,1,1,1,1,1,1,1,1,0},
                            {0,1,1,0,1,1,1,1,0,0,0,0},
                            {0,1,1,0,1,1,1,1,0,0,0,0},
                            {0,0,0,1,1,1,1,1,1,0,0,0},
                            {0,0,0,1,1,1,1,1,1,0,0,0},
                            {0,0,0,1,1,0,0,1,1,0,0,0},
                            {0,0,0,1,1,0,0,1,1,0,0,0} };
                            */
/*{{0,0,1,1,1,0,0},
                            {0,0,1,1,1,0,0},
                            {0,1,1,1,1,1,0},
                            {0,1,1,1,1,1,0},
                            {0,1,1,1,1,1,0},
                            {0,0,1,0,1,0,0},
                            {0,0,1,0,1,0,0}};
                            */

/*{{0,0,0,1,1,1,0,0,0},
                            {0,0,0,1,1,1,0,0,0},
                            {0,0,0,0,1,0,0,0,0},
                            {0,0,1,1,1,1,1,0,0},
                            {0,1,1,1,1,1,1,1,0},                            
                            {1,1,0,1,1,1,0,1,1},
                            {1,0,0,1,1,1,0,0,1},
                            {0,0,1,1,0,1,1,0,0},                            
                            {0,0,1,1,0,1,1,0,0},
                            {0,0,1,1,0,1,1,0,0}};
                            */
/*
/----------------------------------------\
|                                        |
|   ####    ###    ###   #   #           |
|   #   #  #   #  #   #  ## ##           |
|   #   #  #   #  #   #  # # #           |
|   ####    ###    ###   #   # at home   |
|                                        |
|         by Laviero Mancinelli          |
\----------------------------------------/

Press SPACE to begin
*/

vector<pair<vector<int>, vector<vector<char>>>> ui_elements = {
    {{IMAGE_HEIGHT-3, IMAGE_WIDTH-17}, 
       {{'A','M','M','O',':','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0'},
        {'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0'}}}
};
/*vector<vector<char>> = {{'A','M','M','O',':','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0'},
                        {'|','|','|','|','|','|','|','|','|','|','|','|','|','|','|','|'}};
                        */
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

vector<vector<double>> translate_m(double x, double y, double z) {
    return {{1, 0, 0, x},
            {0, 1, 0, y},
            {0, 0, 1, z},
            {0, 0, 0, 1}};
}

vector<vector<double>> translate_m_inv(double x, double y, double z) {
    return {{1, 0, 0, -x},
            {0, 1, 0, -y},
            {0, 0, 1, -z},
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

vector<int> project(const vector<double>& points, const vector<double>& c_pos, const vector<double>& c_rot) { // returns vec of point pairs
    // points is already converted to vec4
    vector<double> projected_points(4, 0);
    vector<int> out;
    // convert to camera space (mult with rot^-1 * trans^-1)
    vector<double> camerad = matrix_mult(matrix_mult(inverse(camera_rot_m(camera_rot)), translate_m_inv(c_pos[0], c_pos[1], c_pos[2])), points);
    projected_points = matrix_mult(projection_m, camerad);
    
    double x = projected_points[0];
    double y = projected_points[1];
    double z = projected_points[2];
    double w = projected_points[3];
    if (w <= 0) return {};
    double xr = double(x)/w, yr = double(y)/w; // x real
    int xs = floor((xr + 1.0) * IMAGE_WIDTH / 2); // x screen
    int ys = floor((-yr + 1.0) * IMAGE_HEIGHT / 2); // x screen
    out = {xs, ys, (int)camerad[2]};
        
    return out;

}

bool offscreen(int x, int y, vector<vector<int>>& canvas) {
    int h = canvas.size(), w = canvas[0].size();
    return x < 0 || x >= w || y < 0 || y >= h;
}

void draw_line(vector<int> a, vector<int> b, vector<vector<int>>& canvas) { // bresenham's line algo
    // -1 is means will make 0 when added
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
    int y = y0;

    int h = canvas.size();
    int w = canvas[0].size();
    for (int x = x0; x <= x1; ++x) {

        int px = steep ? y : x;
        int py = steep ? x : y;
        if (!offscreen(px, py, canvas)) canvas[py][px] = 1;
        
        error -= dy;
        if (error < 0) {
            y += ystep;
            error += dx;
        }
    }

    // clip left and right edges
    int left_top = 0;
    while (canvas[left_top][0] != 1 && left_top < h-1) left_top++;
    int left_bottom = h-1;
    while (canvas[left_bottom][0] != 1 && left_bottom > 0) left_bottom--;
    int right_top = 0;
    while (canvas[right_top][w-1] != 1 && right_top < h-1) right_top++;
    int right_bottom = h-1;
    while (canvas[right_bottom][w-1] != 1 && right_bottom > 0) right_bottom--;
    for (int i = left_top; i <= left_bottom; ++i)
        if (i < h && i > 0) canvas[i][0] = 1;
    for (int i = right_top; i <= right_bottom; ++i)
        if (i < h && i > 0) canvas[i][w-1] = 1;
}

void add_canvas(vector<vector<int>>& out, vector<vector<int>>& part) {// adds part to out
    size_t h = out.size(), w = out[0].size();
    for (size_t i = 0; i < h; ++i) {
        for (size_t j = 0; j < w; ++j) {
            if (part[i][j] == 1)
                out[i][j] = 1;
            else if (part[i][j] == -1)
                out[i][j] = 0;
            else if (part[i][j] > 1) // # on enemy
                out[i][j] = part[i][j];
            //out[i][j] = clamp(out[i][j] + part[i][j], 0, 1); // if part has -1, out will be 0
        }
    }
} 

void fill_poly(vector<vector<int>>& out) { // fills space inbetween 1s with -1s
    size_t h = out.size(), w = out[0].size();
    for (size_t i = 0; i < h; ++i) {
        size_t l = 0, r = 0;
        for (int j = 0; j < w; ++j) {
            if (out[i][j] == 1) {
                while (j < w-1 && out[i][j] == 1) ++j;
                l = j;
                break;
            }
        }
        for (int j = w-1; j > 0; --j) {
            if (out[i][j] == 1) {
                while (j > 0 && out[i][j] == 1) --j;
                r = j;
                break;
            }
        }
        if (l < w && r > 0 && l <= r) {
            for (int j = l; j <= r; ++j) {
                out[i][j] = -1;
            }
        }
    }
}

void render(vector<vector<int>>& image, vector<pair<vector<int>, vector<vector<char>>>> &ui_els) {
    size_t y_size = image.size(), x_size = image[0].size();
    string output = "";
    //image[IMAGE_HEIGHT/2][IMAGE_WIDTH/2-1] = '-';
    //image[IMAGE_HEIGHT/2][IMAGE_WIDTH/2-2] = '-';
    image[IMAGE_HEIGHT/2][IMAGE_WIDTH/2] = -2;
    //image[IMAGE_HEIGHT/2][IMAGE_WIDTH/2+1] = '-';
    //image[IMAGE_HEIGHT/2][IMAGE_WIDTH/2+2] = '-';
    //image[IMAGE_HEIGHT/2-1][IMAGE_WIDTH/2] = '|';
    //image[IMAGE_HEIGHT/2+1][IMAGE_WIDTH/2] = '|';
    for (size_t i = 0; i < y_size; ++i) {
        for (size_t j = 0; j < x_size; ++j) {
            if (image[i][j] < 0) output += '+';
            else if (image[i][j] == 1 || image[i][j] % 2 == 1) {
                //output += image[i][j] + '0';
                output += '#';
                
                /*
                if (i == 0 || i == y_size-1) output += '-';
                else if (j == 0 || j == x_size-1) output += '|';
                else if (image[i-1][j] == 1 && image[i+1][j] == 1) output += '|';
                else if (image[i-1][j+1] == 1 && image[i+1][j-1] == 1) output += '/';
                else if (image[i-1][j-1] == 1 && image[i+1][j+1] == 1) output += '\\\\';
                else  output += '-';
                */
            }
            else if (image[i][j] == 0 || image[i][j] % 2 == 0) output += ' ';
            //else output += image[i][j];
        }
        output += '\n';
    }
    size_t ui_el_count = ui_els.size();
    for (size_t u = 0; u < ui_el_count; ++u) {
        size_t uh = ui_els[u].second.size(), uw = ui_els[u].second[0].size();
        size_t uy = ui_els[u].first[0], ux = ui_els[u].first[1];
        for (size_t i = 0; i < uh; ++i) {
            for (size_t j = 0; j < uw; ++j) {
                output[((IMAGE_WIDTH + 1) * (uy)) + (ux) // ux and uy
                     + ((IMAGE_WIDTH + 1) * (i)) + (j) - 1] = ui_els[u].second[i][j];
            }
        }
    }
    for (int i = 0; i < ammo; ++i) {
        output[((IMAGE_WIDTH + 1) * (IMAGE_HEIGHT-2)) + (IMAGE_WIDTH-17) // ux and uy
                     + (i) - 1] = '|';
    }
    if (reloading) {
        for (int i = 0; i < AMMO_MAX; ++i) {
            output[((IMAGE_WIDTH + 1) * (IMAGE_HEIGHT-2)) + (IMAGE_WIDTH-17) // ux and uy
                     + (i) - 1] = '-';
        }
    }
    cout << "\x1b[H" << output;
}

bool isKeyDown(int k) {
    return GetAsyncKeyState(k) & 0b1000000000000000;
}

void draw_sprite(vector<int> base, vector<vector<int>>& sprite, vector<vector<int>>& canvas, int sprite_id) {

    int x = base[0], y = base[1], z = base[2];

    if (z <= 0) return;
    double scale = 0.25 / ((double)z / 12); // 5 constant
    //if (scale < 0.1) return;

    // add scaled h and w
    size_t sprite_h = sprite.size(), sprite_w = sprite[0].size();
    
    size_t scaled_h = sprite_h * scale, scaled_w = sprite_w * scale;
    int key = (sprite_id)*2;
    for (size_t i = 0; i < scaled_h; ++i) {
        for (size_t j = 0; j < scaled_w; ++j) {
            int og_i = i / scale, og_j = j / scale;
            int py = y - scaled_h / 2 + i, px = x - scaled_w / 2 + j;
            if (!offscreen(px,py,canvas)) {
                canvas[py][px] = sprite[og_i][og_j] == 0 ? 0 : sprite[og_i][og_j] + key;
            }
        }
    }
}

class Plane { // atm quads with vertical sides
private:
    vector<vector<double>> vec;
    size_t s = 0;
    vector<vector<int>> sprite_image;
public:
    bool is_sprite;
    int sprite_id;
    bool invis;
    Plane(const vector<vector<double>>& v, bool is_sprite=false, vector<vector<int>> sprite_image={{}}, int sprite_id=sprites_created, bool invis=false) : vec{v}, s{vec.size()}, is_sprite{is_sprite}, sprite_image{sprite_image}, sprite_id{sprite_id}, invis{invis} {sprites_created += is_sprite;};
    void draw() {
        vector<vector<int>> plane_points;
        vector<vector<int>> plane_canvas {IMAGE_HEIGHT, vector<int>(IMAGE_WIDTH)};
        
        if(is_sprite) {
            if (!invis) {
                draw_sprite(project(three_to_four(vec[0]), camera_pos, camera_rot), sprite_image, plane_canvas, sprite_id);
                add_canvas(image, plane_canvas);
            }
            return;
        }

        for (size_t i = 0; i < s; ++i) {
            plane_points.push_back(project(vec[i], camera_pos, camera_rot));
        }
        for (size_t i = 0; i < s-1; ++i) {
            if (plane_points[i].empty() || plane_points[i+1].empty())
                continue; // skip this edge
            draw_line(plane_points[i], plane_points[i+1], plane_canvas);
        }
        if (!plane_points[s-1].empty() && !plane_points[0].empty())
            draw_line(plane_points[s-1], plane_points[0], plane_canvas);

        // seems like fill_poly is being applied within each poly but those polys are not occluding on the image canvas
        fill_poly(plane_canvas);    
        add_canvas(image, plane_canvas);
    }
    
    int compare(const vector<double>& point) { // 1 -> point is on + side of plane, -1 on - side, 0 on plane
        // assume length >= 3
        if (is_sprite) return 1;
        vector<double> a = vec[0], b = vec[1], c = vec[2]; 
        vector<double> ab = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        vector<double> ac = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        vector<double> cross = {ab[1]*ac[2]-ab[2]*ac[1], ab[2]*ac[0]-ab[0]*ac[2], ab[0]*ac[1]-ab[1]*ac[0]};
    

        double sum = cross[0]*(point[0]-a[0]) + cross[1]*(point[1]-a[1]) + cross[2]*(point[2]-a[2]);
        const double EPS = 1e-9;
        return  (sum > EPS) ? 1 : 
                (sum < EPS) ? -1 :
                0;
    }

    int compare_camera(const vector<double>& c_pos, const vector<double>& c_rot) { // 1 -> point is on + side of plane, -1 on - side, 0 on plane
        if (is_sprite) {
            return 1;
            //vector<double> cam = matrix_mult(matrix_mult(inverse(camera_rot_m(camera_rot)), translate_m_inv(c_pos[0], c_pos[1], c_pos[2])), vec[0]);
            //double z = cam[2];
            //return (z < 0 ? 1 : -1);
        }
        vector<double> a = vec[0], b = vec[1], c = vec[2]; 
        vector<double> ab = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        vector<double> ac = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        vector<double> cross = {ab[1]*ac[2]-ab[2]*ac[1], ab[2]*ac[0]-ab[0]*ac[2], ab[0]*ac[1]-ab[1]*ac[0]};
        vector<double> camv = {c_pos[0] - a[0], c_pos[1] - a[1], c_pos[2] - a[2]};
        double dot = cross[0]*camv[0] + cross[1]*camv[1] + cross[2]*camv[2];

        const double EPS = 1e-9;

        return (dot > EPS) ? 1 :
            (dot < -EPS ? -1 : 0);
    }

    vector<double> point() {
        double x=0, y=0, z=0;
        for (auto& p : vec) {
            x += p[0];
            y += p[1];
            z += p[2];
        }
        double s = vec.size();
        return {x/s, y/s, z/s};
        //return vec[0];
    }

    bool operator==(const Plane& other) const {
        bool same = true;
        if (s != other.s) return false;
        for (size_t i = 0; i < s; ++i) {
            if (vec[i] != other.vec[i]) return false;
        }
        return same;
    }
    
};

vector<vector<double>> translate(vector<vector<double>>& points, vector<double> trans) {
    size_t s = points.size();
    vector<vector<double>> out;
    for (size_t i = 0; i < s; ++i) {
        vector<double> ttf = three_to_four(points[i]);
        out.push_back(matrix_mult(translate_m(trans[0], trans[1], trans[2]), ttf));
    }
    return out;
}

struct BSP_node {
    Plane* val;
    BSP_node* lc;
    BSP_node* rc;
    BSP_node(Plane* p, BSP_node* left, BSP_node* right) : val{p}, lc{left}, rc{right} {};
};

// each call returns head (random pick out of list), and list of planes still need to be added
BSP_node* create_BSP_tree(vector<Plane*>& planes) {
    size_t s = planes.size();
    if (s == 0) return nullptr;

    Plane* root = planes[0];
    // plane chosen as partition is planes[0]
    vector<Plane*> left_tree, right_tree;
    for (size_t i = 1; i < s; ++i) {
        if (root->compare(planes[i]->point()) < 0)
            left_tree.push_back(planes[i]);
        else
            right_tree.push_back(planes[i]);
    }
    

    return new BSP_node(root, create_BSP_tree(left_tree), create_BSP_tree(right_tree)); // choose a plane from which to divide
}

void draw_painter(BSP_node* plane, const vector<double>& c_pos, const vector<double>& c_rot) {
    if (plane == nullptr) return;
    
    
    
    if (plane->val->compare(c_pos) > 0) {
        draw_painter(plane->lc, c_pos, c_rot);
        //if (plane->val.compare_camera(c_pos, c_rot) > 0)
        plane->val->draw();
        draw_painter(plane->rc, c_pos, c_rot);
    } else {
        //if (plane->val->is_sprite) {
        //    plane->val->draw();
        //return;
        //}
        draw_painter(plane->rc, c_pos, c_rot);
        //if (plane->val.compare_camera(c_pos, c_rot) > 0)
        plane->val->draw();
        draw_painter(plane->lc, c_pos, c_rot);
    }        
}

void waitOnMenu() {
    while(1) {
        
        if (isKeyDown(VK_SPACE))
            return;
    }
}

int main() {    
    cout << "\x1b[2J";
    cout << "\x1b[?25l";

    vector<vector<int>> cube_points;
    camera_pos = {0.0, 0.0, 0.0};
    camera_rot = {0.0, 0.0, 0.0};


    vector<Plane*> planes;

    vector<vector<double>> cube_translated = translate(cube, {-4.0, 0.0, 3.0});
        
    /*
    planes.push_back(Plane({cube_translated[3], cube_translated[2], cube_translated[1], cube_translated[0]}));
    planes.push_back(Plane({cube_translated[4], cube_translated[7], cube_translated[3], cube_translated[0]}));
    planes.push_back(Plane({cube_translated[5], cube_translated[4], cube_translated[0], cube_translated[1]}));
    planes.push_back(Plane({cube_translated[6], cube_translated[5], cube_translated[1], cube_translated[2]}));
    planes.push_back(Plane({cube_translated[7], cube_translated[6], cube_translated[2], cube_translated[3]}));
    planes.push_back(Plane({cube_translated[4], cube_translated[5], cube_translated[6], cube_translated[7]}));
    */
    
    
    planes.push_back(new Plane({cube_translated[0], cube_translated[1], cube_translated[2], cube_translated[3]}));
    planes.push_back(new Plane({cube_translated[4], cube_translated[7], cube_translated[6], cube_translated[5]}));
    planes.push_back(new Plane({cube_translated[0], cube_translated[3], cube_translated[7], cube_translated[4]}));
    planes.push_back(new Plane({cube_translated[0], cube_translated[4], cube_translated[5], cube_translated[1]}));
    planes.push_back(new Plane({cube_translated[1], cube_translated[2], cube_translated[6], cube_translated[5]}));
    planes.push_back(new Plane({cube_translated[3], cube_translated[7], cube_translated[6], cube_translated[2]}));
    
    
    vector<vector<double>> cube1_translated = translate(cube, {4.0, 0.0, 3.0});
    planes.push_back(new Plane({cube1_translated[0], cube1_translated[1], cube1_translated[2], cube1_translated[3]}));
    planes.push_back(new Plane({cube1_translated[4], cube1_translated[7], cube1_translated[6], cube1_translated[5]}));
    planes.push_back(new Plane({cube1_translated[0], cube1_translated[3], cube1_translated[7], cube1_translated[4]}));
    planes.push_back(new Plane({cube1_translated[0], cube1_translated[4], cube1_translated[5], cube1_translated[1]}));
    planes.push_back(new Plane({cube1_translated[1], cube1_translated[2], cube1_translated[6], cube1_translated[5]}));
    planes.push_back(new Plane({cube1_translated[3], cube1_translated[7], cube1_translated[6], cube1_translated[2]}));
    

    planes.push_back(new Plane({{-1.0, 0.0, 2.8}}, true, enemy, 0, false));
    planes.push_back(new Plane({{0.0, 0.0, 2.8}}, true, enemy, 1, false));
    planes.push_back(new Plane({{1.0, 0.0, 2.8}}, true, enemy, 2, false));

    BSP_node* planes_painter = create_BSP_tree(planes);

    int s = 0;
    for (int i = 0; i < 1000000; ++i) {
        //int t = i % 6;
        //if (i % 6 == 0) s = (s+1) % 4;

        if (isKeyDown(VK_LEFT))
            camera_rot[1] += 0.03;
        if (isKeyDown(VK_RIGHT))
            camera_rot[1] -= 0.03;
        if (isKeyDown(0x41)) { // a
            camera_pos[0] -= cos(camera_rot[1])*0.05;
            camera_pos[2] -= sin(camera_rot[1])*0.05;
        }
        if (isKeyDown(0x44)) { // d
            camera_pos[0] += cos(camera_rot[1])*0.05;
            camera_pos[2] += sin(camera_rot[1])*0.05;
        }
            
        if (isKeyDown(0x57)) { // w
            camera_pos[0] -= sin(camera_rot[1])*0.05;
            camera_pos[2] += cos(camera_rot[1])*0.05;
        }
        if (isKeyDown(0x53)) {// s
            camera_pos[0] += sin(camera_rot[1])*0.05;
            camera_pos[2] -= cos(camera_rot[1])*0.05;
        }
        
        
        if (isKeyDown(VK_SPACE)) {
            if (ammo > 0) {
                int hit_pixel = image[IMAGE_HEIGHT/2-1][IMAGE_WIDTH/2];
                if (hit_pixel > 1) {
                    size_t s = planes.size();
                    //cout << hit_pixel << ", " << sprites_created << endl;
                    for (size_t i = 0; i < s; ++i) {
                        if (planes[i]->is_sprite && planes[i]->sprite_id == hit_pixel / 2 - 1) {
                            //cout << i << ", " << planes[i]->invis << endl;
                            planes[i]->invis = true;
                        }
                    }
                }
                --ammo;
            }
            if (ammo <= 0) {
                reloading = true;
                reload_counter = RELOAD_DUR;
            }
        }

        if (isKeyDown(0x52)) { // r
            reloading = true;
            ammo = 0;
            reload_counter = RELOAD_DUR;
        }

        if (reloading) {
            --reload_counter;
            if (reload_counter <= 0) {
                ammo = AMMO_MAX;
                reloading = false;
            }
        }
        
        for (auto &row : image) fill(row.begin(), row.end(), 0);
        draw_painter(planes_painter, camera_pos, camera_rot);
        render(image, ui_elements);
        Sleep(11.11); // 90 hz max
    }
        
    return 0;
}