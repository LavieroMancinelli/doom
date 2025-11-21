#include <iostream>
#include <vector>
#include <string>
using namespace std; 

vector<vector<int>> image = {   {0,0,1,0,0},
                                {0,1,0,1,0},
                                {1,0,0,0,1},
                                {0,1,0,1,0},
                                {0,0,1,0,0}
                            };

void render(vector<vector<int>> image) {
    size_t x_size = image.size(), y_size = image[0].size();
    string output = "";
    for (size_t i = 0; i < y_size; ++i) {
        for (size_t j = 0; j < x_size; ++j) {
            if (image[i][j] == 1) {
                if (i == 0 || i == y_size-1) output += '-';
                else if (j == 0 || j == x_size-1) output += '|';
                else if (image[i-1][j+1] == 1 && image[i+1][j-1] == 1) output += '/';
                else if (image[i-1][j-1] == 1 && image[i+1][j+1] == 1) output += '\\\\';
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


    render(image);
    return 0;
}