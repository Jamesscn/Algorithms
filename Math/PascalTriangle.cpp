#include <vector>

using namespace std;

/**
 * Constructs a pascal triangle of a given size.
 * @param size The maximum size of the triangle.
 * @return A 2D vector of long long ints with each of the coefficients of the triangle.
 */
vector<vector<long long int>> PascalTriangle(int size) {
	vector<vector<long long int>> triangle = {{1}};
	for(int y = 1; y < size; y++) {
		vector<long long int> row;
		row.push_back(1);
		for(int x = 1; x < y; x++) {
			row.push_back(triangle[y - 1][x - 1] + triangle[y - 1][x]);
		}
		row.push_back(1);
		triangle.push_back(row);
	}
	return triangle;
}